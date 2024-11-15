/**
 * @file        Cmaes_MPI_alltoall.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with CMA-ES \n
 *              using fully connected all-to-all island population model.
 * 
 * @version     0.1
 * 
 * @date        2021-01-31 (created) \n
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>

#include <random>
#include <vector>
#include <chrono>
#include <functional>

#include <mpi.h>

#include "LaunchParams.h"
#include "PrintFunctions.h"

#include "CMA-ESpp/cma-es/cmaes.h"
#include "FitnessFunc/FitnessFunction.h"
#include "LoggingFuncs/LoggingHooks_CMAESpp.h"

#ifndef HIFU_PROBLEM
  #define HIFU_PROBLEM  
#endif


/**************************************************************************************
 * GLOBALS
 * Launch parameters
***************************************************************************************/
constexpr int   MAX_SECONDS       = 28800;    /// Maximum seconds to run
constexpr int   MAX_FITEVALS      = 150000;   /// Maximum allowed fitness evaluations
constexpr int   MAX_GEN           = 50000;    /// Maximum allowed generations
constexpr float OPT_SCORE         = 0.0f;     /// Optimal score 
constexpr float EPS               = 0.001f;   /// Acceptable leeway from optimal score
constexpr int   SON_NUM           = 6;        /// The amount of sonications to perform
constexpr int   MIG_NUM           = 2;        /// The amount of immigrants to send to broadcast from one node
constexpr int   MIG_INT           = 10;       /// How often (in generatiosn) should the migration occur
constexpr int   TRNMNT_SIZE       = 3;        /// The size of the tournament to select accepted immigrants
// END OF GLOBALS DECLARATION *********************************************************


/**
 * @brief Calculates the overall score of a given individual.
 *
 * @tparam T                      - data type of individual genomes and also a return type.
 * @param [in, out] x             - the individual to evaluate in an array.
 * @param [in, out] p             - the general optimization and multicpu distribution \n
 *                                  parameters.
 * @return T                      - the score given to the evaluated individual.
 */
template<typename T> 
T scoreFunction(T *x,
                LaunchParams & p)
{

  std::vector<float> convertedValues(p.problemDimension);
  int penalize = 0;

  //covert the normalized varibale to the original axis range
  for (int i = 0; i < p.problemDimension; i++)
  {
    if (x[i] < 0.0 || x[i] > 10.0)
    {
      penalize++;
    }

    float scaleFactor = (1.0 - std::cos(M_PI * x[i] * 0.1f)) * 0.5f;
    convertedValues[i] = p.lowerBounds[i] + (p.upperBounds[i] - p.lowerBounds[i]) * scaleFactor;
  }

  T score = fitnessFunction(p.problemDimension, convertedValues.data());
  return score + ((float)penalize/(float)p.problemDimension * score);
}
// END OF scoreFunction **************************************************************

/**
 * @brief A roulette wheel selection to choose the surviving individual.
 *
 * @tparam T                      - data type of individual genomes and also a return type.
 * @param [in, out] arFunvals     - a vector of fitness values of competing individuals.
 * @param [in, out] p             - the general optimization and multicpu distribution \n
 *                                  parameters.
 * @return int                    - index of the selected individual
 */
template <class T> 
int rouletteWheelSelectIndex(std::vector<T> & arFunvals,
                             LaunchParams & p)
{
  const int rouletteSize = arFunvals.size();
  int victorIndex = 0;
  T totalFitness = std::accumulate(arFunvals.begin(), arFunvals.end(), T(0.0));

  std::vector<float> rouletteProbs(rouletteSize);
  std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
  for (int i = 0; i < rouletteSize; i++)
  {
    rouletteProbs[i] = 1.0f - arFunvals[i]/totalFitness;
  }

  float rouletteTick = distribution(p.rngEngine);
  float offset = 0.0f;
  for (int i = 0; i < rouletteSize; i++)
  {
    offset += rouletteProbs[i];
    if (rouletteTick < offset)
    {
      return i;
    }
  }
  return 0;
}
// END OF rouletteWheelSelectIndex **************************************************


/**
 * @brief A tournament selection to choose the surviving individual.
 *
 * @tparam T                      - data type of individual genomes.
 * @param [in, out] arFunvals     - a vector of fitness values of competing individuals.
 * @param [in, out] p             - the general optimization and multicpu distribution \n
 *                                  parameters.
 * @return int                    - index of the selected individual
 */
template <class T>
int tournamentSelectIndex(std::vector<T> & arFunvals,
                          LaunchParams & p,
                          std::function<bool(T, T)> cmpOp)
{
  assert(p.tournamentSize > 1 && arFunvals.size() > p.tournamentSize && arFunvals.size() >= 2 && "Tournament selection error.");
  std::uniform_int_distribution<int> distribution(0, arFunvals.size()-1);
  int victorIndex = distribution(p.rngEngine);
  T minFitness = arFunvals[victorIndex];

  for (int i = 1; i < p.tournamentSize; i++)
  {
    int parentIndex = distribution(p.rngEngine);
    while (parentIndex == victorIndex)
    {
      parentIndex = distribution(p.rngEngine);
    }
    if (cmpOp(arFunvals[parentIndex], minFitness))
    {
      victorIndex = parentIndex;
      minFitness = arFunvals[parentIndex];
    }
  }
  return victorIndex;
}
// END OF tournamentSelectIndex **************************************************

/**
 * @brief Packs the n best individuals in population into the Emmigrants data structure.
 *
 * @tparam T                        - data type of individual genomes.
 * @param [in, out] population      - an array of the local population.
 * @param [in, out] arFunvals       - a vector of fitness values of the population.
 * @param [in] sortedFunvalsIndexes - a vector of indexes sorted by their individual fitness. \n
 *                                  - 0 is the best individual.
 * @param [in, out] p               - the general optimization and multicpu distribution \n
 *                                    parameters.
 * @return Emmigrants<T>            - a object containing the local emmigrating individuals
 */
template <class T> 
Emigrants<T>* prepareDataForEmigration(T *const * population,
                                      std::vector<T> & arFunvals,
                                      const std::vector<size_t> & sortedFunvalsIndexes,
                                      LaunchParams & p)
{
  Emigrants<T>* leavingPop = new Emigrants<T>(p.numberOfImmigrants, p.problemDimension);

  for (int i = 0; i < p.numberOfImmigrants; i++)
  {
    size_t minIndex = sortedFunvalsIndexes[i];
    //copy n-1 random tournament selected individuals
    //int victorIndex = tournamentSelectIndex<T>(arFunvals, p, [](T a, T b) -> bool {return a < b;});
    leavingPop->fitness[i] = arFunvals[minIndex];
    std::copy(population[minIndex],
              population[minIndex] + p.problemDimension,
              leavingPop->population.begin() + i * p.problemDimension);
  }
  return leavingPop;
}
//END of prepareDataForEmigration

/**
 * @brief Broadcasts local emmigrants to all other islands and recieves immigrants from others \n
 *        in the Immigrants data structure. Send and recieve are done in a non-blocking way, however \n
 *        we wait for data to arrive. Thus the islands are synchronized on every migrations.
 *
 * @tparam T                        - data type of individual genomes.
 * @param [in, out] population      - an array of the local population.
 * @param [in, out] arFunvals       - a vector of fitness values of the population.
 * @param [in] sortedFunvalsIndexes - a vector of indexes sorted by their individual fitness. \n
 *                                  - 0 is the best individual.
 * @param [in, out] p               - the general optimization and multicpu distribution \n
 *                                    parameters.
 * @return Immigrants<T>            - a object containing the immigrating individuals from all islands
 */
template <class T> 
Immigrants<T>* broadcastEmigrateAndImmigrate(T *const * population,
                                            std::vector<T> & arFunvals,
                                            const std::vector<size_t> & sortedFunvalsIndexes,
                                            LaunchParams & p)
{
  constexpr int requestsNumber = 2;
  MPI_Request commRequests[requestsNumber];

  Emigrants<T>*  leavingPop = prepareDataForEmigration(population, arFunvals, sortedFunvalsIndexes, p);
  Immigrants<T>* incomingPop = new Immigrants<T>(p.numberOfImmigrants * p.islandCount, p.problemDimension);

  MPI_Iallgather(leavingPop->fitness.data(), 
                 p.numberOfImmigrants,
                 MPI_FLOAT,
                 incomingPop->fitness.data(),
                 p.numberOfImmigrants,
                 MPI_FLOAT,
                 MPI_COMM_WORLD,
                 &commRequests[0]);

  MPI_Iallgather(leavingPop->population.data(), 
                 p.numberOfImmigrants * p.problemDimension,
                 MPI_FLOAT,
                 incomingPop->population.data(),
                 p.numberOfImmigrants * p.problemDimension,
                 MPI_FLOAT,
                 MPI_COMM_WORLD,
                 &commRequests[1]);

  MPI_Status status[requestsNumber];
  MPI_Waitall(requestsNumber, commRequests, status);
  delete leavingPop;
  return incomingPop;
}
//END of broadcastEmigrateAndImmigrate

/**
 * @brief Decides whether or not to accept the incoming immigrants based on the fitness of local \n
 *        population. For HIFU problems, decides based on a prob. table. For others, the decision \n
 *        is based on ration between current best and worst of all bests.
 *
 * @tparam T                        - data type of individual genomes.
 * @param [in, out] bestFit         - best fitness value in current population.
 * @param [in, out] worstBestFit    - worst of the best fitnesses accross all generations. Usually \n
 *                                    the best fitness in the first generation.
 * @param [in, out] p               - the general optimization and multicpu distribution \n
 *                                    parameters.
 * @return bool                     - whether or not to accept the immigrants into local population.
 */
template <class T> 
bool immigrantAccept(const T & bestFit,
                     const T & worstBestFit,
                     LaunchParams & p)
{
  #ifdef HIFU_PROBLEM
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    if (bestFit < T(30))
    {
      return false;
    }
    else if (bestFit < T(50))
    {
      return dist(p.rngEngine) < 0.05f;
    }
    else if (bestFit < T(100))
    {
      return dist(p.rngEngine) < 0.1f;
    }
    else if (bestFit < T(200))
    {
      return dist(p.rngEngine) < 0.4f;
    }
    else
    {
      return dist(p.rngEngine) < 0.75f;
    }
  #else
    return float(bestFit)/float(worstBestFit);
  #endif
}
//END of immigrantAccept


/**
 * @brief Returns sorted indexes based on the fitness values of their individualas.
 *
 * @tparam T                        - data type of individual genomes.
 * @param [in, out] values          - fitness values of the population.
 * @return vector<T>                - vector of sorted indexes.
 */
template <class T> 
std::vector<size_t> indexSort(const std::vector<T> & values)
{
  std::vector<size_t> sortedIndexes(values.size());
  std::iota(sortedIndexes.begin(), sortedIndexes.end(), 0);
  std::stable_sort(sortedIndexes.begin(), sortedIndexes.end(), [&values](size_t index1, size_t index2) { return values[index1] < values[index2]; });
  return sortedIndexes;
}
//END of indexSort

/**
 * @brief Runs the CMAES evolution in fully connected all-to-all islands model.
 *
 * @tparam T                      - data type of genomes of the individuals.
 * @param [in, out]    evoParams  - the cmaes parameters
 * @param [in, out]    arFunvals  - vector to contain fitness values
 * @param [in, out]    p          - the general optimization and multicpu distribution parameters
 * @return CMAES<T>               - all results of the optimization are saved and encoded in the \n
 *                                  final position of particles in the CMAES class
 */
template <class T> 
CMAES<T>* runCMAES(Parameters<T> & evoParams,
                   std::vector<T> & arFunvals,
                   LaunchParams & p)
{
  T *const *pop         = nullptr;
  CMAES<T>* evo         = nullptr;
  T startSigma          = 0;

  int finished          = 0;

  T* bestEverChromo = new T[p.problemDimension];
  T bestEverFitness = std::numeric_limits<float>::max();
  T worstBestFitness = std::numeric_limits<float>::max();

  std::vector<size_t> sortedFunvalsIndexes;
  
  evo = new CMAES<T>();
  evo->init(evoParams, p.seed); // freed by cmaes destructor
  startSigma = evo->get(CMAES<T>::Sigma);
  
  //p.popSize = int(evo->get(CMAES<T>::Lambda));
  
  // Iterate until stop criteria holds
  while (p.generationsNeeded < p.generations && finished == 0)
  {
    p.generationsNeeded++;
    // Generate lambda new search points, sample population
    pop = evo->samplePopulation();
    
    for (int i = 0; i < p.popSize; i++)
    {
      arFunvals[i] = scoreFunction<T>(pop[i], p);
      sortedFunvalsIndexes = indexSort<T>(arFunvals);
      p.fitnessCalls++;
    }
    if (worstBestFitness == std::numeric_limits<float>::max())
    {
      worstBestFitness = arFunvals[sortedFunvalsIndexes[0]];
    }

    // check if migration is about to occur.
    if (p.generationsNeeded % p.migrationInterval == 0 && p.numberOfImmigrants > 0)
    {
      //no need to sync as long as the size of population is the same
      //if it isnt, we are in for it in a bad way
      Immigrants<T> *outsiders = broadcastEmigrateAndImmigrate<T>(pop, arFunvals, sortedFunvalsIndexes, p);

      if (immigrantAccept<T>(arFunvals[sortedFunvalsIndexes[0]], worstBestFitness, p))
      {

        for (int i = 0; i < p.numberOfImmigrants; i++)
        {
          size_t indexToReplace = sortedFunvalsIndexes[sortedFunvalsIndexes.size() - 1];
          int victorIndex = tournamentSelectIndex<T>(outsiders->fitness, p, [](T val1, T val2){ return val1 < val2; });
          arFunvals[indexToReplace] = outsiders->fitness[victorIndex];
          evo->replaceIndividual(indexToReplace, outsiders->population.data() + victorIndex * p.problemDimension);
        }
      }
      delete outsiders;
    }
    evo->updateDistribution(arFunvals.data());
    T newSigma = evo->get(CMAES<T>::Sigma);

    //the immigrants broke the cma, we need to restart. Should not be happening anymore but just in case.
    if (!std::isfinite(newSigma) || newSigma > 2 * startSigma)  
    {
      // save some state variables before reset
      T* mean         = evo->getNew(CMAES<T>::XMean);
      T bestSoFar     = evo->get(CMAES<T>::FBestEver);
      if (bestSoFar < bestEverFitness)
      {
        bestEverFitness = bestSoFar;
        evo->getInto(CMAES<T>::XBestEver, bestEverChromo);
      }
      
      // create a new instance
      delete evo;
      evo = new CMAES<T>();

      // evo instance contains reference to our Parameters instance, dtor destroyed it
      evoParams = Parameters<float>();  // need to construct and reinitilize
      evoParams.init(p.problemDimension, p.xStart.data(), p.stdDevs.data());
      evoParams.stopMaxIter = p.generations;

      //reinit evo
      evo->init(evoParams, p.seed);
      evo->setMean(mean); //set the mean to the last one so we can continue where we were before reset
      evo->setBestEver(bestEverFitness, bestEverChromo, p.fitnessCalls);

      delete[] mean;
    }

    if (p.verboseOutput)
    {
      printProgress<T>(*evo, p); fflush(p.logFile);
    }

    bool limit      = reachedLimit(p.fitnessCalls, p.fitEvalsLimit);
    bool converged  = isConverged(p.optScore, evo->get(CMAES<T>::Fitness), p.eps);
    bool timesUp    = int(MPI_Wtime() - p.startTime.d)  > p.maxSeconds; //seconds since start
    bool optEnded   = p.generationsNeeded > p.generations || limit || converged || timesUp;

    if (optEnded) // if we ended set the variable
    {
      finished = 1;
    }

    // reduction to check whether any island converged or ended 
    int finishedReduced = 0;
    MPI_Allreduce(
      &finished,
      &finishedReduced,
      1,
      MPI_INT,
      MPI_MAX,
      MPI_COMM_WORLD
    );

    if (finishedReduced > 0)
    {
      break;
    }
  }
  delete[] bestEverChromo;
  return evo;
}
// END OF runCMAES *************************************************************


/**
 * @brief the main funcion of CMAES island all-to-all populations model \n
 *        for HIFU problems.
 * 
 * @param [in]      argc      - a number of arguments passed to optimizer
 * @param [in]      argv      - argumenst passed to optimizer
 * @return int                - the return value of the program
 */
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  if (argc != 3)
  {
    fprintf(stderr, "Invalid arguments.\n");
    printUsage();
    MPI_Finalize();
    return 1;
  }

  int masterSeed = strtol(argv[1], nullptr, 10);
  std::string outDir(argv[2]);
  if (outDir.back() != '/')
  {
    outDir += "/";
  }

  std::default_random_engine seedGenerator(masterSeed);

  int rank;
  int numProcs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  std::vector<int> procSeeds(numProcs);
  for (int i = 0; i < numProcs; i++)
  {
    procSeeds[i] = seedGenerator();
  }
  
  LaunchParams p(procSeeds[rank]);
   
  p.islandId            = rank;
  p.islandCount         = numProcs;

  p.numberOfImmigrants  = MIG_NUM;
  p.migrationInterval   = MIG_INT;

  p.optScore            = OPT_SCORE;
  p.fitEvalsLimit       = MAX_FITEVALS;
  p.maxSeconds          = MAX_SECONDS;
  p.eps                 = EPS;
  p.generations         = MAX_GEN;
  p.problemDimension    = 4 * SON_NUM;
  p.tournamentSize      = TRNMNT_SIZE;

  p.fitnessCalls        = 0;
  p.generationsNeeded   = 0;

  p.verboseOutput       = true;

  std::string outFileName = outDir + std::string("node") + std::to_string((int)p.islandId) + std::string("_log.optOut");
  p.logFile = fopen(outFileName.c_str(), "w+");
  if (p.logFile == nullptr)
  {
    fprintf(stderr, "failed to open file\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  p.xStart      = std::vector<float>(p.problemDimension);
  p.stdDevs     = std::vector<float>(p.problemDimension);
  p.upperBounds = std::vector<float>(p.problemDimension);
  p.lowerBounds = std::vector<float>(p.problemDimension);

  getConstraints(p.problemDimension, p.upperBounds, p.lowerBounds);

  //p.rngEngine(p.seed);
  std::uniform_real_distribution<float> distribution(0.1f, 10.0f);

  /* Here you may resample each solution point pop[i] until it
    becomes feasible, e.g. for box constraints (variable
    boundaries).
    Assumptions: the feasible domain is convex, the optimum is
    not on (or very close to) the domain boundary, initialX is
    feasible and initialStandardDeviations are sufficiently small
    to prevent quasi-infinite looping.
  */

  //normalizing the range of every axis to 0 to 10
  //later conversion by using the formula:
  //lowerBound[i] + (upperBound[i] - lowerBound[i]) * (1.0 - std::cos(M_PI * x[i]/10.0))/2.0;
  for(int i = 0; i < p.problemDimension; i++)
  {
    p.xStart[i] = distribution(p.rngEngine);              // start at random position between bounds
    p.stdDevs[i] = 10.0f/3.0f;
  }

  Parameters<float> params;
  params.init(p.problemDimension, p.xStart.data(), p.stdDevs.data());
  params.stopMaxIter = p.generations;

  p.popSize = params.lambda;
  printHeader<float>(p, masterSeed); 
  fprintf(p.logFile, "\n@@@\n"); fflush(p.logFile); // separator of launch and run section

  std::vector<float> arFunVals(params.lambda);
  p.startTime.d  = MPI_Wtime();
  CMAES<float>* evo = runCMAES<float>(params, arFunVals, p);
  double endTime    = MPI_Wtime();
  fprintf(p.logFile, "\n@@@\n"); fflush(p.logFile); // separator of run and results section
  
  double elapsed = (endTime - p.startTime.d) * 1000000;
  float best     = evo->get(CMAES<float>::FBestEver);
  float* chromo  = evo->getNew(CMAES<float>::XBestEver);

  std::vector<float> convertedValues(p.problemDimension);
  //covert the normalized varibale to the original axis range
  for (int i = 0; i < p.problemDimension; i++)
  {
    float scaleFactor = (1.0 - std::cos(M_PI * chromo[i] * 0.1f)) * 0.5f;
    convertedValues[i] = p.lowerBounds[i] + (p.upperBounds[i] - p.lowerBounds[i]) * scaleFactor;
  }
  printFooter<float>(p, best, convertedValues.data(), (unsigned long long)elapsed); fflush(p.logFile);
  fclose(p.logFile);

  delete evo;

  MPI_Finalize();
  return 0;
}
//END of main