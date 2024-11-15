/**
 * @file        Cmaes_MPI_PP.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for solving real value minimization problems with CMA-ES \n
 *              using Farmer-Worker parallelization model
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


/**************************************************************************************
 * GLOBALS
 * Launch parameters
***************************************************************************************/
constexpr int   PP_ISLAND_PARENT  = 0;          /// Farmer island rank
constexpr int   MAX_SECONDS       = 28800;      /// Maximum seconds to run
constexpr int   MAX_FITEVALS      = 150000;     /// Maximum allowed fitness evaluations
constexpr int   MAX_GEN           = 50000;      /// Maximum allowed generations
constexpr float OPT_SCORE         = 0.0f;       /// Optimal score 
constexpr float EPS               = 0.001f;     /// Acceptable leeway from optimal score
constexpr int   SON_NUM           = 6;          /// The amount of sonications to perform
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
  return score + ((float)penalize/(float)p.problemDimension) * score);
}
// END OF scoreFunction **************************************************************


/**
 * @brief Gets part of the population to evaluate. Called by both farmer and workers. \n
 *        Distribution carried out by blocking scatter function call.
 *
 * @tparam T                      - data type of individual genomes and also a return type.
 * @param [in] population         - the population to distribute across all workers. \n
 *                                  Supplied by the farmer, workers are free to pass  \n
 *                                  nullptr.
 * @param [in, out] p             - the general optimization and multicpu distribution \n
 *                                  parameters.
 * @return vector<T>              - the part of the population this node needs to evaluate.
 */
template <class T>
std::vector<T> getPopulationPart(T *const *population, 
                                 LaunchParams & p)
{
  const int chunkSize = p.popSizePerIsland * p.problemDimension; 
  std::vector<T> flatPopulation;
  if (p.islandId == PP_ISLAND_PARENT)
  {
    flatPopulation = std::vector<T>(p.popSize * p.problemDimension);
    for (int i = 0; i < p.popSize; i++)
    {
      std::copy(population[i], population[i] + p.problemDimension, flatPopulation.begin() + i * p.problemDimension);
    }
  }

  std::vector<T> flatPopulationPart(chunkSize);

  MPI_Scatter(
    flatPopulation.data(),
    chunkSize,
    MPI_FLOAT,
    flatPopulationPart.data(),
    chunkSize,
    MPI_FLOAT,
    PP_ISLAND_PARENT,
    MPI_COMM_WORLD
  );

  return flatPopulationPart;
}
// END OF getPopulationPart **************************************************************


/**
 * @brief A wrapper around scoreFunction call using vector of multiple individuals.
 *
 * @tparam T                      - data type of individual genomes and also a return type.
 * @param [in] flatPop            - multiple individuals to evaluate in a flattened vector.
 * @param [in, out] p             - the general optimization and multicpu distribution \n
 *                                  parameters.
 * @return vector<T>              - the scored given to the individuals. Index to this equates \n
 *                                  to the order of individuals in flatPop.
 */
template <class T>
std::vector<T> scorePopulation(std::vector<T> flatPop, 
                               LaunchParams & p)
{
  std::vector<T> results(p.popSizePerIsland);
  for (int i = 0; i < p.popSizePerIsland; i++)
  {
    results[i] = scoreFunction<T>(flatPop.data() + i*p.problemDimension, p);
  }
  return results;
}
// END OF scorePopulation **************************************************************


/**
 * @brief Sends evaluated results back to farmer using a blocking gather function call.
 *
 * @tparam T                      - data type of individual genomes and also a return type.
 * @param [in] results            - a vector containing fitness values of individuals \n
 *                                  evaluated by this node.
 * @param [in] p                  - the general optimization and multicpu distribution \n
 *                                  parameters.
 * @return vector<T>              - fitness values of the entire population across all nodes. \n
 *                                  Only farmer gets the data.
 */
template <class T>
std::vector<T> sendResultsToRoot(std::vector<T> results, 
                                 LaunchParams & p)
{
  std::vector<T> fitVals(p.popSize);
  MPI_Gather(
    results.data(),
    p.popSizePerIsland,
    MPI_FLOAT,
    fitVals.data(),
    p.popSizePerIsland,
    MPI_FLOAT,
    PP_ISLAND_PARENT,
    MPI_COMM_WORLD
  );
  return fitVals;
}
// END OF sendResultsToRoot **************************************************************

/**
 * @brief Broadcasts the actual population size of the entire population across all nodes.
 *
 * @param [in, out] p                  - the general optimization and multicpu distribution parameters \n
 *                                       The size is stored in this structure.
 * @return int                         - the population size in question.
 */
int distributePopSize(LaunchParams & p)
{
  MPI_Bcast(&p.popSize, 1, MPI_INT, PP_ISLAND_PARENT, MPI_COMM_WORLD);
  p.popSizePerIsland = p.popSize / p.islandCount;
  return p.popSize;
}
// END OF distributePopSize **************************************************************

/**
 * @brief Broadcasts information whether the optimization should continue or not. Only \n
 *        the farmer decides. 
 *
 * @param [in, out] shouldEnd           - wheter or not the optimization should end. \n
 *                                        This information is then broadcasted workers.                                    
 * @return bool                         - should the optimization end.
 */
bool endOpt(bool shouldEnd)
{
  int tmp = (int)shouldEnd;
  MPI_Bcast(&tmp, 1, MPI_INT, PP_ISLAND_PARENT, MPI_COMM_WORLD);
  //bool res = tmp > 0;
  return tmp > 0;
}
// END OF endOpt **************************************************************

/**
 * @brief Wrapper around endOpt for workers. Returns whether the optimization \n
 *        is over or not.
 *                                 
 * @return bool                         - should the optimization end.
 */
bool parentEnded()
{
  return endOpt(false);
}
// END OF parentEnded **********************************************************


/**
 * @brief Runs the CMAES evolution. Distributes work to workers and gahters results. \n
 *        Called by the farmer.
 *
 * @tparam T                      - data type of genomes of the individuals.
 * @param [in, out]    evoParams  - the cmaes parameters
 * @param [in, out]    p          - the general optimization and multicpu distribution parameters
 * @return CMAES<T>               - all results of the optimization are saved and encoded in the \n
 *                                  final position of particles in the CMAES class
 */
template <class T> 
CMAES<T>* runCMAES(Parameters<T> & evoParams,
                   LaunchParams & p)
{
  T *const *pop = nullptr;
  CMAES<T>* evo = nullptr;
  
  std::vector<T> flatPopulationPart;
  std::vector<T> populationResults;
  std::vector<T> popFitVals;
  
  double timeStart, evalTimeAcc = 0, printTimeAcc = 0;

  evo = new CMAES<T>();
  evo->init(evoParams, p.seed); // freed by cmaes destructor

  //#pragma pomp inst begin(mpiSecRoot)
  // Iterate until stop criteria holds
  while (true)
  {
    p.generationsNeeded++;
    
    // Generate lambda new search points, sample population
    pop = evo->samplePopulation();

    
    flatPopulationPart  = getPopulationPart<T>(pop, p);
    
    //timeStart = MPI_Wtime();
    populationResults   = scorePopulation<T>(flatPopulationPart, p);
    //evalTimeAcc += (MPI_Wtime() - timeStart);

    popFitVals          = sendResultsToRoot(populationResults, p);
    

    p.fitnessCalls += p.popSize;

    evo->updateDistribution(popFitVals.data());
    if (p.verboseOutput)
    {
      timeStart = MPI_Wtime();
      printProgress<T>(*evo, p);
      fflush(p.logFile);
      printTimeAcc += (MPI_Wtime() - timeStart);
    }

    bool limit      = reachedLimit(p.fitnessCalls, p.fitEvalsLimit);
    bool converged  = isConverged(p.optScore, evo->get(CMAES<T>::Fitness), p.eps); 
    bool timesUp    = int(MPI_Wtime() - p.startTime.d)  > p.maxSeconds; //seconds since start
    bool optEnded   = p.generationsNeeded > p.generations || limit || converged || timesUp;

    if(endOpt(optEnded))
    {
      //#pragma pomp inst altend(mpiSecRoot)
      break;
    }
    
  }
  //#pragma pomp inst end(mpiSecRoot)
  //fprintf(p.logFile, "%g %g\n", evalTimeAcc, printTimeAcc);
  
  return evo;
}
// END OF runCMAES *************************************************************

/**
 * @brief Runs the farmer/parent part of the farmer-worker distribution model. \n
 *        Drives forward the optimization process. 
 *
 * @tparam T                        - data type of genomes of the individuals.
 * @param [in, out]    p            - the general optimization and multicpu distribution \n
 *                                    parameters.
 * @param [in, out]    outFileName  - name of the file to output progress into.
 * @return int                      - the return code. 0 is success.
 */
template <class T>
int runParent(LaunchParams & p, 
              std::string outFileName)
{
  p.optScore            = OPT_SCORE;
  p.fitEvalsLimit       = MAX_FITEVALS * p.islandCount;
  p.maxSeconds          = MAX_SECONDS;
  p.eps                 = EPS;
  p.generations         = MAX_GEN;
  p.fitnessCalls        = 0;
  p.generationsNeeded   = 0;
  p.verboseOutput       = true;
  p.xStart              = std::vector<T>(p.problemDimension);
  p.stdDevs             = std::vector<T>(p.problemDimension);
  p.logFile             = fopen(outFileName.c_str(), "w+");
  if (p.logFile == nullptr)
  {
    fprintf(stderr, "failed to open file\n");
    return 1;
  }

  //normalizing the range of every axis to 0 to 10
  //later conversion by using the formula:
  //lowerBound[i] + (upperBound[i] - lowerBound[i]) * (1.0 - std::cos(M_PI * x[i]/10.0))/2.0;
  std::uniform_real_distribution<T> distribution(0.1f, 10.0f);
  for(int i = 0; i < p.problemDimension; i++)
  {
    p.xStart[i]   = distribution(p.rngEngine);              // start at random position between bounds
    p.stdDevs[i]  = 10.0f/3.0f;
  }

  Parameters<T> params;
  // lambda MUST be set before calling init
  params.lambda = (4 + (int) (3.0f * std::log((float) p.problemDimension))) * p.islandCount;
  params.init(p.problemDimension, p.xStart.data(), p.stdDevs.data());
  params.stopMaxIter = p.generations;

  p.popSize = params.lambda;
  printHeader<T>(p, p.seed); fprintf(p.logFile, "\n@@@\n"); fflush(p.logFile); // separator of launch and run section
  
  distributePopSize(p);
  p.startTime.d     = MPI_Wtime();
  CMAES<T>* evo     = runCMAES<T>(params, p);
  double endTime    = MPI_Wtime();
  fprintf(p.logFile, "\n@@@\n"); fflush(p.logFile); // separator of run and results section
  
  double elapsed = (endTime - p.startTime.d) * 1000000;
  float best     = evo->get(CMAES<float>::FBestEver);
  float* chromo  = evo->getNew(CMAES<float>::XBestEver);

  std::vector<T> convertedValues(p.problemDimension);
  for (int i = 0; i < p.problemDimension; i++)
  {
    float scaleFactor   = (1.0 - std::cos(M_PI * chromo[i] * 0.1f)) * 0.5f;   
    convertedValues[i]  = p.lowerBounds[i] + (p.upperBounds[i] - p.lowerBounds[i]) * scaleFactor; //covert the normalized variable to the original axis range
  }
  printFooter<T>(p, best, convertedValues.data(), (unsigned long long)elapsed); fflush(p.logFile);
  
  fclose(p.logFile);
  delete chromo;
  delete evo;

  return 0;
}
//END of runParent

/**
 * @brief Runs the worker/child part of the farmer-worker distribution model. \n
 *        Recieves data, evaluate individuals and sends back the results.
 *
 * @tparam T                        - data type of genomes of the individuals.
 * @param [in, out]    p            - the general optimization and multicpu distribution \n
 *                                    parameters.
 * @param [in, out]    outFileName  - name of the file to output progress into.
 * @return int                      - the return code. 0 is success.
 */
template <class T>
void runChild(LaunchParams & p)
{
  distributePopSize(p);

  std::vector<T> flatPopulationPart;
  std::vector<T> populationResults;

  while(true)
  {
    flatPopulationPart  = getPopulationPart<T>(nullptr, p);
    populationResults   = scorePopulation<T>(flatPopulationPart, p);
    sendResultsToRoot(populationResults, p);
    if(parentEnded())
    { 
      break;
    }
  }
}
//END of runChild

/**
 * @brief the main funcion of CMAES Pan Population optimization model for \n
 *        HIFU problems.
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
  std::string outFileName = outDir + std::string("run.optOut");

  std::default_random_engine seedGenerator(masterSeed);

  int rank;
  int numProcs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  
  LaunchParams p(masterSeed);
   
  p.islandId            = rank;
  p.islandCount         = numProcs;

  p.problemDimension    = 4 * SON_NUM;
  p.upperBounds         = std::vector<float>(p.problemDimension);
  p.lowerBounds         = std::vector<float>(p.problemDimension);

  getConstraints(p.problemDimension, p.upperBounds, p.lowerBounds);

  int retCode;
  if(p.islandId == PP_ISLAND_PARENT)
  {
    retCode = runParent<float>(p, outFileName);
  }
  else
  {
    runChild<float>(p);
  }

  if (retCode != 0)
  {
    MPI_Abort(MPI_COMM_WORLD, 1);
    return retCode;
  }

  MPI_Finalize();
  return 0;
}
//END of main