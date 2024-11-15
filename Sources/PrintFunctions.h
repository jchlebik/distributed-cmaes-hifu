#pragma once

#include <string>
#include <cstdio>
#include <chrono>
#include <cstdlib>

#include "CMA-ESpp/cma-es/cmaes.h"

#include "LaunchParams.h"
#include "LoggingFuncs/LoggingHooks_CMAESpp.h"


/**
 * @brief the hook function called with the start of evaluation of every generation.
 * 
 * @param [in]  evo    - a class containing all cmaes data for the optimalization run
 * @param [in]  pop    - a current population
 * @return boolean     - should the optimization continue
 */
template <class T>
bool printProgress(CMAES<T> & evo, 
                   LaunchParams & p)
{
  T* popFitness = (T*)evo.getPtr(CMAES<T>::FVals);
  std::string textToLog = createPopulationLogInfo(evo, popFitness);
  std::string popDump   = createPopulationDump<T>(popFitness, evo.get(CMAES<T>::PopSize));

  if (!textToLog.empty() && !popDump.empty())
  {
    fprintf(p.logFile, "%s", textToLog.c_str());
    fprintf(p.logFile, "$fitnessEvaluations:%d", p.fitnessCalls);
    fprintf(p.logFile, "%s\n", popDump.c_str());
  }
  else
  {
    fprintf(p.logFile, "population died out.\n");
    return false;
  }
  return true;

}
// END of printProgress *************************************************************



/**
 * @brief Prints help for the program.
 * 
 * @return Nothing.
 */
void printUsage()
{
  printf("USAGE: all arguments required\n");
  printf("./... [rngSeed] [outputsDirectoryPath]\n");
  printf("starting positions and standard deviations are calculated as mean and stdev of normal distribution from the specified bounds in fitness function.\n");
}

template<class T>
void printArray(T* arr, 
                const size_t arrLen, 
                FILE* logFile, 
                std::string format = "%g",
                std::string separator = ",")
{
  for (int i = 0; i < arrLen; i++)
  {
    fprintf(logFile, format.c_str(), arr[i]);
    if (i + 1 < arrLen)
    {
      fprintf(logFile, separator.c_str());
    }
  }
}

template <class T>
void printHeader(LaunchParams &p, const int masterSeed)
{
  fprintf(p.logFile, "$lambda:%d$genesInChromosome:%d$iterations:%d", 
          p.popSize, p.problemDimension, p.generations);
  
  fprintf(p.logFile, "$lowerBounds:");
  printArray<T>(p.lowerBounds.data(), p.problemDimension, p.logFile);
  
  fprintf(p.logFile, "$upperBounds:");
  printArray<T>(p.upperBounds.data(), p.problemDimension, p.logFile);

  fprintf(p.logFile, "$optimalScore:%g$masterSeed:%u$rngSeed:%u",p.optScore, masterSeed, p.seed);
  
  std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::string timeString(std::ctime(&now));
  fprintf(p.logFile, "$date:%s$", timeString.substr(0, timeString.length() - 1).c_str());
}

template <class T>
void printFooter(LaunchParams & p, 
                 T bestFit, 
                 T* bestChromo, 
                 unsigned long long elapsed)
{
  fprintf(p.logFile, "$fitness:%g$time[uS]:%llu", bestFit, elapsed);
  fprintf(p.logFile, "$fitnessEvaluations:%d$generationsNeeded:%d", p.fitnessCalls, p.generationsNeeded);
  fprintf(p.logFile, "$best:");
  
  printArray<T>(bestChromo, p.problemDimension, p.logFile);
  fprintf(p.logFile, "$\n");
}