/**
 * @file        LoggingHooks.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file of logging functions used to log progression of optimizers.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */
#pragma once
#ifndef LOGGING_HOOKS_H_INCLUDED
#define LOGGING_HOOKS_H_INCLUDED

#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>

#include "stats_t.h"

/**
 * @brief gets a statistical data from the current population 
 * 
 * @param [in] pop           - fitness values of population
 * @param [in] N             - size of the population
 * @return stats*            - structure containing all relevant statistics
 */

template<class T>
stats* getStats(T const* pop, const unsigned N)
{
  double max, min, mean, median, q1, q3, irq;
  std::vector<T> sortedPop(pop, pop + N);
  std::sort(sortedPop.begin(), sortedPop.end(), [](T first, T second){
    return first < second;
  });
  //const entity *best = ga_get_entity_from_rank(pop, 0);
  if (N < 4)
  {
    return new stats {q1 = NAN, q3 = NAN, median = NAN, mean = NAN, min = NAN, max = NAN, irq = NAN};
  }

  bool isOddLength = N & 1;
  int half = N / 2;
  int quarter = N / 4;
  int threeFourths = N * 3.0 / 4.0;
  if (isOddLength)
  {
    median = sortedPop[half];
    q1 = (sortedPop[quarter] + sortedPop[quarter + 1]) / 2.0;
    q3 = (sortedPop[threeFourths] + sortedPop[threeFourths + 1]) / 2.0;
  }
  else
  {
    median = (sortedPop[half] + sortedPop[half + 1]) / 2.0;
    q1 = sortedPop[quarter];
    q3 = sortedPop[threeFourths];
  }
  std::for_each(sortedPop.begin(), sortedPop.end(), [&] (T a) {
    mean += a;
  });
  irq = q3 - q1;
  mean = mean / (float)N;
  min = sortedPop[0];
  max = sortedPop[N-1];

  stats* result = new stats();
  result->q1 = q1;
  result->q3 = q3;
  result->median = median;
  result->mean = mean;
  result->min = min;
  result->max = max;
  result->irq = irq;
  return result;
}
// END OF getStats ******************************************************************
/**
 * @brief creates a string representation of all fitness values in given population. 
 * 
 * @param [in] pop           - fitness values of population
 * @param [in] N             - size of the population
 * @return char*             - allocated string of information to print out
 */
template<class T>
std::string createPopulationDump(T const* pop, const unsigned N)
{
  int i = 0; 
  std::string s = "$populationDump:";
  char buffer[32];

  for(int i = 0; i < N; i++)
  {
    snprintf(buffer, sizeof(buffer), "%g", pop[i]);
    s += std::string(buffer);
    if (i + 1 != N)
    {
      s += ",";
    }
  }
  return s;
/*
  char* cstr = (char*)malloc(s.size()+1*sizeof(char));
  if (cstr == nullptr)
  {
    return nullptr;
  }
  strcpy(cstr, s.c_str());
  return cstr;
*/
}
// END OF createPopulationDump ******************************************************************

/**
 * @brief creates a string from chromosome data and text to be logged so far.
 * 
 * @param [in] textToLog          - prefix for the chromosome data
 * @param [in, out] chromo        - the data of the chromosome
 * @param [in, out] chromoLen     - the amount of fields in chromosome
 * @param [in, out] format        - format to print a chromosome into the resulting string
 * @return string                 - string of the chromosome data
 */
template <class T>
std::string appendChromoData(std::string textToLog, 
                            T* chromo,
                            const int chromoLen,
                            const char* format)
{
  std::string resultCopy(textToLog);
  char buffer[32];
  for (int i = 0; i < chromoLen; i++)
  {
    sprintf(buffer, format, chromo[i]);
    //std::string tmp = std::to_string(chromo[i]);
    std::string tmp(buffer);
    resultCopy.append(tmp);
    if (i + 1 != chromoLen)
    {
      resultCopy.append(",");
    }
  }
  return resultCopy;
}
#endif //LOGGING_HOOKS_H_INCLUDED