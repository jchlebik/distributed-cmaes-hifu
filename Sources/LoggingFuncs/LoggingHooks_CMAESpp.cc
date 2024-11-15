/**
 * @file        LoggingHooks_CMAESpp.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file of logging functions used to log progression of CMA-ESpp library.
 * 
 * @version     0.3
 * 
 * @date        2019-10-20 (created) \n
 *              2020-02-16 (revised)
 */

#include "LoggingHooks_CMAESpp.h"

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] evo           - cmaes class containing the current progress of the optimization
 * @param [in] pop           - current population
 * @return char*             - allocated string of information to print out
 */

std::string createPopulationLogInfo(CMAES<float> &evo,
                                    float const* pop)
{
  unsigned int generation = evo.get(CMAES<float>::Generation);
  unsigned int N = evo.get(CMAES<float>::PopSize);
  stats* st = getStats<float>(pop, N);
  if (!st)
  {
    //char* res = (char*)malloc(sizeof("-1"));
    //sprintf(res, "-1");
    return {};
  }
  const float bestEver = evo.get(CMAES<float>::FBestEver);
  const float bestCurrent = evo.get(CMAES<float>::Fitness);

  const float* bestChromo = evo.getPtr(CMAES<float>::XBest);
  int requiredBytes = snprintf(nullptr, 0,
                               "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
                               generation, st->min, st->max, bestCurrent, st->mean, st->median, st->q1, st->q3, st->irq);
  char* buffer = new(std::nothrow) char[requiredBytes+1];
  if (!buffer)
  {   delete st;
      return {};
  }
  sprintf(buffer, "$gen:%d$best:%g$worst:%g$current:%g$mean:%g$median:%g$q1:%g$q3:%g$irq:%g$chromosome:",
          generation, st->min, st->max, bestCurrent, st->mean, st->median, st->q1, st->q3, st->irq);

  std::string textToLog = appendChromoData<float>(
    std::string(buffer), 
    (float*)bestChromo, 
    (int)evo.get(CMAES<float>::Dimension), 
    "%g"
  );

  delete st;
  return textToLog;
  /*
  char* cstr = (char*)new char[textToLog.size()+1];
  if (cstr == nullptr)
  {
    delete[] buffer;
    return nullptr;
  }
  delete[] buffer;
  strcpy(cstr, textToLog.c_str());
  return cstr;
  */
}