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
#pragma once
#ifndef LOGGING_HOOKS_CMAESPP_H_INCLUDED
#define LOGGING_HOOKS_CMAESPP_H_INCLUDED


#include "../CMA-ESpp/cma-es/cmaes.h"
#include "LoggingHooks.h"

/**
 * @brief creates a string of information to log for the current generation of given population. 
 *        Used with real value population based optimizers.
 * 
 * @param [in] evo           - cmaes class containing the current progress of the optimization
 * @param [in] pop           - current population
 * @return char*             - allocated string of information to print out
 */


std::string createPopulationLogInfo(CMAES<float> &evo,
                                    float const* pop);

#endif // LOGGING_HOOKS_CMAESPP_H_INCLUDED