/**
 * @file        calculateHeating.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for calculating the heating equations
 * 
 * @version     0.1
 * 
 * @date        2020-03-20 (created) \n
 */
#pragma once
#ifndef CALCULATE_HEATING_H
#define CALCULATE_HEATING_H

#include <vector>
#include <cmath>

/*
#include "utils.h"
#include "medium.h"
*/

class Medium;
class Sonication;

#include "KWaveDiffusionSolver/kWaveDiffusionSolver.h"

/**
 * @brief Compute heating diffusion for simple trajectory optimisation problem.
 * 
 * @param [in, out] medium              - the structure representing the medium where the heat is being diffused \n
 *                                        contains all relevant matrices like temperature and so on. the results \n
 *                                        of the simulation are stored in their respective matrices in this class
 * @param [in]      sonicationParams    - parameters of the simulation, the amount of steps and step size
 */
 void calculateHeating(Medium & medium, 
                       Sonication sonicationParams);
#endif