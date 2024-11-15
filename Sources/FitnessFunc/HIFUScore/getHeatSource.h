/**
 * @file        getHeatSource.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for creating a Gaussian heating source
 * 
 * @version     0.5
 * 
 * @date        2020-03-02 (created) \n
 *              2020-11-25 (update)
 */
#pragma once
#ifndef GET_HEAT_SOURCE_H
#define GET_HEAT_SOURCE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "utils.h"

/**
 * @brief Create a Gaussian heating source centered at the position x, y \n
 *        discretised using the given grid parameters.
 * 
 * @param [in] nx     - the amount of points in grid on the x axis
 * @param [in] ny     - the amount of points in grid on the y axis
 * @param [in] x      - x coordinate of the heating source position
 * @param [in] y      - y coordinate of the heating source position
 * @param [in] dx     - spacing of actual samples in the grid (in meters)
 * @return std::vector<float> 
 */
float* getHeatSource(const unsigned nx,
                    const unsigned ny,
                    const unsigned x,
                    const unsigned y,
                    const float dx,
                    const float fullWidthHalfMax);
#endif