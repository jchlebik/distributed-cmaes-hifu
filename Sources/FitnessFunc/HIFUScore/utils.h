/**
 * @file        utils.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Header file containing all declarations for auxilliar classes used during the simulation.
 * 
 * @version     0.3
 * 
 * @date        2020-03-02 (created) \n
 *              2020-11-24 (update) \n
 */

#pragma once
#ifndef UTILS_STRUCTS_H
#define UTILS_STRUCTS_H

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

#ifndef SEP
  #ifdef _WIN32
    #define SEP '\\';
  #else 
    #define SEP '/'
  #endif
#endif

#ifndef __INTEL_COMPILER
  #define _mm_malloc(a, b) aligned_alloc((b), (a))
  #define _mm_free(p) free(p)
  #define __assume_aligned(p, n); 
#endif

#define _ALIGN_ 64

/**
 * @brief A structure carrying information about the discrete grid \n
 *        representing the search space
 * 
 */
class GridInfo
{
public:
  /**
   * @brief Size of the grid on the x axis
   * 
   */
  unsigned nx;

  /**
   * @brief Size of the grid on the y axis
   * 
   */
  unsigned ny;

  /**
   * @brief The distance between the point on the grid in meters
   * 
   */
  float dx;
protected:
  /**
   * @brief Construct a new auxilliary Grid object
   * 
   * @param [in] xSize    - size of the grid on the x axis
   * @param [in] ySize    - size of the grid on the y axis
   * @param [in] dX       - the distance between the point on the grid in meters
   */
  GridInfo(int xSize,
           int ySize,
           float dX);

};
//END of Grid struct

/**
 * @brief An auxilliary sonication object carrying all the neccessary data \n
 *        to start a sonication simulation given a medium.
 * 
 */
struct Sonication
{
  /**
   * @brief The x coordinate for the sonication
   * 
   */
  const float xPos;

  /**
   * @brief The y coordinate for the sonication
   * 
   */
  const float yPos;

  /**
   * @brief For how long the sonication runs
   * 
   */
  const float timeOn;

  /**
   * @brief For how long to cool down after the sonication
   * 
   */
  const float timeOff;

  /**
   * @brief The step size for the simulation
   * 
   */
  const float dT;

  /**
   * @brief Constructs a new auxilliary Sonication object
   * 
   * @param [in] xP     - the x coordinate for the sonication
   * @param [in] yP     - the y coordinate for the sonication
   * @param [in] tOn    - for how long the sonication runs
   * @param [in] tOff   - for how long to cool down after the sonication
   * @param [in] dt     - the step size for the simulation
   */
  Sonication(float xP,
             float yP,
             float tOn,
             float tOff,
             float dt);

  /**
   * @brief Get the discrete amount of steps for how long the sonication runs
   * 
   * @return int    - the number of steps
   */            
  unsigned getOnSteps() const;

  /**
   * @brief Get the discrete amount of steps for how long the medium cools down after sonication
   * 
   * @return int    - the number of steps
   */
  unsigned getOffSteps() const;

  /**
   * @brief Get the discrete x coordinate for the sonication in the matrix
   * 
   * @return int    - the x coordinate
   */
  unsigned getXPos() const;

  /**
   * @brief Get the discrete y coordinate for the sonication in the matrix
   * 
   * @return int    - the y coordinate
   */
  unsigned getYPos() const;

};
//END of Sonication struct


#endif