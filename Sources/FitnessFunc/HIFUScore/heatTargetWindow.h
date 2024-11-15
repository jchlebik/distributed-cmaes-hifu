/**
 * @file        heatTargetWindow.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Header file containing declaration of HeatTargetWindow class. This class \n
 *              represents a "window" that is used to crop the whole medium to a subview \n
 *              based on either dimension size of the generated heat energy or a target or \n
 *              penalization map.
 * 
 * @version     0.1
 * 
 * @date        2020-11-24 (created) \n
 */

#pragma once

#ifndef HEAT_TARGET_WINDOW_CLASS_H
#define HEAT_TARGET_WINDOW_CLASS_H

#include <cmath>
#include <algorithm>

class DiscreteMap;

/**
 * @brief A class carrying information about a "window" that \n
 *        is used to crop the whole medium to a subview based \n
 *        on dimension size of generated heat. Index zero-zero \n
 *        is assumed to be the top-leftmost element. Window is \n
 *        assumed to be a square one to ease the computations.
 * 
 */
class HeatTargetWindow
{
  /**
   * @brief Smallest index along the x-axis. Giving the left side of the window.
   * 
   */
  unsigned xLow;

  /**
   * @brief Highest index in the original medium along the x-axis. Giving the right side of the window.
   * 
   */
  unsigned xHigh;

  /**
   * @brief Smallest index in the original medium along the y-axis. Giving the top side of the window.
   * 
   */
  unsigned yLow;

  /**
   * @brief Highest index along the y-axis. Giving the bottom side of the window.
   * 
   */
  unsigned yHigh;

  /**
   * @brief Size of the window x-axis. Short for xHigh - xLow.
   * 
   */
  unsigned xAxisLen;

  /**
   * @brief Size of the window y-axis. Short for xHigh - xLow.
   * 
   */
  unsigned yAxisLen;

  /**
   * @brief Resolution of the data considered inside the window. Not used right now.
   * 
   */ 
  float    resolution;

  /**
   * @brief Size of the original x-axis.
   * 
   */
  unsigned originalXAxisSize;

  /**
   * @brief Size of the original y-axis.
   * 
   */
  unsigned originalYAxisSize;

public:

  /**
   * @brief Constructs a new HeatTargetWindow object from the discrete map. Used for \n
   *        creating a constant window based on penalize or target map data.
   * 
   * @param [in] dataMap    - The discrete map describing the penalization or target map. 
   */
  HeatTargetWindow(DiscreteMap & dataMap);


  /**
   * @brief Constructs a new HeatTargetWindow object.
   * 
   * @param [in] xStart         - Smallest index in the original medium along the x-axis for the window.
   * @param [in] yStart         - Smallest index in the original medium along the y-axis for the window.
   * @param [in] xLenWindow     - Size of the window x-axis.
   * @param [in] yLenWindow     - Size of the window y-axis.
   * @param [in] xLenOriginal   - Size of the x-axis of the original medium.
   * @param [in] yLenOriginal   - Size of the y-axis of the original medium.
   * @param [in] dx             - Resolution inside the window. Not used right now.
   */
  HeatTargetWindow(const unsigned xStart,
                   const unsigned yStart,
                   const unsigned xLenWindow,
                   const unsigned yLenWindow,
                   const unsigned xLenOriginal,
                   const unsigned yLenOriginal,
                   const float dx);

  /**
   * @brief Gets the resolution inside the window. 
   * 
   */
  float getResolution() const;

  /**
   * @brief Gets the length of the x axis of the window. 
   * 
   */
  unsigned getXAxisLen() const;

  /**
   * @brief Gets the length of the y axis of the window. 
   * 
   */
  unsigned getYAxisLen() const;

  /**
   * @brief Gets the starting index for the x dimension (columns) of the window \n 
   *        in context of the original medium map. 
   * 
   */
  unsigned getXStart() const;

  /**
   * @brief Gets the starting index for the y dimension (rows) of the window \n
   *        in context of the original medium map. 
   * 
   */
  unsigned getYStart() const;

  /**
   * @brief Gets the length of the x axis of the original medium. 
   * 
   */
  unsigned getXAxisOriginalLen() const;

  /**
   * @brief Gets the length of the y axis of the original medium. 
   * 
   */
  unsigned getYAxisOriginalLen() const;
};
//END of HeatTargetWindow class

#endif