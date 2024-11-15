/**
 * @file        discreteMap.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Header file containing all declarations for DiscreteMap class and its data structure.
 * 
 * @version     0.2
 * 
 * @date        2020-03-12 (created) \n
 *              2020-03-18 (update) \n
 */
#pragma once
#ifndef TARGET_MAP_CLASS_H
#define TARGET_MAP_CLASS_H

#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

#include "utils.h"

class Medium;
class HeatTargetWindow;

/*
#include "utils.h"
#include "medium.h"
*/

/**
 * @brief An auxilliary data structure to wrap the definition of a circle in a target map.
 * 
 */
struct CircleInfo
{
  /**
   * @brief An x coordinate in the grid of the circle center.
   * 
   */
  const unsigned x;
  /**
   * @brief An y coordinate in the grid of the circle center.
   * 
   */
  const unsigned y;
  /**
   * @brief A radius of the circle.
   * 
   */
  const unsigned radius;
  /**
   * @brief A value that the circle is going to be filled with.
   * 
   */
  const int      circleValue;
};

/**
 * @brief A class representing the discrete map of values. 
 * 
 */
class DiscreteMap : public GridInfo
{
  /**
   * @brief A vector of values containing the disrete map of values.
   * 
   */
  std::vector<int> data;

public:
 /**
  * @brief Construct a new 2D TargetMap object of given Grid size with circles specified in coordinates. \n
  *        Edge cases are wrapped around its opposite edge calculated by modulus operation. Identical to the 
  *        Matlab version used for fitness maps in k-Wave EUD.
  * 
  * @param [in] g              - A GridInfo struct carrying relevant information about the map (xn, yn, dx).
  * @param [in] data           - The flattened data values of the map.
  */  
  DiscreteMap(GridInfo g, 
              std::vector<int>& mapData);

  /**
   * @brief Construct a new 2D TargetMap object of given Grid size with circles specified in coordinates. \n
   *        Edge cases are wrapped around its opposite edge calculated by modulus operation. Identical to the 
   *        Matlab version used for fitness maps in k-Wave EUD.
   * 
   * @param [in] xDimSize       - The size of the x axis of the discrete medium space.
   * @param [in] yDimSize       - The size of the y axis of the discrete medium space.
   * @param [in] dx             - The distance between two points in the discrete space in meters.
   * @param [in] data           - The flattened data values of the map..
   */  
  DiscreteMap(const unsigned xDimSize,
              const unsigned yDimSize,
              const float dx,
              std::vector<int>& mapData);

  /**
   * @brief Sums the data contained in the map.
   * 
   * @return int - Sum of the data in the map.
   */
  int sum();

  
  DiscreteMap getWindowedMap(HeatTargetWindow & targetWindow);


  /**
   * @brief Creates a lesion map specific to the current state of the given medium.
   * 
   * @param [in] m             - Medium instance to create the lesion map for.
   * @return DiscreteMap       - A binary DiscreteMap representing the lesion map.
   */
  static DiscreteMap createLesionMap(Medium& m);

/**
 * @brief Creates a DiscreteMap object of the required size from the circles defining file
 * 
 * @param [in] nx                     - size of the map on x dimension
 * @param [in] ny                     - size of the map on y dimension
 * @param [in] dx                     - distance between the neighbours in the real medium [m]
 * @param [in] circlesMapFilePath     - path to the file containing the definitions of circles that creates the \n
 *                                      discrete map
 * @return DiscreteMap                - a nx x ny map filled with circles given in the file
 */
  static DiscreteMap getCirclesMap(const unsigned nx,
                                   const unsigned ny,
                                   const float dx,
                                   std::string circlesMapFilePath,
                                   std::function<int(int, int)> ftor);

  /**
   * @brief Creates a DiscreteMap object of the required size from the circles defining file. \n 
   *        Instead of iterating over the entire matrix and calculating distances to every \n
   *        given center, we create a "box" around each circle and  iterate over fields contained \n
   *        in that box.
   * 
   * @param [in] nx                     - Size of the map on x dimension.
   * @param [in] ny                     - Size of the map on y dimension.
   * @param [in] dx                     - Distance between the neighbours in the real medium [m].
   * @param [in] circles                - A vector container of the CircleInfo instances representing the circles on the \n
   *                                      discrete map.
   * @return DiscreteMap                - A nx x ny map filled with circles given in container.
   */
  static DiscreteMap createTargetsMap(const unsigned nx, 
                                      const unsigned ny, 
                                      const float dx, 
                                      std::vector<CircleInfo>& circles,
                                      std::function<int(int, int)> ftor);


  /**
   * @brief Debuggin tool. Prints the map to the standard error stream.
   * 
   * @param [in] sep    - Separator to print between values. Defaults to space.
   */
  void printMap(std::string sep = " ");

  /**
   * @brief Definition of the substration operation for the discrete maps. \n
   *        Equivalent to substraction of matrices. Requires L an R to \n
   *        be of the same size
   * 
   * @param [in] L          - L side argument of the substraction
   * @param [in] R          - R side argument of the substraction
   * @return DiscreteMap    - L - R
   */
  friend DiscreteMap operator-(DiscreteMap const &, DiscreteMap const &);

  /**
   * @brief Definition of the multiplication operation of the discrete maps. \n
   *        Equivalent to scalar product of matrices. Requires L an R to \n
   *        be of the same size
   * 
   * @param [in] L          - L side argument of the multiplication
   * @param [in] R          - R side argument of the multiplication
   * @return DiscreteMap    - L * R
   */
  friend DiscreteMap operator*(DiscreteMap const &, DiscreteMap const &);

  /**
   * @brief Definition of the addition operation of the discrete maps. \n
   *        Equivalent to dot addition of matrices. Requires L an R to \n
   *        be of the same size
   * 
   * @param [in] L          - L side argument of the addition
   * @param [in] R          - R side argument of the addition
   * @return DiscreteMap    - L .+ R
   */
  friend DiscreteMap operator+(DiscreteMap const &, DiscreteMap const &);

  /**
   * @brief Definition of the subscription operation of the discrete maps.
   * 
   * @param [in] index      - index of desired data
   * @return int            - data of discrete map on given index 
   */
  int operator[](int index);

};
#endif