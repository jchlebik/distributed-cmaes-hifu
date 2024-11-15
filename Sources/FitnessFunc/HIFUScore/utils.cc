/**
 * @file        utils.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Source code file containing all definition for auxilliar classes used during the simulation.
 * 
 * @version     0.3
 * 
 * @date        2020-03-02 (created) \n
 *              2020-03-06 (update) \n
 */

#include "utils.h"

/**
 * @brief Construct a new auxilliary GridInfo struct
 * 
 * @param [in] xSize    - size of the grid on the x axis
 * @param [in] ySize    - size of the grid on the y axis
 * @param [in] dX       - the distance between the point on the grid in meters
 */
GridInfo::GridInfo(int xSize, 
                   int ySize, 
                   float dX): 
  nx(xSize), 
  ny(ySize), 
  dx(dX) {}
// END of Grid::ctor

/**
 * @brief Construct a new auxilliary Sonication object
 * 
 * @param [in] xP     - the x coordinate for the sonication
 * @param [in] yP     - the y coordinate for the sonication
 * @param [in] tOn    - how long the sonication runs in seconds
 * @param [in] tOff   - for how long to cool down after the sonication
 * @param [in] dt     - the step size for the simulation
 */
Sonication::Sonication(float xP, 
                       float yP, 
                       float tOn, 
                       float tOff, 
                       float dt):
  xPos(xP),
  yPos(yP),
  timeOn(tOn),
  timeOff(tOff),
  dT(dt) {}
// END of Sonication::ctor

/**
 * @brief Get the discrete amount of steps for how long the sonication runs
 * 
 * @return int    - the number of steps
 */
unsigned Sonication::getOnSteps() const
{
  return unsigned(std::round(timeOn/dT));
}
// END of Sonication::getOnSteps

/**
 * @brief Get the discrete amount of steps for how long the medium cools down after sonication
 * 
 * @return int    - the number of steps
 */
unsigned Sonication::getOffSteps() const
{
  return unsigned(std::round(timeOff/dT));
}
// END of Sonication::getOffSteps

/**
 * @brief Get the discrete x coordinate for the sonication in the matrix
 * 
 * @return int    - the x coordinate
 */
unsigned Sonication::getXPos() const
{
  return unsigned(std::round(xPos));
}
// END of Sonication::getXPos

/**
 * @brief Get the discrete y coordinate for the sonication in the matrix
 * 
 * @return int    - the y coordinate
 */
unsigned Sonication::getYPos() const
{
  return unsigned(std::round(yPos));
}
// END of Sonication::getYPos
