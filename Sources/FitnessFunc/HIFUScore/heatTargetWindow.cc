 /**
 * @file        heatTargetWindow.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Source file containing definitions of HeatTargetWindow class. This class \n
 *              represents a "window" that is used to crop the whole medium to a subview \n
 *              based on either dimension size of the generated heat energy or a target or \n
 *              penalization map.
 * 
 * @version     0.1
 * 
 * @date        2020-11-24 (created) \n
 */

#include "heatTargetWindow.h"
#include "discreteMap.h"
 
/**
 * @brief Constructs a new HeatTargetWindow object from the discrete map. Used for \n
 *        creating a constant window based on penalize or target map data.
 * 
 * @param [in] dataMap    - The discrete map describing the penalization or target map. 
 */
HeatTargetWindow::HeatTargetWindow(DiscreteMap & dataMap)
{
  originalXAxisSize = dataMap.nx;
  originalYAxisSize = dataMap.ny;

  yLow = dataMap.ny + 1;
  xLow = dataMap.nx + 1;

  yHigh = 0;
  xHigh = 0;

  bool rowRelevant = false;

  for (unsigned row = 0; row < dataMap.ny; row++)
  {
    float value             = dataMap[row];
    unsigned oneDIndexStart = row * dataMap.nx;
    rowRelevant = false;

    for (unsigned col = 0; col < dataMap.nx; col++)
    {
      if (dataMap[oneDIndexStart + col] > 0.0f)
      {
        rowRelevant = true;
        xLow = std::min(xLow, col);
        xHigh = std::max(xHigh, col);
      }
    }

    if (rowRelevant)
    {
      yLow = std::min(yLow, row);
      yHigh = row;
    }
  }

  //if (xLow > 0) xLow--;
  //if (yLow > 0) yLow--;
  //if (xHigh < originalXAxisSize-1) xHigh++;
  //if (yHigh < originalYAxisSize-1) yHigh++;

  xAxisLen = xHigh - xLow + 1;
  yAxisLen = yHigh - yLow + 1;
}
//END of HeatTargetWindow ctor

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
HeatTargetWindow::HeatTargetWindow(const unsigned xStart,
                                   const unsigned yStart,
                                   const unsigned xLenWindow,
                                   const unsigned yLenWindow,
                                   const unsigned xLenOriginal,
                                   const unsigned yLenOriginal,
                                   const float dx):
  xLow(xStart),
  yLow(yStart),
  xAxisLen(xLenWindow),
  yAxisLen(yLenWindow),
  originalXAxisSize(xLenOriginal),
  originalYAxisSize(yLenOriginal),
  resolution(dx)
{
  xHigh = xLow + xAxisLen - 1;
  yHigh = yLow + yAxisLen - 1;

  //if (xLow > 0) xLow--;
  //if (yLow > 0) yLow--;
  //if (xHigh < originalXAxisSize-1) xHigh++;
  //if (yHigh < originalYAxisSize-1) yHigh++;
}
//END of HeatTargetWindow ctor

/**
 * @brief Gets the resolution inside the window. 
 * 
 */
float HeatTargetWindow::getResolution() const
{
  return resolution;
}
//END of getResolution()

/**
 * @brief Gets the length of the x axis of the window. 
 * 
 */
unsigned HeatTargetWindow::getXAxisLen() const
{
  return xAxisLen;
}
//END of getXAxisLen()

/**
 * @brief Gets the length of the y axis of the window. 
 * 
 */
unsigned HeatTargetWindow::getYAxisLen() const
{
  return yAxisLen;
}
//END of getYAxisLen()

/**
 * @brief Gets the starting index for the x dimension (columns) of the window \n 
 *        in context of the original medium map. 
 * 
 */
unsigned HeatTargetWindow::getXStart() const
{
  return xLow;
}
//END of getXStart()

/**
 * @brief Gets the starting index for the y dimension (rows) of the window \n
 *        in context of the original medium map. 
 * 
 */
unsigned HeatTargetWindow::getYStart() const
{
  return yLow;
}
//END of getYStart()

/**
 * @brief Gets the length of the x axis of the original medium. 
 * 
 */
unsigned HeatTargetWindow::getXAxisOriginalLen() const
{
  return originalXAxisSize;
}
//END of getXAxisOriginalLen

/**
 * @brief Gets the length of the y axis of the original medium. 
 * 
 */
unsigned HeatTargetWindow::getYAxisOriginalLen() const
{
  return originalYAxisSize;
}
//END of getYAxisOriginalLen