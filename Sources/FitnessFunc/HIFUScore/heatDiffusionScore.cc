/**
 * @file        heatDiffusionScore.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The scoreing function for heat diffusion simulation fitness.
 * 
 * @version     0.2
 * 
 * @date        2020-11-11 (created) \n
 */

#include <vector>
#include <string>
#include <cmath>

#include "../FitnessFunction.h"

#include "utils.h"
#include "discreteMap.h"
#include "heatTargetWindow.h"
#include "medium.h"
#include "getHeatSource.h"
#include "calculateHeating.h"
#include "KWaveDiffusionSolver/kWaveDiffusionSolver.h"

#ifndef HIFU_PROBLEM
  #define HIFU_PROBLEM
#endif

#ifndef BLOB 
  #define BLOB 1
#endif
#ifndef FLOWER
  #define FLOWER 0
#endif

#ifndef PROBLEM_TYPE
  #define PROBLEM_TYPE BLOB
#endif

#ifndef SEP
  #ifdef _WIN32
    #define SEP '\\';
  #else 
    #define SEP '/'
  #endif
#endif

float fitnessFunction(const unsigned n, 
                      float* chromosome)
{
  const int sonicationsAmount = n / 4;
  const float dt = 0.1f;

  float xPos, yPos, timeOn, timeOff;

  std::vector<Sonication> sonications;

  for (int i = 0; i < n; i+=4)
  {
    xPos      = chromosome[i]  ;
    yPos      = chromosome[i+1];
    timeOn    = chromosome[i+2];
    timeOff   = chromosome[i+3];




    sonications.push_back(
      Sonication{  
        std::round(xPos),  
        std::round(yPos),  
        std::round(timeOn),  
        std::round(timeOff), 
        dt 
      }
    );
  }

  // static variables to save computation time on multiple fitness executions
  static const std::string dataFolderPath          = std::string("FitnessFunc") + SEP + "HIFUScore" + SEP + "Data";
  static const std::string austinWomanFolderPath   = dataFolderPath + SEP + "AustinWoman";
  static const std::string targetMapsFolderPath    = dataFolderPath + SEP + "Targets";

  #if PROBLEM_TYPE == FLOWER
    static const std::string problemType = std::string("Flower");
  #else
    static const std::string problemType = std::string("Blob");
  #endif

  static const std::string targetMapCirclesFilePath      = targetMapsFolderPath + SEP + problemType + SEP + "targetMapCircles.dat";
  static const std::string prohibitedMapCirclesFilePath  = targetMapsFolderPath + SEP + problemType + SEP + "prohibitedMapCircles.dat";
    
  static Medium mediumBeforeTreatment(austinWomanFolderPath);

  const unsigned nx = mediumBeforeTreatment.nx;
  const unsigned ny = mediumBeforeTreatment.ny;
  const float    dx = mediumBeforeTreatment.dx;

  static DiscreteMap targetMap = DiscreteMap::getCirclesMap(nx,
                                                            ny,
                                                            dx,
                                                            targetMapCirclesFilePath,
                                                            [](int a, int b) -> int { return a > b ? a + b : b;});

  static DiscreteMap penalizeMap = DiscreteMap::getCirclesMap(nx,
                                                              ny,
                                                              dx,
                                                              prohibitedMapCirclesFilePath,
                                                              [](int a, int b) -> int {return a + b;});

  //save some computation cost between many runs by using statics. TODO: rework into class based approach
  //static HeatTargetWindow targetWindow(penalizeMap);
  //int xS = targetWindow.getXStart();
  //int yS = targetWindow.getYStart();
  //int xL = targetWindow.getXAxisLen();
  //int yL = targetWindow.getYAxisLen();
  //printf("%d %d %d %d\n", xS, yS, xL, yL);
  /*
  static DiscreteMap targetMapWindow    = targetMap.getWindowedMap(targetWindow);
  static DiscreteMap penalizeMapWindow  = penalizeMap.getWindowedMap(targetWindow);

  static Medium windowedMediumBeforeTreatment = mediumBeforeTreatment.getWindowedMedium(targetWindow);

  Medium workingWindowedMedium(windowedMediumBeforeTreatment);
  */
  const int magnitude = 1;
  std::vector<float*> heatSources(sonicationsAmount);

  for (int i = 0; i < sonicationsAmount; i++)
  {
    heatSources[i] = getHeatSource(mediumBeforeTreatment.nx, 
                                   mediumBeforeTreatment.ny, 
                                   sonications[i].xPos, 
                                   sonications[i].yPos, 
                                   mediumBeforeTreatment.dx, 
                                   magnitude);
  }


  Medium mediumCopy(mediumBeforeTreatment);
  for (int i = 0; i < sonicationsAmount; i++)
  {
    // we set the heating energy matrix into the medium instance
    //workingWindowedMedium.setHeatDepositEnergyOnWindow(targetWindow, heatSources[i]);
    // and calculate heating on this window
    //calculateHeating(workingWindowedMedium, sonications[i]);
    //fprintf(stderr, "sonication %d: %0.1f %0.1f %0.1f %0.1f\n", i, sonications[i].xPos, sonications[i].yPos, sonications[i].timeOn, sonications[i].timeOff);
    mediumCopy.setHeatDeposit(heatSources[i]);
    calculateHeating(mediumCopy, sonications[i]);
    _mm_free(heatSources[i]);
  }

  DiscreteMap lesionMap = DiscreteMap::createLesionMap(mediumCopy);
  int stillAlive        = (targetMap - (targetMap * lesionMap)).sum();
  int mistreated        = (penalizeMap * lesionMap).sum();
/*
  fprintf(stderr, "%d %d\n", targetMap.sum(), penalizeMap.sum());
  fprintf(stderr, "%d %d\n", stillAlive, mistreated);
*/
  /*


  DiscreteMap lesionMap = DiscreteMap::createLesionMap(workingWindowedMedium);

  int stillAlive = (targetMapWindow - (targetMapWindow * lesionMap)).sum();
  int mistreated = (penalizeMapWindow * lesionMap).sum();



  fprintf(stderr, "%d %d\n", targetMapWindow.sum(), penalizeMapWindow.sum());
  fprintf(stderr, "%d %d\n", stillAlive, mistreated);
*/
  return stillAlive + mistreated;
}

bool reachedLimit(const unsigned fitnessEvaluations,
                  const unsigned evalLimit)
{
  return fitnessEvaluations >= evalLimit;
}

bool isConverged(const float optimum, 
                 const float currentBest, 
                 const float eps)
{
  return std::abs(optimum - currentBest) < eps;
}

int _isInConstraints_blob(std::vector<float> chromozome)
{
  const int chromoLen = chromozome.size();
  std::vector<float> upperBounds(4);
  std::vector<float> lowerBounds(4);
  getConstraints(4, upperBounds, lowerBounds);

  float x, y, t_on, t_off;
  int retValue = -1;
  for (int i = 0; i < chromoLen; i+=4)
  {
    x     = chromozome[i];
    y     = chromozome[i+1];
    t_on  = chromozome[i+2];
    t_off = chromozome[i+3];

    if (x > upperBounds[0] || x < lowerBounds[0])
    {
      retValue = i;
      break;
    }
    if (y > upperBounds[1] || y < lowerBounds[1])
    {
      retValue = i + 1;
      break;
    }
    if (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      retValue = i + 2;
      break;
    }
    if (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      retValue = i + 3;
      break;
    }
  }
  return retValue;
}

int _isInConstraints_flower(std::vector<float> chromozome)
{
  const int chromoLen = chromozome.size();
  std::vector<float> upperBounds(4);
  std::vector<float> lowerBounds(4);
  getConstraints(4, upperBounds, lowerBounds);

  const float innerCircleRadius    = 10;
  const float innerCircleCenter[]  = {247.5, 247.5};
  const float outerCircleRadius    = 40;

  float pointDistanceSquared;
  float radiusSquared;

  //std::cerr << "isInConstraints ";
  float x, y, t_on, t_off;
  int retValue = -1;
  for (int i = 0; i < chromoLen; i+=4)
  {
    x     = chromozome[i];
    y     = chromozome[i+1];
    t_on  = chromozome[i+2];
    t_off = chromozome[i+3];

	//cirlce with cirlce inside
    /*
    float pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + 
                                   (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    float radiusSquared        = (outerCircleRadius * outerCircleRadius);
    if (pointDistanceSquared > radiusSquared)
    {
      const float distX = std::abs(innerCircleCenter[0] - x);
      const float distY = std::abs(innerCircleCenter[1] - y);

      if (distX < distY)
      {
        //std::cerr << "y" << "\n";
        retValue = i + 1;
        break;
      }
      else
      {
        //std::cerr << "x" << "\n";
        retValue = i;
        break;
      }
    }*/

    if (x > upperBounds[0] || x < lowerBounds[0])
    {
      //std::cerr << "x" << "\n";
      retValue = i;
      break;
    }

    if (y > upperBounds[1] || y < lowerBounds[1])
    {
      //std::cerr << "y" << "\n";
      retValue = i + 1;
      break;
    }

    if (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      //std::cerr << "tOn" << "\n";
      retValue = i + 2;
      break;
    }
    if (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      //std::cerr << "tOff" << "\n";
      retValue = i + 3;
      break;
    }
    
	// is inside the middle circle, if so, return the axis closest to the middle
    pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + 
                                   (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    radiusSquared        = (innerCircleRadius * innerCircleRadius);
    if (pointDistanceSquared < radiusSquared)
    {
      const float distX = std::abs(innerCircleCenter[0] - x);
      const float distY = std::abs(innerCircleCenter[1] - y);

      if (distX < distY)
      {
        retValue = i;

      }
      else
      {
        retValue = i + 1;
      }
      break;
    }
    
  }
  return retValue;
}

/**
 * @brief User defined function to checks if the given chromozome is inside the constraints.
 * 
 * @param chromoLen     - length of the chromozome to check
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int isInConstraints(std::vector<float> chromozome)
{
  #if PROBLEM_TYPE == FLOWER
    return _isInConstraints_flower(chromozome);
  #else
    return _isInConstraints_blob(chromozome);
  #endif
}
// END of isInConstrains

void _applyConstraints_blob(std::vector<float> &data)
{
  const int n = data.size();
  std::vector<float> upperBounds(4);
  std::vector<float> lowerBounds(4);
  getConstraints(4, upperBounds, lowerBounds);

  float ranges[] = { upperBounds[0] - lowerBounds[0], 
                      upperBounds[1] - lowerBounds[1], 
                      upperBounds[2] - lowerBounds[2], 
                      upperBounds[3] - lowerBounds[3] };

  float x, y, t_on, t_off;
  for (int i = 0; i < n; i += 4)
  {
    x     = data[i];
    y     = data[i+1];
    t_on  = data[i+2];
    t_off = data[i+3];

    while (x > upperBounds[0] || x < lowerBounds[0])
    {
      if (x > upperBounds[0])
      { 
        data[i] = upperBounds[0] - ((int)std::round(std::abs(x - upperBounds[0])) % (int)ranges[0]);
      }
      else if (x < lowerBounds[0])
      {
        data[i] = lowerBounds[0] + ((int)std::round(std::abs(lowerBounds[0] - x)) % (int)ranges[0]);
      }
      x = data[i];
    }
    
    while (y > upperBounds[1] || y < lowerBounds[1])
    {
      if (y > upperBounds[1])
      {
        data[i+1] = upperBounds[1] - ((int)std::round(std::abs(y - upperBounds[1])) % (int)ranges[1]);
      }
      else if (y < lowerBounds[1])
      {
        data[i+1] = lowerBounds[1] + ((int)std::round(std::abs(lowerBounds[1] - y)) % (int)ranges[1]);
      }
      y = data[i+1];
    }

    while (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      if (t_on > upperBounds[2])
      {
        data[i+2] = upperBounds[2] - ((int)std::round(std::abs(t_on - upperBounds[2])) % (int)ranges[2]);
      }
      else if (t_on < lowerBounds[2])
      {
        data[i+2] = lowerBounds[2] + ((int)std::round(std::abs(t_on - lowerBounds[2])) % (int)ranges[2]);
      }
      t_on  = data[i+2];
    }

    while (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      if (t_off > upperBounds[3])
      {
        data[i+3] = upperBounds[3] - ((int)std::round(std::abs(t_off - upperBounds[3])) % (int)ranges[3]);
      }
      else if (t_off < lowerBounds[3])
      {
        data[i+3] = lowerBounds[3] + ((int)std::round(std::abs(t_off - lowerBounds[3])) % (int)ranges[3]);
      }
      t_off = data[i+3];
    }
  }
}

void _applyConstraints_flower(std::vector<float>& data)
{

  const int n = data.size();
  std::vector<float> upperBounds(4);
  std::vector<float> lowerBounds(4);
  getConstraints(4, upperBounds, lowerBounds);

  //const float upperBounds[] = {287.5, 287.5,  5.0, 20.0};
  //const float lowerBounds[] = {207.5, 207.5,  0.0,  1.0};
  float ranges[] = { upperBounds[0] - lowerBounds[0], 
                      upperBounds[1] - lowerBounds[1], 
                      upperBounds[2] - lowerBounds[2], 
                      upperBounds[3] - lowerBounds[3] };

  const float innerCircleRadius    = 10;
  const float innerCircleCenter[]  = {247.5, 247.5};
  
  const float outerCircleRadius    = 40;
  //const float outerCircleCenter[]  = {247.5, 247.5};

  float x, y, t_on, t_off;
  for (int i = 0; i < n; i += 4)
  {
    x     = data[i];
    y     = data[i+1];
    t_on  = data[i+2];
    t_off = data[i+3];
  
  /*
    // pushes back inside the circle in the direction created by the line joining the point and the middle 
	// using the distance as a scale.
    float pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + 
                                   (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    float radiusSquared        = (outerCircleRadius * outerCircleRadius);
    if (pointDistanceSquared > radiusSquared)
    {
      //std::cerr << " = ";
      float scale = (std::sqrt(pointDistanceSquared)/outerCircleRadius) - 1.0;
      float shift = outerCircleRadius * scale;

      float theta = std::atan2(y - innerCircleCenter[1], x - innerCircleCenter[0]);
      float intersection[] = { 
          innerCircleCenter[0] + (outerCircleRadius - shift) * std::cos(theta), 
          innerCircleCenter[1] + (outerCircleRadius - shift) * std::sin(theta)
      };
      x = data[i]    = intersection[0];
      y = data[i+1]  = intersection[1];
      //std::cerr << x << " " << y << "\n";
    }
  */

    while (x > upperBounds[0] || x < lowerBounds[0])
    {
      //std::cerr << "tOn apply constraints " << t_on << "\n";
      if (x > upperBounds[0])
      {
        data[i] = upperBounds[0] - ((int)std::round(std::abs(x - upperBounds[0])) % (int)ranges[0]);
      }
      else if (x < lowerBounds[0])
      {
        data[i] = lowerBounds[0] + ((int)std::round(std::abs(x - lowerBounds[0])) % (int)ranges[0]);
      }
      x  = data[i];
    }

    while (y > upperBounds[1] || y < lowerBounds[1])
    {
      //std::cerr << "tOff apply constraints " << t_off << "\n";
      if (y > upperBounds[1])
      {
        data[i+1] = upperBounds[1] - ((int)std::round(std::abs(y - upperBounds[1])) % (int)ranges[1]);
      }
      else if (y < lowerBounds[1])
      {
        data[i+1] = lowerBounds[1] + ((int)std::round(std::abs(y - lowerBounds[1])) % (int)ranges[1]);
      }
      y = data[i+1];
    }

    while (t_on > upperBounds[2] || t_on < lowerBounds[2])
    {
      //std::cerr << "tOn apply constraints " << t_on << "\n";
      if (t_on > upperBounds[2])
      {
        data[i+2] = upperBounds[2] - ((int)std::round(std::abs(t_on - upperBounds[2])) % (int)ranges[2]);
      }
      else if (t_on < lowerBounds[2])
      {
        data[i+2] = lowerBounds[2] + ((int)std::round(std::abs(t_on - lowerBounds[2])) % (int)ranges[2]);
      }
      t_on  = data[i+2];
    }

    while (t_off > upperBounds[3] || t_off < lowerBounds[3])
    {
      //std::cerr << "tOff apply constraints " << t_off << "\n";
      if (t_off > upperBounds[3])
      {
        data[i+3] = upperBounds[3] - ((int)std::round(std::abs(t_off - upperBounds[3])) % (int)ranges[3]);
      }
      else if (t_off < lowerBounds[3])
      {
        data[i+3] = lowerBounds[3] + ((int)std::round(std::abs(t_off - lowerBounds[3])) % (int)ranges[3]);
      }
      t_off = data[i+3];
    }
  
  
  /*
    x     = data[i];
    y     = data[i+1];

    float pointDistanceSquared = ((x - innerCircleCenter[0]) * (x - innerCircleCenter[0]) + (y - innerCircleCenter[1]) * (y - innerCircleCenter[1]));
    float radiusSquared        = (innerCircleRadius * innerCircleRadius);
    if (pointDistanceSquared < radiusSquared)
    {
      //std::cerr << "applyConstraints " << x << " " << y << " = ";
      float scale = std::sqrt(pointDistanceSquared)/innerCircleRadius;
      float shift = (outerCircleRadius - innerCircleRadius) * scale;

      float theta = std::atan2(y - innerCircleCenter[1], x - innerCircleCenter[0]);
      float intersection[] = { 
          innerCircleCenter[0] + (innerCircleRadius + shift) * std::cos(theta), 
          innerCircleCenter[1] + (innerCircleRadius + shift) * std::sin(theta)
      };
      x = data[i]    = intersection[0];
      y = data[i+1]  = intersection[1];
      //std::cerr << "" << x << " " << y << "\n";
    }
  */
  
  }
}

/**
 * @brief user defined way of how to handle the current optimized variables that are outside of constraints
 * 
 * @param [in] n          - size of the data (for problem dimension 'n' and population size 'p' its 'n * p')
 * @param [in, out] data  - current data of points in search space (variables) 
 */
void applyConstraints(std::vector<float>& data)
{
  #if PROBLEM_TYPE == FLOWER
    _applyConstraints_flower(data);
  #else
    _applyConstraints_blob(data);
  #endif
}
// END of applyConstraints

void _getConstraints_flower(const unsigned n, 
                            std::vector<float> &uB,
                            std::vector<float> &lB)
{
  const float upperBounds[] = {287.5, 287.5, 20.0,  150.0};
  const float lowerBounds[] = {207.5, 207.5,  0.0,    0.0};

  for (int i = 0; i < n; i += 4)
  {
    uB[i] = upperBounds[0];
    uB[i+1] = upperBounds[1];
    uB[i+2] = upperBounds[2];
    uB[i+3] = upperBounds[3];

    lB[i] = lowerBounds[0];
    lB[i+1] = lowerBounds[1];
    lB[i+2] = lowerBounds[2];
    lB[i+3] = lowerBounds[3];
  }
}

void _getConstraints_blob(const unsigned n, 
                          std::vector<float> &uB,
                          std::vector<float> &lB)
{
  const float upperBounds[] = {342.0, 293.0, 20.0, 20.0};
  const float lowerBounds[] = {272.0, 233.0,  0.0,  0.0};

  for (int i = 0; i < n; i += 4)
  {
    uB[i] = upperBounds[0];
    lB[i] = lowerBounds[0];

    uB[i+1] = upperBounds[1];
    lB[i+1] = lowerBounds[1];

    uB[i+2] = upperBounds[2];
    lB[i+2] = lowerBounds[2];

    uB[i+3] = upperBounds[3];
    lB[i+3] = lowerBounds[3];
  }
}

/**
 * @brief way for the user to set outer bounds on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void getConstraints(const unsigned n, 
                    std::vector<float> &uB,
                    std::vector<float> &lB)
{
  #if PROBLEM_TYPE == FLOWER
    _getConstraints_flower(n, uB, lB);
  #else
    _getConstraints_blob(n, uB, lB);
  #endif
}
//END of getConstraints 