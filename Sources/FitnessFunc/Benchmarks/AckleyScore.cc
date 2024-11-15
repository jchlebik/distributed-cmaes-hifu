/**
 * @file        AcleyScore.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file of Acleys function fitness score real value optimization.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */

#include "../FitnessFunction.h"
#include <cmath>
#include <vector>

//# define M_PI 3.14159265358979323846 // pi
//# define M_E 2.7182818284590452354   // e

#ifndef MAX_FIT_EVALS
  #define MAX_FIT_EVALS 2000
#endif

/**
 * @brief calculates the fitness score of a given chromosome by Acleys function
 *        optimum is at x_i = 0 forall genes, domain defined for all x_i in [-15, 30]
 * 
 * @param [in]          n                - size of the chromosome / dimension of the problem
 * @param [in, out]     chromosome       - array of size n of genes to evaluate
 * @return float                         - fitness score of given chromosome
 */
float fitnessFunction(const unsigned n, 
                      float* chromosome)
{
    const float gamma = 2.0 * M_PI;
    float squaredInputsSum = 0.0;
    float cosSum = 0.0;

    for (unsigned i = 0; i < n; i++)
    {
        float value = chromosome[i];
        squaredInputsSum += std::pow(value, 2.0);
        cosSum += std::cos(gamma * value);
    }

    float inputAvgCoeff = 1.0f/float(n);
    return (20.0f - 20.0f * std::exp(-0.2f *std::sqrt(inputAvgCoeff * squaredInputsSum)) + M_E - std::exp(inputAvgCoeff * cosSum));
}
// end of fitnessFunction **********************************************************

bool reachedLimit(const unsigned fitnessEvaluations,
                  const unsigned evalLimit)
{
  return fitnessEvaluations > evalLimit;
}

bool isConverged(const float optimum, 
                 const float currentBest, 
                 const float eps)
{
  return std::abs(optimum - currentBest) < eps;
}

/**
 * @brief user defined way of how to handle the current optimized variables that are outside of bounds
 * 
 * @param [in, out] data  - current data of all points in search space (variables) 
 */
void applyConstraints(std::vector<float>& data)
{
  float value = 0.0;
  for (int i = 0; i < data.size(); i++)
  {
    value = data[i];
    if ( value < -20.0)
    {
      data[i] = -20.0;
    }
    else if (value > 20.0)
    {
      data[i] = 20.0;
    }
  }
}
//END of applyConstraints

/**
 * @brief way for the user to set bounds on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void getConstraints(const unsigned n, 
                    std::vector<float>& upperBounds,  // out var
                    std::vector<float>& lowerBounds)  // out var
  
{
  for (int i = 0; i < n; i++)
  {
    upperBounds[i] = 20.0;
    lowerBounds[i] = -20.0;
  }
}
//END getConstraints

/**
 * @brief Checks if the given chromozome is inside the constraints.
 * 
 * @param chromoLen     - length of the chromozome to check
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int isInConstraints(std::vector<float> chromozome)
{
  std::vector<float> upperBounds;
  std::vector<float> lowerBounds;
  getConstraints(chromozome.size(), upperBounds, lowerBounds);

  int retVal = -1;
  for (int i = 0; i < chromozome.size(); i)
  {
    if (chromozome[i] > upperBounds[i] || chromozome[i] < lowerBounds[i])
    {
      retVal = i;
      break;
    }
  }
  return retVal;
}
//END isInConstraints