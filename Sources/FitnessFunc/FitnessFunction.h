/**
 * @file        FitnessFunction.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The header file for defining fitness functions to optimize by optimizers.
 * 
 * @version     0.2
 * 
 * @date        2019-10-20 (created) \n
 *              2020-01-11 (revised)
 */
#pragma once
#ifndef FITNESS_FUNCTION_H_INCLUDED
#define FITNESS_FUNCTION_H_INCLUDED

#include <vector>

/**
 * @brief calculates the fitness score of a given chromosome
 * 
 * @param [in]          n                - size of the chromosome / dimension of the problem
 * @param [in, out]     chromosome       - array of size n of genes to evaluate
 * @return float*                        - fitness score of given chromosome
 */
float fitnessFunction(const unsigned n, 
                      float* chromosome);

bool reachedLimit(const unsigned fitnessEvaluations,
                  const unsigned evalLimit);

bool isConverged(const float optimum, 
                 const float currentBest, 
                 const float eps);

/**
 * @brief Checks if the given chromozome is inside the constraints.
 * 
 * @param chromozome    - the chromozome data 
 * @return int          - index of the first gene that is out of constraints. Returns -1 if none is.
 */
int isInConstraints(std::vector<float> chromozome);


/**
 * @brief user defined way of how to handle the current optimized variables that are outside of bounds
 * 
 * @param [in, out] data  - current data of all points in search space (variables) 
 */
void applyConstraints(std::vector<float>& data);

/**
 * @brief way for the user to set bounds on searched variables (domain constraints)
 * 
 * @param [in] n                  - dimension of the problem (number of searched variables)
 * @param [in, out] upperBounds   - upper constraints on search space for each variable
 * @param [in, out] lowerBounds   - lower constraints on search space for each variable
 */
void getConstraints(const unsigned n, 
                    std::vector<float>& upperBounds,
                    std::vector<float>& lowerBounds);
#endif //FITNESS_FUNCTION_H_INCLUDED