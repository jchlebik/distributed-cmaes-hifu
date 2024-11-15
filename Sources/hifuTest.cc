/**
 * @file        hifuTest.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       A small program capable of running one heat diffusion simulation \n
 *				and printing result to the standard error output
 * 
 * @version     0.1
 * 
 * @date        2020-04-06 (created) \n
 */

#include <iostream>
#include "FitnessFunc/FitnessFunction.h"
#include <ctime>
#include <random>
#include <algorithm>
#include <time.h>

int main()
{

  struct timespec   start, end;
  /*float arr[] = {219.612, 276.122,  9.565,    89.7943, 
                   217.881, 219.77,   8.39193,  87.1744,
                   228.734, 256.377,  7.65081,  83.1919,
                   256.724, 276.753,  9.24244,  96.2361,
                   239.166, 281.089,  7.07247,  77.2956,
                   218.265, 216.67,   7.45452,  48.7604,
                   217.265, 248.73,   8.93117,  59.4931,
                   271.995, 255.494,  9.14906,  91.7417,
                   238.4,   208.947,  7.02642,  73.5941,
                   269.383, 214.727,  9.09386,  60.7772,
                 };
  */
  //float arr[] = {229.246,244.562,5.11265,85.5674,262.621,257.156,5.62753,68.6044,258.611,276.499,9.53613,80.233,276.383,235.531,6.56051,44.5602,274.428,257.47,9.60898,54.3651,270.419,220.252,5.24716,83.8084,236.183,214.068,8.62188,83.3931,253.56,208.707,6.27754,64.092,222.2,235.625,5.29283,56.2962,212.067,249.64,5.68405,74.3387,221.608,231.062,4.32317,53.6417,225.382,272.342,4.38868,90.0175,274.144,220.433,2.08807,85.5653,272.36,214.298,3.15047,84.4163,236.256,282.762,9.95524,73.2463};

  
  //float arr[] = {219.612,276.122,9.565,89.7943};
  // float arr[] = { 327.304, 274.763,  1.70586, 20.1175,
  //                 287.353, 248.787,  5.89634, 32.5261,
  //                 311.56,  269.152, 18.0208,  37.942,
  //                 322.58,  276.502,  4.74801, 31.684,
  //                 320.678, 274.029,  3.98928, 46.7107,
  //                 334.495, 275.434,  2.22895, 54.3376,
  //                 282.69,  263.369, 12.6796,  27.3497,
  //                 297.011, 265.132,  3.85737, 49.3633,
  //                 292.861, 267.037,  5.15391,  1.77128,
  //                 281.605, 263.245,  4.17929, 16.7277 };

       

  


  std::default_random_engine rngEngine(time(NULL));


  std::vector<float> arr(6*4);

  std::vector<float> times(15);
  std::vector<float> uB(4);
  std::vector<float> lB(4);
  getConstraints(4, uB, lB);

  std::uniform_real_distribution<float> distribution(0.0f, 1.0f);

  int n = 6*4;
  std::cerr << n << "\n";

  for (int i = 0; i < 30; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      for (int k = 0; k < 4; k++)
      {
        arr[j*4 + k] = lB[k] + distribution(rngEngine) * (uB[k] - lB[k]);
        //arr[j*4 + k] = uB[k];
        std::cerr << arr[j*4 + k] << " ";
      }
    }
    std::cerr << "\n";
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    fitnessFunction(n, arr.data());
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    times.push_back(elapsedTime/1000000.0f);
  }
  std::sort(times.begin(), times.end());


  const float median = times[times.size()/2];
  std::cerr << "Run time: " << median << "[S]\n";

	/*
  double* upperBound = (double*)malloc(n*sizeof(double));
  double* lowerBound = (double*)malloc(n*sizeof(double));
  getConstraints(n, &upperBound, &lowerBound); 
  applyConstraints(n, &arr);
	

  std::cerr << "Starting ...\n" << "Performing " << n/4 << " sonications\n";
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  int res = fitnessFunction(n, arr);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  unsigned long long elapsedTime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
  std::cerr << "Run time: " << elapsedTime/1000000.0f<< "[S]\n";
  std::cerr << "Resulting fitness: " << res << "\n";
  return 0;
  */

}
