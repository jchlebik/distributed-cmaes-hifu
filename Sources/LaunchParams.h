#pragma once

#include <vector>
#include <random>

struct LaunchParams
{
  //global
  int                 popSize;
  int                 iterations;
  int                 generations;
  int                 problemDimension;
  float               optScore;
  bool                verboseOutput;
  std::vector<float>  upperBounds;
  std::vector<float>  lowerBounds;
  int                 seed;
  int                 fitnessCalls;
  int                 generationsNeeded;
  std::mt19937        rngEngine;
  FILE*               logFile;
  int                 fitEvalsLimit;
  float               eps;
  unsigned long       maxSeconds;
  union sT
  {
    double d;
    unsigned long long l;
  } startTime;

  //arcipelego parameters
  int                 migrationInterval;
  int                 numberOfImmigrants;
  int                 islandId;
  int                 islandCount;
  int                 targetIsland;
  int                 sourceIsland;

  //pan population
  int                 popSizePerIsland;


  //evo
  std::vector<float>  xStart;
  std::vector<float>  stdDevs;
  int                 tournamentSize;
  //...

  LaunchParams(int _seed): rngEngine(_seed)
  {
    seed = _seed;
  }

};

template <class T>
struct Emigrants
{
  std::vector<T> fitness;
  std::vector<T> population;  //flattened

  Emigrants(const size_t numberOfImmigrants, 
            const size_t problemDimension) : 
    fitness(numberOfImmigrants), 
    population(numberOfImmigrants * problemDimension)
  {
  }
};

template <class T>
struct Immigrants
{
  std::vector<T> fitness;
  std::vector<T> population;  //flattened

  Immigrants(const size_t numberOfImmigrants, 
             const size_t problemDimension):
    fitness(numberOfImmigrants),
    population(numberOfImmigrants * problemDimension)
  {}
};

