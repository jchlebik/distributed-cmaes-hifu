/**
 * @file        discreteMap.c
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       A source code file containing all definitions for DiscreteMap class.
 * 
 * @version     0.3
 * 
 * @date        2020-03-12 (created) \n
 *              2020-11-24 (update) \n
 */

#include "utils.h"
#include "medium.h"
#include "discreteMap.h"
#include "heatTargetWindow.h"

/**
 * @brief Construct a new 2D TargetMap object of xDimSize to yDimSize with circles specified in coordinates. \n
  *        Edge cases are wrapped around its opposite edge calculated by modulus operation. Identical to the 
  *        Matlab version used for fitness maps in k-Wave EUD.
  * 
  * @param [in] g              - a Grid class carrying relevant information about the map (xn, yn, dx)
  * @param [in] data           - the flattened data values of the map.
  */  
DiscreteMap::DiscreteMap(GridInfo g,
                         std::vector<int>& mapData):
GridInfo(g), data(mapData)
{}
//END of TargetMap(Grid) ctor

/**
 * @brief Construct a new 2D TargetMap object of xDimSize to yDimSize with circles specified in coordinates. \n
  *       Edge cases are wrapped around its opposite edge calculated by modulus operation. Result is  \n 
  *       dentical to the Matlab version used for fitness maps in k-Wave HIFU optimization.
  * 
  * @param [in] xDimSize       - the size of the x axis of the discrete medium space
  * @param [in] yDimSize       - the size of the y axis of the discrete medium space
  * @param [in] dx             - the distance between two points in the discrete space in meters
  * @param [in] data           - the flattened data values of the map.
  */
DiscreteMap::DiscreteMap(const unsigned xDimSize, 
                         const unsigned yDimSize, 
                         const float dx,
                         std::vector<int>& mapData): 
GridInfo(int(xDimSize), int(yDimSize), dx), data(mapData)
{}
// END of DiscreteMap(unsigned, unsigned, float) ctor

/**
 * @brief sums the data contained in the map
 * 
 * @return int - sum of the data in the map
 */
int DiscreteMap::sum()
{
  int sumOfElements = 0;

  #pragma omp parallel for reduction(+:sumOfElements)
  for (int i = 0; i < data.size(); i++)
  {
    sumOfElements += data[i];
  }
  return sumOfElements;
}
//END of sum

  DiscreteMap DiscreteMap::getWindowedMap(HeatTargetWindow & targetWindow)
  {
    const unsigned xRange = targetWindow.getXAxisLen();
    const unsigned yRange = targetWindow.getYAxisLen();
    const unsigned xStart = targetWindow.getXStart();
    const unsigned yStart = targetWindow.getYStart();

    const unsigned xRangeOriginal = targetWindow.getXAxisOriginalLen();

    std::vector<int> windowData(xRange * yRange);

    for (int y = 0; y < yRange; y++)
    {
      for (int x = 0; x < xRange; x++)
      {
        windowData[y * xRange + x]  = data[(y + yStart) * xRangeOriginal + xStart + x];
      }
    }
    return DiscreteMap(xRange, yRange, dx, windowData);
  }

/**
 * @brief Definition of the multiplication operation of the discrete maps. \n
 *        Equivalent to scalar product of the matrices. Requires L an R to \n
 *        be of the same size
 * 
 * @param [in] L          - L side argument of the multiplication
 * @param [in] R          - R side argument of the multiplication
 * @return DiscreteMap    - L * R
 */
DiscreteMap operator*(DiscreteMap const &L, DiscreteMap const &R)
{
  if (L.data.size() != R.data.size()) 
  {     
    throw std::logic_error("Incompatible map sizes for operator *."); 
  }

  std::vector<int> res(L.data.size());

  #pragma omp parallel for
  for (int i = 0; i < res.size(); i++)
  {
    res[i] = L.data[i] * R.data[i];
  }
  return DiscreteMap(L.nx, L.ny, L.dx, res);
}
//END of operator*()

/**
 * @brief Definition of the substration operation for the discrete maps. \n
 *        Equivalent to substraction of the matrices. Requires L an R to \n
 *        be of the same size
 * 
 * @param [in] L          - L side argument of the substraction
 * @param [in] R          - R side argument of the substraction
 * @return DiscreteMap    - L - R
 */
DiscreteMap operator-(DiscreteMap const &L, DiscreteMap const &R)
{
  if (L.data.size() != R.data.size()) 
  { 
    throw std::logic_error("Incompatible map sizes for operator -.");
  }

  std::vector<int> res(L.data.size());
  
  #pragma omp parallel for
  for (int i = 0; i < res.size(); i++)
  {
    res[i] = L.data[i] - R.data[i];
  }
  return DiscreteMap(L.nx, L.ny, L.dx, res);
}
//END of operator-()

/**
 * @brief Definition of the addition operation of the discrete maps. \n
 *        Equivalent to dot addition of matrices. Requires L an R to \n
 *        be of the same size
 * 
 * @param [in] L          - L side argument of the addition
 * @param [in] R          - R side argument of the addition
 * @return DiscreteMap    - L .+ R
 */
DiscreteMap operator+(DiscreteMap const &L, DiscreteMap const &R)
{
  if (L.data.size() != R.data.size()) 
  { 
    throw std::logic_error("Incompatible map sizes for operator +.");
  }

  std::vector<int> res(L.data.size());

  #pragma omp parallel for
  for (int i = 0; i < res.size(); i++)
  {
    res[i] = L.data[i] + R.data[i];
  }
  return DiscreteMap(L.nx, L.ny, L.dx, res);
}
//END of operator+()

/**
 * @brief Definition of the subscription operation of the discrete maps.
 * 
 * @param [in] index      - index of desired data
 * @return int            - data of discrete map on given index 
 */
int DiscreteMap::operator[](int index)
{
  return data[index];
}
//END of operator[]()


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
DiscreteMap DiscreteMap::getCirclesMap(const unsigned nx,
                                       const unsigned ny,
                                       const float dx,
                                       std::string circlesMapFilePath,
                                       std::function<int(int, int)> ftor)
{
  std::vector<CircleInfo> circlesCollection;
  std::string circleLine;
  std::stringstream lineStream;

  std::fstream mapCirclesStream(circlesMapFilePath, std::ios_base::in);
  if (mapCirclesStream.bad())
  {
    throw std::runtime_error("File " + circlesMapFilePath + " with circles specification does not exist");
  }

  std::vector<DiscreteMap> maps;

  std::string sx, sy, radius, val;

  try
  {
    while (std::getline(mapCirclesStream, circleLine))
    {
      lineStream.str(circleLine);
      if (lineStream.peek() != '#')
      {
        if (lineStream.peek() == '-')
        {
          maps.push_back(DiscreteMap::createTargetsMap(nx, ny, dx, circlesCollection, ftor));
          circlesCollection.clear();
        }
        else
        {
          std::getline(lineStream,  sx,      ',');
          std::getline(lineStream,  sy,      ',');
          std::getline(lineStream,  radius,  ',');
          std::getline(lineStream,  val,     ',');

          auto n = val.find("#");
          if (n != std::string::npos)
          {
            val = val.substr(0, n);
          }

          circlesCollection.push_back(
            CircleInfo { 
              unsigned(std::round(std::strtof(sx.c_str(), nullptr))),
              unsigned(std::round(std::strtof(sy.c_str(), nullptr))),
              unsigned(std::round(std::strtof(radius.c_str(), nullptr))),
              int(std::round(std::strtof(val.c_str(), nullptr))),
            }
          );
        }
      }
      lineStream.clear();
    }
  }
  catch(...)
  {
    throw std::runtime_error("File " + circlesMapFilePath + " contains invalid characters.");
  }
  
  if (!circlesCollection.empty())
  {
    maps.push_back(DiscreteMap::createTargetsMap(nx, ny, dx, circlesCollection, ftor));
    circlesCollection.clear();
  }

  for (int i = maps.size()-1; i > 0; i--)
  {
    maps[i-1] = maps[i] + maps[i-1];
  }
  return maps[0];
}
//END of GetCirclesMap

/**
 * @brief Creates a lesion map specific to the current state of the given medium
 * 
 * @param [in] m             - Medium instance to create the lesion map for
 * @return DiscreteMap       - a binary DiscreteMap representing the lesion map
 */
DiscreteMap DiscreteMap::createLesionMap(Medium& m)
{
  std::vector<int> lesionMap(m.matSize);

  #pragma omp parallel for
  for (int i = 0; i < m.matSize; i++)
  {
    lesionMap[i] = m.cem43[i] >= 240 ? 1 : 0;
  }

  return DiscreteMap(m.nx, m.ny, m.dx, lesionMap);
}
// END of CreateLesionMap

/**
 * @brief Creates a DiscreteMap object of the required size from the circles defining file. \n 
 *        Instead of iterating over the entire matrix and calculating distances to every \n
 *        given center, we create a "box" around each circle and  iterate over fields contained \n
 *        in that box.
 * 
 * @param [in] nx                     - size of the map on x dimension
 * @param [in] ny                     - size of the map on y dimension
 * @param [in] dx                     - distance between the neighbours in the real medium [m]
 * @param [in] circles                - a vector container of the CircleInfo instances representing the circles on the \n
 *                                      discrete map
 * @param [in] ftor                   - a function - how apply the value to its position (e.g. by addition or replacing)
 * @return DiscreteMap                - a nx x ny map filled with circles given in container
 */
DiscreteMap DiscreteMap::createTargetsMap(const unsigned nx, 
                                          const unsigned ny, 
                                          const float dx, 
                                          std::vector<CircleInfo>& circles,
                                          std::function<int(int, int)> ftor)
{
  std::vector<int> res(nx*ny, 0);

  // negative modulo lambda function
  auto mod = [](int a, int b) -> int{
    if (b == 0)
    {
      return (unsigned)!((int)0);
    }
    int ret = a % b;
    while(ret < 0)
    {
      ret += b;
    }
    return ret; 
  };

  // instead of iterating over the entire matrix and calculating distances to every 
  // given center, we create a "box" around each circle and 
  // iterate over fields contained in that box
  for(auto &coord : circles)
  {
    // find box bounds
    // with edge cases, the bounds are allowed to be negative for smooth iterations (extending the marix)
    int lowBound    = coord.x - coord.radius;
    int highBound   = coord.x + coord.radius;
    int leftBound   = coord.y - coord.radius;
    int rightBound  = coord.y + coord.radius;

    // iterate over the fields in box
    //#pragma omp parallel for
    for (int row = lowBound; row < highBound; row++)
    {
      int moduloRow = mod(row, int(nx)); // cmath modulus does not calculate negative modulos correctly
      //#pragma omp simd
      for (int col = leftBound; col < rightBound; col++)
      {
        int moduloCol = mod(col, int(ny));
        int index = moduloRow * nx + moduloCol;

        // no need to round the distance because we are using <. we can just cut the floating part off using integer
        // faster than calculating the square root of distance and power operation is congruent on Inequality split
        int distanceSquared = (row - (int(coord.x) - 1)) * (row - (int(coord.x) - 1))  + (col - (int(coord.y) - 1)) * (col - (int(coord.y) - 1));
        int radiusSquared   = coord.radius * coord.radius;  

        if (distanceSquared < radiusSquared)
        {
          res[index] = ftor(res[index], coord.circleValue);
        }
      }
    }
  }
  return DiscreteMap(nx, ny, dx, res);
}
//END of CreateTargetsMap

/**
 * @brief Debuggin tool. Prints the map to the standard error stream.
 * 
 * @param [in] sep    - Separator to print between values. Defaults to space.
 */
void DiscreteMap::printMap(std::string sep)
{
  for (int i = 0; i < ny; i++)
  {
    for (int j = 0; j < nx; j++)
    {
      std::cout << data[i * nx + j] << sep;
    }
    std::cout << std::endl;
  }
}
//END of PrintMap
