/**
 * @file        medium.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Source file containing all definitions for Medium class used during the simulation. \n
 *              This class represents and contains all the data about the medium we are using for \n
 *              heat deposit and tissue ablation simulations.
 * 
 * @version     0.1
 * 
 * @date        2020-11-24 (created) \n
 */


#include "utils.h"    //GridInfo
#include "heatTargetWindow.h"
#include "medium.h"

/**
 * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
 *        Achieved by multiple matrices, each giving a value to a point on the grid.
 * 
 * @param [in] folderPath   - the path to the folder containing density.txt, heat.txt and conductivity.txt
 */
Medium::Medium(std::string folderPath):
  GridInfo(0, 0, 0.0f)
{
  init(folderPath);
}
// END of Medium::ctor

Medium::~Medium()
{
  _mm_free(q);
  _mm_free(t);
  _mm_free(cem43);
  _mm_free(rho);
  _mm_free(c);
  _mm_free(lambda);
}

/**
 * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
 *        Achieved by multiple matrices, each giving a value to a point on the grid and cropped by \n
 *        a given window represented by HeatTargetWindow class.
 * 
 * @param [in] targetWindow   - The instance of the HeatTargetWindow
 */
Medium::Medium(HeatTargetWindow   targetWindow, 
               float* heat, 
               float* temperature, 
               float* density, 
               float* capacity,
               float* cem,
               float* conductivity) :
  GridInfo(targetWindow.getXAxisLen(), targetWindow.getYAxisLen(), targetWindow.getResolution())
{
  matSize = nx * ny;
  q       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  t       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  cem43   = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  rho     = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );  
  c       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ ); 
  lambda  = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ ); 


  for (int i = 0; i < matSize; i++)
  {
    q[i] = heat[i];
    t[i] = temperature[i];
    cem43[i] = cem[i];
    rho[i] = density[i];
    c[i] = capacity[i];
    lambda[i] = conductivity[i];
  }
}
//END of Medium::ctor


Medium::Medium(const Medium & cpy) :
  GridInfo(cpy.nx, cpy.ny, cpy.dx)
{
  matSize = nx * ny;
  q       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  t       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  cem43   = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  rho     = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );  
  c       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ ); 
  lambda  = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ ); 

  for (int i = 0; i < matSize; i++)
  {
    q[i] = cpy.q[i];
    t[i] = cpy.t[i];
    cem43[i] = cpy.cem43[i];
    rho[i] = cpy.rho[i];
    c[i] = cpy.c[i];
    lambda[i] = cpy.lambda[i];
  }
}


/**
 * @brief Initializes the medium to a begining state. Prepares all data for simulation.
 * 
 * @param [in] folderPath   - the path to the folder containing mediumSpecs.dat, density.dat, heat.dat and conductivity.dat
 */
void Medium::init(std::string folderPath)
{
  if (folderPath.back() != SEP)
  {
    folderPath += SEP;
  }

  std::fstream densityFile(       folderPath + "density.dat",        std::ios_base::in);
  std::fstream specificHeatFile(  folderPath + "heat.dat",           std::ios_base::in);
  std::fstream conductivityFile(  folderPath + "conductivity.dat",   std::ios_base::in);
  std::fstream mediumSpecsFile(   folderPath + "mediumSpecs.dat",    std::ios_base::in);

  if (densityFile.bad() || specificHeatFile.bad() || conductivityFile.bad() || mediumSpecsFile.bad())
  {
    throw std::runtime_error("One or more of the needed medium specifying files could not be opened.");
  }
  
  float              density,          specificHeat,     thermalConductivity;
  std::string        xDimSizeString,   yDimSizeString,   dxSizeString;          // for parsing the mediumSpecs file
  std::string        densityLine,      heatLine,         conductivityLine;      // for parsing lines of data files
  std::string        densityVal,       heatVal,          conductivityVal;       // for single values on each line
  std::stringstream  densityStream,    heatStream,       conductivityStream;    // to convert lines of files into streams

  std::getline(mediumSpecsFile, xDimSizeString, ',');   // nx
  std::getline(mediumSpecsFile, yDimSizeString, ',');   // ny
  std::getline(mediumSpecsFile, dxSizeString,   ',');   // dx

  nx = int(std::strtof(xDimSizeString.c_str(), nullptr));
  ny = int(std::strtof(yDimSizeString.c_str(), nullptr));
  dx = std::strtof(dxSizeString.c_str(), nullptr);
  mediumSpecsFile.close();

  matSize = nx * ny;
  q       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  t       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  cem43   = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );
  rho     = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ );  
  c       = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ ); 
  lambda  = (float*)_mm_malloc(sizeof(float) * matSize, _ALIGN_ ); 

  int i = 0;
  // comma operator to perform every reading and avoid shortcircuiting on first
  while(std::getline(densityFile,      densityLine),
        std::getline(specificHeatFile, heatLine),
        std::getline(conductivityFile, conductivityLine))
  {
    densityStream.clear();
    heatStream.clear();
    conductivityStream.clear();

    densityStream.str(densityLine);
    heatStream.str(heatLine);
    conductivityStream.str(conductivityLine);

    densityVal      = "";
    heatVal         = "";
    conductivityVal = "";

    while(std::getline(densityStream,      densityVal,      ','),
          std::getline(heatStream,         heatVal,         ','),
          std::getline(conductivityStream, conductivityVal, ','))
    {
      density             = strtof(densityVal.c_str(), nullptr);
      specificHeat        = strtof(heatVal.c_str(), nullptr);
      thermalConductivity = strtof(conductivityVal.c_str(), nullptr);
      
      q[i]      = 0.0f;                          // volume rate of heat deposition [W/m^3]
      t[i]      = 37.0f;                         // temperature distribution [degC]
      cem43[i]  = 0.0f;                      // thermal dose given in cumulative equivalent minutes (cem) relative to T = 43 degC
      rho[i]    = density;                     // tissue mass density [kg/m^3]
      c[i]      = specificHeat;                  // tissue specific heat capacity [J/(kg.K)]
      lambda[i] = thermalConductivity;      // tissue thermal conductivity [W/(m.K)]
      i++;
    }
  }

  if (!(densityFile.eof() && specificHeatFile.eof() && conductivityFile.eof()))
  {
    throw std::runtime_error("One or more of the needed medium files are not equal in their dimensions.");
  }
  if (densityFile.bad() || specificHeatFile.bad() || conductivityFile.bad())
  {
    throw std::runtime_error("Error reading medium defining files.");
  }
}
// END of Medium::init



void Medium::setHeatDeposit(float* heatAdded)
{
  HeatTargetWindow wholeMedium = HeatTargetWindow(0,0,nx,ny,nx,ny,dx);
  setHeatDepositEnergyOnWindow(wholeMedium, heatAdded);
}

/**
 * @brief Sets the heat deposit energy into the medium after cropping it by the provided window.
 * 
 * @param [in] heatTargetWindow       - Instance of HeatTargetWindow used to crop the heat doposition \n
 *                                      vector to the appropriate window.
 * @param [in] heatAdded              - Vector of data containing the values for the added heat energy \n
 *                                      onto the medium. Expected to be of the original size and to be \n
 *                                      cropped by the provided window.
 */
void Medium::setHeatDepositEnergyOnWindow(HeatTargetWindow & heatTargetWindow, float* heatAdded)
{
  const unsigned xRange = heatTargetWindow.getXAxisLen();
  const unsigned yRange = heatTargetWindow.getYAxisLen();
  const unsigned xStart = heatTargetWindow.getXStart();
  const unsigned yStart = heatTargetWindow.getYStart();

  const unsigned xRangeOriginal = heatTargetWindow.getXAxisOriginalLen();

  for (int row = 0; row < yRange; row++)
  {
    const int rowStartIndex = (yStart + row) * xRangeOriginal + xStart;
    for(int col = 0; col < xRange; col++)
    {
      q[row * xRange + col] = heatAdded[rowStartIndex + col];
    }
  }
}
//END of setHeatDepositEnergyOnWindow

/**
 * @brief Static method to create an original sized medium with updated values \n
 *        at locations in the medium cropped by window. Does not modify its \n
 *        parameters, instead opting to create a new instance.
 * 
 * @param [in] heatTargetWindow   - A window that was used to create the croppedMedium \n
 *                                  instance.
 * @param [in] originalMedium     - Original sized Medium containing the data from before \n 
 *                                  the sonications. Copy of this medium with updated values \n
 *                                  for t and cem43 vector is returned.
 * @param [in] croppedMedium      - A medium that was cropped by the provided window carrying \n
 *                                  all relevant information after the sonications.
 * 
 * @returns Medium                - A copy of the originalMedium parameter with updated values \n
 *                                  of t and cem43 vectors from the cropped medium.
 * 
 */
Medium Medium::updateOriginalByCroppedMedium(HeatTargetWindow & heatTargetWindow,
                                             Medium originalMedium,
                                             Medium & croppedMedium)
{
  const int xRange = heatTargetWindow.getXAxisLen();
  const int yRange = heatTargetWindow.getYAxisLen();
  const int xStart = heatTargetWindow.getXStart();
  const int yStart = heatTargetWindow.getYStart();

  Medium m(originalMedium);

  for (int row = 0; row < yRange; row++)
  {
    const int rowStartIndex = (yStart + row) * originalMedium.nx + xStart;
    for (int col = 0; col < xRange; col++)
    {
      m.t[rowStartIndex + col] = croppedMedium.t[row * xRange + col];
      m.cem43[rowStartIndex + col] = croppedMedium.cem43[row * xRange + col];
    }
  }

  return m;
}
//END of updateOriginalByCroppedMedium
/**
 * @brief   Provides a new instance of the Medium class that contains the data \n
 *          that fit inside the provided HeatTargetWindow. The new instance is \n
 *          of smaller size (the size of the window).
 * 
 * @param [in] heatTargetWindow     - A window with which to crop the medium.
 * 
 * @returns Medium                  - A new smaller instance of the Medium class \n
 *                                    with data only inside the window.
 * 
 */
Medium Medium::getWindowedMedium(HeatTargetWindow & heatTargetWindow)
{
  const int xRange = heatTargetWindow.getXAxisLen();
  const int yRange = heatTargetWindow.getYAxisLen();
  const int xStart = heatTargetWindow.getXStart();
  const int yStart = heatTargetWindow.getYStart();

  const unsigned xRangeOriginal = heatTargetWindow.getXAxisOriginalLen();

  std::vector<float> heat(xRange * yRange);
  std::vector<float> temperature(xRange * yRange);
  std::vector<float> density(xRange * yRange);
  std::vector<float> capacity(xRange * yRange);
  std::vector<float> cem(xRange * yRange);
  std::vector<float> conductivity(xRange * yRange);
  
  
  std::vector<std::vector<float>*> croppedDatasVector
  {
    &heat, &temperature, &density, &capacity, &cem, &conductivity
  };

  std::vector<float**> originalDatasVector
  {
    &q, &t, &rho, &c, &cem43, &lambda
  };

  for (int i = 0; i < originalDatasVector.size(); i++)
  {
    for (int row = 0; row < yRange; row++)
    {
      const int rowStartIndex = (yStart + row) * xRangeOriginal + xStart;
      for(int col = 0; col < xRange; col++)
      {
        (*croppedDatasVector[i])[row * xRange + col] = (*originalDatasVector[i])[rowStartIndex + col];
      }
    }
  }

  return Medium(heatTargetWindow,
                heat.data(),
                temperature.data(),
                density.data(),
                capacity.data(),
                cem.data(),
                conductivity.data());
}
//END of getWindowedMedium
