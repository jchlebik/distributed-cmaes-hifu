/**
 * @file        medium.h
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       Header file containing all declarations for Medium class used during the simulation. \n
 *              This class represents and contains all the data about the medium we are using for \n
 *              heat deposit and tissue ablation simulations.
 * 
 * @version     0.1
 * 
 * @date        2020-11-24 (created) \n
 */

#pragma once
#ifndef MEDIUM_CLASS_H_INCLUDED
#define MEDIUM_CLASS_H_INCLUDED

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

class GridInfo;
class HeatTargetWindow;
/*
#include "utils.h"    //GridInfo
#include "heatTargetWindow.h"
*/

/**
 * @brief A class representing the medium on which the heat diffusion simulation is performed \n
 *        Inherits from Grid class that represents the informations about the size of the \n
 *        discrete grid.
 */
class Medium : public GridInfo
{
private:
  /**
   * @brief Initializes the medium to a begining state. Prepares all data for simulation.
   * 
   * @param [in] folderPath   - the path to the folder containing density.txt, heat.txt and conductivity.txt
   */
  void init(std::string folderPath);

public:
  /**
   * @brief volume rate of heat deposition [W/m^3]
   * 
   */
  float* q;

  /**
   * @brief temperature distribution [degC]
   * 
   */
  float* t;

  /**
   * @brief thermal dose given in cumulative equivalent minutes (cem) relative to T = 43 degC
   * 
   */
  float* cem43;

  /**
   * @brief tissue mass density [kg/m^3]
   * 
   */
  float* rho;

  /**
   * @brief tissue specific heat capacity [J/(kg.K)]
   * 
   */
  float* c;

  /**
   * @brief tissue thermal conductivity [W/(m.K)]
   * 
   */
  float* lambda;

  size_t matSize;
  
  Medium() = delete;

  /**
   * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
   *        Achieved by multiple matrices, each giving a value to a point on the grid. Complete medium \n
   *        as defined in the dat files without any windowing.
   * 
   * @param [in] folderPath   - the path to the folder containing density.txt, heat.txt and conductivity.txt
   */
  Medium(std::string folderPath);

  /**
   * @brief Construct a new Medium object representing the medium on which the sonication is done. \n
   *        Achieved by multiple matrices, each giving a value to a point on the grid and cropped by \n
   *        a given window represented by HeatTargetWindow class.
   * 
   * @param [in] targetWindow   - The instance of the HeatTargetWindow
   */
  Medium(HeatTargetWindow   targetWindow, 
         float* heat, 
         float* temperature, 
         float* density, 
         float* capacity,
         float* cem,
         float* conductivity);


  Medium(const Medium& cpy);

  ~Medium();


  void setHeatDeposit(float* heatAdded);

  /**
   * @brief Sets the heat deposit energy into the medium after cropping it by the provided window.
   * 
   * @param [in] heatTargetWindow       - Instance of HeatTargetWindow used to crop the heat doposition \n
   *                                      vector to the appropriate window.
   * @param [in] heatAdded              - Vector of data containing the values for the added heat energy \n
   *                                      onto the medium. Expected to be of the original size and to be \n
   *                                      cropped by the provided window.
   */
  void setHeatDepositEnergyOnWindow(HeatTargetWindow & heatTargetWindow, 
                                    float* heatAdded);

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
  static Medium updateOriginalByCroppedMedium(HeatTargetWindow & heatTargetWindow,
                                              Medium originalMedium,
                                              Medium & croppedMedium);

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
  Medium getWindowedMedium(HeatTargetWindow & heatTargetWindow);
};
//END of Medium class
#endif