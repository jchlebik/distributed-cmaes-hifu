/**
 * @file      kWaveDiffusionSolver.h
 *
 * @author    Filip Kuklis \n
 *            Faculty of Information Technology\n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing the main class of the project responsible for the entire simulation.
 *
 * @version   KSpacePlaneRecon
 *
 * @date      05 November      2018, 09:45 (created) \n
 *            12 November      2018, 15:00 (revised)
 *
 * @copyright Copyright (C) 2018 Filip Kuklis and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */

#ifndef KWAVE_DIFFUSION__2D_SOLVER_H
#define KWAVE_DIFFUSION__2D_SOLVER_H

#include <fftw3.h>
#include <omp.h>
#include <sstream>
#include <complex>
//#include "../Parameters/CommandLineParameters.h"


/**
 * @class   kWaveDiffusionSolver
 * @brief   Class responsible for running the k-Wave diffusion solver 2D.
 * @details Class responsible for running the k-Wave diffusion solver 2D. This class maintain
 *          the whole simulation.
 *
 */
class kWaveDiffusionSolver2D {

public:
  /// Constructor.
  kWaveDiffusionSolver2D(int NX, 
                         int NY, 
                         const float* mediumDensity, 
                         const float* mediumSpecificHeat, 
                         const float* mediumThermalConductivity);

  /// Copy constructor not allowed for public.
  kWaveDiffusionSolver2D(const kWaveDiffusionSolver2D &) = delete;

  /// Destructor.
  ~kWaveDiffusionSolver2D();

  /// operator= not allowed for public.
  kWaveDiffusionSolver2D &operator=(const kWaveDiffusionSolver2D &) = delete;

  /// Memory allocation.
  void allocateMemory();

  /// Sets the data for sonication
  void setup(float* addedQ,
             float* mediumTemperature,
             float* objCem43);

  /// Memory deallocation.
  void freeMemory();

  /**
   * @brief Load simulation data.
   *
   * Copy data from matlab to local allocated memory (mex input array should not be changed)
   */
   void loadInputData();

  /**
   * @brief Create diffusion terms
   *
   * Compute diffusion terms from medium parameters
   */
   void createDiffusionTerms();

  /**
   * @brief This method computes k-Wave Diffusion simulation.
   *
   * This method computes k-Wave Diffusion 2D simulation. It launches calculation on a given
   * dataset going through FFT initialization, pre-processing, main loop and post-processing phases.
   */
  void takeTimeStep(int nt, bool useQTerm, float dt);



    /**
   * @brief  Get memory usage in MB on the host side.
   * @return Memory consumed on the host side in MB.
   */
   size_t getMemoryUsage() const;

  /**
   * @brief  Get code name - release code version.
   * @return Release code version.
   */
  std::string getCodeName() const;

  /// Print the code name and license.
  void printFullCodeNameAndLicense() const;


protected:

  /// Initialize FFTW plans.
  void InitializeFftwPlans();

  /**
   * @brief Compute pre-processing phase.
   *
   * Initialize all indices, pre-compute constants such as c^2, rho0Sgx * dt  and create kappa,
   * absorbEta, absorbTau, absorbNabla1, absorbNabla2  matrices.
   *
   */
  void preProcessing();

  /// make Wavenumber vectors
  template<bool reduced>
  void makeKVec(float* k_vec, float  d, size_t N);

  /// make cartesian spatial different operators
  void makeDerivMatrix();

  /// make K wavenumbers matrix squared
  void makeK2();

  /// make k-space term kappa
  void makeKappa();


  /// Compute PTerm
  void createPTerm();

  /// Compute QTerm using scaling
  void createQTermScale();

  /// Compute QTerm using diffusion
  void createQTermDiff();

  /// Compute the main time loop of the kWaveDiffusionSolver.
  void computeMainLoop();

  /// Compute the dterm for homogenous simulation
  void computeDTermHomogenous();

  /// Compute the dterm for heterogenous simulation
  void computeDTermHeterogenous();

  /// Post processing, and closing the output streams.
  void postProcessing();



private:
  //CommandLineParameters &mParameters;
  ///--------- Pointers to matrices created in MATLAB! ----------///
  // mObjQ is not const because it is used in fftw execute (but it should not be changed as well)
  // if not used with Matlab, output matrices are not necessery, temperature and cem43 matrices could be used as input/output (if you created it on your own outside the object)
  float* mObjQ;
  float* mObjTemperature;
  float* mObjCem43;
  const float* mMediumDensity;
  const float* mMediumSpecificHeat;
  const float* mMediumThermalConductivity;

  ///--------- Pointers to matrices created and allocated in C++ code ----------///
  /// medium temperature
  float* mTemperature;
  /// Q source term
  float* mQTerm;
  /// perfusion term
  float* mPTerm;
  /// diffusion terms
  float* mDiffusionP1;
  float* mDiffusionP2;
  /// cumulative equivalent minutes at temperature 43°C
  float* mCem43;
  /// k-space term kappa
  float* mKappa;
  /// wavenumbers K matrix squared
  float* mK2;
  /// wavenumbers Kx vector
  float* mKx_vec;
  /// wavenumbers Ky vector
  float* mKy_vec;


  /// temporary real matrices
  float* mTempRealArray1;
  float* mTempRealArray2;

  /// Cartesian spatial different operators
  std::complex<float>* mDerivX;
  std::complex<float>* mDerivY;

  /// temporary complex matrices
  std::complex<float>* mTempFftwArray1;
  std::complex<float>* mTempFftwArray2;

  /// FFTW plans
  fftwf_plan mFftwPlanForward;
  fftwf_plan mFftwPlanBackward;

  /// reference diffusion coefficient
  float mDdiffusionCoeffRef;
  float mQScaleFac;
  /// Dimensions
  int mNx;
  int mNxr;
  int mNy;
  int mNt;

  float mDt;

  /// flags
  bool mUseQterm;
};// end of kWaveDiffusionSolver
//----------------------------------------------------------------------------------------------------------------------

#endif	/* KWAVE_DIFFUSION__2D_SOLVER_H */
