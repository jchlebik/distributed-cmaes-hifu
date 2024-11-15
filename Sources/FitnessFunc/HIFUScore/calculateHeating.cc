/**
 * @file        calculateHeating.cc
 * 
 * @author      Jakub Chlebik \n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              xchleb07@stud.fit.vutbr.com
 * 
 * @brief       The source file for calculating the heating equations
 * 
 * @version     0.1
 * 
 * @date        2020-03-20 (created) \n
 */

#include "utils.h"
#include "medium.h"
#include "calculateHeating.h"

/**
 * @brief Compute heating diffusion for simple trajectory optimisation problem.
 * 
 * @param [in, out] medium              - the structure representing the medium where the heat is being diffused \n
 *                                        contains all relevant matrices like temperature and so on. the results \n
 *                                        of the simulation are stored in their respective matrices in this class
 * @param [in]      sonicationParams    - parameters of the simulation, the amount of steps and step size
 */
 void calculateHeating(Medium & medium, 
                       Sonication sonicationParams)
 {
   // create output containers
   //std::vector<float> output_cem43(medium.cem43.size());
   //std::vector<float> output_temp(medium.t.size());

  // create the diffusion solver as static to avoid uneccessary computation
  static kWaveDiffusionSolver2D kWaveDiffSolver(
      medium.nx,                    // medium grid size on x axis
      medium.ny,                    // medium grid size on y axis
      medium.rho,                   // medium density matrix - tissue mass density [kg/m^3]
      medium.c,                     // medium specific heat capacity matrix - tissue specific heat capacity [J/(kg.K)]
      medium.lambda                // medium thermal conductivity - tissue thermal conductivity [W/(m.K)]
    );

  kWaveDiffSolver.setup(medium.q, medium.t, medium.cem43);

  // perform the sonication
  int onSteps = sonicationParams.getOnSteps();
  if (onSteps != 0)
    kWaveDiffSolver.takeTimeStep(onSteps,  true,  sonicationParams.dT); // perform heating

  //medium.q = std::vector<float>(medium.q.size(), 0.0f);
  std::fill(medium.q, medium.q + medium.matSize, 0.0f);
  int offSteps = sonicationParams.getOffSteps();
  if (offSteps != 0)
    kWaveDiffSolver.takeTimeStep(offSteps, false, sonicationParams.dT); // cool down

  // move semantics to store the output and update the medium state
  //medium.t      = std::move(output_temp); 
  //medium.cem43  = std::move(output_cem43);
 }
// END of calculateHeating