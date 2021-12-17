/* This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
 *
 * GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
 * based on the alternating direction method of multipliers (ADMM).
 *
 * Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
 * All rights reserved.
 *
 * GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */

#pragma once

#include "grampcd/util/types.hpp"

namespace grampcd
{

    /**
     * @brief Information about the optimization.
     */
    struct OptimizationInfo
    {
    public:
        OptimizationInfo() {}

        /*Optimization horizon*/
		typeRNum COMMON_Thor_ = 1;
        /*Step size*/
        typeRNum COMMON_dt_ = 0.05;
        /*Number of discretization points*/
        unsigned int COMMON_Nhor_ = 21;
        /*Flag for shift control*/
        bool COMMON_ShiftControl_ = true;
        /*Used integration method*/
        std::string COMMON_Integrator_ = "heun";

        /*Maximum number of gradient iterations*/
        unsigned int GRAMPC_MaxGradIter_ = 0;
        /*Maximum number of multiplier iterations*/
        unsigned int GRAMPC_MaxMultIter_ = 0;
        /*Minimum value for penalty parameter*/
        typeRNum GRAMPC_PenaltyMin_ = 0;
        /*Maximum value for penalty parameter*/
        typeRNum GRAMPC_PenaltyMax_ = 0;
        /*Threshold for the maximum relative gradient of the inner minimization problem*/
        typeRNum GRAMPC_AugLagUpdateGradientRelTol_ = 0;
        /*Used integration method*/
        std::string GRAMPC_Integrator_ = "";
        /*Method to solve the line search problem*/
        std::string GRAMPC_LineSearchType_ = "";
        /*Penalty increase factor*/
        typeRNum GRAMPC_PenaltyIncreaseFactor_ = 0;
        /*Penalty decrease factor*/
        typeRNum GRAMPC_PenaltyDecreaseFactor_ = 0;
        /*Maximum value for line search*/
        typeRNum GRAMPC_LineSearchMax_ = 0;
        /*Minimum value for line search*/
        typeRNum GRAMPC_LineSearchMin_ = 0;
        /*Activates the convergence check*/
        std::string GRAMPC_ConvergenceCheck_ = "";
        /*Threshold for the convergence check*/
        typeRNum GRAMPC_ConvergenceGradientRelTol_ = 0;
        /*Threshold for convergence of constraints*/
        std::vector<typeRNum> GRAMPC_ConstraintsAbsTol_ = std::vector<typeRNum>();
        /*Threshold to increase penalty*/
        typeRNum GRAMPC_PenaltyIncreaseThreshold_ = 0;
        /*Initial value for line search problem*/
        typeRNum GRAMPC_LineSearchInit_ = 0;

        /*Maximum number of ADMM iterations*/
	    unsigned int ADMM_maxIterations_ = 20;
        /*Maximum number of inner iterations*/
	    unsigned int ADMM_innerIterations_ = 1;
        /*Convergence tolerance*/
        typeRNum ADMM_ConvergenceTolerance_ = 0.02;
        /*Penalty increase factor*/
        typeRNum ADMM_PenaltyIncreaseFactor_ = 1.5;
        /*Penalty decrease factor*/
        typeRNum ADMM_PenaltyDecreaseFactor_ = 0.75;
        /*Minimum value for penalty*/
        typeRNum ADMM_PenaltyMin_ = 1e-4;
        /*Maximum value for penalty*/
        typeRNum ADMM_PenaltyMax_ = 1e4;
        /*Initial value for penalty*/
        typeRNum ADMM_PenaltyInit_ = 1;
        /*Activate penalty adaption*/
        bool ADMM_AdaptPenaltyParameter_ = true;
        /*Activate debug cost*/
        bool ADMM_DebugCost_ = false;

        /*Activate approximation of neighbors cost*/
        bool APPROX_ApproximateCost_ = false;
        /*Activate approximation of neighbors constraints*/
        bool APPROX_ApproximateConstraints_ = false;
        /*Activate approximation of neighbors dynamics*/
		bool APPROX_ApproximateDynamics_ = false;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(COMMON_Thor_, COMMON_dt_, COMMON_Nhor_, COMMON_ShiftControl_, COMMON_Integrator_,
                GRAMPC_MaxGradIter_, GRAMPC_MaxMultIter_, GRAMPC_PenaltyMin_, GRAMPC_PenaltyMax_, 
                GRAMPC_AugLagUpdateGradientRelTol_, GRAMPC_Integrator_, GRAMPC_LineSearchType_, 
                GRAMPC_PenaltyIncreaseFactor_, GRAMPC_PenaltyDecreaseFactor_, GRAMPC_LineSearchMax_, 
                GRAMPC_LineSearchMin_, GRAMPC_ConvergenceCheck_, GRAMPC_ConvergenceGradientRelTol_, 
                GRAMPC_ConstraintsAbsTol_, GRAMPC_PenaltyIncreaseThreshold_, GRAMPC_LineSearchInit_,
                ADMM_maxIterations_, ADMM_innerIterations_, ADMM_ConvergenceTolerance_, 
                ADMM_PenaltyIncreaseFactor_, ADMM_PenaltyDecreaseFactor_, ADMM_PenaltyMin_, 
                ADMM_PenaltyMax_, ADMM_PenaltyInit_, ADMM_AdaptPenaltyParameter_, ADMM_DebugCost_,
                APPROX_ApproximateCost_, APPROX_ApproximateConstraints_, APPROX_ApproximateDynamics_);
		}
    };

}