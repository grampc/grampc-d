/* This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
 *
 * GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
 * 
 *
 * Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
 * All rights reserved.
 *
 * GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */

#include "grampcd/optim/optim_util.hpp"

#include "grampcd/state/agent_state.hpp"
#include "grampcd/state/coupling_state.hpp"
#include "grampcd/state/multiplier_state.hpp"
#include "grampcd/state/penalty_state.hpp"
#include "grampcd/state/sensi_state.hpp"
#include "grampcd/state/constraint_state.hpp"

#include "grampcd/info/optimization_info.hpp"

#include <algorithm>

namespace grampcd
{

    void interpolateState(const AgentState& state, double t, AgentState& out)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        const unsigned int Nx = static_cast<unsigned int>(state.x_.size() / Nt);
        const unsigned int Nu = static_cast<unsigned int>(state.u_.size() / Nt);
        const unsigned int Nv = static_cast<unsigned int>(state.v_.size() / Nt);

        out.i_ = state.i_;
        out.t_ = state.t_;

        // interpolate x
        if(Nx > 0)
        {
            out.x_.resize(Nx);
            out.lambda_.resize(Nx);
            interplin(&out.x_[0], &state.t_[0], &state.x_[0], t, Nx, Nt, 1);
            interplin(&out.lambda_[0], &state.t_[0], &state.lambda_[0], t, Nx, Nt, 1);
        }

        // interpolate v
        if(Nv > 0)
        {
            out.v_.resize(Nv);
            interplin(&out.v_[0], &state.t_[0], &state.v_[0], t, Nv, Nt, 1);
        }

        // interpolate u
        if(Nu > 0)
        {
            out.u_.resize(Nu);
            interplin(&out.u_[0], &state.t_[0], &state.u_[0], t, Nu, Nt, 1);
        }
    }

    void interpolateState(const CouplingState& state, double t, CouplingState& out)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        const unsigned int Nx = static_cast<unsigned int>(state.z_x_.size() / Nt);
        const unsigned int Nu = static_cast<unsigned int>(state.z_u_.size() / Nt);
        const unsigned int Nv = static_cast<unsigned int>(state.z_v_.size() / Nt);

        out.i_ = state.i_;
        out.t_ = state.t_;

        // interpolate x
        if(Nx > 0)
        {
            out.z_x_.resize(Nx);
            interplin(&out.z_x_[0], &state.t_[0], &state.z_x_[0], t, Nx, Nt, 1);
        }

        // interpolate v
        if(Nv > 0)
        {
            out.z_v_.resize(Nv);
            interplin(&out.z_v_[0], &state.t_[0], &state.z_v_[0], t, Nv, Nt, 1);
        }

        // interpolate u
        if(Nu > 0)
        {
            out.z_u_.resize(Nu);
            interplin(&out.z_u_[0], &state.t_[0], &state.z_u_[0], t, Nu, Nt, 1);
        }
    }

    void interpolateState(const MultiplierState& state, double t, MultiplierState& out)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        const unsigned int Nx = static_cast<unsigned int>(state.mu_x_.size() / Nt);
        const unsigned int Nu = static_cast<unsigned int>(state.mu_u_.size() / Nt);
        const unsigned int Nv = static_cast<unsigned int>(state.mu_v_.size() / Nt);

        out.i_ = state.i_;
        out.t_ = state.t_;

        // interpolate x
        if(Nx > 0)
        {
            out.mu_x_.resize(Nx);
            interplin(&out.mu_x_[0], &state.t_[0], &state.mu_x_[0], t, Nx, Nt, 1);
        }

        // interpolate v
        if(Nv > 0)
        {
            out.mu_v_.resize(Nv);
            interplin(&out.mu_v_[0], &state.t_[0], &state.mu_v_[0], t, Nv, Nt, 1);
        }

        // interpolate u
        if(Nu > 0)
        {
            out.mu_u_.resize(Nu);
            interplin(&out.mu_u_[0], &state.t_[0], &state.mu_u_[0], t, Nu, Nt, 1);
        }
    }

    void interpolateState(const PenaltyState& state, double t, PenaltyState& out)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        const unsigned int Nx = static_cast<unsigned int>(state.rho_x_.size() / Nt);
        const unsigned int Nu = static_cast<unsigned int>(state.rho_u_.size() / Nt);
        const unsigned int Nv = static_cast<unsigned int>(state.rho_v_.size() / Nt);

        out.i_ = state.i_;
        out.t_ = state.t_;

        // interpolate x
        if(Nx > 0)
        {
            out.rho_x_.resize(Nx);
            interplin(&out.rho_x_[0], &state.t_[0], &state.rho_x_[0], t, Nx, Nt, 1);
        }

        // interpolate v
        if(Nv > 0)
        {
            out.rho_v_.resize(Nv);
            interplin(&out.rho_v_[0], &state.t_[0], &state.rho_v_[0], t, Nv, Nt, 1);
        }

        // interpolate u
        if(Nu > 0)
        {
            out.rho_u_.resize(Nu);
            interplin(&out.rho_u_[0], &state.t_[0], &state.rho_u_[0], t, Nu, Nt, 1);
        }
    }

    void interpolateState(const SensiState& state, double t, SensiState& out)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        const unsigned int Nx = static_cast<unsigned int>(state.psi_x_.size() / Nt);
        const unsigned int Nu = static_cast<unsigned int>(state.psi_u_.size() / Nt);

        const unsigned int Nxx = static_cast<unsigned int>(state.psi_xx_.size() / Nt);
        const unsigned int Nuu = static_cast<unsigned int>(state.psi_uu_.size() / Nt);
        const unsigned int Nxu = static_cast<unsigned int>(state.psi_xu_.size() / Nt);

        out.i_ = state.i_;
        out.t_ = state.t_;

        // interpolate psi_x ( no interpolation needed for psi_V)
        if (Nx > 0)
        {
            out.psi_x_.resize(Nx);
            interplin(&out.psi_x_[0], &state.t_[0], &state.psi_x_[0], t, Nx, Nt, 1);
        }

        // interpolate psi_u
        if (Nu > 0)
        {
            out.psi_u_.resize(Nu);
            interplin(&out.psi_u_[0], &state.t_[0], &state.psi_u_[0], t, Nu, Nt, 1);
        }
        // interpolate psi_u
        if (Nxx > 0)
        {
            out.psi_xx_.resize(Nxx);
            interplin(&out.psi_xx_[0], &state.t_[0], &state.psi_xx_[0], t, Nxx, Nt, 1);
        }

        // interpolate psi_u
        if (Nuu > 0)
        {
            out.psi_uu_.resize(Nuu);
            interplin(&out.psi_uu_[0], &state.t_[0], &state.psi_uu_[0], t, Nuu, Nt, 1);
        }

        // interpolate psi_u
        if (Nxu > 0)
        {
            out.psi_xu_.resize(Nxu);
            interplin(&out.psi_xu_[0], &state.t_[0], &state.psi_xu_[0], t, Nxu, Nt, 1);
        }
    }

    void interpolateState(const ConstraintState& state, double t, ConstraintState& out)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        const unsigned int Ng = static_cast<unsigned int>(state.mu_g_.size() / Nt);
        const unsigned int Nh = static_cast<unsigned int>(state.mu_h_.size() / Nt);

        out.i_ = state.i_;
        out.t_ = state.t_;

        // interpolate psi_x ( no interpolation needed for psi_V)
        if (Ng > 0)
        {
            out.mu_g_.resize(Ng);
            interplin(&out.mu_g_[0], &state.t_[0], &state.mu_g_[0], t, Ng, Nt, 1);
            interplin(&out.c_g_[0], &state.t_[0], &state.mu_g_[0], t, Ng, Nt, 1);
        }

        // interpolate psi_u
        if (Nh > 0)
        {
            out.mu_h_.resize(Nh);
            interplin(&out.mu_h_[0], &state.t_[0], &state.mu_h_[0], t, Nh, Nt, 1);
            interplin(&out.c_h_[0], &state.t_[0], &state.mu_g_[0], t, Ng, Nt, 1);
        }
    }

    void shiftState(AgentState& state, typeRNum dt, typeRNum t0)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        if( Nt > 0 )
        {
            const unsigned int Nx = static_cast<unsigned int>(state.x_.size() / Nt);
            const unsigned int Nu = static_cast<unsigned int>(state.u_.size() / Nt);
            const unsigned int Nv = static_cast<unsigned int>(state.v_.size() / Nt);

            // shift x
            if (Nx > 0)
            {
                shiftTrajectory(&state.x_[0], Nt, Nx, Nx, dt, &state.t_[0]);
                shiftTrajectory(&state.lambda_[0], Nt, Nx, Nx, dt, &state.t_[0]);
            }

            // shift v
            if(Nv > 0)
                shiftTrajectory(&state.v_[0], Nt, Nv, Nv, dt, &state.t_[0]);

            // shift u
            if(Nu > 0)
                shiftTrajectory(&state.u_[0], Nt, Nu, Nu, dt, &state.t_[0]);
	    }

	    // This maximum ensures that all agents have the same time base 
	    // even in case of plug-and-play scenarios
	    state.t0_ = std::max(state.t0_ + dt, t0 + dt);
    }

    void shiftState(CouplingState& state, typeRNum dt, typeRNum t0)
    {
	    // This maximum ensures that all agents have the same time base 
	    // even in case of plug-and-play scenarios
        state.t0_ = std::max(state.t0_ + dt, t0);

        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        if( Nt > 0 )
        {
            const unsigned int Nx = static_cast<unsigned int>(state.z_x_.size() / Nt);
            const unsigned int Nu = static_cast<unsigned int>(state.z_u_.size() / Nt);
            const unsigned int Nv = static_cast<unsigned int>(state.z_v_.size() / Nt);

            // shift x
            if(Nx > 0)
                shiftTrajectory(&state.z_x_[0], Nt, Nx, Nx, dt, &state.t_[0]);

            // shift v
            if(Nv > 0)
                shiftTrajectory(&state.z_v_[0], Nt, Nv, Nv, dt, &state.t_[0]);

            // shift u
            if(Nu > 0)
                shiftTrajectory(&state.z_u_[0], Nt, Nu, Nu, dt, &state.t_[0]);
        }
    }

    void shiftState(MultiplierState& state, typeRNum dt, typeRNum t0)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        if( Nt > 0 )
        {
            const unsigned int Nx = static_cast<unsigned int>(state.mu_x_.size() / Nt);
            const unsigned int Nu = static_cast<unsigned int>(state.mu_u_.size() / Nt);
            const unsigned int Nv = static_cast<unsigned int>(state.mu_v_.size() / Nt);

            // shift x
            if(Nx > 0)
                shiftTrajectory(&state.mu_x_[0], Nt, Nx, Nx, dt, &state.t_[0]);

            // shift v
            if(Nv > 0)
                shiftTrajectory(&state.mu_v_[0], Nt, Nv, Nv, dt, &state.t_[0]);

            // shift u
            if(Nu > 0)
                shiftTrajectory(&state.mu_u_[0], Nt, Nu, Nu, dt, &state.t_[0]);
        }

	    // This maximum ensures that all agents have the same time base 
	    // even in case of plug-and-play scenarios
	    state.t0_ = std::max(state.t0_ + dt, t0);
    }

    void shiftState(PenaltyState& state, typeRNum dt, typeRNum t0)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        if( Nt > 0 )
        {
            const unsigned int Nx = static_cast<unsigned int>(state.rho_x_.size() / Nt);
            const unsigned int Nu = static_cast<unsigned int>(state.rho_u_.size() / Nt);
            const unsigned int Nv = static_cast<unsigned int>(state.rho_v_.size() / Nt);

            // shift x
            if(Nx > 0)
                shiftTrajectory(&state.rho_x_[0], Nt, Nx, Nx, dt, &state.t_[0]);

            // shift v
            if(Nv > 0)
                shiftTrajectory(&state.rho_v_[0], Nt, Nv, Nv, dt, &state.t_[0]);

            // shift u
            if(Nu > 0)
                shiftTrajectory(&state.rho_u_[0], Nt, Nu, Nu, dt, &state.t_[0]);
        }

	    // This maximum ensures that all agents have the same time base 
	    // even in case of plug-and-play scenarios
	    state.t0_ = std::max(state.t0_ + dt, t0);
    }

    void shiftState(SensiState& state, typeRNum dt, typeRNum t0)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        if (Nt > 0)
        {
            const unsigned int Nx = static_cast<unsigned int>(state.psi_x_.size() / Nt);
            const unsigned int Nu = static_cast<unsigned int>(state.psi_u_.size() / Nt);

            const unsigned int Nxx = static_cast<unsigned int>(state.psi_xx_.size() / Nt);
            const unsigned int Nuu = static_cast<unsigned int>(state.psi_uu_.size() / Nt);
            const unsigned int Nxu = static_cast<unsigned int>(state.psi_xu_.size() / Nt);

            // shift x
            if (Nx > 0)
            {
                shiftTrajectory(&state.psi_x_[0], Nt, Nx, Nx, dt, &state.t_[0]);
            }

            // shift u
            if (Nu > 0)
                shiftTrajectory(&state.psi_u_[0], Nt, Nu, Nu, dt, &state.t_[0]);

            // shift psi_xx
            if (Nxx > 0)
            {
                shiftTrajectory(&state.psi_xx_[0], Nt, Nxx, Nxx, dt, &state.t_[0]);
            }

            // interpolate psi_u
            if (Nuu > 0)
            {
                shiftTrajectory(&state.psi_uu_[0], Nt, Nuu, Nuu, dt, &state.t_[0]);
            }

            // interpolate psi_u
            if (Nxu > 0)
            {
                shiftTrajectory(&state.psi_xu_[0], Nt, Nxu, Nxu, dt, &state.t_[0]);
            }
        }

        // This maximum ensures that all agents have the same time base 
        // even in case of plug-and-play scenarios
        state.t0_ = std::max(state.t0_ + dt, t0 + dt);
    }

    void shiftState(ConstraintState& state, typeRNum dt, typeRNum t0)
    {
        const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
        if (Nt > 0)
        {
            const unsigned int Ng = static_cast<unsigned int>(state.mu_g_.size() / Nt);
            const unsigned int Nh = static_cast<unsigned int>(state.mu_h_.size() / Nt);

            // shift x
            if (Ng > 0)
            {
                shiftTrajectory(&state.mu_g_[0], Nt, Ng, Ng, dt, &state.t_[0]);
                shiftTrajectory(&state.c_g_[0], Nt, Ng, Ng, dt, &state.t_[0]);
            }

            // shift u
            if (Nh > 0)
            {
                shiftTrajectory(&state.mu_h_[0], Nt, Ng, Ng, dt, &state.t_[0]);
                shiftTrajectory(&state.c_h_[0], Nt, Ng, Ng, dt, &state.t_[0]);
            }
        }

        // This maximum ensures that all agents have the same time base 
        // even in case of plug-and-play scenarios
        state.t0_ = std::max(state.t0_ + dt, t0 + dt);
    }

    void configureSolver(std::shared_ptr< grampc::Grampc>& solver, const OptimizationInfo& info)
    {
        solver->setparam_real("Thor", info.COMMON_Thor_);
        solver->setparam_real("dt", info.COMMON_dt_);

        solver->setopt_string("ShiftControl", info.COMMON_ShiftControl_ ? "on" : "off");

        if(info.COMMON_Nhor_ > 0)
            solver->setopt_int("Nhor", info.COMMON_Nhor_);

        if(info.GRAMPC_MaxGradIter_ > 0)
            solver->setopt_int("MaxGradIter", info.GRAMPC_MaxGradIter_);

        if(info.GRAMPC_MaxMultIter_ > 0)
            solver->setopt_int("MaxMultIter", info.GRAMPC_MaxMultIter_);

        if( info.GRAMPC_PenaltyMin_ > 0 )
            solver->setopt_real("PenaltyMin", info.GRAMPC_PenaltyMin_);

        if( info.GRAMPC_PenaltyMax_ > 0 )
            solver->setopt_real("PenaltyMax", info.GRAMPC_PenaltyMax_);

        if( info.GRAMPC_AugLagUpdateGradientRelTol_ > 0 )
            solver->setopt_real("AugLagUpdateGradientRelTol", info.GRAMPC_AugLagUpdateGradientRelTol_);

        if( info.GRAMPC_LineSearchMax_ > 0 )
            solver->setopt_real("LineSearchMax", info.GRAMPC_LineSearchMax_);

        if( info.GRAMPC_LineSearchMin_ > 0 )
            solver->setopt_real("LineSearchMin", info.GRAMPC_LineSearchMin_);

        if( info.GRAMPC_PenaltyIncreaseFactor_ > 0 )
            solver->setopt_real("PenaltyIncreaseFactor", info.GRAMPC_PenaltyIncreaseFactor_);

        if( info.GRAMPC_PenaltyDecreaseFactor_ > 0 )
            solver->setopt_real("PenaltyDecreaseFactor", info.GRAMPC_PenaltyDecreaseFactor_);

        if( info.GRAMPC_Integrator_.size() > 0 )
            solver->setopt_string("Integrator", info.GRAMPC_Integrator_.c_str());

        if( info.GRAMPC_LineSearchType_.size() > 0 )
            solver->setopt_string("LineSearchType", info.GRAMPC_LineSearchType_.c_str());

        if( info.GRAMPC_ConvergenceCheck_.size() > 0 )
            solver->setopt_string("ConvergenceCheck", info.GRAMPC_ConvergenceCheck_.c_str());

        if( info.GRAMPC_ConstraintsAbsTol_.size() > 0 )
            solver->setopt_real_vector("ConstraintsAbsTol", info.GRAMPC_ConstraintsAbsTol_.data());

        if( info.GRAMPC_ConvergenceGradientRelTol_ > 0 )
            solver->setopt_real("ConvergenceGradientRelTol", info.GRAMPC_ConvergenceGradientRelTol_);

        if( info.GRAMPC_PenaltyIncreaseThreshold_ > 0 )
            solver->setopt_real("PenaltyIncreaseThreshold", info.GRAMPC_PenaltyIncreaseThreshold_);

        if( info.GRAMPC_LineSearchInit_ > 0 )
            solver->setopt_real("LineSearchInit", info.GRAMPC_LineSearchInit_);
    }

    void resetState(AgentState& state, int i, std::vector<typeRNum> t)
    {
        state.i_ = i;
        state.t_ = t;
        state.t0_ = 0;
        state.u_.clear();
        state.v_.clear();
        state.x_.clear();
        state.lambda_.clear();
    }

    void resetState(CouplingState& state, int i, std::vector<typeRNum> t)
    {
        state.i_ = i;
        state.t_ = t;
        state.z_u_.clear();
        state.z_v_.clear();
        state.z_x_.clear();
    }

    void resetState(MultiplierState& state, int i, std::vector<typeRNum> t)
    {
        state.i_ = i;
        state.t_ = t;
        state.mu_u_.clear();
        state.mu_v_.clear();
        state.mu_x_.clear();
    }

    void resetState(PenaltyState& state, int i, std::vector<typeRNum> t)
    {
        state.i_ = i;
        state.t_ = t;
        state.rho_u_.clear();
        state.rho_v_.clear();
        state.rho_x_.clear();
    }

    void resetState(SensiState& state, int i, std::vector<typeRNum> t)
    {
        state.i_ = i;
        state.t_ = t;
        state.t0_ = 0;
        state.psi_u_.clear();
        state.psi_x_.clear();
        state.psi_V_.clear();
        state.psi_xx_.clear();
        state.psi_uu_.clear();
        state.psi_xu_.clear();
        state.psi_VV_.clear();
    }

    void resetState(ConstraintState& state, int i, std::vector<typeRNum> t)
    {
        state.i_ = i;
        state.t_ = t;
        state.t0_ = 0;
        state.mu_g_.clear();
        state.mu_h_.clear();
        state.c_g_.clear();
        state.c_h_.clear();
    }

    const bool compare_stateDimensions( const AgentState& state_1, const AgentState& state_2 )
    {
        return state_1.x_.size() == state_2.x_.size() && state_1.u_.size() == state_2.u_.size()
                && state_1.v_.size() == state_2.v_.size() && state_1.t_.size() == state_2.t_.size() 
                && state_1.lambda_.size() == state_2.lambda_.size();
    }

    const bool compare_stateDimensions( const CouplingState& state_1, const CouplingState& state_2 )
    {
        return state_1.z_x_.size() == state_2.z_x_.size() && state_1.z_u_.size() == state_2.z_u_.size()
                && state_1.z_v_.size() == state_2.z_v_.size() && state_1.t_.size() == state_2.t_.size();
    }

    const bool compare_stateDimensions( const MultiplierState& state_1, const MultiplierState& state_2 )
    {
        return state_1.mu_x_.size() == state_2.mu_x_.size() && state_1.mu_u_.size() == state_2.mu_u_.size()
                && state_1.mu_v_.size() == state_2.mu_v_.size() && state_1.t_.size() == state_2.t_.size();
    }

    const bool compare_stateDimensions( const PenaltyState& state_1, const PenaltyState& state_2 )
    {
        return state_1.rho_x_.size() == state_2.rho_x_.size() && state_1.rho_u_.size() == state_2.rho_u_.size()
                && state_1.rho_v_.size() == state_2.rho_v_.size() && state_1.t_.size() == state_2.t_.size();
    }

    const bool compare_stateDimensions(const SensiState& state_1, const SensiState& state_2)
    {
        return state_1.psi_x_.size() == state_2.psi_x_.size() && state_1.psi_u_.size() == state_2.psi_u_.size()
            && state_1.psi_V_.size() == state_2.psi_V_.size() && state_1.psi_xx_.size() == state_2.psi_xx_.size()
            && state_1.psi_uu_.size() == state_2.psi_uu_.size() && state_1.psi_xu_.size() == state_2.psi_xu_.size()
            && state_1.psi_VV_.size() == state_2.psi_VV_.size() && state_1.t_.size() == state_2.t_.size();
    }

    const bool compare_stateDimensions(const ConstraintState& state_1, const ConstraintState& state_2)
    {
        return state_1.mu_g_.size() == state_2.mu_g_.size() && state_1.mu_h_.size() == state_2.mu_h_.size()
            && state_1.c_g_.size() == state_2.c_g_.size() && state_1.c_h_.size() == state_2.c_h_.size() 
            && state_1.t_.size() == state_2.t_.size();
    }

}
