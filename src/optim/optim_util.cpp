/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#include "dmpc/optim/optim_util.hpp"
#include <algorithm>

namespace dmpc
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
        interplin(&out.x_[0], &state.t_[0], &state.x_[0], t, Nx, Nt, 1);
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

void shiftState(AgentState& state, typeRNum dt, typeRNum t0)
{
    const unsigned int Nt = static_cast<unsigned int>(state.t_.size());
    if( Nt > 0 )
    {
        const unsigned int Nx = static_cast<unsigned int>(state.x_.size() / Nt);
        const unsigned int Nu = static_cast<unsigned int>(state.u_.size() / Nt);
        const unsigned int Nv = static_cast<unsigned int>(state.v_.size() / Nt);

        // shift x
        if(Nx > 0)
            shiftTrajectory(&state.x_[0], Nt, Nx, Nx, dt, &state.t_[0]);

        // shift v
        if(Nv > 0)
            shiftTrajectory(&state.v_[0], Nt, Nv, Nv, dt, &state.t_[0]);

        // shift u
        if(Nu > 0)
            shiftTrajectory(&state.u_[0], Nt, Nu, Nu, dt, &state.t_[0]);
	}

	// This maximum ensures that all agents have the same time base 
	// even in case of plug-and-play scenarios
	state.t0_ = std::max(state.t0_ + dt, t0);
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

const bool compare_stateDimensions( const AgentState& state_1, const AgentState& state_2 )
{
    return state_1.x_.size() == state_2.x_.size() && state_1.u_.size() == state_2.u_.size()
            && state_1.v_.size() == state_2.v_.size() && state_1.t_.size() == state_2.t_.size();
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

}
