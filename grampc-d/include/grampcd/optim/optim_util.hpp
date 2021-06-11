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

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{

    //*********************************************
    // interpolate state
    //*********************************************

    /*Compute interpolated states x(t) and u(t)*/
    void interpolateState(const AgentState& state, double t, AgentState& out);
    /*Compute interpolated states z_x(t) and z_u(t)*/
    void interpolateState(const CouplingState& state, double t, CouplingState& out);

    /*Compute interpolated states mu_x(t) and mu_u(t)*/
    void interpolateState(const MultiplierState& state, double t, MultiplierState& out);

    /*Compute interpolated states mu_x(t) and mu_u(t)*/
    void interpolateState(const PenaltyState& state, double t, PenaltyState& out);

    //*********************************************
    // shift states
    //*********************************************

    /*Shift agent states.*/
    void shiftState(AgentState& state, typeRNum dt, typeRNum t0);
    /*Shift coupling states.*/
    void shiftState(CouplingState& state, typeRNum dt, typeRNum t0);
    /*Shift multiplier states.*/
    void shiftState(MultiplierState& state, typeRNum dt, typeRNum t0);
    /*Shift penalty states.*/
    void shiftState(PenaltyState& state, typeRNum dt, typeRNum t0);

    //*********************************************
    // reset states
    //*********************************************

    /*Reset agent states.*/
    void resetState(AgentState& state, int i, std::vector<typeRNum> t);
    /*reset coupling states.*/
    void resetState(CouplingState& state, int i, std::vector<typeRNum> t);
    /*Reset multiplier states.*/
    void resetState(MultiplierState& state, int i, std::vector<typeRNum> t);
    /*Reset penalty states.*/
    void resetState(PenaltyState& state, int i, std::vector<typeRNum> t);

    //*********************************************
    // compare states
    //*********************************************

    /*Compare agent states.*/
    const bool compare_stateDimensions( const AgentState& state_1, const AgentState& state_2 );
    /*Compare coupling states.*/
    const bool compare_stateDimensions( const CouplingState& state_1, const CouplingState& state_2 );
    /*Compare multiplier states.*/
    const bool compare_stateDimensions( const MultiplierState& state_1, const MultiplierState& state_2 );
    /*Compare penalty states.*/
    const bool compare_stateDimensions( const PenaltyState& state_1, const PenaltyState& state_2 );

    /** @brief Definition of steps for the alternating direction method of multipliers */
    enum class ADMMStep
    {
        UPDATE_AGENT_STATE,
        SEND_AGENT_STATE,
        UPDATE_COUPLING_STATE,
        SEND_COUPLING_STATE,
        UPDATE_MULTIPLIER_STATE,
        SEND_MULTIPLIER_STATE,
        SEND_CONVERGENCE_FLAG,
        INITIALIZE,
        SEND_TRUE_STATE,
        PRINT
    };

    /*Configure the solver regarding the optimization info.*/
    void configureSolver(std::shared_ptr< grampc::Grampc>& solver, const OptimizationInfo& info);

}