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

    /*Compute interpolated states psi_x(t) and psi_u(t)*/
    void interpolateState(const SensiState& state, double t, SensiState& out);

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
    /*Shift Sensi states.*/
    void shiftState(SensiState& state, typeRNum dt, typeRNum t0);
    /*Shift Constraint states.*/
    void shiftState(ConstraintState& state, typeRNum dt, typeRNum t0);

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
    /*Reset Sensi states.*/
    void resetState(SensiState& state, int i, std::vector<typeRNum> t);
    /*Reset Constraint states.*/
    void resetState(ConstraintState& state, int i, std::vector<typeRNum> t);

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
    /*Compare Sensi states.*/
    const bool compare_stateDimensions(const SensiState& state_1, const SensiState& state_2);
    /*Compare Constraint states.*/
    const bool compare_stateDimensions(const ConstraintState& state_1, const ConstraintState& state_2);

    /** @brief Definition of steps for the alternating direction method of multipliers and the sensitivity-based Algorithm*/
    enum class AlgStep
    {
        ADMM_UPDATE_AGENT_STATE,
        ADMM_SEND_AGENT_STATE,
        ADMM_UPDATE_COUPLING_STATE,
        ADMM_SEND_COUPLING_STATE,
        ADMM_UPDATE_MULTIPLIER_STATE,
        ADMM_SEND_MULTIPLIER_STATE,
        ADMM_INITIALIZE,
        ADMM_SEND_TRUE_STATE,
        ADMM_START_ASYNC_ADMM,
        SENSI_INITIALIZE,
        SENSI_UPDATE_AGENT_STATE,
        SENSI_SEND_AGENT_STATE,
        SENSI_UPDATE_SENSI_STATE,
        SENSI_SEND_SENSI_STATE,
        SENSI_START_ASYNC_SENSI,
        GEN_PRINT,
        GEN_SEND_CONVERGENCE_FLAG,
    };

    /*Configure the solver regarding the optimization info.*/
    void configureSolver(std::shared_ptr< grampc::Grampc>& solver, const OptimizationInfo& info);

}