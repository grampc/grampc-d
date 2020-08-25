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

#ifndef NEIGHBOR_HPP
#define NEIGHBOR_HPP

#include "dmpc/model/agent_model.hpp"
#include "dmpc/model/coupling_model.hpp"

#include "dmpc/util/logging.hpp"

#include "dmpc/optim/approximate_neighbor.hpp"

namespace dmpc
{

/**
 * @brief Coupled agents are implemented in form of a neighbor.
 */
class Neighbor
{
public:

    Neighbor(const int id, const AgentModelPtr& agent_model, const LoggingPtr& log);

    //*********************************************
    // get and set functions
    //*********************************************

    /*Returns the neighbors id.*/
    const int get_id() const;

    /*Returns the number of controls of the agent.*/
    const unsigned int get_Nui() const;
    /*Returns the number of states of the agent.*/
    const unsigned int get_Nxi() const;
    /*Returns the number of controls of the neighbor.*/
    const unsigned int get_Nuj() const;
    /*Returns the number of states of the neighbor.*/
    const unsigned int get_Nxj() const;

    /*Set number of neighbors.*/
    void set_numberOfNeighbors( const unsigned int number );
    /*Returns number of neighbors.*/
    const unsigned int get_numberOfNeighbors() const;

    /*Returns the neighbors agent model.*/
    const AgentModelPtr& get_agentModel() const;
    /*Returns the neighbors coupling model.*/
    const CouplingModelPtr& get_couplingModel() const;
    /*Returns the copied coupling model for neighbor approximation.*/
    const CouplingModelPtr& get_copied_couplingModel() const;
    /*Set the neighbors coupling model.*/
    void set_couplingModel( const CouplingModelPtr& model );
    /*Set the copied coupling model for neighbor approximation.*/
    void set_copied_couplingModel( const CouplingModelPtr& model );

    //*********************************************
    // get states
    //*********************************************

    /*Returns the local copies for the neighbor.*/
    const AgentState& get_localCopies() const;
    /*Returns the multiplier regarding the coupling with the neighbor.*/
    const MultiplierState& get_coupled_multiplierState() const;
    /*Returns the penalties regarding the coupling with the neighbor.*/
    const PenaltyState& get_coupled_penaltyState() const;
    /*Returns the neighbors local copies.*/
    const AgentState& get_neighbors_localCopies() const;
    /*Returns the neighbors coupling state.*/
    const CouplingState& get_neighbors_couplingState() const;
    /*Returns the multipliers regarding the neighbors coupling.*/
    const MultiplierState& get_neighbors_coupled_multiplierState() const;
    /*Returns the penalties regarding the neighbors coupling.*/
    const PenaltyState& get_neighbors_coupled_penaltyState() const;
    /*Returns the neighbors desired agent state.*/
    const AgentState& get_neighbors_desiredAgentState() const;
    /*Returns the external influence agent state.*/
    const AgentState& get_externalInfluence_agentState() const;
    /*Returns the external influence coupling state.*/
    const CouplingState& get_externalInfluence_couplingState() const;
    /*Returns the external influence multiplier state.*/
    const MultiplierState& get_externalInfluence_multiplierState() const;
    /*Returns the external influence penalty state.*/
    const PenaltyState& get_externalInfluence_penaltyState() const;

    /*Returns the neighbors external influence coupling state.*/
    const CouplingState& get_neighbors_externalInfluence_couplingState() const;
    /*Returns the neighbors external influence multiplier state.*/
    const MultiplierState& get_neighbors_externalInfluence_multiplierState() const;
    /*Returns the neighbors external influence penalty state.*/
    const PenaltyState& get_neighbors_externalInfluence_penaltyState() const;

    /*Returns the previous external influence multiplier state.*/
    const MultiplierState& get_previous_externalInfluence_multiplierState() const;
    /*Returns the previous external influence coupling state.*/
    const CouplingState& get_previous_externalInfluence_couplingState() const;
    /*Returns the previous coupling state.*/
    const CouplingState& get_previous_couplingState() const;
    /*Returns the previous multiplier state.*/
    const MultiplierState& get_previous_multiplierState() const;
    /*Returns the neighbors previous external influence coupling state.*/
    const CouplingState& get_previous_neighbors_externalInfluence_couplingState() const;
    /*Returns the neighbors previous coupling state.*/
    const CouplingState& get_previous_neighbors_couplingState() const;

    /*Returns the neighbors desired agent state.*/
    const AgentState& get_desiredAgentState() const;

    /*Returns description for neighbor approximation.*/
    const ApproximateNeighborPtr get_neighborApproximation() const;

    //*********************************************
    // set states
    //*********************************************

    /*Set the local copies.*/
    void set_localCopies( const AgentState& state );
    /*Set multiplier states regarding the coupling with the neighbor.*/
    void set_coupled_multiplierState( const MultiplierState& multiplier );
    /*Set penalty states regarding the coupling with the neighbor.*/
    void set_coupled_penaltyState( const PenaltyState& penalty );

    /*Set neighbors local copies.*/
    void set_neighbors_localCopies( const AgentState& state );
    /*Set neighbors coupling state.*/
    void set_neighbors_couplingState( const CouplingState& coupling );
    /*Set multiplier states regarding the neighbors coupling.*/
    void set_neighbors_coupled_multiplierState( const MultiplierState& multiplier );
    /*Set penalty states regarding the neighbors coupling.*/
    void set_neighbors_coupled_penaltyState( const PenaltyState& penalty );
    /*Set neighbors desired agent state.*/
    void set_neighbors_desiredAgentState( const AgentState& state );

    /*Set external influence agent state.*/
    void set_externalInfluence_agentState( const AgentState& state );
    /*Set external influence coupling state.*/
    void set_externalInfluence_couplingState( const CouplingState& coupling );
    /*Set external influence multiplier state.*/
    void set_externalInfluence_multiplierState( const MultiplierState& multiplier );
    /*Set external influence penalty state.*/
    void set_externInfluence_penaltyState( const PenaltyState& penalty );

    /*Set neighbors external influence coupling state.*/
    void set_neighbors_externalInfluence_couplingState( const CouplingState& coupling );
    /*Set neighbors external influence multiplier state.*/
    void set_neighbors_externalInfluence_multiplierState( const MultiplierState& multiplier );
    /*Set neighbors external influence penalty state.*/
    void set_neighbors_externalInfluence_penaltyState( const PenaltyState& penalty );

    //*********************************************
    // basic functions
    //*********************************************

    /*Reset all of the neighbors states.*/
    void reset_neighborStates();

    /*Initialize the neighbor.*/
    void initialize_neighbor(const std::vector<typeRNum>& t, const OptimizationInfo& optimization_info, Agent* agent);

    /*Return true if neighbor is sending neighbor.*/
    const bool is_sendingNeighbor() const;
    /*Return true if neighbor is receiving neighbor.*/
    const bool is_receivingNeighbor() const;
    /*Return true if neighbor is using neighbor approximation.*/
    const bool is_approximating() const;
    /*Return true if neighbor is approximation cost.*/
    const bool is_approximatingCost() const;
    /*Return true if neighbor is approximating constraints.*/
    const bool is_approximatingConstraints() const;
    /*Return true if neighbor is approximating dynamics.*/
    const bool is_approximatingDynamics() const;

    /*Define neighbor as sending neighbor.*/
    void defineAs_sendingNeighbor( const bool isSending );
    /*Define neighbor as receiving neighbor.*/
    void defineAs_receivingNeighbor( const bool isReceiving );

    /*Shift the neighbors states.*/
    void shift_states(const typeRNum dt, const typeRNum t0);

private:
    LoggingPtr log_;

    /************************
     Neighbor parameters
    ************************/

    int id_;
    int agent_id_;
    bool sending_neighbor_;
    bool receiving_neighbor_;

    unsigned int numberOfNeighbors_ = 0;
    AgentModelPtr agentModel_;
    CouplingModelPtr couplingModel_;
    CouplingModelPtr copied_couplingModel_;

    /************************
     States for neighbor
    ************************/

    AgentState local_copies_;
    MultiplierState coupled_multiplierState_;
    PenaltyState coupled_penaltyState_;

    AgentState neighbors_localCopies_;
    CouplingState neighbors_couplingState_;
    MultiplierState neighbors_coupled_multiplierState_;
    PenaltyState neighbors_coupled_penaltyState_;

    AgentState externalInfluence_agentState_;
    CouplingState externalInfluence_couplingState_;
    MultiplierState externalInfluence_multiplierState_;
    PenaltyState externalInfluence_penaltyState_;

    CouplingState neighbors_externalInfluence_couplingState_;
    MultiplierState neighbors_externalInfluence_multiplierState_;
    PenaltyState neighbors_externalInfluence_penaltyState_;

    CouplingState previous_couplingState_;
    MultiplierState previous_multiplierState_;

    CouplingState previous_externalInfluence_couplingState_;
    MultiplierState previous_externalInfluence_multiplierState_;
    CouplingState previous_neighbors_externalInfluence_couplingState_;
    CouplingState previous_neighbors_couplingState_;

    /************************
     States for cost approximation
    ************************/

    AgentState neighbors_desiredAgentState_;

    /************************
     Neighbor approximation
    ************************/

    ApproximateNeighborPtr approximate_neighbor_;

    bool is_approximating_ = false;
    bool is_approximatingCost_ = false;
    bool is_approximatingConstraints_ = false;
    bool is_approximatingDynamics_ = false;

};

typedef std::shared_ptr<Neighbor> NeighborPtr;

}

#endif // NEIGHBOR_HPP
