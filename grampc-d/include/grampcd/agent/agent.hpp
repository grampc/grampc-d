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

#include "grampcd/info/agent_info.hpp"
#include "grampcd/info/coupling_info.hpp"
#include "grampcd/info/optimization_info.hpp"

#include "grampcd/state/agent_state.hpp"
#include "grampcd/state/coupling_state.hpp"
#include "grampcd/state/multiplier_state.hpp"
#include "grampcd/state/penalty_state.hpp"

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{  
    /**
     * @brief The agent is a single node in the network.
     * All communication with other agents is handled over the communication interface.
     */
    class Agent
    {
    public:

        Agent(const CommunicationInterfacePtr& communication_interface,
              const ModelFactoryPtr& model_factory,
              const AgentInfo& agent_info, 
            const LoggingPtr& log);

        /*Returns the agents id.*/
        const int get_id() const;

        /*Returns the agents number of states.*/
        const unsigned int get_Nxi() const;
        /* Returns the agents number of controls.*/
        const unsigned int get_Nui() const;

        /* Returns true if the agent is using neighbor approximation.*/
        const bool is_approximating() const;
        /*Returns true if the agent is approximating cost.*/
        const bool is_approximatingCost() const;
	    /*Returns true if the agent is approximating constraints.*/
        const bool is_approximatingConstraints() const;
	    /*Returns true if the agent is approximating dynamics.*/
        const bool is_approximatingDynamics() const;

        /*Returns the agent info.*/
        const AgentInfo& get_agentInfo() const;
        /*Return the agent model.*/
        const AgentModelPtr& get_agentModel() const;
        /*Returns the optimization info.*/
        const OptimizationInfo& get_optimizationInfo() const;
        /*returns the step selector*/
        const StepSelectorPtr& get_stepSelector() const;

        /*************************************************************************
         neighbor functions
         *************************************************************************/

        /*Add a neighbor.*/
        const bool add_neighbor(const CouplingInfo& coupling_info, const AgentInfo& neighbor_info);
        /*Add a receiving neighbor.*/
        const bool add_receivingNeighbor(const CouplingInfo& coupling_info, const AgentInfo& neighbor_info);
        /*Add a sending neighbor.*/
        const bool add_sendingNeighbor(const CouplingInfo& coupling_info, const AgentInfo& neighbor_info);
        /*Remove a neighbor.*/
        void remove_neighbor( const CouplingInfo& coupling_info );

        /*************************************************************************
         state functions
         *************************************************************************/

        /*Initialize the agent.*/
        void initialize(const OptimizationInfo& optimization_info);

        /*Shift the agents states.*/
        void shift_states(const typeRNum dt, const typeRNum t0);

        /*Returns the agents desired agent state.*/
        const AgentState& get_desiredAgentState() const;

        /*Sets the agents desired agent state.*/
	    void set_desiredAgentState(const AgentState& state);
	    /*Sets the agents desired agent state.*/
	    void set_desiredAgentState(const std::vector<typeRNum>& x_des);
	    /*Sets the agents desired agent state.*/
        void set_desiredAgentState(const std::vector<typeRNum>& x_des, const std::vector<typeRNum>& u_des);

        /*Sets the initial state.*/
        void set_initialState(const std::vector<typeRNum>& x_init);
        void set_initialState(const std::vector<typeRNum>& x_init, const std::vector<typeRNum>& u_init);

        /*Return the agents agent state.*/
	    const AgentState& get_agentState() const;
        /*Returns the agents previous agent state*/
        const AgentState&  get_previous_agentState() const;
	    /*Return the agents coupling state.*/
	    const CouplingState& get_couplingState() const;
	    /*Return the agents previous coupling state.*/
	    const CouplingState& get_previous_couplingState() const;
	    /*Return the agents multiplier state.*/
	    const MultiplierState& get_multiplierState() const;
	    /*Return the agents previous multiplier state.*/
	    const MultiplierState& get_previous_multiplierState() const;
	    /*Return the agents penalty state.*/
	    const PenaltyState& get_penaltyState() const;

        /*Set the agents agent state.*/
	    void set_agentState(const AgentState& state);
	    /*Set the agents coupling state.*/
	    void set_couplingState(const CouplingState& state);
	    /*Set the agents multiplier state.*/
	    void set_multiplierState(const MultiplierState& state);
	    /*Set the agents penalty state.*/
        void set_penaltyState( const PenaltyState& penalty );

	    /*Returns all coupling models.*/
        const std::shared_ptr< std::map<int, CouplingModelPtr> > get_couplingModels() const;

        /*Sets the simulated state.*/
        void set_updatedState(const std::vector<typeRNum>& new_state, const typeRNum dt, const typeRNum t0, const typeRNum cost);

        /*Returns the agents neighbors.*/
        const std::vector<NeighborPtr>& get_neighbors() const;
        /*Returns the agents receiving neighbors.*/
        const std::vector<NeighborPtr>& get_receivingNeighbors() const;
        /*Returns the agents sending neighbors.*/
        const std::vector<NeighborPtr>& get_sendingNeighbors() const;

        /*Returns the predicted cost.*/
        const typeRNum get_predicted_cost() const;

        /*************************************************************************
         communication functions
         *************************************************************************/

        /*This function is called if a message is received to register a coupling.*/
        void fromCommunication_registered_coupling(const CouplingInfo& coupling_info, const AgentInfo& neighbor_info);
        /*This function is called if a message is received to de-register a coupling.*/
        void fromCommunication_deregistered_coupling(const CouplingInfo& coupling_info);
        /*This function is called if a message is received including the number of a neighbors neighbors..*/
        void fromCommunication_received_numberOfNeighbors( const unsigned int number, const int from );
        /*This function is called if a message is received including the neighbors local copies.*/
        void fromCommunication_received_localCopies(const AgentState& state, int from);
        /*This function is called if a message is received including the neighbors agent state.*/
        void fromCommunication_received_agentState(const AgentState& state, const ConstraintState& constr_state, int from);
        /*This function is called if a message is received including the neighbors desired agent state.*/
        void fromCommunication_received_desiredAgentState(const AgentState& state, int from);
        /*This function is called if a message is received including the neighbors coupling state.*/
	    void fromCommunication_received_couplingState(const CouplingState& state, int from);
	    /*This function is called if a message is received including the neighbors coupling states.*/
        void fromCommunication_received_couplingState(const CouplingState& state, const CouplingState& state2, int from);
        /*This function is called if a message is received including the neighbors multiplier state.*/
        void fromCommunication_received_multiplierState(const MultiplierState& state, PenaltyState penalty, int from);
        /*This function is called if a message is received including the neighbors sensi states*/
        void fromCommunication_received_sensiState(const SensiState& state, int from);
        /*This function is called if a message is received including the optimization info.*/
        void fromCommunication_configured_optimization(const OptimizationInfo& info);
        /*This function is called if a message is received to trigger an ADMM step.*/
        void fromCommunication_trigger_step(const AlgStep& step);
        /*this function is called if a message is received to stop the ADMM algorithm*/
        void fromCommunication_received_flagToStopAdmm(const bool flag);
        /*this function is called if a message is received containing the ADMM iterations of the neighbors*/
        void fromCommunication_received_flagStoppedAdmm(const bool flag, int from);

        /*************************************************************************
        end of communication functions
        *************************************************************************/

        /*Set initial states for neighbors.*/
        void set_neighbors_initial_states();
        /*Return the current solution.*/
        const SolutionPtr& get_solution() const;
        /*Sets a new solution.*/
        void set_solution(const SolutionPtr& solution);
        /*Reset the current solution.*/
        void reset_solution();
        /*writes predicted debug cost to solution*/
        void print_debugCost();

        /************************************************************************
        asynchronous functions/
        *************************************************************************/

        /*increase all the corresponding delays of the neighbors*/
        void increase_all_delays(const AlgStep& step) const;
        /*initialize the delays of the neighbors*/
        void initialize_allNeighborDelays() const ;
        /*returns the delay of the corresponding neighbor*/
        const int get_delay_sending_neighbors(const AlgStep& step) const;
        /*returns the delay of the corresponding neighbor*/
        const int get_delay_receiving_neighbors(const AlgStep& step) const;
        /*resets the flag stopAdmmflag*/
        void reset_stopAdmmflag_of_neighbors();

        /************************************************************************
        end of asynchronous functions/
       *************************************************************************/

       /*************************************************************************
        sensi functions
        *************************************************************************/


        /*************************************************************************
        end of sensi functions
        *************************************************************************/
       

    private:
        //*********************************************
        // log
        //*********************************************
        const LoggingPtr log_;

        //*********************************************
        // agent parameters
        //*********************************************
        CommunicationInterfacePtr communication_interface_;
        const ModelFactoryPtr model_factory_;
        const AgentModelPtr model_;
        const AgentInfo info_;
        std::vector<typeRNum> x_init_;
	    std::vector<typeRNum> u_init_;
	    std::vector<typeRNum> x_des_;
	    std::vector<typeRNum> u_des_;

        //*********************************************
        // vector for neighbors
        //*********************************************

        std::vector<NeighborPtr> neighbors_;
        std::vector<NeighborPtr> receiving_neighbors_;
        std::vector<NeighborPtr> sending_neighbors_;
        std::vector<NeighborPtr> stopped_neighbors_;


        //*********************************************
        // agent states
        //*********************************************

        AgentState desired_agentState_;
        AgentState agentState_;
        AgentState previous_agentState_;

        CouplingState couplingState_;
        CouplingState previous_couplingState_;

        MultiplierState multiplierState_;
        MultiplierState previous_multiplierState_;

        PenaltyState penaltyState_;
        typeRNum initial_penalty_;

        //*********************************************
        // optimization parameters
        //*********************************************

        SolverLocalPtr local_solver_;
        OptimizationInfo optimizationInfo_;
        SolutionPtr solution_;
        StepSelectorPtr step_selector_;

        //*********************************************
        // neighbor approximation
        //*********************************************

        bool is_approximating_ = false;
        bool is_approximatingCost_ = false;
        bool is_approximatingConstraints_ = false;
        bool is_approximatingDynamics_ = false;
    };

}
