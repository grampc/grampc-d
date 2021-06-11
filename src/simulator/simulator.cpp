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

#include "dmpc/simulator/simulator.hpp"

#include "dmpc/agent/agent.hpp"

#include "dmpc/model/agent_model.hpp"
#include "dmpc/model/coupling_model.hpp"

#include "dmpc/optim/optim_util.hpp"

#include "dmpc/optim/solver_central.hpp"

#include "dmpc/comm/communication_interface.hpp"

#include "dmpc/util/logging.hpp"

#include <algorithm>

namespace dmpc
{

    Simulator::Simulator(const CommunicationInterfacePtr& communication_interface, const LoggingPtr& log)
        : communication_interface_(communication_interface),
        log_(log)
    {}

    void Simulator::distributed_simulation(const std::string& Integrator, typeRNum dt)
    {
        // integration method
        Integrator_ = Integrator;

        // time step
        dt_ = dt;

        // get map with all agent infos
        const auto agents = communication_interface_->get_agentInfo_from_coordinator();
        if (agents == nullptr)
            return;

        // collect AgentStates
        agentStates_.clear();
        for(const auto& [id, info] : *agents )
        {
            // get agentState
			const auto agentStatePtr = communication_interface_->get_agentState_from_agent(id);
			const auto desired_agentStatePtr = communication_interface_->get_desiredAgentState_from_agent(id);

            // stop simulation if pointer is nullptr
            if (agentStatePtr == nullptr || desired_agentStatePtr == nullptr)
            {
                log_->print(DebugType::Warning) << "[Simulator::distributed_simulation] Agent " << id
                    << " did not answer. Simulation is interrupted." << std::endl;
                return;
            }

			agentStates_.insert(std::make_pair(id, agentStatePtr));
			desired_agentStates_.insert(std::make_pair(id, desired_agentStatePtr));
        }

        //start simulation
        simulate();
    }

    void Simulator::centralized_simulation(const SolverCentral *solver, const std::string& Integrator, typeRNum dt)
    {
        // integration method
        Integrator_ = Integrator;

        // time step
        dt_ = dt;

        // get all agents
        const auto& agents = solver->get_agents();

        // collect AgentStates
        agentStates_.clear();
        for (const auto& agent : agents)
        {
            agentStates_.insert(std::make_pair(agent->get_id(), std::make_shared<AgentState>(agent->get_agentState())));
            desired_agentStates_.insert(std::make_pair(agent->get_id(), std::make_shared<AgentState>(agent->get_desiredAgentState())));
        }

	    //start simulation
        simulate();
    }

    void Simulator::simulate()
    {
        // evaluate t0
        // The initial time step may differ, e.g. due to plug-and-play scenarios.
        // Hence, they are synchronized here.
        for (const auto& [id, state] : agentStates_)
            t0_ = std::max(t0_, state->t0_);

        if( Integrator_ == "euler" )
        {
            for(const auto& [id, state] : agentStates_)
            {
                // get agent model
                const auto agent_model = communication_interface_->get_agentModel( id );

                if (agent_model == nullptr)
                {
                    log_->print(DebugType::Warning) << "[Simulator::simulate] Agent " << id
                        << " did not send its agent model. Simulation is interrupted." << std::endl;
                    return;
                }

                // get all coupling model from agent
                const auto couplingModels = communication_interface_->get_couplingModels_from_agent( id );

                if (couplingModels == nullptr)
                {
                    log_->print(DebugType::Warning) << "[Simulator::simulate] Agent " << id
                        << " did not send its coupling models. Simulation is interrupted." << std::endl;
                    return;
                }

                // get Nxi
                const auto Nxi = agent_model->get_Nxi();

                // initialize vector for new states
                std::vector<typeRNum> x_next( Nxi, 0.0 );

                // evaluate f_i( x_i, u_i )
                agent_model->ffct( &x_next[0], state->t_[0], &state->x_[0], &state->u_[0] );

                // consider each sending neighbor
                for( const auto& [neighbor_id, neighbor_model] : *couplingModels )
                {
                    // get neighbor state
                    const auto& neighbor_state = agentStates_.find(neighbor_id)->second;

                    // evaluate f_{ij}( xi, ui, x_j, u_j )
                    neighbor_model->ffct( &x_next[0], state->t_[0], &state->x_[0], &state->u_[0], &neighbor_state->x_[0], &neighbor_state->u_[0] );
                }

                // x_{k+1} = x_k + dt*x'_k
                for( unsigned int i = 0; i < Nxi; ++i )
                    x_next[i] = state->x_[i] + dt_ * x_next[i];

                const auto cost = evaluate_cost(id, agent_model, couplingModels);

                // send simulated state to agent
                communication_interface_->set_simulatedState_for_agent( id, x_next, dt_, t0_, cost);
            }
        }
        else if( Integrator_ == "heun" )
        {
            std::map<unsigned int, std::shared_ptr< std::vector< typeRNum > > > out0;
            std::map<unsigned int, std::shared_ptr< std::vector< typeRNum > > > out1;

            for(const auto& [id, state] : agentStates_)
            {
                // get agent model
                const auto agent_model = communication_interface_->get_agentModel( id );

			    if (agent_model == nullptr)
			    {
				    log_->print(DebugType::Warning) << "[Simulator::simulate] Agent " << id
					    << " did not send its agent model. Simulation is interrupted." << std::endl;
				    return;
			    }

                // get all coupling model from agent
                auto couplingModels = communication_interface_->get_couplingModels_from_agent( id );

			    if (couplingModels == nullptr)
			    {
				    log_->print(DebugType::Warning) << "[Simulator::simulate] Agent " << id
					    << " did not send its coupling models. Simulation is interrupted." << std::endl;
				    return;
			    }

			    // initialize vector
                std::shared_ptr< std::vector< typeRNum > > ptr_to_data(new std::vector<typeRNum>(agent_model->get_Nxi(), 0.0));

                out0.insert(std::make_pair(id, ptr_to_data));

                // evaluate f_i( x_i, u_i )
                agent_model->ffct( &((*ptr_to_data)[0]), state->t_[0], &state->x_[0], &state->u_[0] );

                // consider each sending neighbor
                for(const auto& [neighbor_id, neighbor_model] : *couplingModels)
                {
                    // get neighbor state
                    const auto& neighbor_state = agentStates_.find(neighbor_id)->second;

                    // evaluate f_{ij}( xi, ui, x_j, u_j )
                    neighbor_model->ffct( &((*ptr_to_data)[0]), state->t_[0], &state->x_[0], &state->u_[0], &neighbor_state->x_[0], &neighbor_state->u_[0] );
                }

                // compute euler
                // x^*_{k+1} = x_k + dt*x'_k
                for( unsigned int i = 0; i < agent_model->get_Nxi(); ++i )
                    (*ptr_to_data)[i] = state->x_[i] + dt_ * (*ptr_to_data)[i];
            }

            for(const auto& [id, state] : agentStates_)
            {
                // get agent model
                const auto agent_model = communication_interface_->get_agentModel( id );

			    if (agent_model == nullptr)
			    {
				    log_->print(DebugType::Warning) << "[Simulator::simulate] Agent " << id
					    << " did not send its agent model. Simulation is interrupted." << std::endl;
				    return;
			    }

                // get all coupling model from agent
                auto couplingModels = communication_interface_->get_couplingModels_from_agent( id );

			    if (couplingModels == nullptr)
			    {
				    log_->print(DebugType::Warning) << "[Simulator::simulate] Agent " << id
					    << " did not send its coupling models. Simulation is interrupted." << std::endl;
				    return;
			    }

			    // initialize vector
			    auto ptr_to_data0 = out0.find(id)->second;
                std::shared_ptr<std::vector<typeRNum>> ptr_to_data1(new std::vector<typeRNum>(agent_model->get_Nxi(), 0.0));
			    out1.insert(std::make_pair(id, ptr_to_data1));

                // evaluate f_i( x_i, u_i )
                agent_model->ffct( &((*ptr_to_data1)[0]), state->t_[0], &((*ptr_to_data0)[0]), &state->u_[0] );

                // consider each sending neighbor
                for(const auto& [neighbor_id, neighbor_model] : *couplingModels)
                {
                    // get neighbor state
                    const auto& neighbor_state = agentStates_.find(neighbor_id)->second;
                    const auto& ptr_to_data_of_neighbor = out0.find(neighbor_id)->second;

                    // evaluate f_{ij}( xi, ui, x_j, u_j )
                    neighbor_model->ffct( &((*ptr_to_data1)[0]), state->t_[0], &((*ptr_to_data0)[0]), &state->u_[0], &((*ptr_to_data_of_neighbor)[0]), &neighbor_state->u_[0] );
                }

                // compute heun
                std::vector<typeRNum> x_next_heun(agent_model->get_Nxi(), 0.0);
                for(unsigned int i = 0; i < agent_model->get_Nxi(); ++i)
                    x_next_heun[i] = 0.5 * state->x_[i] + 0.5 * ((*ptr_to_data0)[i]) + 0.5 * dt_ * ((*ptr_to_data1)[i]);


				const auto cost = evaluate_cost(id, agent_model, couplingModels);

                communication_interface_->set_simulatedState_for_agent( id, x_next_heun, dt_, t0_, cost);
            }
        }
        else
            log_->print(DebugType::Error) << "[Simulator::simulate]: Unknown Integrator." << std::endl;
    }

    void Simulator::set_t0(typeRNum t0)
    {
        t0_ = t0;
    }

    const typeRNum Simulator::evaluate_cost
	(
        const unsigned int agent_id,
        const AgentModelPtr& agent_model,
        const std::shared_ptr< std::map<int, CouplingModelPtr> >& coupling_models
    ) const
    {
        typeRNum cost = 0;

		const auto& agent_states = agentStates_.find(agent_id)->second;
		const auto& desired_agent_state = desired_agentStates_.find(agent_id)->second;
        const auto& t = agent_states->t_;
        const auto Nhor = t.size();

        const auto Nxi = agent_model->get_Nxi();
        const auto Nui = agent_model->get_Nui();

        for (unsigned int i = 0; i < Nhor; ++i)
        {
            agent_model->lfct(&cost, t0_ + t[i], &agent_states->x_[i*Nxi], &agent_states->u_[i*Nui], &desired_agent_state->x_[i*Nxi]);

            for (const auto& [neighbor_id, coupling_model] : *coupling_models)
            {
				const auto Nxj = coupling_model->get_Nxj();
                const auto Nuj = coupling_model->get_Nuj();
                const auto& neighbor_states = agentStates_.find(neighbor_id)->second;

                coupling_model->lfct
                (
                    &cost,
                    t0_ + t[i],
                    &agent_states->x_[i * Nxi],
                    &agent_states->u_[i * Nui],
                    &neighbor_states->x_[i * Nxj],
                    &neighbor_states->u_[i * Nuj]
                );
            }
        }

        cost *= t[1] - t[0];

        agent_model->Vfct(&cost, t0_ + t.back(), &agent_states->x_[(Nhor - 1) * Nxi], &desired_agent_state->x_[(Nhor - 1) * Nxi]);

		for (const auto& [neighbor_id, coupling_model] : *coupling_models)
		{
			const auto Nxj = coupling_model->get_Nxj();
			const auto Nuj = coupling_model->get_Nuj();
			const auto& neighbor_states = agentStates_.find(neighbor_id)->second;

			coupling_model->Vfct
			(
				&cost,
				t0_ + t.back(),
				&agent_states->x_[(Nhor - 1) * Nxi],
				&neighbor_states->x_[(Nhor - 1) * Nxj]
			);
		}

        return cost;
    }
}
