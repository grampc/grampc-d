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

#include "grampcd/interface/dmpc_interface.hpp"

#include "general_model_factory.hpp"

#include "grampcd/comm/communication_interface_local.hpp"
#include "grampcd/comm/communication_interface_central.hpp"

#include "grampcd/agent/agent.hpp"

#include "grampcd/coord/coordinator.hpp"

#include "grampcd/simulator/simulator.hpp"

#include "grampcd/optim/solution.hpp"

#include "grampcd/util/logging.hpp"
#include "grampcd/util/data_conversion.hpp"
#include "grampcd/util/auto_tune.hpp"

#include "grampcd/info/tuning_info.hpp"

#include "chrono"
#include <fstream>
#include <algorithm>

namespace grampcd
{
	DmpcInterface::DmpcInterface() :
		log_(std::make_shared<Logging>())
	{
	}

	void DmpcInterface::run_MPC
	(
		const std::vector<AgentPtr>& agents,
		const SimulatorPtr& simulator,
		const OptimizationInfo& oi,
		double Tsim,
		double t_0
	)
	{
		SolverCentral solver(agents, oi);
		const unsigned int maxSimIter = static_cast<unsigned int>(Tsim / oi.COMMON_dt_);
		std::chrono::milliseconds CPUtime(0);

		// initialize agents
		for (const auto& agent : agents)
			agent->initialize(oi);

		// main loop for centralized solution
		log_->print(DebugType::Base) << "MPC running ..." << std::endl;

		for (unsigned int iMPC = 0; iMPC <= maxSimIter; ++iMPC)
		{
			// optimize
			const auto tstart = std::chrono::steady_clock::now();
			solver.solve();
			const auto tend = std::chrono::steady_clock::now();
			CPUtime += std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);

			// update state and time
			simulator->centralized_simulation(&solver, oi.COMMON_Integrator_, oi.COMMON_dt_);

			// update progress bar
			print_progressbar(iMPC, maxSimIter);
		}

		log_->print(DebugType::Progressbar) << std::endl;

		log_->print(DebugType::Base) << "MPC finished. Average computation time: "
			<< static_cast<typeRNum>(static_cast<typeRNum>(CPUtime.count()) / static_cast<typeRNum>(maxSimIter + 1)) << " ms." << std::endl << std::endl;
	}

	void DmpcInterface::run_MPC(const std::vector<AgentPtr>& agents, const SimulatorPtr& simulator, const OptimizationInfo& oi)
	{
		SolverCentral solver(agents, oi);

		// initialize agents
		for (const auto& agent : agents)
			agent->initialize(oi);

		unsigned int CPUtime_iteration(0);
		const unsigned int dt_in_ms = static_cast<unsigned int>(oi.COMMON_dt_ * 1000);
		log_->print(DebugType::Base) << "MPC running in endless mode." << std::endl;

		// main loop for centralized solution
		while (true)
		{
			// optimize
			const auto tstart = std::chrono::steady_clock::now();
			solver.solve();
			const auto tend = std::chrono::steady_clock::now();
			CPUtime_iteration = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count());

			// simulate realtime
			if (realtime_ && static_cast<int>(dt_in_ms - CPUtime_iteration) > 0)
				std::this_thread::sleep_for(std::chrono::milliseconds(dt_in_ms - CPUtime_iteration));

			// update state and time
			typeRNum simulate_timestep;
			if (realtime_)
				simulate_timestep = static_cast<typeRNum>(std::max(dt_in_ms, CPUtime_iteration)) / 1000.0;
        else
				simulate_timestep = oi.COMMON_dt_;

			simulator->centralized_simulation(&solver, oi.COMMON_Integrator_, simulate_timestep);
		}
	}

	void DmpcInterface::run_DMPC
	(
		const SimulatorPtr& simulator,
		const OptimizationInfo& oi, 
		const typeRNum Tsim,
		const typeRNum t_0
	)
	{
		const unsigned int maxSimIter = static_cast<unsigned int>(Tsim / oi.COMMON_dt_);
		std::chrono::milliseconds CPUtime(0);

		// main loop for centralized solution
		coordinator_->initialize_ADMM(oi);
		log_->print(DebugType::Base) << "DMPC running ..." << std::endl;

		simulator_->set_t0(t_0);

		std::chrono::milliseconds::rep CPUtime_max = 0;

		for (unsigned int iMPC = 0; iMPC <= maxSimIter; ++iMPC)
		{
			// optimize
			const auto tstart = std::chrono::steady_clock::now();
			coordinator_->solve_ADMM(oi.ADMM_maxIterations_, oi.ADMM_innerIterations_);
			const auto tend = std::chrono::steady_clock::now();
			CPUtime += std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);
			CPUtime_max = std::max(CPUtime_max, std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count());

			// update state and time
			simulator->distributed_simulation(oi.COMMON_Integrator_, oi.COMMON_dt_);

			// update progress bar
			print_progressbar(iMPC, maxSimIter);
		}

		log_->print(DebugType::Progressbar) << std::endl;

		const auto number_of_agents = communication_interface_->get_numberOfAgents();
		log_->print(DebugType::Base) << "DMPC finished." << std::endl
			<< "Maximum computation time : "
			<< CPUtime_max << " ms in total or " << CPUtime_max / number_of_agents << " ms per agent." << std::endl
			<< "Average computation time : "
			<< CPUtime.count() / static_cast<typeRNum>(maxSimIter + 1) << " ms in total or "
			<< CPUtime.count() / static_cast<typeRNum>(maxSimIter + 1) / number_of_agents << " ms per agent." << std::endl << std::endl;
	}

	void DmpcInterface::run_DMPC(const SimulatorPtr& simulator, const OptimizationInfo& oi)
	{
		// main loop for centralized solution
		coordinator_->initialize_ADMM(oi);

		unsigned int CPUtime_iteration(0);
		const unsigned int dt_in_ms = static_cast<unsigned int>(oi.COMMON_dt_ * 1000);

		log_->print(DebugType::Base) << "DMPC running in endless mode..." << std::endl;

		while (true)
		{
			// optimize
			const auto tstart = std::chrono::steady_clock::now();
			coordinator_->solve_ADMM(oi.ADMM_maxIterations_, oi.ADMM_innerIterations_);
			const auto tend = std::chrono::steady_clock::now();
			CPUtime_iteration = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count());

			// simulate real time
			if (realtime_ && static_cast<int>(dt_in_ms - CPUtime_iteration) > 0)
				std::this_thread::sleep_for(std::chrono::milliseconds(dt_in_ms - CPUtime_iteration));

			// update state and time
			typeRNum simulate_timestep;
			if (realtime_)
				simulate_timestep = static_cast<typeRNum>(std::max(dt_in_ms, CPUtime_iteration)) / 1000.0;
			else
				simulate_timestep = oi.COMMON_dt_;

			// update state and time
			simulator->distributed_simulation(oi.COMMON_Integrator_, simulate_timestep);
		}
	}

	void DmpcInterface::print_progressbar(const int current_iter, const int max_iter)
	{
		if (max_iter == 0)
			return;

		const int progress = static_cast<int>(std::round(static_cast<double>(current_iter) / max_iter * 100));

		log_->print(DebugType::Progressbar) << "\r" << progress << "% completed: " << std::string(progress / 2, '|');
		log_->print(DebugType::Progressbar).flush();

		if (current_iter == max_iter)
			log_->print(DebugType::Progressbar) << std::endl;
	}

	void DmpcInterface::initialize_central_communicationInterface(const int number_of_threads)
	{
		// create communication interface
		CommunicationInterfacePtr communication_interface(new CommunicationInterfaceCentral(log_, number_of_threads));
		communication_interface_ = communication_interface;

		// create coordinator
		coordinator_ = CoordinatorPtr(new Coordinator(communication_interface_, true, log_));

		// set coordinator
		std::static_pointer_cast<CommunicationInterfaceCentral>(communication_interface_)->set_coordinator(coordinator_);

		// create simulator
		simulator_ = SimulatorPtr(new Simulator(communication_interface_, log_));

		// set simulator
		std::static_pointer_cast<CommunicationInterfaceCentral>(communication_interface_)->set_simulator(simulator_);

		// create factory
		factory_ = ModelFactoryPtr(new GeneralModelFactory(log_));
	}

	void DmpcInterface::initialize_local_communicationInterface_as_agent(const CommunicationInfo& adress_coordinator)
	{
		// create communication interface
		communication_interface_ = CommunicationInterfacePtr(new CommunicationInterfaceLocal(log_, adress_coordinator));

		// create factory
		factory_ = ModelFactoryPtr(new GeneralModelFactory(log_));
	}

	void DmpcInterface::initialize_local_communicationInterface_as_coordinator(const unsigned short port)
	{
		// create communication interface
		communication_interface_ = CommunicationInterfacePtr(new CommunicationInterfaceLocal(log_, port));

		// create coordinator
		coordinator_ = CoordinatorPtr(new Coordinator(communication_interface_, true, log_));

		// set coordinator
		std::static_pointer_cast<CommunicationInterfaceLocal>(communication_interface_)->set_coordinator(coordinator_);

		// create simulator
		simulator_ = SimulatorPtr(new Simulator(communication_interface_, log_));

		// set simulator
		std::static_pointer_cast<CommunicationInterfaceLocal>(communication_interface_)->set_simulator(simulator_);

		// create factory
		factory_ = ModelFactoryPtr(new GeneralModelFactory(log_));
	}

	void DmpcInterface::register_agent
	(
		const AgentInfo& info,
		const std::vector<typeRNum>& x_init,
		const std::vector<typeRNum>& u_init,
		const std::vector<typeRNum>& x_des,
		const std::vector<typeRNum>& u_des
	)
	{
		if (info.id_ < 0)
		{
			log_->print(DebugType::Error) << "[DmpcInterface::register_agent] "
				<< "Registration of agent rejected since negative agent ids are not allowed." << std::endl;

			return;
		}

		if (x_init.size() != x_des.size())
		{
			log_->print(DebugType::Error) << "[DmpcInterface::register_agent] "
				<< "Registration of agent rejected since negative dimensions of x_init and x_des do not fit." << std::endl;

			return;
		}

		if (u_init.size() != u_des.size())
		{
			log_->print(DebugType::Error) << "[DmpcInterface::register_agent] "
				<< "Registration of agent rejected since negative dimensions of u_init and u_des do not fit." << std::endl;

			return;
		}

		// create new agent
		AgentPtr agent = std::make_shared<Agent>(Agent(communication_interface_, factory_, info, log_));

		// register agent
		communication_interface_->register_agent(agent);

		// initialize
		agent->set_initialState(x_init, u_init);
		agent->initialize(optimizationInfo_);
		agent->set_desiredAgentState(x_des, u_des);

		// safe in list
		agents_.push_back(agent);
	}

	void DmpcInterface::deregister_agent(const AgentInfo& info)
	{
		if (communication_interface_->deregister_agent(info))
			DataConversion::erase_element_from_vector(agents_, info);
	}

	void DmpcInterface::set_desiredAgentState(const int agent_id, const std::vector<typeRNum>& x_des, const std::vector<typeRNum>& u_des)
	{
		for (const auto& agent : agents_)
		{
			if (agent->get_id() == agent_id)
			{
				agent->set_desiredAgentState(x_des, u_des);
				break;
			}
		}
	}

	void DmpcInterface::register_coupling(const CouplingInfo& info)
	{
		communication_interface_->register_coupling(info);
	}

	void DmpcInterface::run_MPC(ctypeRNum t_0, ctypeRNum Tsim)
	{
		run_MPC(agents_, simulator_, optimizationInfo_, Tsim, t_0);
	}

	void DmpcInterface::run_MPC()
	{
		run_MPC(agents_, simulator_, optimizationInfo_);
	}

	void DmpcInterface::run_DMPC(ctypeRNum t_0, ctypeRNum Tsim)
	{
		run_DMPC(simulator_, optimizationInfo_, Tsim, t_0);
	}

	void DmpcInterface::run_DMPC()
	{
		run_DMPC(simulator_, optimizationInfo_);
	}

	SolutionPtr DmpcInterface::get_solution(const unsigned int agent_id) const
	{
		return communication_interface_->get_solution(agent_id);
	}

	std::vector< SolutionPtr > DmpcInterface::get_solution(const std::string& agents) const
	{
		return communication_interface_->get_solution(agents);
	}

	void DmpcInterface::reset_solution(const unsigned int agent_id)
	{
		communication_interface_->reset_solution(agent_id);
	}

	void DmpcInterface::reset_solution(const std::string& agents)
	{
		communication_interface_->reset_solution(agents);
	}

	OptimizationInfoPtr DmpcInterface::get_optimizationInfo() const
	{
		return std::make_shared<OptimizationInfo>(optimizationInfo_);
	}

	void DmpcInterface::set_optimizationInfo(const OptimizationInfo& optimization_info)
	{
		optimizationInfo_ = optimization_info;
	}

	void DmpcInterface::wait_for_connections(int agents, int couplings)
	{
		communication_interface_->waitFor_connection(agents, couplings);
	}

	void DmpcInterface::set_passive()
	{
		communication_interface_->set_passive();
	}

	void DmpcInterface::send_flag_to_agents(const int agent_id) const
	{
		communication_interface_->send_flag_to_agents(agent_id);
	}

	void DmpcInterface::send_flag_to_agents(const std::vector<int>& agent_ids) const
	{
		communication_interface_->send_flag_to_agents(agent_ids);
	}

	void DmpcInterface::send_flag_to_agents(const std::string& agents) const
	{
		communication_interface_->send_flag_to_agents(agents);
	}

	void DmpcInterface::waitFor_flag_from_coordinator()
	{
		communication_interface_->waitFor_flag_from_coordinator();
	}

	void DmpcInterface::deregister_coupling(const CouplingInfo& info)
	{
		communication_interface_->deregister_coupling(info);
	}

	void DmpcInterface::wait_blocking_s(const unsigned int s)
	{
		std::this_thread::sleep_for(std::chrono::seconds(s));
	}

	void DmpcInterface::simulate_realtime(const bool realtime)
	{
		realtime_ = realtime;
	}

	void DmpcInterface::cap_stored_data(const unsigned int data_points)
	{
		communication_interface_->cap_stored_data(data_points);
	}

	void DmpcInterface::set_print_base(const bool print)
	{
		log_->print_base_ = print;
	}

	void DmpcInterface::set_print_error(const bool print)
	{
		log_->print_error_ = print;
	}

	void DmpcInterface::set_print_message(const bool print)
	{
		log_->print_message_ = print;
	}

	void DmpcInterface::set_print_warning(const bool print)
	{
		log_->print_warning_ = print;
	}

	void DmpcInterface::set_print_progressbar(const bool print)
	{
		log_->print_progressbar_ = print;
	}

	void DmpcInterface::print_solution_to_file(const unsigned int agent_id, const std::string& prefix) const
	{
		std::string filename = prefix + std::to_string(agent_id);

		std::ofstream outputFile(filename);
		if (outputFile)
		{
			// get solution
			const auto solutions = communication_interface_->get_solution(agent_id);

			// print solution
			outputFile << *solutions;
		}
		else
			log_->print(DebugType::Error) << "[DmpcInterface::print_solution_to_file] "
			<< "Failed to open file." << std::endl;
	}

	void DmpcInterface::print_solution_to_file(const std::string& agents, const std::string& prefix) const
	{
		// get solutions
		const auto solutions = communication_interface_->get_solution(agents);

		for (unsigned int k = 0; k < solutions.size(); ++k)
		{
			std::string filename = prefix + std::to_string(solutions[k]->agentState_.i_) + ".txt";

			// open file
			std::ofstream outputFile(filename);

			// print solutions
			if (outputFile)
				outputFile << *solutions[k];
			else
				log_->print(DebugType::Error) << "[DmpcInterface::print_solution_to_file] "
				<< "Failed to open file." << std::endl;
		}
	}

	AgentInfo DmpcInterface::agent_info() const
	{
		return AgentInfo();
	}

	CouplingInfo DmpcInterface::coupling_info() const
	{
		return CouplingInfo();
	}

	OptimizationInfo DmpcInterface::optimization_info() const
	{
		return OptimizationInfo();
	}

	CommunicationInfo DmpcInterface::communication_info() const
	{
		return CommunicationInfo();
	}

	TuningInfo DmpcInterface::tuning_info() const
	{
		return TuningInfo();
	}

	void DmpcInterface::set_initialState(const unsigned int agent_id, const std::vector<typeRNum>& x_init)
	{
		// search for agent
		for (const auto& agent : agents_)
		{
			if (agent->get_id() == agent_id)
			{
				agent->set_initialState(x_init);
				return;
			}
		}

		// agent was not found
		log_->print(DebugType::Error) << "[DmpcInterface::set_initialState]: "
			<< "Agent " << agent_id << " is not registered. " << std::endl;
	}

	const OptimizationInfo DmpcInterface::auto_tune_parameters
	(
		const TuningInfo& tuning_info,
		ctypeRNum convergence_tolerance,
		const std::string& type,
		const int size_of_population,
		const int number_of_generations
	)
	{
		Autotune auto_tune(log_, optimizationInfo_, agents_, coordinator_ );

		return auto_tune.auto_tune_parameters
		(
			tuning_info, 
			convergence_tolerance, 
			type, size_of_population, 
			number_of_generations
		);
	}
}
