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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dmpc/interface/python_interface.hpp>
#include "dmpc/interface/dmpc_interface.hpp"

#include "dmpc/state/agent_state.hpp"

#include "dmpc/optim/solution.hpp"

namespace dmpc
{
	namespace py = pybind11;

	PythonInterface::PythonInterface()
		: dmpc_interface_(std::make_shared<DmpcInterface>())
	{}

	void PythonInterface::initialize_central_communicationInterface(int number_of_threads)
	{
		dmpc_interface_->initialize_central_communicationInterface(number_of_threads);
	}

	void PythonInterface::initialize_local_communicationInterface_as_agent(CommunicationInfo adress_coordinator)
	{
		dmpc_interface_->initialize_local_communicationInterface_as_agent(adress_coordinator);
	}

	void PythonInterface::initialize_local_communicationInterface_as_coordinator(unsigned short port)
	{
		dmpc_interface_->initialize_local_communicationInterface_as_coordinator(port);
	}

	void PythonInterface::register_agent(AgentInfo info, std::vector<typeRNum> x_init, std::vector<typeRNum> u_init)
	{
		dmpc_interface_->register_agent(info, x_init, u_init);
	}

	void PythonInterface::set_desiredAgentState(int agent_id, std::vector<typeRNum> x_des, std::vector<typeRNum> u_des)
	{
		dmpc_interface_->set_desiredAgentState(agent_id, x_des, u_des);
	}

	void PythonInterface::register_coupling(CouplingInfo info)
	{
		dmpc_interface_->register_coupling(info);
	}

	void PythonInterface::run_MPC(typeRNum Tsim, typeRNum t0)
	{
		dmpc_interface_->run_MPC(Tsim, t0);
	}

	void PythonInterface::run_MPC()
	{
		dmpc_interface_->run_MPC();
	}

	void PythonInterface::run_DMPC(typeRNum Tsim, typeRNum t0)
	{
		dmpc_interface_->run_DMPC(Tsim, t0);
	}

	void PythonInterface::run_DMPC()
	{
		dmpc_interface_->run_DMPC();
	}

	SolutionPtr PythonInterface::get_solution(unsigned int agent_id) const
	{
		return dmpc_interface_->get_solution(agent_id);
	}

	std::vector< SolutionPtr > PythonInterface::get_solution(std::string agents) const
	{
		return dmpc_interface_->get_solution(agents);
	}

	void PythonInterface::reset_solution(unsigned int agent_id)
	{
		dmpc_interface_->reset_solution(agent_id);
	}

	void PythonInterface::reset_solution(std::string agents)
	{
		dmpc_interface_->reset_solution(agents);
	}

	OptimizationInfoPtr PythonInterface::get_optimizationInfo() const
	{
		return dmpc_interface_->get_optimizationInfo();
	}

	void PythonInterface::set_optimizationInfo(OptimizationInfo optimization_info)
	{
		dmpc_interface_->set_optimizationInfo(optimization_info);
	}

	void PythonInterface::wait_for_connections(int agents, int couplings)
	{
		dmpc_interface_->wait_for_connections(agents, couplings);
	}

	void PythonInterface::set_passive()
	{
		dmpc_interface_->set_passive();
	}

	void PythonInterface::send_flag_to_agents(const int agent_id) const
	{
		dmpc_interface_->send_flag_to_agents(agent_id);
	}

	void PythonInterface::send_flag_to_agents(const std::vector<int> agent_ids) const
	{
		dmpc_interface_->send_flag_to_agents(agent_ids);
	}

	void PythonInterface::send_flag_to_agents(const std::string agents) const
	{
		dmpc_interface_->send_flag_to_agents(agents);
	}

	void PythonInterface::waitFor_flag_from_coordinator()
	{
		dmpc_interface_->waitFor_flag_from_coordinator();
	}

	void PythonInterface::deregister_coupling(CouplingInfo info)
	{
		dmpc_interface_->deregister_coupling(info);
	}

	void PythonInterface::wait_blocking_s(unsigned int s)
	{
		dmpc_interface_->wait_blocking_s(s);
	}

	void PythonInterface::simulate_realtime(bool realtime)
	{
		dmpc_interface_->simulate_realtime(realtime);
	}

	void PythonInterface::deregister_agent(AgentInfo info)
	{
		dmpc_interface_->deregister_agent(info);
	}

	void PythonInterface::cap_stored_data(unsigned int data_points)
	{
		dmpc_interface_->cap_stored_data(data_points);
	}

	void PythonInterface::set_print_base(bool print)
	{
		dmpc_interface_->set_print_base(print);
	}

	void PythonInterface::set_print_error(bool print)
	{
		dmpc_interface_->set_print_error(print);
	}

	void PythonInterface::set_print_message(bool print)
	{
		dmpc_interface_->set_print_message(print);
	}

	void PythonInterface::set_print_warning(bool print)
	{
		dmpc_interface_->set_print_warning(print);
	}

	void PythonInterface::print_solution_to_file(const unsigned int agent_id, const std::string prefix) const
	{
		dmpc_interface_->print_solution_to_file(agent_id, prefix);
	}

	void PythonInterface::print_solution_to_file(const std::string agents, const std::string prefix) const
	{
		dmpc_interface_->print_solution_to_file(agents, prefix);
	}

	void PythonInterface::set_initialState(const unsigned int agent_id, const std::vector<typeRNum>& x_init)
	{
		dmpc_interface_->set_initialState(agent_id, x_init);
	}

	PYBIND11_MODULE(grampcd_interface, m) {

		py::class_<AgentState>(m, "AgentState")
			.def(py::init<>())
			.def_readonly("i_", &AgentState::i_)
			.def_readonly("t0_", &AgentState::t0_)
			.def_readonly("t_", &AgentState::t_)
			.def_readonly("u_", &AgentState::u_)
			.def_readonly("v_", &AgentState::v_)
			.def_readonly("x_", &AgentState::x_);

		py::class_<Solution, std::shared_ptr<Solution>>(m, "Solution")
			.def(py::init<>())
			.def_readonly("agentState_", &Solution::agentState_)
			.def_readonly("predicted_agentState_", &Solution::predicted_agentState_)
			.def_readonly("cost_", &Solution::cost_)
			.def_readonly("predicted_cost_", &Solution::predicted_cost_)
			.def_readonly("debug_cost_", &Solution::debug_cost_);

		py::class_<PythonInterface>(m, "interface")
			.def(py::init<>())
			.def("initialize_central_communicationInterface", &PythonInterface::initialize_central_communicationInterface, py::arg("number_of_threads") = 0)
			.def("initialize_local_communicationInterface_as_agent", &PythonInterface::initialize_local_communicationInterface_as_agent)
			.def("initialize_local_communicationInterface_as_coordinator", &PythonInterface::initialize_local_communicationInterface_as_coordinator)
			.def("register_agent", &PythonInterface::register_agent)
			.def("deregister_agent", &PythonInterface::deregister_agent)
			.def("set_desiredAgentState", &PythonInterface::set_desiredAgentState)
			.def("set_initialState", &PythonInterface::set_initialState)
			.def("register_coupling", &PythonInterface::register_coupling)
			.def("deregister_coupling", &PythonInterface::deregister_coupling)
			.def("run_MPC", (void (PythonInterface::*)(void)) & PythonInterface::run_MPC)
			.def("run_MPC", (void (PythonInterface::*)(typeRNum, typeRNum)) & PythonInterface::run_MPC)
			.def("run_DMPC", (void (PythonInterface::*)(void)) & PythonInterface::run_DMPC)
			.def("run_DMPC", (void (PythonInterface::*)(typeRNum, typeRNum)) & PythonInterface::run_DMPC)
			.def("set_optimizationInfo", &PythonInterface::set_optimizationInfo)
			.def("get_optimizationInfo", &PythonInterface::get_optimizationInfo)
			.def("wait_for_connections", &PythonInterface::wait_for_connections)
			.def("send_flag_to_agents", (void (PythonInterface::*)(int) const) & PythonInterface::send_flag_to_agents)
			.def("send_flag_to_agents", (void (PythonInterface::*)(std::vector<int>) const) & PythonInterface::send_flag_to_agents)
			.def("send_flag_to_agents", (void (PythonInterface::*)(std::string) const) & PythonInterface::send_flag_to_agents)
			.def("wait_blocking_s", &PythonInterface::wait_blocking_s)
			.def("waitFor_flag_from_coordinator", &PythonInterface::waitFor_flag_from_coordinator)
			.def("set_passive", &PythonInterface::set_passive)
			.def("get_solution", (SolutionPtr(PythonInterface::*)(unsigned int) const) & PythonInterface::get_solution)
			.def("get_solution", (std::vector< SolutionPtr >(PythonInterface::*)(std::string) const) & PythonInterface::get_solution)
			.def("reset_solution", (void(PythonInterface::*)(unsigned int)) & PythonInterface::reset_solution)
			.def("reset_solution", (void(PythonInterface::*)(std::string)) & PythonInterface::reset_solution)
			.def("print_solution_to_file", (void(PythonInterface::*)(const unsigned int, const std::string) const) & PythonInterface::print_solution_to_file, py::arg("agent_id"), py::arg("prefix") = "Solution_agent")
			.def("print_solution_to_file", (void(PythonInterface::*)(const std::string, const std::string) const) & PythonInterface::print_solution_to_file, py::arg("agents"), py::arg("prefix") = "Solution_agent")
			.def("simulate_realtime", &PythonInterface::simulate_realtime)
			.def("cap_stored_data", &PythonInterface::cap_stored_data)
			.def("set_print_base", &PythonInterface::set_print_base)
			.def("set_print_error", &PythonInterface::set_print_error)
			.def("set_print_message", &PythonInterface::set_print_message)
			.def("set_print_warning", &PythonInterface::set_print_warning);

		py::class_<AgentInfo>(m, "AgentInfo")
			.def(py::init<>())
			.def_readwrite("id_", &AgentInfo::id_)
			.def_readwrite("model_name_", &AgentInfo::model_name_)
			.def_readwrite("model_parameters_", &AgentInfo::model_parameters_)
			.def_readwrite("cost_parameters_", &AgentInfo::cost_parameters_);

		py::class_<CouplingInfo>(m, "CouplingInfo")
			.def(py::init<>())
			.def_readwrite("agent_id_", &CouplingInfo::agent_id_)
			.def_readwrite("neighbor_id_", &CouplingInfo::neighbor_id_)
			.def_readwrite("model_name_", &CouplingInfo::model_name_)
			.def_readwrite("model_parameters_", &CouplingInfo::model_parameters_);

		py::class_<CommunicationInfo>(m, "CommunicationInfo")
			.def(py::init<>())
			.def_readwrite("agent_type_", &CommunicationInfo::agent_type_)
			.def_readwrite("id_", &CommunicationInfo::id_)
			.def_readwrite("ip_", &CommunicationInfo::ip_)
			.def_readwrite("port_", &CommunicationInfo::port_);

		py::class_<OptimizationInfo>(m, "OptimizationInfo")
			.def(py::init<>())
			// Common parameters
			.def_readwrite("COMMON_Thor_", &OptimizationInfo::COMMON_Thor_)
			.def_readwrite("COMMON_dt_", &OptimizationInfo::COMMON_dt_)
			.def_readwrite("COMMON_Nhor_", &OptimizationInfo::COMMON_Nhor_)
			.def_readwrite("COMMON_ShiftControl_", &OptimizationInfo::COMMON_ShiftControl_)
			.def_readwrite("COMMON_Integrator_", &OptimizationInfo::COMMON_Integrator_)

			// parameters for GRAMPC
			.def_readwrite("GRAMPC_MaxGradIter_", &OptimizationInfo::GRAMPC_MaxGradIter_)
			.def_readwrite("GRAMPC_MaxMultIter_", &OptimizationInfo::GRAMPC_MaxMultIter_)
			.def_readwrite("GRAMPC_PenaltyMin_", &OptimizationInfo::GRAMPC_PenaltyMin_)
			.def_readwrite("GRAMPC_PenaltyMax_", &OptimizationInfo::GRAMPC_PenaltyMax_)
			.def_readwrite("GRAMPC_AugLagUpdateGradientRelTol_", &OptimizationInfo::GRAMPC_AugLagUpdateGradientRelTol_)
			.def_readwrite("GRAMPC_Integrator_", &OptimizationInfo::GRAMPC_Integrator_)
			.def_readwrite("GRAMPC_LineSearchType_", &OptimizationInfo::GRAMPC_LineSearchType_)
			.def_readwrite("GRAMPC_PenaltyIncreaseFactor_", &OptimizationInfo::GRAMPC_PenaltyIncreaseFactor_)
			.def_readwrite("GRAMPC_PenaltyDecreaseFactor_", &OptimizationInfo::GRAMPC_PenaltyDecreaseFactor_)
			.def_readwrite("GRAMPC_LineSearchMax_", &OptimizationInfo::GRAMPC_LineSearchMax_)
			.def_readwrite("GRAMPC_LineSearchMin_", &OptimizationInfo::GRAMPC_LineSearchMin_)
			.def_readwrite("GRAMPC_ConvergenceCheck_", &OptimizationInfo::GRAMPC_ConvergenceCheck_)
			.def_readwrite("GRAMPC_ConvergenceGradientRelTol_", &OptimizationInfo::GRAMPC_ConvergenceGradientRelTol_)
			.def_readwrite("GRAMPC_ConstraintsAbsTol_", &OptimizationInfo::GRAMPC_ConstraintsAbsTol_)
			.def_readwrite("GRAMPC_PenaltyIncreaseThreshold_", &OptimizationInfo::GRAMPC_PenaltyIncreaseThreshold_)
			.def_readwrite("GRAMPC_LineSearchInit_", &OptimizationInfo::GRAMPC_LineSearchInit_)

			// parameters for ADMM
			.def_readwrite("ADMM_maxIterations_", &OptimizationInfo::ADMM_maxIterations_)
			.def_readwrite("ADMM_ConvergenceTolerance_", &OptimizationInfo::ADMM_ConvergenceTolerance_)
			.def_readwrite("ADMM_PenaltyIncreaseFactor_", &OptimizationInfo::ADMM_PenaltyIncreaseFactor_)
			.def_readwrite("ADMM_PenaltyDecreaseFactor_", &OptimizationInfo::ADMM_PenaltyDecreaseFactor_)
			.def_readwrite("ADMM_PenaltyMin_", &OptimizationInfo::ADMM_PenaltyMin_)
			.def_readwrite("ADMM_PenaltyMax_", &OptimizationInfo::ADMM_PenaltyMax_)
			.def_readwrite("ADMM_PenaltyInit_", &OptimizationInfo::ADMM_PenaltyInit_)
			.def_readwrite("ADMM_AdaptPenaltyParameter_", &OptimizationInfo::ADMM_AdaptPenaltyParameter_)
			.def_readwrite("ADMM_innerIterations_", &OptimizationInfo::ADMM_innerIterations_)
			.def_readwrite("ADMM_DebugCost_", &OptimizationInfo::ADMM_DebugCost_)

			// parameters for neighbor approximation
			.def_readwrite("APPROX_ApproximateCost_", &OptimizationInfo::APPROX_ApproximateCost_)
			.def_readwrite("APPROX_ApproximateConstraints_", &OptimizationInfo::APPROX_ApproximateConstraints_)
			.def_readwrite("APPROX_ApproximateDynamics_", &OptimizationInfo::APPROX_ApproximateDynamics_);
	}
}
