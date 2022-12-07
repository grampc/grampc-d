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

#include "grampcd/interface/dmpc_interface.hpp"


int main(int argc, char** argv)
{
    // create interface
    const auto interface = std::make_shared<grampcd::DmpcInterface>();

    // initialize communication interface
    interface->initialize_central_communicationInterface();

    // set optimization info
    auto optimization_info = interface->optimization_info();
    optimization_info.COMMON_Nhor_ = 21;
    optimization_info.COMMON_Thor_ = 2;
    optimization_info.COMMON_dt_ = 0.02;
    optimization_info.GRAMPC_MaxGradIter_ = 15;
    optimization_info.GRAMPC_MaxMultIter_ = 1;
    optimization_info.COMMON_Solver_ = "ADMM";
    optimization_info.ADMM_maxIterations_ = 10;
    optimization_info.ADMM_ConvergenceTolerance_ = 0.001;

    optimization_info.ASYNC_Active_ = false;
    optimization_info.ASYNC_Delay_ = 0;

    interface->set_optimizationInfo(optimization_info);

    typeRNum Tsim = 4;

    // parameters for cost function
    typeRNum P_x = 1; typeRNum P_vx = 1; typeRNum P_y = 1; typeRNum P_vy = 1;
    typeRNum Q_x = 5; typeRNum Q_vx = 2; typeRNum Q_y = 5; typeRNum Q_vy = 2;
    typeRNum R_ux = 0.01; typeRNum R_uy = 0.01;

    // model parameters
    typeRNum m_agent = 7.5; typeRNum c = 0.5;

    // number of agents
    const unsigned int n_agents_x = 5;
    const unsigned int n_agents_y = 5;
    const unsigned int n_agents = n_agents_x * n_agents_y;

    // register agents
    auto agentInfo = interface->agent_info();
    agentInfo.model_name_ = "ssms2d_agentModel";
    agentInfo.model_parameters_ = { m_agent, 1 };
    agentInfo.cost_parameters_ = { P_x, P_vx, P_y, P_vy, Q_x, Q_vx, Q_y, Q_vy, R_ux, R_uy };

    for (unsigned int i = 0; i < n_agents_x; ++i)
    {
        for (unsigned int j = 0; j < n_agents_y; ++j)
        {
            agentInfo.id_ = i * n_agents_x + j;

            //define offset in x and y
            ctypeRNum offset_x = i % 2 ? -0.2 : 0.2;
            ctypeRNum offset_y = j % 2 ? -0.2 : 0.2;

            const std::vector<typeRNum> x_init = { static_cast<typeRNum>(i + offset_x), 0.0, static_cast<typeRNum>(j + offset_y), 0.0 };
            const std::vector<typeRNum> x_des = { static_cast<typeRNum>(i), 0, static_cast<typeRNum>(j), 0 };
            interface->register_agent(agentInfo, x_init, { 0, 0 }, x_des, { 0, 0 });
        }
    }

    // register couplings
    auto coupling_info = interface->coupling_info();
    coupling_info.model_name_ = "ssms2d_couplingModel";
    coupling_info.model_parameters_ = { m_agent, c };

    unsigned int idx = 0;
    for (unsigned int i = 0; i < n_agents_y; ++i)
    {
        for (unsigned j = 0; j < n_agents_x; ++j)
        {
            coupling_info.agent_id_ = idx;

            // coupling with neighbor on the left
            if (j > 0)
            {
                coupling_info.neighbor_id_ = idx - 1;
                interface->register_coupling(coupling_info);
            }

            // coupling with neighbor on the right
            if (j < n_agents_x - 1)
            {
                coupling_info.neighbor_id_ = idx + 1;
                interface->register_coupling(coupling_info);
            }

            // coupling with neighbor above
            if (i > 0)
            {
                coupling_info.neighbor_id_ = idx - n_agents_x;
                interface->register_coupling(coupling_info);
            }

            // coupling with neighbor below
            if (i < n_agents_y - 1)
            {
                coupling_info.neighbor_id_ = idx + n_agents_x;
                interface->register_coupling(coupling_info);
            }

            idx = idx + 1;
        }
    }

    //activate printing of progressbar
    interface->set_print_progressbar(true);

    // run MPC
    interface->run_MPC(0, Tsim);

    // print solution
    interface->print_solution_to_file("all", "MPC_");

    // deregister agents
    for (unsigned int i = 0; i < n_agents_x; ++i)
    {
        for (unsigned int j = 0; j < n_agents_y; ++j)
        {
            agentInfo.id_ = i * n_agents_x + j;
            interface->deregister_agent(agentInfo);
        }
    }

    // register agents
    for (unsigned int i = 0; i < n_agents_x; ++i)
    {
        for (unsigned int j = 0; j < n_agents_y; ++j)
        {
            agentInfo.id_ = i * n_agents_x + j;

            //define offset in x and y
            ctypeRNum offset_x = i % 2 ? -0.4 : 0.4;
            ctypeRNum offset_y = j % 2 ? -0.3 : 0.3;
            const std::vector<typeRNum> x_init = { static_cast<typeRNum>(i + offset_x), 0.0, static_cast<typeRNum>(j + offset_y), 0.0 };
            const std::vector<typeRNum> x_des = { static_cast<typeRNum>(i), 0, static_cast<typeRNum>(j), 0 };
            interface->register_agent(agentInfo, x_init, { 0, 0 }, x_des, { 0, 0 });
        }
    }

    idx = 0;
    for (unsigned int i = 0; i < n_agents_y; ++i)
    {
        for (unsigned j = 0; j < n_agents_x; ++j)
        {
            coupling_info.agent_id_ = idx;

            // coupling with neighbor on the left
            if (j > 0)
            {
                coupling_info.neighbor_id_ = idx - 1;
                interface->register_coupling(coupling_info);
            }

            // coupling with neighbor on the right
            if (j < n_agents_x - 1)
            {
                coupling_info.neighbor_id_ = idx + 1;
                interface->register_coupling(coupling_info);
            }

            // coupling with neighbor above
            if (i > 0)
            {
                coupling_info.neighbor_id_ = idx - n_agents_x;
                interface->register_coupling(coupling_info);
            }

            // coupling with neighbor below
            if (i < n_agents_y - 1)
            {
                coupling_info.neighbor_id_ = idx + n_agents_x;
                interface->register_coupling(coupling_info);
            }

            idx = idx + 1;
        }
    }

    // run DMPC
    interface->run_DMPC(0, Tsim);

    // print solution
    interface->print_solution_to_file("all", "DMPC_");

    return 0;
}
