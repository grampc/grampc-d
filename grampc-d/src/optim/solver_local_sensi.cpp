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

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

#include "grampcd/optim/solver_local_sensi.hpp"
#include "grampcd/optim/optim_util.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

#include "grampcd/util/logging.hpp"

#include "grampcd/agent/step_selector.hpp"
#include "grampcd/agent/async_step_selector.hpp"

#include "grampcd/optim/solution.hpp"
#include <cmath>
#include <algorithm>

namespace grampcd
{

    SolverLocalSensi::SolverLocalSensi(Agent* agent, const OptimizationInfo& info, const LoggingPtr& log, const CommunicationInterfacePtr& communication_interface)
        : agent_(agent),
        sensi_problem_description_(agent),
        solver_(new grampc::Grampc(&sensi_problem_description_)),
        info_(info),
        log_(log),
        communication_interface_(communication_interface)
    {

        configureSolver(solver_, info);
    }

    void SolverLocalSensi::initialize_Sensi()
    {
        // steps needed for asynchronous execution 
        send_numberofNeighbors();
        agent_->initialize_allNeighborDelays();
        agent_->reset_stopAdmmflag_of_neighbors();
    }

    void SolverLocalSensi::update_agentStates()
    {
        const unsigned int Nx = solver_->getParameters()->Nx;
        const unsigned int Nu = solver_->getParameters()->Nu;
        const unsigned int Ng = solver_->getParameters()->Ng;
        const unsigned int Nh = solver_->getParameters()->Nh;
        const unsigned int Nxi = agent_->get_Nxi();
        const unsigned int Nui = agent_->get_Nui();
        const unsigned int Nhor = solver_->getOptions()->Nhor;

        // set initial state x0 and initial control
        {
            const AgentState& state = agent_->get_agentState();

            // set initial time
            solver_->setparam_real("t0", state.t0_);

            // set initial state
            for (unsigned int k = 0; k < Nxi; ++k)
                solver_->getParameters()->x0[k] = state.x_[k];

            // set initial control trajectory
            for (unsigned int i = 0; i < Nhor; ++i)
            {
                for (unsigned int k = 0; k < Nui; ++k)
                    solver_->getWorkspace()->u[i * Nu + k] = state.u_[i * Nui + k];
            }
        }

        // construct control limits
        std::vector<typeRNum> umin(solver_->getParameters()->Nu, -INF);
        std::vector<typeRNum> umax(solver_->getParameters()->Nu, INF);

        agent_->get_agentModel()->get_controlLimits(umin, umax);

        solver_->setparam_real_vector("umin", &umin[0]);
        solver_->setparam_real_vector("umax", &umax[0]);

        // run GRAMPC
        solver_->run();

        // update x_i and u_i
        AgentState state = agent_->get_agentState();
        for (unsigned int i = 0; i < Nhor; ++i)
        {
            // update state trajectory x_i
            for (unsigned int k = 0; k < Nxi; ++k)
            {
                state.x_[i * Nxi + k] = solver_->getWorkspace()->x[i * Nx + k];
            }

            // update control trajectory u_i
            for (unsigned int k = 0; k < Nui; ++k)
            {
                state.u_[i * Nui + k] = solver_->getWorkspace()->u[i * Nu + k];

                // consider simple convex sum 
                if (agent_->get_optimizationInfo().SENSI_ConvexSum_)
                {
                    typeRNum alpha = agent_->get_optimizationInfo().SENSI_ConvexSumAlpha_;
                    state.u_[i * Nui + k] = alpha * solver_->getWorkspace()->u[i * Nu + k] + (1 - alpha) * agent_->get_previous_agentState().u_[i * Nu + k];
                }
            }

            // update adjoint states
            for (unsigned int k = 0; k < Nxi; ++k)
            {
                state.lambda_[i * Nxi + k] = solver_->getWorkspace()->adj[i * Nx + k];
            }
        }
        agent_->set_agentState(state);

        // update Lagrange multipliers of coupled constraints 
        int idx_eq = agent_->get_agentModel()->get_Ngi();
        int idx_ineq = agent_->get_agentModel()->get_Nhi();

        for (const auto neighbor : agent_->get_sendingNeighbors())
        {
            ConstraintState constraintState = neighbor->get_coupled_constraintState();

            unsigned int N_gij = neighbor->get_couplingModel()->get_Ngij();
            unsigned int N_hij = neighbor->get_couplingModel()->get_Nhij();

            for (unsigned int i = 0; i < Nhor; i++)
            {
                for (unsigned int k = 0; k < N_gij; ++k)
                {
                    // get equality constraint info
                    constraintState.mu_g_[i * N_gij + k] = solver_->getWorkspace()->mult[i * (Ng + Nh) + idx_eq + k];
                    constraintState.c_g_[i * N_gij + k] = solver_->getWorkspace()->pen[i * (Ng + Nh) + idx_eq + k];                  
                }
                for (unsigned int k = 0; k < N_hij; ++k)
                {
                    // get inequality constraint info
                    constraintState.mu_h_[i * N_hij + k] = solver_->getWorkspace()->mult[i * (Ng + Nh) + idx_ineq + k + Ng]; 
                    constraintState.c_h_[i * N_hij + k] = solver_->getWorkspace()->pen[i * (Ng + Nh) + idx_ineq + k + Ng];
                }             
            }
            // increase index
            idx_eq += N_gij;
            idx_ineq += N_hij;

            // set states
            neighbor->set_coupled_constraintState(constraintState);
        }
    }
       

    void SolverLocalSensi::update_sensiStates()
    {
        // in neighbor-affine form the agent can directly calculate its needed sensitivites 
        // dJ_j/du_i, dJ_j/dx_i
            const unsigned int Nxi = agent_->get_Nxi();
            const unsigned int Nui = agent_->get_Nui();
            unsigned int Nhor = solver_->getOptions()->Nhor;

            AgentState agentState = agent_->get_agentState();

            for (const NeighborPtr& neighbor : agent_->get_receivingNeighbors())
            {
               const unsigned int Nxj = neighbor->get_Nxj();
               const unsigned int Nuj = neighbor->get_Nuj();
               SensiState sensiState = neighbor->get_sensiState();

               const AgentState& neighbors_agentState = neighbor->get_neighbors_agentState();
               const ConstraintState& neighbors_constraintState = neighbor->get_neighbors_coupled_constraintState();

               // We need the coupling model from the perspective of the receiving neighbor f_ji(x_j,u_j,x_i,u_i,lambda_j), l_ji(x_j,u_j,x_i,u_i), Vji(x_j,x_i)
               // gji(x_j,u_j,x_i,u_i,mu_ji) and hji(x_j,u_j,x_i,u_i,mu_ji)
               const CouplingModelPtr& copied_coupling_model = neighbor->get_copied_couplingModel();
               const unsigned int Ngji = copied_coupling_model->get_Ngij();
               const unsigned int Nhji = copied_coupling_model->get_Nhij();

               // reset sensi state 
               resetState(sensiState, neighbor->get_id(),sensiState.t_);
               sensiState.psi_x_.resize(Nxi * Nhor, 0.0);
               sensiState.psi_u_.resize(Nui * Nhor, 0.0);
               sensiState.psi_V_.resize(Nxi, 0.0);
               if (info_.SENSI_higherOrder_)
               {
                   sensiState.psi_xx_.resize(Nxi * Nxi * Nhor, 0.0);
                   sensiState.psi_uu_.resize(Nui * Nui * Nhor, 0.0);
                   sensiState.psi_xu_.resize(Nxi * Nui * Nhor, 0.0);
                   sensiState.psi_VV_.resize(Nxi * Nxi, 0.0);
               }

               // calculate sensitivites
               for (unsigned int i = 0; i < Nhor;++i)
               {   
                   // first-order sensitivities
                   // sensi w.r.t dynamics and costs
                   copied_coupling_model->dfdxj_vec(&sensiState.psi_x_[i * Nxi], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_agentState.lambda_.data() + i * Nxj);
                   copied_coupling_model->dldxj(&sensiState.psi_x_[i * Nxi], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                   copied_coupling_model->dfduj_vec(&sensiState.psi_u_[i * Nui], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_agentState.lambda_.data() + i * Nxj);
                   copied_coupling_model->dlduj(&sensiState.psi_u_[i * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                   // sensi w.r.t equality constraints 
                   copied_coupling_model->dgdxj_vec(&sensiState.psi_x_[i * Nxi], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_g_.data() + i * Ngji);
                   copied_coupling_model->dgduj_vec(&sensiState.psi_u_[i * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_g_.data() + i * Ngji);
                   
                   // sensi w.r.t transformed inequality constraints 
                   for (unsigned int k = 0; k < Nhji; ++k)
                   {
                      typeRNum cfct = 0;
                      copied_coupling_model->hfct(&cfct, agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                      if (cfct > - neighbors_constraintState.mu_h_[i * k] / neighbors_constraintState.c_h_[i * k])
                      {
                          copied_coupling_model->dhdxj_vec(&sensiState.psi_x_[i * Nxi], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_h_.data() + i * Nhji);
                          copied_coupling_model->dhduj_vec(&sensiState.psi_u_[i * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_h_.data() + i * Nhji);
                      }
                   }

                   // second-order sensitivities 
                   if (info_.SENSI_higherOrder_)
                   {                    
                       // sensi w.r.t dynamics and costs
                       copied_coupling_model->dfdxjdxj_vec(&sensiState.psi_xx_[i * Nxi * Nxi], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_agentState.lambda_.data() + i * Nxj);
                       copied_coupling_model->dldxjdxj(&sensiState.psi_xx_[i * Nxi * Nxi], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                       copied_coupling_model->dfdujduj_vec(&sensiState.psi_uu_[i * Nui * Nui], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_agentState.lambda_.data() + i * Nxj);
                       copied_coupling_model->dldujduj(&sensiState.psi_uu_[i * Nui * Nui], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                       copied_coupling_model->dfdxjduj_vec(&sensiState.psi_xu_[i * Nxi * Nui], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_agentState.lambda_.data() + i * Nxj);
                       copied_coupling_model->dldxjduj(&sensiState.psi_xu_[i * Nxi * Nui], agentState.t0_ + agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                       // sensi w.r.t equality constraints 
                       copied_coupling_model->dgdxjdxj_vec(&sensiState.psi_xx_[i * Nxi * Nxi], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_g_.data() + i * Ngji);
                       copied_coupling_model->dgdujduj_vec(&sensiState.psi_uu_[i * Nui * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_g_.data() + i * Ngji);
                       copied_coupling_model->dgdxjduj_vec(&sensiState.psi_xu_[i * Nxi * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_g_.data() + i * Ngji);
                      
                       // sensi w.r.t transformed inequality constraints 
                       for (unsigned int k = 0; k < Nhji; ++k)
                       {
                           typeRNum cfct = 0;
                           copied_coupling_model->hfct(&cfct, agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui);

                           if (cfct > -neighbors_constraintState.mu_h_[i * k] / neighbors_constraintState.c_h_[i * k])
                           {
                               copied_coupling_model->dhdxjdxj_vec(&sensiState.psi_xx_[i * Nxi * Nxi], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_h_.data() + i * Nhji);
                               copied_coupling_model->dhdujduj_vec(&sensiState.psi_uu_[i * Nui * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_h_.data() + i * Nhji);
                               copied_coupling_model->dhdxjduj_vec(&sensiState.psi_xu_[i * Nxi * Nui], agentState.t_[i], neighbors_agentState.x_.data() + i * Nxj, neighbors_agentState.u_.data() + i * Nuj, agentState.x_.data() + i * Nxi, agentState.u_.data() + i * Nui, neighbors_constraintState.mu_h_.data() + i * Nhji);
                           }
                       }
                   }
               }
               // calculate sensitivity of terminal cost 
               copied_coupling_model->dVdxj(&sensiState.psi_V_[0], agentState.t0_ + agentState.t_.back(), agentState.x_.data() + ((Nhor - 1) * Nxi), neighbors_agentState.x_.data() + ((Nhor - 1) * Nxj));

               // calculate second order sensitivity of terminal cost 
               if (info_.SENSI_higherOrder_)
               {
                   copied_coupling_model->dVdxjdxj(&sensiState.psi_VV_[0], agentState.t0_ + agentState.t_.back(), agentState.x_.data() + ((Nhor - 1) * Nxi), neighbors_agentState.x_.data() + ((Nhor - 1) * Nxj));
               }

               // set sensi State
               neighbor->set_sensiState(sensiState);
            }
    }

    void SolverLocalSensi::send_agentStates()
    {
        const auto state = agent_->get_agentState();

        // send agent State to all neighbors
        for (const auto& neighbor : agent_->get_neighbors())
        {
            communication_interface_->send_agentState(state, neighbor->get_coupled_constraintState(), agent_->get_id(), neighbor->get_id());
        }
    }

    const bool SolverLocalSensi::is_converged() 
    {   
       previous_cost_ = cost_;
       cost_ = agent_->get_predicted_cost();
       return info_.SENSI_ConvergenceTolerance_ > std::abs(cost_ - previous_cost_);
    }

    void SolverLocalSensi::print_debugCost()
    {
        agent_->get_solution()->update_debug_cost(agent_->get_predicted_cost());
    }

    void SolverLocalSensi::send_convergenceFlag()
    {
        communication_interface_->send_convergenceFlag(is_converged(), agent_->get_id());
    }

    void SolverLocalSensi::send_stoppedAlgFlag()
    {
        // get and cast step selector 
        AsyncStepSelectorPtr stepselector = std::static_pointer_cast<AsyncStepSelector>(agent_->get_stepSelector());

        // send flag to neighbors
        for (const auto& neighbor : agent_->get_neighbors())
        {
            communication_interface_->send_stoppedAlgFlag(stepselector->get_flagStopAlg(), agent_->get_id(), neighbor->get_id());
        }
        // send flag to Coordinator 
        communication_interface_->send_stoppedAlgFlag(stepselector->get_flagStopAlg(), agent_->get_id());
    }
}
