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

#include "grampcd/optim/solver_local.hpp"


namespace grampcd
{

    SolverLocal::~SolverLocal()
    {}

    void SolverLocal::initialize_ADMM()
    {
        // implemented in derived classes
    }

    void SolverLocal::update_couplingStates()
    {
        // implemented in derived classes
    }

    void SolverLocal::send_couplingStates()
    {
        // implemented in derived classes
    }

    void SolverLocal::update_multiplierStates()
    {
        // implemented in derived classes
    }

    void SolverLocal::send_multiplierStates()
    {
        // implemented in derived classes
    }

    void SolverLocal::send_numberofNeighbors()
    {
        // implemented in derived classes
    }

    void SolverLocal::initialize_Sensi()
    {
        // implemented in derived classes
    }

    void SolverLocal::update_sensiStates()
    {
        // implemented in derived classes
    }
}
