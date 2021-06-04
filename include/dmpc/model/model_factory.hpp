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

#include "dmpc/util/class_forwarding.hpp"

namespace dmpc
{

    /* @brief The model factory creates AgentModels and CouplingModels given AgentInfo and CouplingInfo.*/
    class ModelFactory
    {
    public:

	    virtual ~ModelFactory();

        /*Creates and returns an agent model regarding the agent info.*/
        virtual AgentModelPtr create_agentModel(const AgentInfo& info) const = 0;

	    /*Creates and returns a coupling model regarding the coupling info.*/
        virtual CouplingModelPtr create_couplingModel(const CouplingInfo& info) const = 0;
    };

}