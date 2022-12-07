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