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

#include "dmpc/util/types.hpp"

namespace dmpc
{

    /*@brief Info describing the coupling between agents.*/
    struct CouplingInfo
    {
    public:
        /*Id of the agent.*/
        int agent_id_;
        /*Id of the neighbor.*/
        int neighbor_id_;
        /*Name of the coupling model.*/
        std::string model_name_;
        /*Parameters for the coupling model.*/
        std::vector<typeRNum> model_parameters_;

        const bool operator== (const CouplingInfo &info) const
        {
            return ( agent_id_ == info.agent_id_ ) && ( neighbor_id_ == info.neighbor_id_ ) && ( model_name_ == info.model_name_ );
        }
    };

}