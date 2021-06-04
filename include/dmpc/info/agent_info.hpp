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

    /**
     * @brief Information about a single agent.
     */
    struct AgentInfo
    {
    public:
        AgentInfo() {}

        int id_ = -1;
        std::string model_name_ = "";
        std::vector<typeRNum> model_parameters_ = std::vector<typeRNum>(0, 0);
        std::vector<typeRNum> cost_parameters_ = std::vector<typeRNum>(0, 0);

        const bool operator== (const AgentInfo& info) const { return (id_ == info.id_) && (model_name_ == info.model_name_); }
    };

}
