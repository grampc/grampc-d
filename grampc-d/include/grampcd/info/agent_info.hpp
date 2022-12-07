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

#include "grampcd/util/types.hpp"

namespace grampcd
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

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(id_, model_name_, model_parameters_, cost_parameters_);
		}
    };

}
