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

    /*@brief Class that provides informations for a TCP connection.*/
    struct CommunicationInfo
    {
    public:
        std::string agent_type_ = "";
        int id_ = -1;
        std::string ip_ = "";
		std::string port_ = "";

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_type_, id_, ip_, port_);
		}
    };

}