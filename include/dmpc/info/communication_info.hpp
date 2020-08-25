/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#ifndef COMMUNICATION_INFO_HPP
#define COMMUNICATION_INFO_HPP

#include "dmpc/util/types.hpp"

namespace dmpc
{

/*@brief Class that provides informations for a TCP connection.*/
struct CommunicationInfo
{
public:
    CommunicationInfo() {}

    std::string agent_type_ = "";
    int id_ = -1;
    std::string ip_ = "";
    std::string port_ = "";
};

typedef std::shared_ptr<CommunicationInfo> CommunicationInfoPtr;

}

#endif // COMMUNICATION_INFO_HPP
