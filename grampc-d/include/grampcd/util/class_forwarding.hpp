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

#include "grampcd/util/types.hpp"

namespace grampcd
{
	enum class ADMMStep;

	DMPC_CLASS_FORWARD(CommunicationInterfaceLocal);
	DMPC_CLASS_FORWARD(CommunicationInterface);

	DMPC_CLASS_FORWARD(SolverLocal);
	DMPC_CLASS_FORWARD(Solution);

	DMPC_CLASS_FORWARD(Agent);
	DMPC_CLASS_FORWARD(Neighbor);
	DMPC_CLASS_FORWARD(Coordinator);

	DMPC_CLASS_FORWARD(DmpcInterface);

	DMPC_CLASS_FORWARD(Simulator);

	DMPC_CLASS_FORWARD(CommunicationData);

	DMPC_CLASS_FORWARD(ModelFactory);
	DMPC_CLASS_FORWARD(AgentModel);
	DMPC_CLASS_FORWARD(CouplingModel);

	DMPC_CLASS_FORWARD(Logging);

	DMPC_CLASS_FORWARD(ApproximateNeighbor);
	DMPC_CLASS_FORWARD(SolverCentral);
	DMPC_CLASS_FORWARD(SolverLocal);
	DMPC_CLASS_FORWARD(Simulator);

	DMPC_STRUCT_FORWARD(AgentInfo);
	DMPC_STRUCT_FORWARD(CouplingInfo);

	DMPC_STRUCT_FORWARD(AgentState);
	DMPC_STRUCT_FORWARD(CouplingState);
	DMPC_STRUCT_FORWARD(MultiplierState);
	DMPC_STRUCT_FORWARD(PenaltyState);

	DMPC_STRUCT_FORWARD(AgentInfo);
	DMPC_STRUCT_FORWARD(CouplingInfo);
	DMPC_STRUCT_FORWARD(OptimizationInfo);
	DMPC_STRUCT_FORWARD(CommunicationInfo);

	DMPC_STRUCT_FORWARD(Message);
	DMPC_STRUCT_FORWARD(MessageHandler);
}