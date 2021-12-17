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

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{
	struct MessageHandler
	{
	public:
		MessageHandler(CommunicationInterfaceLocal* communication_interface, const LoggingPtr& log);

		void handle_message(const CommunicationDataPtr& comm_data, const MessagePtr& message);

		CommunicationInterfaceLocal* communication_interface_;
		const LoggingPtr log_;
	};
}