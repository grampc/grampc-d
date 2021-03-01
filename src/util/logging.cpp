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

#include "dmpc/util/logging.hpp"

std::ostream& Logging::print_debug(MessageType type) const
{
	std::cout.clear();

	if (type == Logging::Error && set_print_error_)
		return std::cerr;
	else if (type == Logging::Warning && set_print_warning_)
		return std::cerr;
	else if (type == Logging::Message && set_print_message_)
		return std::cout;
	else if (type == Logging::Base && set_print_base_)
		return std::cout;
	else
	{
		std::cout.setstate(std::ios_base::badbit);
		return std::cout;
	}
}