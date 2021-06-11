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

#include "grampcd/util/logging.hpp"

namespace grampcd
{
	std::ostream& Logging::print(const DebugType type) const
	{
		std::cout.clear();

		switch (type)
		{
		case DebugType::Error:
			if(set_print_error_)
				return std::cerr;
			break;
		case DebugType::Warning:
			if (set_print_warning_)
				return std::cerr;
			break;
		case DebugType::Message:
			if(set_print_message_)
				return std::cout;
			break;
		case DebugType::Base:
			if(set_print_base_)
				return std::cout;
			break;
		default:
			break;
		}

		std::cout.setstate(std::ios_base::badbit);
		return std::cout;
	}
}