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

#include "grampcd/util/logging.hpp"

namespace grampcd
{
	std::ostream& Logging::print(const DebugType type) const
	{
		std::cout.clear();

		switch (type)
		{
		case DebugType::Error:
			if(print_error_)
				return std::cerr;
			break;
		case DebugType::Warning:
			if (print_warning_)
				return std::cerr;
			break;
		case DebugType::Message:
			if(print_message_)
				return std::cout;
			break;
		case DebugType::Base:
			if(print_base_)
				return std::cout;
			break;
		case DebugType::Progressbar:
			if (print_progressbar_)
				return std::cout;
			break;
		default:
			break;
		}

		std::cout.setstate(std::ios_base::badbit);
		return std::cout;
	}
}