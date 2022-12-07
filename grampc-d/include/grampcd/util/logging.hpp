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

#include <iostream>

#include <memory>

namespace grampcd
{
	enum class DebugType
	{
		Error,
		Warning,
		Message,
		Base,
		Progressbar
	};

	class Logging
	{
	public:
		std::ostream& print(const DebugType type) const;

		bool print_base_ = true;
		bool print_message_ = false;
		bool print_warning_ = false;
		bool print_error_ = false;
		bool print_progressbar_ = false;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(print_base_, print_message_, print_warning_, print_error_, print_progressbar_);
		}
	};
}