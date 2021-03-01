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

#ifndef Logging_HPP
#define Logging_HPP

#include <iostream>

#include <memory>

class Logging
{
public:
	enum MessageType
	{
		Error,
		Warning,
		Message,
		Base
	};

	/*Activate printing messages of type base*/
	void set_print_base(bool print) { set_print_base_ = print; }
	/*Activate printing messages of type message*/
	void set_print_message(bool print) { set_print_message_ = print; }
	/*Activate printing messages of type warning*/
	void set_print_warning(bool print) { set_print_warning_ = print; }
	/*Activate printing messages of type error*/
	void set_print_error(bool print) { set_print_error_ = print; }

	std::ostream& print_debug(MessageType type) const;

private:
	bool set_print_base_ = true;
	bool set_print_message_ = false;
	bool set_print_warning_ = false;
	bool set_print_error_ = false;
};

typedef std::shared_ptr<Logging> LoggingPtr;

#endif // Logging_HPP

