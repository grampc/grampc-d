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

#include <string>
#include <vector>
#include <map>
#include <memory>
// Note: shared_ptr is included in memory, but not detected by code completion of QtCreator
//#include <bits/shared_ptr.h>
#include <grampc.hpp>

namespace grampcd
{

    // Macro that defines shared_ptr to type
    #define DMPC_DECLARE_PTR(Name, Type) \
        typedef std::shared_ptr<Type> Name##Ptr; \
        typedef std::shared_ptr<const Type> Name##ConstPtr;

    // Macro that forward declares a class
    #define DMPC_CLASS_FORWARD(C) \
        class C; \
        DMPC_DECLARE_PTR(C, C);

    // Macro that forward declares a struct
    #define DMPC_STRUCT_FORWARD(C) \
        struct C; \
        DMPC_DECLARE_PTR(C, C);

    // Define ProblemDescriptionPtr
    DMPC_DECLARE_PTR(ProblemDescription, grampc::ProblemDescription)

    // Define SolverPtr
    DMPC_DECLARE_PTR(Solver, grampc::Grampc)

}