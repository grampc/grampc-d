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

    /**
     * @brief Information about the optimization.
     */
    struct TuningInfo
    {
    public:
        TuningInfo() {}

		std::vector<typeRNum> COMMON_Nhor_ = { 11, 21 };
        std::vector<typeRNum> GRAMPC_MaxGradIter_ = {3, 20};
        std::vector<typeRNum> GRAMPC_MaxMultIter_ = {1, 3};

		std::vector<typeRNum> ADMM_maxIterations_ = { 5, 20 }; 
		std::vector<typeRNum> ADMM_innerIterations_ = { 1, 3 }; 
    };

}