# This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

from pathlib import Path
import matplotlib.pyplot as plt
import sys, os

# generate path to module
path = os.getcwd()
sys.path.append(os.path.join(path, 'bin'))

path = str(Path(path).parents[1])
path = os.path.join(path, 'bin')

# append path to list of folders where python is searching for modules
sys.path.append(path)

#import python interface
import grampcd_interface

# create interface
interface = grampcd_interface.interface()

# show logging
interface.set_print_message(1)
interface.set_print_warning(1)
interface.set_print_error(1)

# initialize communication interface
interface.initialize_local_communicationInterface_as_coordinator(7777)

# set optimization info
optimization_info = grampcd_interface.OptimizationInfo()
optimization_info.COMMON_Nhor_ = 21
optimization_info.COMMON_Thor_ = 5
optimization_info.COMMON_dt_ = 0.1
optimization_info.GRAMPC_MaxGradIter_ = 15
optimization_info.GRAMPC_MaxMultIter_ = 1
optimization_info.ADMM_maxIterations_ = 20
optimization_info.ADMM_ConvergenceTolerance_ = 0.05

approx = 1
optimization_info.APPROX_ApproximateCost_ = approx
optimization_info.APPROX_ApproximateConstraints_ = approx
optimization_info.APPROX_ApproximateDynamics_ = approx

interface.set_optimizationInfo(optimization_info)

interface.simulate_realtime(1)

# run distributed controller
interface.run_DMPC()

