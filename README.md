# GRAMPC-D
GRAMPC-D is a framework for nonlinear distributed model predictive control based on [GRAMPC](https://github.com/grampc/grampc). The implementation is in C++ allowing for sampling times in the millisecond range. A manual of GRAMPC-D is provided in the folder [documentation](doc/manual.pdf).

More details about the algorithm and its performance can be found in the corresponding article published in Optimization and Engeneering. The article is available online under open access at: https://doi.org/10.1007/s11081-021-09605-3. Please cite the paper when you are using results obtained with GRAMPC-D.

## Features
The most important features are listed below. For a more detailed description of the features, please refer to the [documentation](doc/manual.pdf).
- Real-time solution of large-scale optimal control problems via the alternating direction method of multipliers (ADMM) or sensitivity-based distributed programming (SBDP), where the local problems are solved by [GRAMPC](https://github.com/grampc/grampc)
- Agent-based implementation allowing for a modular problem description 
- Provides a communication interface to enable communication between agents via TCP/IP without additional implementation effort
- Allows for plug-and-play operations during runtime 
- Interface to Python 

## Installation
Clone the repository and all submodules:
```
git clone --recurse-submodules https://github.com/grampc/grampc-d.git
```
To build GRAMPC-D using Linux, run the CMakeLists.txt and call make afterwards:
```
cmake -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
make
```
To build GRAMPC-D using Windows, the editor Visual Studio Code with the [C++ extension](https://code.visualstudio.com/docs/languages/cpp) is recommended. Open the folder "grampc-d" in Visual Studio Code, click the button "Build", and select a compiler from the list.

## License
GRAMPC-D is licensed under BSD 3-clause license, see LICENSE.txt.
