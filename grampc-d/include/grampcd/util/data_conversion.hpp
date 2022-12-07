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

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{
	class DataConversion
	{
	public:
		/*Erase element from list of integers*/
		static void erase_element_from_vector(std::vector<int>& list, const int element);
		/*Erase element from list of coupling infos*/
		static void erase_element_from_vector(std::vector<std::shared_ptr<CouplingInfo>>& list, const CouplingInfo& element);
		/*Erase element from list of agent infos*/
		static void erase_element_from_vector(std::vector< AgentInfo >& list, const AgentInfo& element);
		/*Erase element from list of coupling infos*/
		static void erase_element_from_vector(std::vector< CouplingInfo >& list, const CouplingInfo& element);
		/*Erase element from list of communication data*/
		static void erase_element_from_vector(std::vector< CommunicationDataPtr >& list, const CommunicationDataPtr& element);
		/*Erase element from list of neighbor pointers*/
		static void erase_element_from_vector(std::vector< NeighborPtr >& list, const NeighborPtr& element);
		/*Erase element from list of agents*/
		static void erase_element_from_vector(std::vector< AgentPtr >& list, const AgentInfo& element);

		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector<int>& list, const int element);
		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector<std::shared_ptr<CouplingInfo>>& list, const CouplingInfo& element);
		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector< CouplingInfo >& list, const CouplingInfo& element);
		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector< CommunicationDataPtr >& list, const std::shared_ptr<CommunicationData>& element);
		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector<NeighborPtr>& list, const int neighbor_id);
		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector<NeighborPtr>& list, const NeighborPtr& element);
		/*Returns true if element is in vector*/
		static const bool is_element_in_vector(const std::vector<AgentInfo>& list, const int id);

		/*Returns specific element from vector*/
		static const NeighborPtr get_element_from_vector(const std::vector<NeighborPtr>& list, const int id);
	};
}