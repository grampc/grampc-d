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

#include "grampcd/util/data_conversion.hpp"

#include "grampcd/info/coupling_info.hpp"
#include "grampcd/info/communication_data.hpp"

#include "grampcd/optim/optim_util.hpp"
#include "grampcd/optim/solution.hpp"

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

namespace grampcd
{
	const bool DataConversion::is_element_in_vector(const std::vector<int>& list, const int element)
	{
		for (const auto entry : list)
		{
			if (entry == element)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<NeighborPtr>& list, const int neighbor_id)
	{
		for (const auto& neighbor : list)
		{
			if (neighbor->get_id() == neighbor_id)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<NeighborPtr>& list, const NeighborPtr& element)
	{
		for (const auto& neighbor : list)
		{
			if (neighbor == element)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<std::shared_ptr<CouplingInfo>>& list, const CouplingInfo& element)
	{
		for (const auto& info : list)
		{
			if (*info == element)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector< CouplingInfo >& list, const CouplingInfo& element)
	{
		for (const auto& info : list)
		{
			if (info == element)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector< CommunicationDataPtr >& list, const std::shared_ptr<CommunicationData>& element)
	{
		for (const auto& comm_data : list)
		{
			if (comm_data->communication_info_->id_ == element->communication_info_->id_)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<AgentInfo>& list, const int id)
	{
		for (const auto& info : list)
		{
			if (info.id_ == id)
				return true;
		}
		return false;
	}

	void DataConversion::erase_element_from_vector(std::vector< NeighborPtr >& list, const NeighborPtr& element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (list[i] == element)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	void DataConversion::erase_element_from_vector(std::vector<AgentPtr>& list, const AgentInfo& element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (list[i]->get_id() == element.id_)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	void DataConversion::erase_element_from_vector(std::vector<int>& list, const int element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (list[i] == element)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	void DataConversion::erase_element_from_vector(std::vector< AgentInfo >& list, const AgentInfo& element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (list[i] == element)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	void DataConversion::erase_element_from_vector(std::vector< CouplingInfo >& list, const CouplingInfo& element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (list[i] == element)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	void DataConversion::erase_element_from_vector(std::vector<std::shared_ptr<CouplingInfo>>& list, const CouplingInfo& element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (*list[i] == element)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	void DataConversion::erase_element_from_vector(std::vector< CommunicationDataPtr >& list, const  CommunicationDataPtr& element)
	{
		for (unsigned int i = 0; i < list.size(); ++i)
		{
			if (list[i] == element)
			{
				list.erase(list.begin() + i);
				return;
			}
		}
	}

	const NeighborPtr DataConversion::get_element_from_vector(const std::vector<NeighborPtr>& list, const int id)
	{
		for (const auto& neighbor : list)
		{
			if (neighbor->get_id() == id)
				return neighbor;
		}
		return nullptr;
	}
}