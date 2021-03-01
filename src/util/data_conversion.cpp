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

#include "dmpc/util/data_conversion.hpp"

const int DataConversion::ADMMStep_to_int(dmpc::ADMMStep step)
	{
		switch (step)
		{
		case dmpc::ADMMStep::UPDATE_AGENT_STATE: return 0;
		case dmpc::ADMMStep::SEND_AGENT_STATE: return 1;
		case dmpc::ADMMStep::UPDATE_COUPLING_STATE: return 2;
		case dmpc::ADMMStep::SEND_COUPLING_STATE: return 3;
		case dmpc::ADMMStep::UPDATE_MULTIPLIER_STATE: return 4;
		case dmpc::ADMMStep::SEND_MULTIPLIER_STATE: return 5;
		case dmpc::ADMMStep::SEND_CONVERGENCE_FLAG: return 6;
		case dmpc::ADMMStep::INITIALIZE: return 7;
		case dmpc::ADMMStep::SEND_TRUE_STATE: return 8;
		case dmpc::ADMMStep::PRINT: return 9;
		default: return -1;
		}
	}

const dmpc::ADMMStep DataConversion::Int_to_ADMMStep(int a)
	{
		switch (a)
		{
		case 0: return dmpc::ADMMStep::UPDATE_AGENT_STATE;
		case 1: return dmpc::ADMMStep::SEND_AGENT_STATE;
		case 2: return dmpc::ADMMStep::UPDATE_COUPLING_STATE;
		case 3: return dmpc::ADMMStep::SEND_COUPLING_STATE;
		case 4: return dmpc::ADMMStep::UPDATE_MULTIPLIER_STATE;
		case 5: return dmpc::ADMMStep::SEND_MULTIPLIER_STATE;
		case 6: return dmpc::ADMMStep::SEND_CONVERGENCE_FLAG;
		case 7: return dmpc::ADMMStep::INITIALIZE;
		case 8: return dmpc::ADMMStep::SEND_TRUE_STATE;
		case 9: return dmpc::ADMMStep::PRINT;
		default: return dmpc::ADMMStep::UPDATE_AGENT_STATE;
		}
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const bool a)
	{
		if (a)
			(*data)[pos] = 1;
		else
			(*data)[pos] = 0;
		++pos;
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const unsigned int a)
	{
		const char* ptr = reinterpret_cast<const char*>(&a);
		for (unsigned int i = 0; i < sizeof(int); ++i)
			(*data)[pos + i] = *(ptr + sizeof(int) - 1 - i);
		pos += sizeof(int);
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const int a)
	{
		const char* ptr = reinterpret_cast<const char*>(&a);
		for (unsigned int i = 0; i < sizeof(int); ++i)
			(*data)[pos + i] = *(ptr + sizeof(int) - 1 - i);
		pos += sizeof(int);
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const typeRNum a)
	{
		const char* ptr = reinterpret_cast<const char*>(&a);
		for (unsigned int i = 0; i < sizeof(typeRNum); ++i)
			(*data)[pos + i] = ptr[i];
		pos += sizeof(typeRNum);
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const std::vector<typeRNum>& a)
	{
		// insert length of array
		const unsigned int size_of_array = static_cast<unsigned int>( sizeof(typeRNum) * a.size() );
		insert_into_charArray(data, pos, size_of_array);

		if (size_of_array == 0) 
			return;

		// insert array
		const char* ptr = reinterpret_cast<const char*>(&a[0]);
		for (unsigned int k = 0; k < a.size(); ++k)
			insert_into_charArray(data, pos, a[k]);
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const std::string& a)
	{
		// insert length of string
		const unsigned int size_of_data = (int) a.size();
		insert_into_charArray(data, pos, size_of_data);

		// insert string
		const char* ptr = a.data();

		for (unsigned int k = 0; k < size_of_data; ++k)
			(*data)[pos + k] = *(ptr + k);

		pos += static_cast<unsigned int>( a.size() );
	}

	void DataConversion::insert_into_charArray(const std::shared_ptr<std::vector<char>>& data, unsigned int& pos, const char a)
	{
		(*data)[pos] = a;
		++pos;
	}

	void DataConversion::read_from_charArray(const std::vector<char>& data, unsigned int& pos, int& a)
	{
		char* ptr = reinterpret_cast<char*>(&a);
		for (unsigned int i = 0; i < sizeof(int); ++i)
			*(ptr + sizeof(int) - 1 - i) = data[pos + i];
		pos += sizeof(int);
	}

	void DataConversion::read_from_charArray(const std::vector<char>& data, unsigned int& pos, unsigned int& a)
	{
		char* ptr = reinterpret_cast<char*>(&a);
		for (unsigned int i = 0; i < sizeof(int); ++i)
			*(ptr + sizeof(int) - 1 - i) = data[pos + i];
		pos += sizeof(int);
	}

	void DataConversion::read_from_charArray(const std::vector<char>& data, unsigned int& pos, typeRNum& a)
	{
		char* ptr = reinterpret_cast<char*>(&a);
		for (unsigned int i = 0; i < sizeof(typeRNum); ++i)
			ptr[i] = data[pos + i];
		pos += sizeof(typeRNum);
	}

	void DataConversion::read_from_charArray(const std::vector<char>& data, unsigned int& pos, bool& a)
	{
		a = data[pos]==1;
		++pos;
	}

	void DataConversion::read_from_charArray(const std::vector<char>& data, unsigned int& pos, std::vector<typeRNum>& a)
	{
		// read length of array
		unsigned int size_of_array = 0;

		read_from_charArray(data, pos, size_of_array);
		size_of_array /= sizeof(typeRNum);

		if (size_of_array == 0) return;

		a.resize(size_of_array, 0);

		// read array
		char* ptr = reinterpret_cast<char*>(&a[0]);
		for (unsigned int k = 0; k < a.size(); ++k)
			read_from_charArray(data, pos, a[k]);
	}

	void DataConversion::read_from_charArray(const std::vector<char>& data, unsigned int& pos, std::string& a)
	{
		// read length of string
		unsigned int size_of_data = 0;
		read_from_charArray(data, pos, size_of_data);

		// read string
		a = "";

		for (unsigned int k = 0; k < size_of_data; ++k)
			a += data[pos + k];

		pos += static_cast<unsigned int>(a.size());
	}

	void DataConversion::skip_in_charArray_bool(const std::vector<char>& data, unsigned int& pos)
	{
		++pos;
	}

	void DataConversion::skip_in_charArray_char(const std::vector<char>& data, unsigned int& pos)
	{
		++pos;
	}

	void DataConversion::skip_in_charArray_int(const std::vector<char>& data, unsigned int& pos)
	{
		pos += sizeof(int);
	}

	void DataConversion::skip_in_charArray_typeRNum(const std::vector<char>& data, unsigned int& pos)
	{
		pos += sizeof(typeRNum);
	}

	void DataConversion::skip_in_charArray_vecTypeRNum(const std::vector<char>& data, unsigned int& pos)
	{
		unsigned int size_of_vec = 0;
		read_from_charArray(data, pos, size_of_vec);

		pos += size_of_vec;
	}

	void DataConversion::skip_in_charArray_string(const std::vector<char>& data, unsigned int& pos)
	{
		unsigned int size_of_string = 0;
		read_from_charArray(data, pos, size_of_string);

		pos += size_of_string;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<int>& list, const int element)
	{
		for (const auto entry : list)
		{
			if (entry == element)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<dmpc::NeighborPtr>& list, const int neighbor_id)
	{
		for (const auto& neighbor : list)
		{
			if (neighbor->get_id() == neighbor_id) 
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<dmpc::NeighborPtr>& list, const dmpc::NeighborPtr& element)
	{
		for (const auto& neighbor : list)
		{
			if (neighbor == element)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<std::shared_ptr<dmpc::CouplingInfo>>& list, const dmpc::CouplingInfo& element)
	{
		for (const auto& info : list)
		{
			if (*info == element) 
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector< dmpc::CouplingInfo >& list, const dmpc::CouplingInfo& element)
	{
		for (const auto& info : list)
		{
			if (info == element) 
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector< dmpc::CommunicationDataPtr >& list, const std::shared_ptr<dmpc::CommunicationData>& element)
	{
		for (const auto& comm_data : list)
		{
			if (comm_data->communication_info_->id_ == element->communication_info_->id_)
				return true;
		}
		return false;
	}

	const bool DataConversion::is_element_in_vector(const std::vector<dmpc::AgentInfo>& list, const int id)
	{
		for (const auto& info : list)
		{
			if (info.id_ == id) 
				return true;
		}
		return false;
	}

	void DataConversion::erase_element_from_vector(std::vector< dmpc::NeighborPtr >& list, const dmpc::NeighborPtr& element)
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

	void DataConversion::erase_element_from_vector(std::vector<dmpc::AgentPtr>& list, const dmpc::AgentInfo& element)
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

	void DataConversion::erase_element_from_vector(std::vector< dmpc::AgentInfo >& list, const dmpc::AgentInfo& element)
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

	void DataConversion::erase_element_from_vector(std::vector< dmpc::CouplingInfo >& list, const dmpc::CouplingInfo& element)
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

	void DataConversion::erase_element_from_vector(std::vector<std::shared_ptr<dmpc::CouplingInfo>>& list, const dmpc::CouplingInfo& element)
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

	void DataConversion::erase_element_from_vector(std::vector< dmpc::CommunicationDataPtr >& list, const  dmpc::CommunicationDataPtr& element)
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

	const dmpc::NeighborPtr DataConversion::get_element_from_vector(const std::vector<dmpc::NeighborPtr>& list, const int id)
	{
		for (const auto& neighbor : list)
		{
			if (neighbor->get_id() == id)
				return neighbor;
		}
		return nullptr;
	}