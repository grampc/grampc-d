#pragma once
#include "grampcd/util/types.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>

namespace grampcd {
	DMPC_CLASS_FORWARD(Logging)
}

namespace unit_test
{

	struct SimulationFile_Reference
	{
		std::string filename_ = "";
		std::string hash_code_ = "";

		template <class Archive>
		void serialize(Archive& ar)
		{
			ar(filename_, hash_code_);
		}
	};

	struct SimulationResult
	{
		std::string executable_name_ = "";
		std::vector<std::shared_ptr<SimulationFile_Reference>> simulationFile_References_;
		bool runtest_ = "false";

		template <class Archive>
		void serialize(Archive& ar)
		{
			ar(executable_name_, simulationFile_References_, runtest_);
		}
	};

	struct Results_of_UnitTests
	{
		std::vector<std::shared_ptr<SimulationResult>> simulationResults;

		template <class Archive>
		void serialize(Archive& ar)
		{
			ar(simulationResults);
		}
	};

	class ChecksumHandler {
	private:
		static const std::string PATH_REFERENCE_CHECKSUMS_;
		static const std::shared_ptr<std::string> read_file_to_string(const std::string& filename, grampcd::LoggingPtr log_);
	public:
		static const int REFERENCE_CHECKSUM_FILE_NOT_READABLE_;

		//generates checksum from file
		static const std::string generate_checksum(const std::string& filename, const grampcd::LoggingPtr log_);

		//obtains the reference data of all simulations
		static const std::shared_ptr<Results_of_UnitTests> get_results_of_all_unittests(const grampcd::LoggingPtr log_);

		//add a new unit test
		static void serialize_reference_data(std::shared_ptr<unit_test::Results_of_UnitTests> simulationResult, grampcd::LoggingPtr log_);
	};

}