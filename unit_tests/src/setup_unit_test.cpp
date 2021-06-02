
#include "dmpc/util/logging.hpp"

#include "unit_tests/checksum_handler.hpp"
#include "unit_tests/file_handler.hpp"

int main(int argc, char* argv[]) 
{
	unit_test::ChecksumHandler checksumHandler(LoggingPtr(new Logging()));
	checksumHandler.add_unit_test("coupled_watertanks");
	checksumHandler.add_unit_test("evaluate_neighborApproximation");
	checksumHandler.add_unit_test("high_scaled_system");
	checksumHandler.add_unit_test("plug-and-play");
}