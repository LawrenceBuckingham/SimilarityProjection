#include "KmerRank.hpp"

using namespace KmerQuery;

int main( int argc, char** argv ) {
	try {
		Main::Run( argc, argv );
	}
	catch ( Exception &ex ) {
		(cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl).flush();
	}
	catch ( runtime_error &err ) {
		(cerr << "Unhandled exception:" << endl << err.what() << endl).flush();
	}
	return 0;
}

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
