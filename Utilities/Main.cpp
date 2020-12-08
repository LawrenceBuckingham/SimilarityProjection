#include "Main.hpp"

using namespace Utilities;

namespace QutBio {
mutex QutBio::FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
}

int main( int argc, char** argv ) {
	try {
		Main::Run( argc, argv );
	}
	catch ( Exception &ex ) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
	}
	catch ( runtime_error & err ) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
	}
	return 0;
}
