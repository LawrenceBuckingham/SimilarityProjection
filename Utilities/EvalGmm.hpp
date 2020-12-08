#pragma once

#include "Args.hpp"
#include "GMM1D.hpp"

using namespace QutBio;

namespace Utilities {
	class EvalGmm {
	private:
		struct Parameters {
			string inGmmFile;
			string outFile;
			int index;
			vector<double> xValues;

			Parameters( Args &args ) {
				if ( !args.Get( "inGmmFile", inGmmFile ) ) {
					cerr << "Command line argument 'inGmmFile' is required." << endl;
					exit( 1 );
				}

				if ( !args.Get( "outFile", outFile ) ) {
					cerr << "Command line argument 'outFile' is required." << endl;
					exit( 1 );
				}

				if ( !args.Get( "index", index ) ) {
					cerr << "Command line argument 'index' is required." << endl;
					exit( 1 );
				}

				if ( !args.Get( "xValues", xValues ) ) {
					cerr << "Command line argument 'xValues' is required." << endl;
					exit( 1 );
				}
			}
		};
	public:
		static void Run( Args &args ) {
			Parameters p( args );
			ifstream instream( p.inGmmFile );
			vector<GMM1D> models;
			GMM1D::Parse( instream, models );
			GMM1D & model = models[p.index - 1];

			ofstream out( p.outFile );

			for ( auto x : p.xValues ) {
				out << x << "," << model.Cdf( x ) << "," << model.Pdf( x ) << endl;
			}
		}
	};
}
