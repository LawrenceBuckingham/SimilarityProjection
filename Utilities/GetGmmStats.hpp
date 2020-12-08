#pragma once

#include "Args.hpp"
#include "GMM1D.hpp"

using namespace QutBio;

namespace Utilities {
	class GetGmmStats {
	private:
		struct Parameters {
			string inGmmFile;
			string outFile;
			int index;
			vector<double> pValues;

			Parameters( Args &args) {
				if ( ! args.Get("inGmmFile", inGmmFile) ) {
					cerr << "Command line argument 'inGmmFile' is required." << endl;
					exit(1);
				}

				if ( !args.Get( "outFile", outFile ) ) {
					cerr << "Command line argument 'outFile' is required." << endl;
					exit(1);
				}

				if ( !args.Get("index", index) ) {
					cerr << "Command line argument 'index' is required." << endl;
					exit(1);
				}

				if ( !args.Get("pValues", pValues) ) {
					cerr << "Command line argument 'pValues' is required." << endl;
					exit(1);
				}
			}
		};
	public:
		static void Run(Args &args) {
			Parameters p( args );
			ifstream instream( p.inGmmFile );
			vector<GMM1D> models;
			GMM1D::Parse( instream, models );
			GMM1D & model = models[p.index - 1];

			ofstream out( p.outFile );
			out << "mean," << model.Mean() << endl;
			out << "std," << model.StdDev() << endl;

			for ( auto x: p.pValues ) {
				out << x << "," << model.InverseCdf( x ) << endl;
			}
		}
	};
}
