#pragma once

#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>
#include <random>

#include "QutBio.h"
#include "Homologs.hpp"
#include "OpenHausdorff.hpp"

using namespace QutBio;

namespace Utilities
{

	class PlotMixtureModel
	{
	private:
		struct Parameters
		{
			string gmmFile;
			string outFile;
			double xMin;
			double xMax;
			int steps = 100;

			friend ostream & operator<<( ostream& out, const Parameters & parms )
			{
				out << "Parameters:" << endl;
				out << "--gmmFile " << parms.gmmFile << endl;
				out << "--outFile " << parms.outFile << endl;
				out << "--xMin " << parms.xMin << endl;
				out << "--xMax " << parms.xMax << endl;
				out << "--steps " << parms.steps << endl;
				return out;
			}

			Parameters( Args & arguments )
			{
				if ( !arguments.Get( "gmmFile", gmmFile ) )
				{
					cerr << "Argument 'gmmFile' not defined." << endl;
					exit( 1 );
				}

				if ( !arguments.Get( "outFile", outFile ) )
				{
					cerr << "Argument 'outFile' not defined." << endl;
					exit( 1 );
				}

				if ( !arguments.Get( "xMin", xMin ) )
				{
					cerr << "Argument 'xMin' not valid." << endl;
					exit( 1 );
				}

				if ( !arguments.Get( "xMax", xMax ) )
				{
					cerr << "Argument 'xMax' not valid." << endl;
					exit( 1 );
				}

				if ( arguments.IsDefined( "steps" ) )
				{
					if ( !arguments.Get( "steps", steps ) )
					{
						cerr << "Argument 'steps' not valid." << endl;
						exit( 1 );
					}

					if ( steps < 1 )
					{
						cerr << "Argument 'steps' must be greater than zero." << endl;
						exit( 1 );
					}
				}

				if ( xMax < xMin )
				{
					cerr << "Argument 'xMax' must not be less than argument 'xMin'." << endl;
					exit( 1 );
				}
			}
		};

	public:

		/**
		Processes a ranking file produced by KmerRank to construct intra-class and
		inter-class empirical distributions of distances.
		The format of the input dataset is (laid out in comma-separated columns):

		query_id,ignored_class_id,subject_id,ignored_class_id,distance,rank,isHomolog

		*/
		static void Run( Args & args )
		{
			Parameters parms( args );
			cerr << parms;
			vector<GMM1D> models;
			ifstream gmmFile( parms.gmmFile );
			ReadMixtureModels( gmmFile, models );

			ofstream modelFile( parms.outFile );

			PlotCdf( modelFile, models, parms.xMin, parms.xMax, parms.steps );
			modelFile << "\n\n";
			PlotPdf( modelFile, models, parms.xMin, parms.xMax, parms.steps );
		}

		static void PlotCdf( ostream & out, vector<GMM1D> & model, double xMin, double xMax, int steps )
		{
			out << "AICc,x";

			for ( int i = 0; i <= steps; i++ )
			{
				double x = Interpolate( xMin, xMax, steps, i );
				out << "," << x;
			}

			out << "\n";

			for ( uint m = 0; m < model.size(); m++ )
			{
				out << model[m].AICc() << ",CDF[" << model[m].Size() << "]";

				for ( int i = 0; i <= steps; i++ )
				{
					double x = Interpolate( xMin, xMax, steps, i );
					out << "," << setprecision( 17 ) << model[m].Cdf( x );
				}

				out << "\n";
			}
		}

		static void PlotPdf( ostream & out, vector<GMM1D> & model, double xMin, double xMax, int steps )
		{
			out << "AICc,x";

			for ( int i = 0; i <= steps; i++ )
			{
				double x = Interpolate( xMin, xMax, steps, i );
				out << "," << x;
			}

			out << "\n";

			for ( uint m = 0; m < model.size(); m++ )
			{
				out << model[m].AICc() << ",PDF[" << model[m].Size() << "]";

				for ( int i = 0; i <= steps; i++ )
				{
					double x = Interpolate( xMin, xMax, steps, i );
					out << "," << setprecision( 17 ) << model[m].Pdf( x );
				}

				out << "\n";
			}
		}

		static int Interpolate( double xMin, double xMax, uint steps, uint i )
		{
			return xMin + i * ( xMax - xMin ) / steps;
		}

		static void ReadMixtureModels( istream & gmmFile, vector<GMM1D> & mixtures )
		{
			CsvReader reader( gmmFile );
			vector<vector<string>> records;
			reader.Read( records );
			const int N = (int) records.size();
			int i = 0;

			while ( i < N )
			{
				double aicc = 0;

				while ( i < N && records[i][0] != "alpha" )
				{
					if ( records[i][0].find( "AIC" ) == 0 )
					{
						vector<string> parts = String::Split( records[i][0], " \t" );
						aicc = Double::Parse( parts[1] );
					}

					i++;
				}

				if ( i == N ) break;

				i++;

				vector<double> a;
				vector<double> mu;
				vector<double> sigma;

				while ( i < N && records[i].size() == 3 && records[i][0].length() > 0 )
				{
					a.push_back( Double::Parse( records[i][0] ) );
					mu.push_back( Double::Parse( records[i][1] ) );
					sigma.push_back( Double::Parse( records[i][2] ) );
					i++;
				}

				mixtures.emplace_back( GMM1D( a, mu, sigma, aicc ) );
			}
		}

		static void Tabulate( ostream & out, GMM1D & model, double xMin, double xMax, uint steps )
		{
			out << "distance,CDF,PDF" << endl;

			double stepSize = ( xMax - xMin ) / steps;

			for ( size_t i = 0; i <= steps; i++ )
			{
				double x = xMin + i * stepSize;

				out << model.Cdf( x ) << "," << model.Pdf( x ) << endl;
			}
		}
	};
}
