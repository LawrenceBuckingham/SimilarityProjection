#pragma once

#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>
#include <string>

#include "QutBio.h"
#include "DataLoader.hpp"
#include "Homologs.hpp"
#include "SequenceWrapper.hpp"
#include "Substring.hpp"

using namespace QutBio;

namespace Utilities {

	class GetKmerDistanceDistribution {
	private:
		struct Parameters {
			bool ok = true;

			Alphabet * alphabet;
			SimilarityMatrix * matrix;

			string dbFile;

			string overallAicFile;
			string overallDistributionFile;

			string posAicFile;
			string posDistributionFile;

			string negAicFile;
			string negDistributionFile;

			size_t idIndex = 0;
			int classIndex = -1;
			size_t kmerLength = 0;

			time_t seed = time( 0 );

			size_t sampleSize = 10000;
			size_t maxModelSize = 10;

			string homologFile = "";

			Parameters( Args & arguments ) {
				auto & p = *this;

#define REQ(x,h) arguments.Required(x, #x, h)
#define OPT(x,h) arguments.Optional(x, #x, h)
#define STR(x) STR2(x)
#define STR2(x) #x
				arguments.SetTitle(
					"Monte-Carlo approximation of uniform random pairwise k-mer distance \n"
					"distributions. Separate empirical CDF functions are calculated for sampled \n"
					"k-mers originating in 'related' sequences and 'unrelated' sequences, along \n"
					"with a combined distribution. A Gaussian mixture model is fitted to each ECDF \n"
					"for visual assessment of the distributions, and possible downstream use in \n"
					"assessing statistical significance. Mixture models are fitted via Akaike \n"
					"Information Criterion (AIC).\n"
					"Two sequences are considered to be 'related' if their lists of class labels \n"
					"have a non-empty intersection. Otherwise, they are unrelated." );

				REQ( dbFile,
					"The file path of a FASTA formatted sequence collection." );

				REQ( overallAicFile,
					"The file path of a document that will be populated with the AIC and \n"
					"parameters of a series of Gaussian mixture models fitted to the complete sample." );

				REQ( overallDistributionFile,
					"The file path of a document that will be populated with the ECDF and fitted \n"
					"PDF and CDF of each Gaussian mixture model obtained from the complete sample." );

				REQ( posAicFile,
					"The file path of a document that will be populated with the AIC and \n"
					"parameters of Gaussian mixture models fitted to the distances between \n"
					"k-mers from 'related' sequences." );

				REQ( posDistributionFile,
					"The file path of a document that will be populated with the ECDF and fitted \n"
					"PDF and CDF of each Gaussian mixture model obtained from the 'related' pairs \n"
					"in the sample." );

				REQ( negAicFile,
					"The file path of a document that will be populated with the AIC and \n"
					"parameters of Gaussian mixture models fitted to the distances between \n"
					"k-mers originating in 'unrelated' sequences." );

				REQ( negDistributionFile,
					"The file path of a document that will be populated with the ECDF and fitted \n"
					"PDF and CDF of each Gaussian mixture model obtained from the 'unrelated' pairs \n"
					"in the sample." );

				OPT( idIndex,
					"The zero-origin index of the sequence Id in the pipe-separated FASTA definition line." );

				OPT( classIndex,
					"The zero-origin index of the list of semi-colon separated class labels in the pipe-\n"
					"separated FASTA definition line." );

				REQ( kmerLength, "The word length for k-mer tiling." );

				if ( kmerLength == 0 ) {
					arguments.Fail();
					cerr << "Argument 'kmerLength' must be greater than 0." << endl;
				}

				OPT( sampleSize,
					"The number of random (k1,k2) pairs used to compute the empirical distributions." );

				OPT( homologFile,
					"The name of a homolog file which defines relevant pairs of documents. " );

				if ( !(p.sampleSize > 0) ) {
					cerr << "Argument 'sampleSize' must be greater than zero." << endl;
					arguments.Fail();
				}

				OPT( maxModelSize, "The maximum number of kernels to use in Gaussian mixture models." );

				if ( !(maxModelSize > 0) ) {
					cerr << "Argument 'maxModelSize' must be greater than zero." << endl;
					arguments.Fail();
				}

				arguments.Required( alphabet, matrix );

				OPT( seed, "Seed for random number generator." );

				if ( !arguments.Ok() ) {
					arguments.Help();
					throw Exception( "Invalid arguments supplied.", FileAndLine );
				}

			}
		};
	public:

		static void Run( Args & args ) {
			Parameters p( args );
			UniformRealRandom rand( (int) p.seed );
			SimilarityMatrix * matrix = p.matrix;

			auto dbSeqs = Load::Fasta( p.dbFile, p.idIndex, matrix->Alphabet() );
			cerr << args.ProgName() << ": " << dbSeqs.size() << " sequences loaded.\n";

			struct Seq : public SequenceWrapper {
				Seq( FastaSequence * base ) : SequenceWrapper( base ) {}

				vector<Seq *> homologs;

				bool IsHomolog( Seq * other ) {
					return find( homologs.begin(), homologs.end(), other ) != homologs.end();
				}
			};

			vector<Seq *> db;
			SequenceWrapper::Wrap( dbSeqs, db );

			Index<Seq> idx( db );

			if ( p.homologFile.size() > 0 ) {
				if ( p.homologFile.find_last_of( "qrels" ) == string::npos ) {
					ifstream inStream( p.homologFile );
					Homologs::Parse( inStream, idx, idx );
				}
				else {
					FILE * qrelsFile = fopen( p.homologFile.c_str(), "r" );

					if ( qrelsFile ) {
						Homologs::ParseQrels( qrelsFile, idx, idx );
						fclose( qrelsFile );
					}
				}
			}

			auto addKmerCount = [&p]( long val, Seq * seq ) {
				return val + (int) seq->KmerCount( p.kmerLength );
			};

			size_t totalKmerCount = accumulate( db.begin(), db.end(), 0, addKmerCount );

			Selector selector( rand, p.sampleSize * 2, totalKmerCount );

			vector<Distance> distances, positiveDistances, negativeDistances;

			using P = Pair<Substring *, Seq *>;
			vector<P> sample;

			for ( auto seq : db ) {
				size_t n = seq->Sequence().size() + 1 - p.kmerLength;

				for ( size_t i = 0; i < n; i++ ) {
					if ( selector.SelectThis() ) {
						auto kmer = new Substring( seq->Sequence().data(), i, p.kmerLength );
						P pair( kmer, seq );
						sample.push_back( pair );
					}
				}
			}

			for ( size_t i = 0; i < sample.size(); i++ ) {
				size_t ix1 = size_t( rand() * sample.size() );
				size_t ix2 = size_t( rand() * sample.size() );
				auto pair = sample[ix1];
				sample[ix1] = sample[ix2];
				sample[ix2] = pair;
			}

			for ( size_t i = 0; i < p.sampleSize; i++ ) {
				auto & px = sample[i * 2];
				auto & py = sample[11 * 2 + 1];
				Substring & kx = *(px.item1);
				Substring & ky = *(py.item1);
				auto x = kx.Chars();
				auto y = ky.Chars();
				Distance d = matrix->Difference( x, y, p.kmerLength );
				distances.push_back( d );

				if ( px.item2->IsHomolog( py.item2 ) ) {
					positiveDistances.push_back( d );
				}
				else {
					negativeDistances.push_back( d );
				}
			}

			if ( p.sampleSize != distances.size() ) {
				throw Exception( "Sample size and number of distances computed do not match.", FileAndLine );
			}

			if ( p.sampleSize != positiveDistances.size() + negativeDistances.size() ) {
				throw Exception( "Combined number of distances in positive and negative samples does not equal sample size", FileAndLine );
			}

			if ( positiveDistances.size() == 0 ) {
				throw Exception( "Positive sample size is zero: try a larger sample", FileAndLine );
			}

			if ( negativeDistances.size() == 0 ) {
				throw Exception( "Negative sample size is zero: try a larger sample", FileAndLine );
			}

			//DumpSample( distances, p.overallSampleFile );
			//DumpSample( positiveDistances, p.posSampleFile );
			//DumpSample( negativeDistances, p.negSampleFile );

			ProcessDistribution( p.maxModelSize, p.overallAicFile, p.overallDistributionFile, distances );
			ProcessDistribution( p.maxModelSize, p.posAicFile, p.posDistributionFile, positiveDistances );
			ProcessDistribution( p.maxModelSize, p.negAicFile, p.negDistributionFile, negativeDistances );
		}

		static void DumpSample( vector<Distance> &sample, string & fileName ) {
			ofstream f( fileName );
			bool deja = false;

			for ( auto x : sample ) {
				if ( deja ) f << "\t"; else deja = true;

				f << x;
			}

			f << "\n";
		}

		typedef struct {
			GMM1D model;
			double AICc;
		} FittedModel;

		static void EmitEmpiricalDistribution( ostream & out, const vector<Distance> & distances ) {
			out << "distance";

			for ( size_t i = 0; i < distances.size(); i++ ) {
				if ( i == 0 || distances[i] != distances[i - 1] ) {
					auto dist = distances[i];
					out << "\t" << dist;
				}
			}

			out << "\nECDF";

			for ( size_t i = 0; i < distances.size(); i++ ) {
				if ( i == 0 || distances[i] != distances[i - 1] ) {
					auto dist = distances[i];
					out << "\t" << ((double) (i + 1) / distances.size());
				}
			}

			out << "\n";
		}

		static void EmitDistribution( ostream & out, const vector<Distance> & distances, FittedModel & mixture ) {
			out << "mCDF[" << mixture.model.Size() << "]";

			for ( size_t i = 0; i < distances.size(); i++ ) {
				if ( i == 0 || distances[i] != distances[i - 1] ) {
					auto dist = distances[i];
					out << "\t" << mixture.model.Cdf( dist );
				}
			}

			out << "\n";
			out << "mPDF[" << mixture.model.Size() << "]";

			for ( size_t i = 0; i < distances.size(); i++ ) {
				if ( i == 0 || distances[i] != distances[i - 1] ) {
					auto dist = distances[i];
					out << "\t" << mixture.model.Pdf( dist );
				}
			}
			out << "\n";
		}

		static FittedModel FitMixtureModel( size_t modelSize, const vector<Distance> & distances ) {
			FittedModel m{ GMM1D( modelSize ), 0 };
			m.model.Initialise( distances );
			m.model.Train( distances, 1000, 1e-10, false );
			m.AICc = m.model.AICc( distances );
			cerr << "Gaussian mixture model with " << modelSize << " components: AICc = " << m.AICc << endl;
			return m;
		}

		static void EmitMixtureModel( ofstream & out, FittedModel & mixture ) {
			out << "Model\nAIC\t" << mixture.AICc << "\n";
			out << mixture.model << "\nEndModel\n";
		}

		static void ProcessDistribution(
			size_t maxModelSize,
			const string & aicFileName,
			const string & distributionFileName,
			vector<Distance> & distances
			//
		) {
			ofstream aicFile( aicFileName );
			ofstream distributionFile( distributionFileName );

			vector<Distance> sortedDistances = distances;
			sort( sortedDistances.begin(), sortedDistances.end() );
			EmitEmpiricalDistribution( distributionFile, sortedDistances );

			for ( uint i = 1; i <= maxModelSize; i++ ) {
				FittedModel mixture = FitMixtureModel( i, distances );
				EmitMixtureModel( aicFile, mixture );
				EmitDistribution( distributionFile, sortedDistances, mixture );
			}
		}
	};
}
