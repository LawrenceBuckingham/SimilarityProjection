// Complete rewrite of Timothy Chappell's trec_eval program.
//
//	This program reads a .homologs file (in place of .qrels)
//	and a "compact rankings" file as emitted by AAD2, KmerRank, and AAClustSig.
//	Then, rather than computing the interpolated precision out to 100%, it
//	computes the precision at a designated number of specified rankings.

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Array.hpp"

#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace QutBio;

	struct TrecEval {

		struct Ranking {
			double score;
			bool isRelevant;
			double recall;
			double precision;

			Ranking(
				double score,
				bool isRelevant
			) : score( score ), isRelevant( isRelevant ), recall( 0 ), precision( 0 ) {}

			friend bool operator<( const Ranking & lhs, const Ranking & rhs ) {
				return lhs.score < rhs.score;
			}
		};

		static int main( int argc, char **argv_ ) {
			if ( argc < 6 ) {
				fprintf( stderr, "Usage: %s homologsFile rankingFile summaryFile ignoreMissing=true|false rank1>=1 [rank2>=1 rank3>=1 ...]\n", argv_[0] );
				exit( 1 );
			}

			const char * qrelsFile = argv_[1];
			const char * rankingFile = argv_[2];
			const char * summaryFile = argv_[3];

			bool ignoreMissing = strcmp( argv_[4], "true" ) == 0;

			vector<uint> ranks;

			for ( int i = 5; i < argc; i++ ) {
				int rank;
				int numConverted = sscanf( argv_[i], "%d", &rank );

				if ( numConverted != 1 ) {
					fprintf( stderr, "Invalid numeric data: '%s'.\n", argv_[i] );
					exit( 1 );
				}

				ranks.push_back(rank);
			}

			fprintf( stderr, "Reading homologFile\n" );

			FILE *fp;
			fp = fopen( qrelsFile, "r" );

			if ( !fp ) {
				fprintf( stderr, "homologFile file %s does not exist.\n", qrelsFile );
				exit( 1 );
			}

			size_t overallRelevant = 0;

			unordered_map<string, size_t> topicIds;
			vector<string> topicNames;

			unordered_map<string, size_t> docIds;
			vector<string> docNames;

			vector<unordered_set<size_t>> qrels;
			vector<size_t> relevantDocumentCount;

			//size_t n = 0;

			for ( ;; ) {
				char topic_[256];
				char delimiter;

				int items = fscanf( fp, "%256s%c", topic_, &delimiter );

				if ( items < 2 ) {
					// fprintf(stderr, "%d items read. Finishing.\n", items);
					break;
				}

				string topic( topic_ );
				size_t topicId = GetTopicId( topic, topicIds, topicNames, qrels, relevantDocumentCount );

				//if ( n++ /* % 1000 == 0 */ ) {
				//	( cerr << "\r" << n << ": Current topic: " << topic << "                    " ).flush();
				//}

				while ( delimiter == ' ' ) {
					char doc_[256];

					items = fscanf( fp, "%256s%c", doc_, &delimiter );

					if ( items < 2 ) break;

					string doc( doc_ );
					size_t docId;

					auto docPos = docIds.find( doc );

					if ( docPos == docIds.end() ) {
						docId = docNames.size();
						docIds.emplace( doc, docId );
						docNames.push_back( doc );
					}
					else {
						docId = docPos->second;
					}

					qrels[topicId].insert( docId );
					relevantDocumentCount[topicId]++;
					overallRelevant++;
				}
			}

			fclose( fp );
			fprintf( stderr, "topicCount: %u\n", (unsigned) topicNames.size() );

			fprintf( stderr, "Reading rankings...\n" );
			FILE * rankingStream = fopen( rankingFile, "r" );

			if ( !rankingStream ) {
				fprintf( stderr, "Rankings file %s does not exist.\n", rankingFile );
				exit( 1 );
			}

			FILE *summaryStream = fopen( summaryFile, "w" "b" );

			if ( !summaryStream ) {
				fprintf( stderr, "Unable to open summary file '%s'\n", summaryFile );
				exit( 1 );
			}

			PrintSummaryHeadings( summaryStream, ranks );

			unordered_map<int,double> averagePrecision;

			for ( auto rank: ranks ) {
				averagePrecision[rank] = 0; 
			}

			size_t topicRetCount = 0;

			unordered_map<size_t, bool> retrievedResultsFor;
			size_t overallReturned = 0;
			size_t overallRelevantReturned = 0;
			string prevTopic = "";

			//n = 0;

			while ( !feof( rankingStream ) ) {
				char topic_[257];
				double score;
				char delimiter;

				int itemsParsed = fscanf( rankingStream, " %256s", topic_ );

				if ( itemsParsed == 0 ) break;

				string topic( topic_ );
				size_t topicId = GetTopicId( topic, topicIds, topicNames, qrels, relevantDocumentCount );

				retrievedResultsFor[topicId] = true;

				//if ( n++ /* % 1000 == 0 */ ) {
				//	( cerr << "\r" << n << ": Current rankings: " << topic_ << "                    " ).flush();
				//}

				vector<Ranking> rankings;
				rankings.reserve( 10000 );

				size_t numReturned = 0;
				size_t numRelevantReturned = 0;

				while ( !feof( rankingStream ) ) {
					char doc_[257];

					itemsParsed = fscanf( rankingStream, " %256s %lf%c", doc_, &score, &delimiter );

					if ( itemsParsed < 3 ) break;

					// cerr << "topic(" << topic_ << ") doc(" << doc_ << ") score(" << score << ") delimiter(" << int(delimiter) << ")\n";


					string doc( doc_ );
					size_t docId;

					auto docPos = docIds.find( doc );

					if ( docPos == docIds.end() ) {
						docId = docNames.size();
						docNames.push_back( doc );
						docIds.emplace( doc, docId );
					}
					else {
						docId = docPos->second;
					}

					bool relevant = qrels[topicId].find( docId ) != qrels[topicId].end();

					rankings.emplace_back( -score, relevant );
					retrievedResultsFor[topicId] = true;

					// Update returned document count for individual query and total.
					numReturned++;
					overallReturned++;

					if ( relevant ) {
						// Update relevant document count for individual query and total.
						numRelevantReturned++;
						overallRelevantReturned++;
					}

					if ( delimiter != ' ' ) {
						if ( delimiter != '\n' && delimiter != '\r' ) {
							( cerr << "\nUnexpected delimiter: ASCII(" << int( delimiter ) << ")\n" ).flush();
						}
						break;
					}
				}

				if ( topic != prevTopic ) {
					ProcessTopic(
						rankings,
						relevantDocumentCount[topicId],
						ranks,
						averagePrecision
					);

					PrintTopic(
						topicNames[topicId].c_str(),
						relevantDocumentCount[topicId],
						numReturned,
						numRelevantReturned,
						ranks,
						rankings,
						summaryStream
					);

					topicRetCount++;
				}

				prevTopic = topic;
			}

			// Emit results for topics where no records returned (if we are not ignoring them).
			if ( !ignoreMissing ) {
				//(cerr << "\nProcessing missing topics\n").flush();

				vector<Ranking> empty;

				for ( size_t topicId = 0; topicId < topicNames.size(); topicId++ ) {
					if ( !retrievedResultsFor[topicId] ) {
						//( cerr << "\rMissing topic: " << topicNames[topicId] << "                    " ).flush();

						PrintTopic(
							topicNames[topicId].c_str(),
							relevantDocumentCount[topicId],
							0,
							0,
							ranks,
							empty,
							summaryStream
						);

						topicRetCount++;
					}
				}
			}

			for ( auto j: ranks ) {
				averagePrecision[j] /= topicRetCount;
			}

			// Emit results for total.
			fprintf( summaryStream, "%s\t%zu\t%zu\t%zu",
				"Overall", overallRelevant, overallRelevantReturned, overallReturned );

			for ( auto j: ranks ) {
				fprintf( summaryStream, "\t%.4f", averagePrecision[j] );
			}

			fprintf( summaryStream, "\n" );

			fclose( summaryStream );

			fprintf( stderr, "Finished.\n" );
			return 0;
		}

		static size_t GetTopicId(
			const string & topic,
			unordered_map<string, size_t> & topicIds,
			vector<string> & topicNames,
			vector<unordered_set<size_t>> & qrels,
			vector<size_t> & relevantDocumentCount
		) {
			auto topicPos = topicIds.find( topic );
			size_t topicId;

			if ( topicPos == topicIds.end() ) {
				topicId = topicNames.size();
				topicNames.push_back( topic );
				topicIds.emplace( topic, topicId );
				relevantDocumentCount.push_back( 0 );
				qrels.resize( qrels.size() + 1 );
			}
			else {
				topicId = topicPos->second;
			}

			return topicId;
		}

		static void PrintTopic(
			const char * topicName,
			size_t relevantDocumentCount,
			size_t returnedDocumentCount,
			size_t returnedRelevantDocumentCount,
			const vector<uint> & ranks,
			const vector<Ranking> & rankings,
			FILE * f
		) {
			fprintf( f, "%s\t%zu\t%zu\t%zu",
				topicName, relevantDocumentCount, returnedRelevantDocumentCount, returnedDocumentCount );

			for ( auto i: ranks ) {
				fprintf( f, "\t%.4f", i <= rankings.size() ? rankings[i-1].precision : 0 );
			}

			fprintf( f, "\n" );
		}

		static void ProcessTopic(
			vector<Ranking> & rankings,
			size_t relevantDocumentCount,
			vector<uint> & ranks,
			unordered_map<int,double> & averagePrecision
		) {
			sort( begin( rankings ), end( rankings ) );

			size_t relFound = 0;
			size_t totalFound = 0;

			for ( auto & t : rankings ) {
				totalFound++;

				if ( t.isRelevant ) relFound++;

				t.recall = double( relFound ) / relevantDocumentCount;
				t.precision = double( relFound ) / totalFound;
			}

			for ( auto rank: ranks ) {
				averagePrecision[rank] += rank <= rankings.size() ? rankings[rank-1].precision : 0;
			}
		}

		static void PrintSummaryHeadings( FILE * f, vector<uint> ranks ) {
			fprintf( f, "Topic\tRelevant\tRelevant Returned\tTotal Returned\tAverage Precision" );

			for ( auto & rank: ranks ) {
				fprintf( f, "\t%d", rank );
			}

			fprintf( f, "\n" );
		}
	};


int main( int argc, char ** argv ) {
	return TrecEval::main( argc, argv );
}
