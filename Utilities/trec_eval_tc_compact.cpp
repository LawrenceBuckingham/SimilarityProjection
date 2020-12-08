// Complete rewrite of Timothy Chappell's trec_eval program.
//
//	This program reads a .homologs file (in place of .qrels)
//	and a "compact rankings" file as emitted by COV-Jacc.
//
//	Not intended to be compatible with trec_eval; however, the
//	interpolated precision/recall curves generated are numerically 
//	equal to hose calculated by Tim's program.

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Array.hpp"

using namespace std;

namespace QutBio {
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
			if ( argc < 5 ) {
				fprintf( stderr, "Usage: %s homologsFile rankingFile summaryFile ignoreMissing=true|false [interpolated_precision_points=11]\n", argv_[0] );
				exit( 1 );
			}

			const char * qrelsFile = argv_[1];
			const char * rankingFile = argv_[2];
			const char * summaryFile = argv_[3];

			bool ignoreMissing = strcmp( argv_[4], "true" ) == 0;

			uint interpolationPoints = 11;

			if ( argc >= 6 ) {
				int numConverted = sscanf( argv_[5], "%d", &interpolationPoints );

				if ( numConverted != 1 ) {
					fprintf( stderr, "Homologs file %s does not exist.\n", qrelsFile );
					exit( 1 );
				}
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

			PrintSummaryHeadings( summaryStream, interpolationPoints );

			vector<double> averageIprec( interpolationPoints, 0.0 );
			size_t topicRetCount = 0;
			double meanAveragePrecision = 0.0;

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
					ProcessAndPrint( 
						interpolationPoints, 
						rankings, 
						relevantDocumentCount, 
						topicId, 
						averageIprec, 
						topicNames, 
						numReturned, 
						numRelevantReturned, 
						summaryStream, 
						meanAveragePrecision, 
						topicRetCount 
					);
				}

				prevTopic = topic;
			}

			// Emit results for topics where no records returned (if we are not ignoring them).
			if ( !ignoreMissing ) {
				//(cerr << "\nProcessing missing topics\n").flush();

				for ( size_t topicId = 0; topicId < topicNames.size(); topicId++ ) {
					if ( !retrievedResultsFor[topicId] ) {
						//( cerr << "\rMissing topic: " << topicNames[topicId] << "                    " ).flush();

						vector<double> emptyGrid( interpolationPoints, 0.0 );

						PrintTopic(
							topicNames[topicId].c_str(),
							relevantDocumentCount[topicId],
							0,
							0,
							0,
							emptyGrid,
							summaryStream
						);

						topicRetCount++;
					}
				}
			}

			meanAveragePrecision /= topicRetCount;

			for ( size_t j = 0; j < interpolationPoints; j++ ) {
				averageIprec[j] /= topicRetCount;
			}

			// Emit results for total.
			fprintf( summaryStream, "%s\t%zu\t%zu\t%zu\t%0.4f",
				"Overall", overallRelevant, overallRelevantReturned, overallReturned, meanAveragePrecision );

			for ( size_t j = 0; j < interpolationPoints; j++ ) {
				fprintf( summaryStream, "\t%.4f", averageIprec[j] );
			}

			fprintf( summaryStream, "\n" );

			fclose( summaryStream );

			fprintf( stderr, "Finished.\n" );
			return 0;
		}

		static void ProcessAndPrint( const QutBio::uint &interpolationPoints, std::vector<QutBio::TrecEval::Ranking> &rankings, std::vector<size_t> &relevantDocumentCount, const size_t &topicId, std::vector<double> &averageIprec, std::vector<string> &topicNames, const size_t &numReturned, const size_t &numRelevantReturned, FILE * summaryStream, double &meanAveragePrecision, size_t &topicRetCount ) {
			double averagePrecision;
			vector<double> interpolatedGrid( interpolationPoints, 0.0 );

			ProcessTopic(
				rankings,
				relevantDocumentCount[topicId],
				averagePrecision,
				averageIprec,
				interpolatedGrid
			);

			PrintTopic(
				topicNames[topicId].c_str(),
				relevantDocumentCount[topicId],
				numReturned,
				numRelevantReturned,
				averagePrecision,
				interpolatedGrid,
				summaryStream
			);

			meanAveragePrecision += averagePrecision;
			topicRetCount++;
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
			double averagePrecision,
			const vector<double> & interpolatedGrid,
			FILE * f
		) {
			fprintf( f, "%s\t%zu\t%zu\t%zu\t%0.4f",
				topicName, relevantDocumentCount, returnedRelevantDocumentCount, returnedDocumentCount, averagePrecision );

			size_t interpolationPoints = interpolatedGrid.size();

			for ( size_t i = 0; i < interpolationPoints; i++ ) {
				fprintf( f, "\t%.4f", interpolatedGrid[i] );
			}

			fprintf( f, "\n" );
		}

		static void ProcessTopic(
			vector<Ranking> & rankings,
			size_t relevantDocumentCount,
			double & averagePrecision,
			vector<double> & averageIprec,
			vector<double> & interpolatedGrid
		) {
			sort( begin( rankings ), end( rankings ) );

			size_t relFound = 0;
			size_t totalFound = 0;

			averagePrecision = 0.0;

			for ( auto & t : rankings ) {
				totalFound++;

				if ( t.isRelevant ) relFound++;

				t.recall = double( relFound ) / relevantDocumentCount;
				t.precision = double( relFound ) / totalFound;

				if ( t.isRelevant ) {
					averagePrecision += t.precision;
				}
			}

			if ( relevantDocumentCount > 0 ) {
				averagePrecision /= relevantDocumentCount;
			}

			double zeroAtRecall = numeric_limits<double>::max();
			size_t numRankings = rankings.size();

			// "Interpolate" the observed precision.
			for ( int i = int( numRankings ) - 1; i > 0; i-- ) {
				double precision = rankings[i].precision;

				if ( rankings[i - 1].precision < precision ) {
					rankings[i - 1].precision = precision;
				}

				if ( precision == 0 ) zeroAtRecall = rankings[i].recall;
			}

			size_t currentRanking = 0;
			size_t interpolationPoints = interpolatedGrid.size();

			// "Interpolate" the synthesised precision at a set of evenly spaced points.
			for ( size_t j = 0; j < interpolationPoints; j++ ) {
				double recall = double( j ) / ( interpolationPoints - 1 );
				double precision = 0;

				if ( currentRanking < numRankings && recall < zeroAtRecall ) {
					while ( currentRanking < numRankings && rankings[currentRanking].recall < recall ) {
						currentRanking++;
					}

					if ( currentRanking < numRankings ) {
						precision = rankings[currentRanking].precision;
					}
				}

				interpolatedGrid[j] = precision;
				averageIprec[j] += precision;
			}
		}

		static void PrintSummaryHeadings( FILE * f, size_t interpolationPoints ) {
			fprintf( f, "Topic\tRelevant\tRelevant Returned\tTotal Returned\tAverage Precision" );

			for ( size_t j = 0; j < interpolationPoints; j++ ) {
				double recall = 1.0 * j / ( interpolationPoints - 1 );
				fprintf( f, "\t%.2f", recall );
			}

			fprintf( f, "\n" );
		}
	};
}

int main( int argc, char ** argv ) {
	return QutBio::TrecEval::main( argc, argv );
}
