#pragma once

#include <FastaSequence.hpp>
#include <KmerClusterPrototype.hpp>
#include <KmerCluster.hpp>

namespace QutBio {
	template<typename DistanceFunction, typename KmerType>
	class KMedoids {
	public:
		using Cluster = KmerCluster<DistanceFunction, KmerType>;

		enum class SortMode {
			SortRandom = 1,
			SortLongestFirst = 2,
			SortShortestFirst = 3
		};

		enum class SelectMode {
			SelectGreedy = 1,
			SelectNearest = 2
		};

		enum class MedoidMode {
			MedoidMin = 1,

			MedoidBruteForce = MedoidMin,
			MedoidMeddit = MedoidBruteForce + 1,
			MedoidNone = MedoidMeddit + 1,

			MedoidMax = MedoidNone
		};

		static void Partition(
			vector<Subsequence> & seqs,
			vector<Kmer *> & clusterProtos,
			vector<Cluster *> & clusters,
			size_t kmerLength,
			Distance threshold,
			int randSeed,
			Alphabet &alphabet,
			DistanceFunction &distance,
			size_t trials = 40,
			size_t iterations = 3,
			SortMode sortMode = SortMode::SortRandom,
			SelectMode selectMode = SelectMode::SelectNearest,
			MedoidMode medoidMode = MedoidMode::MedoidMeddit,
			size_t minMedditSize = 1000
			//
		) {
			KmerIndex kmerIndex( seqs, kmerLength );
			vector<Kmer *> & kmers = kmerIndex.GetKmers();
			const uint N = kmers.size();
			vector<EncodedKmer> kmerCodes( N );
			vector<unsigned long> allocatedDist( N );
			vector<uint> allocatedCount( N );
			vector<vector<uint>> bestKmersPerCluster;
			vector<Kmer *> bestProtos;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
			Subsequence * bestInitialSeq = 0;
#pragma GCC diagnostic pop

			UniformIntRandom<uint> rand( randSeed, 0, seqs.size() - 1 );

			if ( seqs.size() > 1 ) {
				if ( sortMode == SortMode::SortRandom ) {
					for ( uint i = 0; i < seqs.size(); i++ ) {
						uint pos = rand();
						auto t = seqs[i];
						seqs[i] = seqs[pos];
						seqs[pos] = t;
					}
				}
				else if ( sortMode == SortMode::SortLongestFirst ) {
					sort( seqs.begin(), seqs.end(), []( const Subsequence &lhs, const Subsequence &rhs ) {
						return lhs.source->Length() > rhs.source->Length();
					} );
				}
				else {
					sort( seqs.begin(), seqs.end(), []( const Subsequence &lhs, const Subsequence &rhs ) {
						return lhs.source->Length() < rhs.source->Length();
					} );
				}
			}

			uint bestAssignedKmers = 0;

			for ( size_t trial = 0; trial < std::min( seqs.size(), trials); trial++ ) {
				// the shuffle/sort takes care of randomisation if required.
				uint initialSeqPos = trial;

				auto & initialSeq = seqs[initialSeqPos];
				KmerIndex initialKmerIdx( &initialSeq, 1, kmerLength );
				vector<Kmer *> & initialKmers = initialKmerIdx.GetKmers();
				const uint K = initialKmers.size();
				vector<Kmer *> protos( K );
				vector<EncodedKmer> protoCodes( K );
				vector<vector<uint>> kmersPerCluster( K );

				// Use these to estimate sigma 
				vector<ulong> dSumPerCluster( K );
				vector<ulong> dSumSqPerCluster( K );

				for ( uint n = 0; n < N; n++ ) {
					kmerCodes[n] = kmers[n]->PackedEncoding();
				}

				for ( uint k = 0; k < K; k++ ) {
					protos[k] = initialKmers[k];
					protoCodes[k] = protos[k]->PackedEncoding();
					kmersPerCluster[k].reserve( ( N + K - 1 ) / K );
				}

				uint numAssignedKmers;

				for ( uint iter = 0; iter < iterations; iter++ ) {
					numAssignedKmers = 0;

					for ( uint k = 0; k < K; k++ ) {
						kmersPerCluster[k].clear();
						dSumPerCluster[k] = 0;
						dSumSqPerCluster[k] = 0;
					}

					for ( uint n = 0; n < N; n++ ) {
						EncodedKmer currentCode = kmerCodes[n];
						uint nearestProtoPos = numeric_limits<uint>::max();
						Distance smallestDistance = numeric_limits<Distance>::max();

						for ( uint k = 0; k < K; k++ ) {
							if ( !protoCodes[k] ) continue;

							Distance currentDist = distance( currentCode, protoCodes[k], kmerLength );

							if ( selectMode == SelectMode::SelectNearest ) {
								if ( currentDist < smallestDistance ) {
									nearestProtoPos = k;
									smallestDistance = currentDist;
								}
							}
							else {
								// SelectGreedy
								if ( currentDist <= threshold ) {
									nearestProtoPos = k;
									smallestDistance = currentDist;
									break;
								}
							}
						}

						if ( smallestDistance <= threshold ) {
							kmersPerCluster[nearestProtoPos].push_back( n );
							numAssignedKmers += kmers[n]->Instances().size();
							dSumPerCluster[nearestProtoPos] += smallestDistance;
							dSumSqPerCluster[nearestProtoPos] += smallestDistance * smallestDistance;
						}
					}

					if ( medoidMode == MedoidMode::MedoidNone ) {
						for ( uint k = 0; k < K; k++ ) {
							protos[k] = kmerIndex.at( protos[k]->Substr() );
							protoCodes[k] = protos[k]->PackedEncoding();
						}
					}
					else {
						for ( uint k = 0; k < K; k++ ) {
							vector<uint> &clusterAssignment = kmersPerCluster[k];
							Kmer *&protoPtr = protos[k];
							EncodedKmer &protoCode = protoCodes[k];

							//cerr << "clusterAssignment.size() = " << clusterAssignment.size() << "\n";
							//cerr << "parms.minMedditSize = " << parms.minMedditSize << "\n";

							if ( medoidMode == MedoidMode::MedoidBruteForce || clusterAssignment.size() <= minMedditSize ) {
								// Cluster is small or we're not using MEDDIT anyway...
								GetMedoid( clusterAssignment, allocatedDist, allocatedCount, distance, kmerCodes, kmerLength, kmers, protoPtr, protoCode );
								//cerr << "Doing Exact Medoid\n";
							}
							else {
								int n = kmersPerCluster[k].size();
								double mu = double( dSumPerCluster[k] ) / n;
								double sigma = sqrt( double( dSumSqPerCluster[k] ) / n - mu * mu );

								// Use MEDDIT.
								GetMedoid_MEDDIT(
									clusterAssignment,
									allocatedDist,
									allocatedCount,
									distance,
									kmerCodes,
									kmerLength,
									kmers,
									rand,
									sigma,
									protoPtr,
									protoCode
								);
							}
						}

					}
				}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
				if ( numAssignedKmers > bestAssignedKmers ) {
					bestProtos = protos;
					bestKmersPerCluster = kmersPerCluster;
					bestAssignedKmers = numAssignedKmers;
					bestInitialSeq = &seqs[initialSeqPos];
				}

				//cerr << "Trial " << trial << ": " << numAssignedKmers << " kmers assigned to "
				//	<< kmersPerCluster.size()
				//	<< " clusters; best so far: " << bestAssignedKmers << "\n";
#pragma GCC diagnostic pop
			}

			//cerr << "Final kmers assigned to "
			//	<< bestKmersPerCluster.size()
			//	<< " clusters: " << bestAssignedKmers << "\n";
			//cerr << "Seed domain: " << bestInitialSeq->source->Id() << "\n";

			for ( uint k = 0; k < bestProtos.size(); k++ ) {
				Kmer *p = bestProtos[k];

				if ( !p ) continue;

				CreateCluster( p, alphabet, kmerLength, distance, kmers, bestKmersPerCluster[k], clusterProtos, clusters );
			}
		}

		static void GetMedoid(
			vector<uint_least32_t> & clusterAssignment,
			vector<unsigned long> &allocatedDist,
			vector<uint_least32_t> &allocatedCount,
			DistanceFunction &distance,
			vector<QutBio::EncodedKmer> &kmerCodes,
			uint kmerLength,
			vector<QutBio::Kmer *> & kmers,
			Kmer *& protoPtr,
			EncodedKmer & protoCode //
		) {
			uint N_k = clusterAssignment.size();

			if ( N_k > 0 ) {
				for ( uint i = 0; i < N_k; i++ ) {
					uint candidatePos = clusterAssignment[i];
					allocatedDist[candidatePos] = 0;
					allocatedCount[candidatePos] = 0; // Place holder for counters needed by MEDIT.

					for ( uint j = 0; j < N_k; j++ ) {
						// if (j == i) continue;

						uint comparePos = clusterAssignment[j];
						allocatedDist[candidatePos] += distance( kmerCodes[candidatePos], kmerCodes[comparePos], kmerLength );
						allocatedCount[candidatePos] ++;
					}
				}

				uint minimalDistPos = clusterAssignment[0];

				for ( uint i = 0; i < N_k; i++ ) {
					uint candidatePos = clusterAssignment[i];
					if ( allocatedDist[candidatePos] < allocatedDist[minimalDistPos] ) {
						minimalDistPos = candidatePos;
					}
				}

				protoPtr = kmers[minimalDistPos];
				protoCode = kmerCodes[minimalDistPos];
			}
			else {
				protoPtr = nullptr;
				protoCode = nullptr;
			}
		}

		/**
		*	Implement MEDDIT algorithm of Bagaria et al,
		*	I got the formula for the confidence interval
		*	at https://ocw.mit.edu/courses/mathematics/18-s997-high-dimensional-statistics-spring-2015/lecture-notes/MIT18_S997S15_Chapter1.pdf
		*/

		static void GetMedoid_MEDDIT(
			vector<uint> & clusterAssignment,
			vector<ulong> &distSum,
			vector<uint> &distCount,
			DistanceFunction &distance,
			vector<QutBio::EncodedKmer> &kmerCodes,
			uint kmerLength,
			vector<QutBio::Kmer *> & kmers,
			UniformIntRandom<uint> &rand,
			double sigma,
			Kmer *& protoPtr,
			EncodedKmer & protoCode //
		) {
			uint N_k = clusterAssignment.size();

			if ( N_k > 1 ) {
				const double delta = 1e-2;

				auto confidence = [=]( uint n ) { return sigma * sqrt( 2 * log( 2 / delta ) / n ); };

				vector<double> lower( N_k ), upper( N_k );

				// Algorithm Step 1: get initial distSum, distCount = 1.
				for ( uint i = 0; i < N_k; i++ ) {
					uint candidatePos = clusterAssignment[i];
					distSum[candidatePos] = 0;
					distCount[candidatePos] = 0;

					while ( distCount[candidatePos] == 0 ) {
						uint j = rand( 0, N_k - 1 );
						uint comparePos = clusterAssignment[j];

						if ( comparePos == candidatePos ) continue;

						distSum[candidatePos] += distance( kmerCodes[candidatePos], kmerCodes[comparePos], kmerLength );
						distCount[candidatePos] ++;
						double conf = confidence( 1 );
						lower[i] = distSum[candidatePos] - conf;
						upper[i] = distSum[candidatePos] + conf;
						// cerr << "distConf[candidatePos] = " << distConf[candidatePos] << "\n";
					}
				}

				// Step 2-12
				while ( true ) {
					// Keep the best and second-best candidates
					uint turnCandidate = numeric_limits<uint>::max();
					uint turnIndex = 0;
					double turnLimit = numeric_limits<double>::max();

					// Step 3.
					for ( uint i = 0; i < N_k; i++ ) {
						uint candidatePos = clusterAssignment[i];
						double limit = lower[i];

						if ( limit < turnLimit ) {
							turnLimit = limit;
							turnCandidate = candidatePos;
							turnIndex = i;
						}
					}

					// Step 4.
					if ( distCount[turnCandidate] < N_k - 1 ) {
						uint comparePos;
						do {
							uint j = rand( 0, N_k - 1 );
							comparePos = clusterAssignment[j];
						}
						while ( comparePos == turnCandidate );

						distSum[turnCandidate] += distance( kmerCodes[turnCandidate], kmerCodes[comparePos], kmerLength );
						distCount[turnCandidate] ++;
						double conf = confidence( distCount[turnCandidate] );
						double mu = double( distSum[turnCandidate] ) / distCount[turnCandidate];
						lower[turnIndex] = mu - conf;
						upper[turnIndex] = mu + conf;
					}
					else {
						distSum[turnCandidate] = 0;
						distCount[turnCandidate] = N_k - 1;

						for ( uint i = 0; i < N_k; i++ ) {
							uint comparePos = clusterAssignment[i];

							if ( comparePos == turnCandidate ) continue;

							distSum[turnCandidate] += distance( kmerCodes[turnCandidate], kmerCodes[comparePos], kmerLength );
						}

						double mu = double( distSum[turnCandidate] ) / distCount[turnCandidate];
						lower[turnIndex] = mu;
						upper[turnIndex] = mu;
					}

					// Step 9.
					// The only point where anything has changed is the best candidate.
					// Just need to compare it to the second-best candidate to see if
					bool allDone = true;

					for ( uint j = 0; j < N_k; j++ ) {
						if ( j != turnIndex && lower[j] < upper[turnIndex] ) {
							allDone = false;
							break;
						}
					}

					if ( allDone ) {
						protoPtr = kmers[turnCandidate];
						protoCode = kmerCodes[turnCandidate];
						break;
					}
				}
			}
			else if ( N_k == 1 ) {
				protoPtr = kmers[clusterAssignment[0]];
				protoCode = kmerCodes[clusterAssignment[0]];
			}
			else {
				protoPtr = nullptr;
				protoCode = nullptr;
			}
		}

		static void CreateCluster(
			Kmer *protoKmer,
			Alphabet & alphabet,
			uint kmerLength,
			DistanceFunction &distance,
			vector<Kmer *> & kmers,
			vector<uint_least32_t> &clusterAssignment,
			vector<Kmer *> & clusterProtos,
			vector<Cluster *> &clusters //
		) {
			Kmer * prototype = new Kmer( protoKmer->Sequence(), protoKmer->KmerPosition(), protoKmer->Length() );
			clusterProtos.push_back( prototype );

			Cluster *c = new Cluster( prototype, 0, distance );
			clusters.push_back( c );

			uint N_k = clusterAssignment.size();

			for ( uint i = 0; i < N_k; i++ ) {
				c->Add( kmers[clusterAssignment[i]], numeric_limits<Distance>::max() );
			}
		}


	};
}
