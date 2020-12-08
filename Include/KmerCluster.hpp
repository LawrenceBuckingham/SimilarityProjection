#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <functional>
#include <numeric>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <csignal>

#include "FastaSequence.hpp"
#include "Kmer.hpp"
#include "KmerClusterPrototype.hpp"
#include "KmerDistanceCache.hpp"
#include "KmerIndex.hpp"
#include "SimilarityMatrix.hpp"
#include "Selector.hpp"
#include "SimilarityMatrix.hpp"

#if USE_ISSL
#include "IsslKmerIndex.hpp"
#endif

#if defined(SHOW_PROGRESS)
#define PROGRESS(x) x
#else
#define PROGRESS(x)
#endif

#if USE_OMP
#include <omp.h>
#endif

using namespace std;

namespace QutBio {
	/**
		**	<summary>
		**		Template class representing cluster of Kmers with a "central" prototype.
		**	</summary>
		**	<typeparam name=DistanceFunction>
		**		A class which provides a function
		**			Distance GetDistance(EncodedKmer, EncodedKmer, size_t);
		**		which returns the distance between two kmers.
		**	</typeparam>
		*/
	template <typename DistanceFunction, typename KmerType>
	struct KmerCluster {
		// using DistanceFunction = function<Distance (KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength )>;
		using Cluster = KmerCluster<DistanceFunction, KmerType>;

		Kmer * prototype;
		vector<pair<KmerType *, Distance>> kmers;
		size_t expectedSize;
		size_t index;
		const DistanceFunction &distanceFunction;
		vector<vector<pair<KmerType *, Distance>>> kmersPerThread;
		map<string, string> metadata;

#if USE_OMP
		omp_lock_t lock;
#endif

		static bool *interrupted() {
			static bool interrupted = false;
			return &interrupted;
		}

		static void sigintHandler(int arg) {
			auto i = interrupted();
			*i = true;
		}

		KmerCluster(
			Kmer * prototype,
			size_t expectedSize,
			const DistanceFunction &distanceFunction
		) : prototype(prototype),
			expectedSize(expectedSize),
			distanceFunction(distanceFunction) {

#if USE_OMP
			omp_init_lock(&lock);
#endif
		}

		KmerCluster(const KmerCluster & other) = delete;
		KmerCluster & operator=(const KmerCluster & other) = delete;

		//KmerCluster(PrototypeSequenceType * protoSeq, size_t expectedSize, const DistanceFunction & distanceFunction) :
		//	prototype(Substring(protoSeq->Sequence().c_str(), 0, protoSeq->Sequence().length())),
		//	expectedSize(expectedSize),
		//	distanceFunction(distanceFunction)
		//	//
		//{
		//	prototype->Add(protoSeq, 0);
		//}

		virtual ~KmerCluster() {
#if USE_OMP
			omp_destroy_lock(&lock);
#endif
		}

		friend ostream & operator << (ostream &str, const Cluster & cluster) {
			str << "Cluster"
				<< ","
				<< cluster.kmers.size()
				<< ","
				<< *cluster.prototype;

			for (auto & kvp : cluster.metadata) {
				str << kvp.first << ":" << kvp.second << ";";
			}

			str << endl;

			for (auto & kmer : cluster.kmers) {
				if (kmer.second < numeric_limits<Distance>::max()) {
					str << "distance:" << kmer.second << ";";
				}
				str << *kmer.first << '\n';
			}

			return str;
		}

		/**
		 *	Gets the number of kmer instances assigned to the cluster.
		 *	This will generally be greater than the number of kmers, because each
		 *	kmer may appear (identically) in multiple sequences or positions within
		 *	a sequence.
		 */

		size_t InstanceCount() {
			size_t count = std::accumulate(kmers.begin(), kmers.end(), 0,
				[](size_t accumulated, pair<KmerType *, Distance> &kmer) { return accumulated + kmer.first->Instances().size(); });
			return count;
		}

		/*
		**	Appends a kmer to the list of kmers attached to this cluster.
		**	Not OMP thread-safe.
		*/

		void Add(KmerType *kmer, Distance distance) {
			kmers.emplace_back(kmer, distance);
		}

#if 0
		/*
**	Appends a kmer to the list of kmers attached to this cluster.
**	OMP thread-safe.
*/

		void AddParallel(KmerType *kmer) {
#if USE_OMP
			omp_set_lock(&lock);
#endif
			kmers.push_back(kmer);
#if USE_OMP
			omp_unset_lock(&lock);
#endif
		}

#endif // 0

#if 0
		void AddParallel(const vector<Kmer *> &kmers) {
#if USE_OMP
			omp_set_lock(&lock);
#endif
			for (auto kmer : kmers) {
				Add(kmer);
			}
#if USE_OMP
			omp_unset_lock(&lock);
#endif
		}

#endif // 0

		size_t Size() const {
			size_t total = 0;

			for (auto & k : kmers) {
				total += k.first->Instances().size();
			}

			return total;
		}

		void AllocateKmersToThreads(size_t numThreads) {
			kmersPerThread.resize(numThreads);

			for (size_t i = 0; i < kmers.size(); i++) {
				kmersPerThread[i % numThreads].push_back(kmers[i]);
			}
		}

		virtual double DistanceTo(const Kmer *kmer) {
			EncodedKmer thisCode = prototype->PackedEncoding();
			EncodedKmer kmerCode = kmer->PackedEncoding();
			Distance dist = distanceFunction(thisCode, kmerCode, prototype->Substr().Length());
			return dist;
		}

		virtual double DistanceTo(EncodedFastaSequence *seq, size_t kmerPosition) {
			EncodedKmer thisCode = prototype->PackedEncoding();
			EncodedKmer kmerCode = seq->GetEncodedKmer(kmerPosition);
			Distance dist = distanceFunction(thisCode, kmerCode, prototype->Substr().Length());
			return dist;
		}

		virtual double DistanceTo(EncodedKmer encodedKmer) {
			EncodedKmer thisCode = prototype->PackedEncoding();
			Distance dist = distanceFunction(thisCode, encodedKmer, prototype->Substr().Length());
			return dist;
		}

		static void AddKmersToExistingClusterParallel_Good_Slowish(
			vector<KmerType *> & kmers,
			size_t & firstUnalloc,
			DistanceFunction & symbolCodeDist,
			int K,
			double threshold,
			QutBio::KmerCluster<DistanceFunction, KmerType> * newCluster
		) {
			if (firstUnalloc >= kmers.size()) return;

			size_t numRemaining = kmers.size() - firstUnalloc;

			QutBio::Kmer * protoKmer = newCluster->prototype;
			auto prototypeEncoding = protoKmer->PackedEncoding();
			vector<size_t> allAssigned;
			allAssigned.reserve(numRemaining);

#if USE_OMP
#pragma omp parallel
#endif
			{
				vector<size_t> assigned;
				assigned.reserve(numRemaining);
#if USE_OMP
#pragma omp for
#endif
				for (size_t i = firstUnalloc; i < kmers.size(); i++) {
					KmerType *kmer = kmers[i];
					auto kmerEncoding = kmer->PackedEncoding();
					auto dist = symbolCodeDist(kmerEncoding, prototypeEncoding, K);

					//#pragma omp critical
					//					(cerr << "i = " << i << ", dist = " << dist << "\n").flush();

					if (dist <= threshold) {
						kmer->SetDistanceFromPrototype(dist);
						assigned.push_back(i);
					}
				}

#pragma omp critical
				{
					allAssigned.insert(allAssigned.end(), assigned.begin(), assigned.end());
				}
			}

			std::sort(allAssigned.begin(), allAssigned.end());

			for (auto i : allAssigned) {
				newCluster->Add(kmers[i]);

				if (i > firstUnalloc) {
					std::swap(kmers[i], kmers[firstUnalloc]);
				}

				firstUnalloc++;
			}
		}

		static void AddKmersToExistingClusterParallel(
			vector<KmerType *> & kmers,
			size_t & firstUnalloc,
			DistanceFunction & symbolCodeDist,
			int K,
			double threshold,
			QutBio::KmerCluster<DistanceFunction, KmerType> * newCluster
		) {
			if (firstUnalloc >= kmers.size()) return;

			QutBio::Kmer * protoKmer = newCluster->prototype;
			auto prototypeEncoding = protoKmer->PackedEncoding();

#if USE_OMP
#pragma omp parallel
#endif
			{
#if USE_OMP
#pragma omp for
#endif
				for (size_t i = firstUnalloc; i < kmers.size(); i++) {
					KmerType *kmer = kmers[i];
					auto kmerEncoding = kmer->PackedEncoding();
					auto dist = symbolCodeDist(kmerEncoding, prototypeEncoding, K);

					//#pragma omp critical
					//					(cerr << "i = " << i << ", dist = " << dist << "\n").flush();

					if (dist <= threshold) {
#pragma omp critical
						{
							newCluster->Add(kmer, dist);

							if (i > firstUnalloc) {
								std::swap(kmers[i], kmers[firstUnalloc]);
							}

							firstUnalloc++;
						}
					}
				}

			}
		}

		static void AddKmersToExistingClusterParallel(
			vector<KmerType *> & kmers,
			size_t & firstUnalloc,
			DistanceFunction & symbolCodeDist,
			int K,
			double threshold,
			vector<Cluster *> &newClusters
		) {
			const size_t C = newClusters.size();

			if (firstUnalloc >= kmers.size()) return;

			vector<EncodedKmer> prototypeEncoding(C);

			for (size_t c = 0; c < C; c++) {
				auto protoKmer = newClusters[c]->prototype;
				prototypeEncoding[c] = protoKmer->PackedEncoding();
			}

#if USE_OMP
#pragma omp parallel
#endif
			{
#if USE_OMP
#pragma omp for
#endif
				for (size_t i = firstUnalloc; i < kmers.size(); i++) {
					KmerType *kmer = kmers[i];
					auto kmerEncoding = kmer->PackedEncoding();

					for (size_t c = 0; c < C; c++) {
						auto protoEncoding = prototypeEncoding[c];
						auto dist = symbolCodeDist(kmerEncoding, protoEncoding, K);

						//#pragma omp critical
						//					(cerr << "i = " << i << ", dist = " << dist << "\n").flush();

						if (dist <= threshold) {
#if USE_OMP
#pragma omp critical
#endif
							{
								newClusters[c]->Add(kmer, dist);

								if (i > firstUnalloc) {
									std::swap(kmers[i], kmers[firstUnalloc]);
								}

								firstUnalloc++;
							}

							break;
						}
					}
				}

			}
		}

		static vector<Cluster *> DoExhaustiveIncrementalClustering(
			vector<KmerType *> & allKmers,
			int wordLength,
			double threshold,
			size_t alphaSize,
			DistanceFunction &distanceFunction,
			UniformRealRandom &rand,
			size_t increment,
			function<KmerClusterPrototype *(Kmer *kmer)> createPrototype,
			vector<function<void(KmerClusterPrototype *proto, KmerCluster * cluster)>> & process
			//
		) {
			size_t firstUnalloc = 0;
			size_t N = allKmers.size();

			// TODO: select random prototypes, don't shuffle everything. 
			// Shuffle the kmers to avoid pathological selection.
			for (size_t i = 0; i < N; i++) {
				auto newLoc = (size_t)(rand() * N);
				std::swap(allKmers[i], allKmers[newLoc]);
			}

			// TODO: Don't worry about it; we'll get empty clusters anyway
			//  Eliminate all kmers from consideration that are impossible cluster centres.
			//	due to excessive self-match distance.
			for (size_t i = 0; i < allKmers.size(); i++) {
				auto encoding = allKmers[i]->PackedEncoding();
				Distance selfMatchDistance = distanceFunction(encoding, encoding, wordLength);

				if (selfMatchDistance > threshold) {
					std::swap(allKmers[i], allKmers[firstUnalloc]);
					firstUnalloc++;
				}
			}

			vector<Cluster *> newClusters(increment);
			vector<KmerClusterPrototype *> newProtos(increment);

			uint numberAdded = 0;

			while (firstUnalloc < N) {
				(cerr << "\n" << (N - firstUnalloc) << " unassigned kmers.                               ").flush();

				newClusters.clear();
				newProtos.clear();

				for (size_t nc = 0; nc < increment && firstUnalloc + nc < N; nc++) {
					numberAdded ++;
					auto nextKmer = allKmers[firstUnalloc];
					KmerClusterPrototype *proto = createPrototype(nextKmer);
					// (cerr << "\nnumberAdded = " << numberAdded << "\n" << "Proto = " << proto).flush();
					Kmer *protoKmer = proto->SingletonKmer();
					KmerCluster * newCluster = new KmerCluster(protoKmer, 0, distanceFunction);
					auto encoding = nextKmer->PackedEncoding();
					auto dist = distanceFunction(encoding, encoding, wordLength);
					newCluster->Add(nextKmer, dist);
					newClusters.push_back(newCluster);
					newProtos.push_back(proto);
					firstUnalloc++;
				}

				AddKmersToExistingClusterParallel(
					allKmers,
					firstUnalloc,
					distanceFunction,
					wordLength,
					threshold,
					newClusters
				);

				for (size_t i = 0; i < newProtos.size(); i++) {
					auto proto = newProtos[i];
					auto newCluster = newClusters[i];

					for (auto proc : process) {
						proc(proto, newCluster);
					}
				}
			}

			(cerr << "added " << numberAdded << " prototypes.\n").flush();

			return newClusters;
		}

		static void InitialiseClusters(
			vector<KmerClusterPrototype *> & protos,
			size_t wordLength,
			const DistanceFunction &dist,
			vector<KmerCluster *> &clusters
		) {
			for (auto proto : protos) {
				Kmer * protoKmer = proto->SingletonKmer();
				clusters.push_back(new KmerCluster(protoKmer, 0, dist));
			}
		}

		vector<KmerType *> & GetKmers() {
			return kmers;
		}

		void AddMetadata(const string & key, const string & value) {
			metadata[key] = value;
		}

		template <typename T>
		T GetMetadataValue(const string & key) {
			istringstream str(metadata[key]);
			T value;
			str >> value;
			return value;
		}

		template <typename T>
		vector<T> GetMetadataVector(const string & key) {
			vector<T> vect;
			auto parts = String::Split(metadata[key], ',');

			for (auto part : parts) {
				istringstream str();
				T value;
				str >> value;
				vect.push_back(value);
			}

			return vect;
		}

		/**
		 *	Calculate a measure of cluster purity based on the relevance (or otherwise) of the
		 *	sequences that contain the k-mer instances within the cluster. The measure used here is:
		 *	#{k1:(kmers & K(s1)); k2:(kmers & K(s2)) | (s1 related s2) } / #{k1:kmers; k2:kmers}.
		 *
		 *	Result is stored as metadata["purity"].
		 */
		void AssignPurityMeasure(
			map<string, set<string>> & homologs
		) {
			int relevant = 0;
			int total = 0;

			for (auto & k1 : kmers) {
				for (auto & i1 : k1.first->Instances()) {
					auto & id1 = i1.sequence.IdStr();
					auto & s1 = homologs.at(id1);

					for (auto & k2 : kmers) {
						for (auto & i2 : k2.first->Instances()) {
							auto & id2 = i2.sequence.IdStr();
							bool isRelevant = s1.count(id2);

							if (isRelevant) {
								relevant++;
							}

							total++;
						}
					}
				}
			}

			AddMetadata("purity", Double::ToString((double)relevant / total));
		}

		/**
		 *	Calculate a coarse empirical CDF of the pairwise distances between kmer instances
		 *	within the cluster.
		 *	N = #{k1:kmers; k2:kmers}, F = ECDF({k1:kmers; k2:kmers | d(k1, k2)})
		 *	distance bins are {i:(1..bins) | F(N*i/bins-1)}.
		 *
		 *	Result is stored as a comma-separated string in metadata["distance_bins"].
		 */
		void ComputeDistanceDistribution(
			DistanceFunction & distance,
			size_t wordLength,
			size_t bins
		) {
			if (bins == 0) throw Exception("Invalid argument: bins", FileAndLine);

			vector<Distance> distances;

			for (auto & k1 : kmers) {
				auto code1 = k1.first->PackedEncoding();

				for (auto & k2 : kmers) {
					auto code2 = k2.first->PackedEncoding();

					auto dist = distance(code1, code2, wordLength);

					for (size_t i = 0; i < k1.first->Instances().size() * k2.first->Instances().size(); i++) {
						distances.push_back(dist);
					}
				}
			}

			std::sort(distances.begin(), distances.end());

			if (distances.size() == 0) {
				cerr << "No distances for cluster at " << *prototype << "\n";
			}
			else {
				ostringstream buffer;

				for (size_t i = 1; i <= bins; i++) {
					int j = distances.size() * i / bins;
					buffer << distances[j - 1];

					if (i < bins) buffer << ",";
				}

				AddMetadata("distance_bins", buffer.str());
			}
		}
	};
} // namespace QutBio
