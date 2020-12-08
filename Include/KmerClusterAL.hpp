#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include <iterator>
#include <cstdio>
#include <omp.h>

#include "db.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "SimilarityMatrix.hpp"
#include "KmerDistanceCache.hpp"

namespace QutBio {

	/**
	 *  <summary>
	 *	Pairwise symbol distance function mapping 0..(aSize-1) \times 0..(aSize-1) \to DistanceType.
	 *	I use KmerWord as the argument type, which may in some applications may contain 2 or 3
	 *	packed symbols, but in this context, it should contain a single encoded symbol value
	 *	between 0 and (aSize-1).
	 *	</summary>
	 */
	struct SymbolDistanceFunction {
		virtual Distance GetDistance1(KmerWord, KmerWord) const = 0;
	};

	/**
	 *  <summary>
	 *	An Average-linkage Kmer Cluster.
	 *
	 *  The distance of a kmer to the cluster as a whole is the average of the distances from the
	 *  kmer to the individual kmers assigned to the cluster.
	 *
	 *  A position-weight matrix is used to optimise the average calculation (maximum benefit
	 *  is attained for large clusters).
	 *	</summary>
	 */
	template <typename DistanceFunction, const int USE_PWM = 1>
	struct KmerClusterAL {
		typedef KmerClusterAL *P;

		// The length of kmers.
		size_t K;

		/**
		 *	<summary>
		 *	The size of alphabet.
		 *	</summary>
		 */
		size_t alphaSize;

		/**
		 *	<summary>
		 *	Function used to determine symbol distances.
		 *	</summary>
		 */
		const DistanceFunction & symbolCodeDist;

		/**
		 *	<summary>
		 *	The kmers associated with this cluster.
		 *	</summary>
		 */
		vector<Kmer *> kmers;

		/**
		 *	<summary>
		 * A k \times aSize matrix containing the count of each residue
		 * at each of the k positions in the kmer.
		 *
		 * Invariant:
		 * ----------
		 *	Let m = kmers.size().
		 *	For i in 0..(k-1):
		 *		For a in 0..(aSize-1):
		 *			pwm[i][a] == count (j=0..(m-1)) kmers[j][i] == a.
		 *		End
		 *	End
		 *	</summary>
		 */
		FlatMatrix<uint> pwm;

		/**
		 *	<summary>
		 *	Initialise a new, empty, Average-linkage cluster.
		 *	</summary>
		 */
		KmerClusterAL(
			size_t K,
			size_t alphaSize,
			const DistanceFunction & symbolCodeDist
		) :
			K(K),
			alphaSize(alphaSize),
			symbolCodeDist(symbolCodeDist)
			//
		{
			pwm.resize(K, alphaSize);
			pwm.fill(0);
		}

		virtual ~KmerClusterAL() {
			for (auto kmer: kmers) {
				delete kmer;
			}
		}

		double DistanceTo(EncodedFastaSequence * seq, size_t offset) {
			if (kmers.size() == 0) {
				return NAN;
			}

			EncodedKmer code = seq->GetEncodedKmer1(offset);
			return DistanceTo(code);
		}

		double DistanceTo(Kmer & kmer) {
			return DistanceTo(&kmer.Sequence(), kmer.Instances()[0].kmerPosition);
		}

		double DistanceTo(KmerWord code[]) {
			double sum = 0;

			if (USE_PWM) {
				for (uint k = 0; k < K; k++) {
					auto weights = pwm.row(k);

					for (KmerWord j = 0; j < alphaSize; j++) {
						auto d = symbolCodeDist.GetDistance1(j, code[k]);
						sum += d * weights[j];
					}
				}
			}
			else {
				for (auto & kmer : kmers) {
					sum += symbolCodeDist(kmer->PackedEncoding(), code, kmer->Length());
				}
			}

			sum /= kmers.size();

			return sum;
		}

		double DistanceTo(KmerClusterAL & other) {
			double sum = 0;

			for (uint k = 0; k < K; k++) {
				auto w1 = pwm.row(k);
				auto w2 = other.pwm.row(k);

				for (KmerWord i = 0; i < alphaSize; i++) {
					for (KmerWord j = 0; j < alphaSize; j++) {
						auto d = symbolCodeDist.GetDistance1(i, j);
						sum += d * w1[i] * w2[j];
					}
				}
			}

			sum /= (kmers.size() * other.kmers.size());
			return sum;
		}

		bool NearerTo(Kmer & kmer, double & nearest) {
			return NearerTo(&kmer.Sequence(), kmer.Instances()[0].kmerPosition, nearest);
		}

		bool NearerTo(EncodedFastaSequence * seq, size_t offset, double & nearest) {
			return NearerTo(seq->GetEncodedKmer1(offset), nearest);
		}

		bool NearerTo(KmerWord code[], double & nearest) {
			double sum = 0;
			double scaledNearest = nearest * kmers.size();

			if (USE_PWM) {
				for (uint k = 0; k < K; k++) {
					auto weights = pwm.row(k);

					for (KmerWord j = 0; j < alphaSize; j++) {
						auto d = symbolCodeDist.GetDistance1(j, code[k]);
						sum += d * weights[j];
					}
				}
			}
			else {
				for (auto & kmer : kmers) {
					sum += symbolCodeDist(kmer->PackedEncoding(), code, kmer->Length());
				}
			}

			if (sum >= scaledNearest) {
				return false;
			}

			sum /= kmers.size();
			nearest = sum;

			return true;
		}

		bool NearerTo(KmerClusterAL & other, double & nearest) {
			double sum = 0;
			double scaledNearest = nearest * kmers.size() * other.kmers.size();

			for (uint k = 0; k < K; k++) {
				auto w1 = pwm.row(k);
				auto w2 = other.pwm.row(k);

				for (KmerWord i = 0; i < alphaSize; i++) {
					for (KmerWord j = 0; j < alphaSize; j++) {
						auto d = symbolCodeDist.GetDistance1(i, j);
						sum += d * w1[i] * w2[j];
					}
				}
			}

			if (sum > scaledNearest) {
				return false;
			}

			sum /= (kmers.size() * other.kmers.size());
			nearest = sum;
			return true;
		}

		/*
		**	Adds a kmer to the cluster, updating the kmer collection
		**	and the position weight matrix.
		*/

		void Add(EncodedFastaSequence * seq, size_t offset) {
			kmers.push_back(new Kmer(*seq, offset, K));
			EncodedKmer code = kmers[kmers.size() - 1]->UnpackedEncoding();

			for (uint k = 0; k < K; k++) {
				pwm(k, code[k])++;
			}
		}

		/*
		**	Adds a kmer to the cluster, updating the kmer collection
		**	and the position weight matrix.
		*/

		void Add(Kmer & kmer) {
			Add(&kmer.Sequence(), kmer.KmerPosition());
		}

		/*
		**	Adds all kmers from another cluster to this cluster. Equivalent
		**	to calling Add(const Kmer &) once for each kmer in the other
		**	cluster and the position weight matrix.
		*/

		void Add(KmerClusterAL & other) {
			for (auto & kmer : other.kmers) {
				Add(kmer);
			}
		}

		/*
		**	Merges another cluster with this cluster. Equivalent
		**	to calling Add(const Kmer &) once for each kmer in the other
		**	cluster and the position weight matrix, but should be somewhat
		**	faster
		*/

		void Merge(KmerClusterAL & other) {
			kmers.insert(kmers.end(), other.kmers.begin(), other.kmers.end());
			pwm += other.pwm;
		}

		/*
		**	Removes all kmers from the cluster.
		*/

		void Clear() {
			for (auto kmer: kmers) {
				delete kmer;
			}
			kmers.clear();
			pwm.fill(0);
		}

		/*
		**	Gets the average distance from a kmer in the cluster to the
		**	centroid.
		*/

		double AverageDistortion() {
			if (kmers.size() < 2) {
				return NAN;
			}

			double totalDistortion = 0;

			for (auto & kmer : kmers) {
				totalDistortion += DistanceTo(*kmer);
			}

			return totalDistortion / kmers.size();
		}

		using Prototype = pair<Kmer*, double>;

		/*
		**	Gets the (first encountered) stored kmer with the lowest quantisation error.
		*/

		Prototype GetCentralKmer(void) {
			Prototype result;

			if (kmers.size() < 1) {
				result.first = 0;
				result.second = NAN;
			}
			else {
				result.second = numeric_limits<double>::max();

				for (auto & kmer : kmers) {
					auto dist = DistanceTo(kmer);

					if (dist < result.second) {
						result.first = &kmer;
						result.second = dist;
					}
				}
			}

			return result;
		}

		/*
		**	Gets a string of length K made up of the most frequently appearing symbols
		**	at each position in 0..(K-1).
		*/

		string GetConsensus(const string & symbols) {
			string result = "";

			if (kmers.size() >= 1) {
				for (size_t i = 0; i < K; i++) {
					int bestIdx = 0;

					for (size_t j = 1; j < alphaSize; j++) {
						if (pwm(i, j) > pwm(i, bestIdx)) {
							bestIdx = j;
						}
					}

					result.push_back(symbols[bestIdx]);
				}
			}

			return result;
		}

		/*
		**	Gets the average distance from a kmer in the cluster to the
		**	centroid.
		*/

		double AverageDistortionParallel() {
			if (kmers.size() < 2) {
				return NAN;
			}

			int numThreads;

#pragma omp parallel
#pragma omp master
			numThreads = omp_get_num_threads();

			vector<double> totalDistortion(numThreads);
			const int N = (int)kmers.size();

#pragma omp parallel for num_threads(numThreads)
			for (int i = 0; i < N; i++) {
				int threadNum = omp_get_thread_num();
				totalDistortion[threadNum] += DistanceTo(*kmers[i]);
			}

			return sum(totalDistortion.begin(), totalDistortion.end(), 0.0) / kmers.size();
		}

		/*
		**	Gets the (first encountered) stored kmer with the lowest quantisation error.
		**	Uses explicit OpenMP parallel for.
		*/

		Prototype GetCentralKmerParallel() {
			if (kmers.size() < 1) {
				Prototype result;
				result.first = 0;
				result.second = NAN;
				return result;
			}
			else {
				int numThreads;

#pragma omp parallel
#pragma omp master
				numThreads = omp_get_num_threads();

				vector<Prototype> results(numThreads);
				int N = (int)kmers.size();

#pragma omp parallel for
				for (int i = 0; i < numThreads; i++) {
					results[i].second = numeric_limits<double>::max();
				}

#pragma omp parallel for
				for (int i = 0; i < N; i++) {
					int threadNum = omp_get_thread_num();
					auto & prototype = results[threadNum];

					Kmer & kmer = *kmers[i];

					auto dist = DistanceTo(kmer);

					if (dist < prototype.second) {
						prototype.first = &kmer;
						prototype.second = dist;
					}
				}

				auto result = results[0];

				for (int threadNum = 1; threadNum < numThreads; threadNum++) {
					auto & prototype = results[threadNum];

					if (prototype.second < result.second) {
						result = prototype;
					}
				}

				return result;
			}
		}

		/**
		**	Computes the information content of a cluster. This is necessarily
		**	going to be a bit shaky because some kmers belong to sequences that
		**	are members of multiple classes.
		**	I count these once for each distinct class.
		**	TODO: review literature for alternative measure of homogeneity.
		*/

		double Information(void) {
			double info = 0;

			if (kmers.size() > 0) {
				int overallClassCount = EncodedFastaSequence::ClassNumberRegistry().size();
				int * freq = (int *)alloca(overallClassCount * sizeof(int));
				memset(freq, 0, overallClassCount * sizeof(int));
				int count = 0;

				for (auto & kmer : kmers) {
					for (auto & classLabel : kmer->Sequence().classNumbers) {
						freq[classLabel] ++;
						count++;
					}
				}

				for (int i = 0; i < overallClassCount; i++) {
					if (freq[i] > 0) {
						double p = (double)freq[i] / count;
						info += p * log(p) / log(2);
					}
				}
			}

			return -info;
		}

		/*
		**	Serialise the position weight matrix to a stream.
		*/

		void SerializePwm(ostream & out) {
			auto rows = pwm.rows();
			auto cols = pwm.cols();
			auto buffer = pwm.buffer();

			out << "pwm," << rows << ',' << cols;

			for (size_t i = 0; i < rows * cols; i++) {
				out << ',' << buffer[i];
			}

			out << "\n";
		}

		/*
		**	Serialise the kmers, saving the sequence id, position, and length only.
		**	the length may appear redundant, but it might be useful in situations
		**	where kmers of diverse lengths are in simultaneous use in some hypothetical
		**	future codebook design.
		*/

		void SerializeKmers(ostream & out) {
			out << "kmers," << kmers.size() << "\n";

			for (auto & kmer : kmers) {
				out << kmer << "\n";
			}
		}

		/*
		**	Serialise the cluster to a stream in text format.
		*/

		void Serialize(ostream & out) {
			out << "cluster" << "\n";
			SerializePwm(out);
			SerializeKmers(out);
		}

		/*
		**	Serialise a list of clusters to a stream in text format.
		*/

		static void Serialize(ostream & out, const vector<KmerClusterAL *> & clusters) {
			for (auto cluster : clusters) {
				cluster->Serialize(out);
			}
		}

		/*
		**	Detribalise a text-formatted cluster from a stream (if present and all OK)
		**	and return a pointer to the resulting object. Returns null if an error is
		**	encountered.
		*/

		static P Deserialize(
			istream & in,
			const Index<EncodedFastaSequence> & idx,
			const DistanceFunction & symbolCodeDist
		) {
			P result = 0;
			string inputLine;

			getline(in, inputLine);

			if (inputLine.find("cluster") != 0) return 0;

			getline(in, inputLine);

			if (inputLine.find("pwm,") != 0) return 0;

			const char * nextTerm = inputLine.c_str();

			auto skip = [&nextTerm]() {
				while (*nextTerm && *nextTerm != ',') {
					nextTerm++;
				}

				if (*nextTerm == ',') {
					nextTerm++;
				}
			};

			int rows, cols;

			skip(); sscanf(nextTerm, "%d", &rows);
			skip(); sscanf(nextTerm, "%d", &cols);

			result = new KmerClusterAL(rows, cols, symbolCodeDist);

			uint *buffer = (uint *) result->pwm.buffer();

			for (int i = 0; i < rows * cols; i++) {
				skip(); sscanf(nextTerm, "%d", &buffer[i]);
			}

			getline(in, inputLine);

			if (inputLine.find("kmers,") != 0) {
				delete result;
				return 0;
			}

			nextTerm = inputLine.c_str();

			size_t kmerCount;

			skip(); sscanf(nextTerm, "%zu", &kmerCount);

			for (size_t i = 0; i < kmerCount; i++) {
				getline(in, inputLine);
				nextTerm = inputLine.c_str();
				skip();
				size_t idLen = nextTerm - inputLine.c_str() - 1;
				string id = inputLine.substr(0, idLen);

				if (idx.find(id) != idx.end()) {
					auto sequence = idx.at(id);
					size_t pos;
					sscanf(nextTerm, "%zu", &pos);
					result->Add(sequence, pos);
				}
			}

			return result;
		}

		/*
		**	Detribalise a collection of clusters from a stream (if present and all OK)
		**	and save them in vector referenced by parameter 'clusters'.
		*/

		static void Deserialize(
			istream & in,
			const Index<EncodedFastaSequence> & idx,
			const DistanceFunction & symbolCodeDist,
			vector<P> &clusters
		) {
			P cluster;

			while ((cluster = Deserialize(in, idx, symbolCodeDist))) {
				if (cluster->kmers.size() > 0) {
					clusters.push_back(cluster);
				}
				else {
					delete cluster;
				}
			}
		}


		/**
		*	Populates a list of clusters by testing each kmer against all existing clusters;
		*	if distance from kmer to nearest cluster (equivalently, to all existing clusters)
		*	is greater than the threshold, then a new cluster is seeded.
		*
		*	NB: the clusters are dynamically created via new, and must be deleted later
		*	to avoid memory leakage.
		*/

		static void DoIncrementalClusteringParallel(
			vector<Kmer *> & allKmers,
			uint K,
			double threshold,
			size_t alphaSize,
			const DistanceFunction & symbolCodeDist,
			vector<KmerClusterAL<DistanceFunction> *> & clusters
		) {
			using Clusters = vector<KmerClusterAL<DistanceFunction> *>;

			int numThreads = 1;

			#if USE_OMP
#pragma omp parallel
#pragma omp master
			numThreads = omp_get_num_threads();
#endif

			vector<Clusters> newClusters(numThreads);
			size_t origClusterCount = clusters.size();

#if USE_OMP
#pragma omp parallel for     // Create place-holders for any pre-existing clusters. 
#endif
			for (int threadId = 0; threadId < numThreads; threadId++) {
				for (size_t j = 0; j < origClusterCount; j++) {
					newClusters[threadId].push_back(new KmerClusterAL(K, alphaSize, symbolCodeDist));
				}
			}

			vector<int> ctr(numThreads);

			fprintf(stderr, "threadId\tkmers\tclusters\ttime\n");

			auto start = omp_get_wtime();

#if USE_OMP
#pragma omp parallel for
#endif
			for (size_t i = 0; i < allKmers.size(); i++) {
				int threadId = omp_get_thread_num();
				auto & clustersPerThread = newClusters[threadId];
				ctr[threadId] ++;

				Kmer & kmer = * allKmers[i];

				auto nearestClusterIndex = -1;
				auto nearestDistance = numeric_limits<double>::max();

				for (size_t j = 0; j < clustersPerThread.size(); j++) {
					// Try to map kmer to a pre-existing cluster if possible, otherwise search the
					// newly added clusters.
					auto cluster = j < origClusterCount ? clusters[j] : clustersPerThread[j];

					// Always go with the first cluster that falls within the threshold.
					// This way, any well-conserved kmers will be quickly mapped into clusters.
					// Hopefully, these are in the majority.
					auto dist = cluster->DistanceTo(kmer);

					if (dist < nearestDistance) {
						nearestDistance = dist;
					}

					if (nearestDistance <= threshold) {
						nearestClusterIndex = j;
						break;
					}
				}

				if (nearestClusterIndex >= 0) {
					clustersPerThread[nearestClusterIndex]->Add(kmer);
				}
				else {
					auto cluster = new KmerClusterAL<DistanceFunction>(K, alphaSize, symbolCodeDist);
					clustersPerThread.push_back(cluster);
					cluster->Add(kmer);
				}

				if (ctr[threadId] % 1000 == 0) {
					auto now = omp_get_wtime();
#if USE_OMP
#pragma omp critical
#endif
					{
						fprintf(stderr, "%d\t%d\t%zu\t%f\n", threadId, ctr[threadId], clustersPerThread.size(), now - start);
						fflush(stderr);
					}
				}
			}

			//	At this stage we have a list of lists of clusters which collectively 
			//	cover all kmers in the dataset, but need to be merged.

			// Coalesce the items that were in the original set of clusters.
#if USE_OMP
#pragma omp parallel for       
#endif
			for (size_t i = 0; i < origClusterCount; i++) {
				clusters[i]->Clear();

				for (int threadId = 0; threadId < numThreads; threadId++) {
					auto & newItems = newClusters[threadId];
					clusters[i]->Merge(*newItems[i]);
				}
			}

			// Append any new clusters (after the first sweep through, hopefully this 
			// will be "none"). 
			for (int threadId = 0; threadId < numThreads; threadId++) {
				auto & newItems = newClusters[threadId];

				for (size_t j = origClusterCount; j < newItems.size(); j++) {
					clusters.push_back(newItems[j]);
				}
			}

			sort(clusters.begin(), clusters.end(), [](KmerClusterAL<DistanceFunction> * const& a, KmerClusterAL<DistanceFunction> * const& b) {
				return a->kmers.size() > b->kmers.size();
			});

			int i;

			for (i = clusters.size() - 1; i >= 0 && clusters[i]->kmers.size() == 0; i--) {
				delete clusters[i];
				clusters[i] = 0;
			}

			clusters.resize(i + 1);
		}
	};

	/**
	**  This is an equivalence relation. Complexity is linear in the size of the
	**  matrices. Matrices are considered equivalent if their dimensions are equal,
	**  and if corresponding elements compare equal.
	*/
	template <typename T>
	inline bool operator==(const KmerClusterAL<T> & x, const KmerClusterAL<T> & y) {
		return (x.kmers == y.kmers && x.pwm == y.pwm);
	}

	/**
	**  This is an equivalence relation. Complexity is linear in the size of the
	**  matrices. Matrices are considered equivalent if their dimensions are equal,
	**  and if corresponding elements compare equal.
	*/
	template <typename T>
	inline bool operator!=(const KmerClusterAL<T> & x, const KmerClusterAL<T> & y) {
		return x.kmers != y.kmers || x.pwm != y.pwm;
	}

	template <typename T>
	inline ostream & operator << (ostream & stream, const KmerClusterAL<T> & cluster) {
		for (size_t i = 0; i < cluster.kmers.size(); i++) {
			stream << cluster.kmers[i]->Substr() << '\n';
		}

		return stream;
	}
}
