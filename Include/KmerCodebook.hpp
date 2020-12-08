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
#include <math.h>
#include <sstream>
#include <cstring>

#include "Assert.hpp"
#include "Alphabet.hpp"
#include "FastaSequence.hpp"
#include "Kmer.hpp"
#include "KmerClusterPrototype.hpp"
#include "SimilarityMatrix.hpp"
#include "Selector.hpp"
#include "KmerDistanceCache.hpp"
#include "KmerCluster.hpp"
#include "Index.hpp"

using namespace std;

namespace QutBio {
	using SequenceSet = unordered_set<pEncodedFastaSequence, pFastaHash>;

	template <typename DistanceFunction, typename KmerType>
	struct SequenceAssignment {
		using ClusterType = KmerCluster<DistanceFunction, KmerType>;

		ClusterType *cluster;
		SequenceSet sequences;

		SequenceAssignment(ClusterType * cluster) : cluster(cluster) // 
		{}
	};

	/**
		**	<summary>
		**		Template class representing codebook composed of 1-NN clusters.
		**	</summary>
		**	<typeparam name=DistanceFunction>
		**		A class which provides a function
		**			Distance GetDistance(EncodedKmer, EncodedKmer, size_t);
		**		which returns the distance between two kmers.
		**	</typeparam>
		**	<typeparam name=KmerType>
		**		A subclass of Kmer (or more likely, Kmer itself). This was introduced
		**		during a round of experiments in which an attempted hybrid ISSL/Cluster
		**		scheme was contemplated.
		**	</typeparam>
		*/
	template <typename DistanceFunction, typename KmerType>
	class KmerCodebook {
	public:
		enum InitMode {
			Uniform,
			LinearPerCluster,
			LogPerCluster
		};

		// using D = function<Distance ( KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength )>;
		using D = DistanceFunction;
		using K = KmerType;
		using Cluster = KmerCluster<D, K>;

	protected:
		size_t codebook_size;
		Alphabet *alphabet;
		const D &distanceFunction;
		uint charsPerWord;
		uint kmerLength;
		size_t splitClusterThreshold = numeric_limits<size_t>::max();

		// Index used to find sequence by name in the sequence database.
		const Index<EncodedFastaSequence> &dbIndex;

		// Index used to find prototype by name in the prototype database.
		const Index<KmerClusterPrototype> &protoIndex;

		// Private sequences which represent literal kmers, such as cluster
		// centroids which do not necessarily exist in the database, yet
		// appear as prototypes. By convention, I create these with a special
		// ID of the form "__chars:aaa,0" where aaa is the literal kmer
		// character string.
		vector<EncodedFastaSequence *> literalSequences;

		// Index used to find k-mers in the sequence database.
		KmerIndex &kmerIndex;

		bool ignoreInstances;

	public:
		FlatMatrix<KmerWord> kmerData;
		vector<Cluster *> codebook;

	public:
		KmerCodebook(
			//---------------------------------------
			pAlphabet alphabet,
			const D &distanceFunction,
			uint charsPerWord,
			uint kmerLength,
			Index<EncodedFastaSequence> &dbIndex,
			Index<KmerClusterPrototype> &protoIndex,
			KmerIndex &kmerIndex,
			FILE *savedCodebook,
			int splitClusterThreshold = -1,
			bool ignoreInstances = false
			//---------------------------------------
		) : alphabet(alphabet),
			distanceFunction(distanceFunction),
			charsPerWord(charsPerWord),
			kmerLength(kmerLength),
			splitClusterThreshold(splitClusterThreshold),
			dbIndex(dbIndex),
			protoIndex(protoIndex),
			kmerIndex(kmerIndex),
			ignoreInstances(ignoreInstances)
			//---------------------------------------
		{
			pCluster currentCluster = 0;

			const int BufferSize = 1000000;
			vector<char> _buffer(BufferSize);
			char *buffer = _buffer.data();
			size_t nextReadLoc = 0;
			size_t availSpace = BufferSize - 1;
			size_t bytesRead;

			vector<char *> lines;
			lines.reserve(BufferSize / 10);

			vector<pCluster> pendingClusters;
			vector<int> clusterFirstLine;
			vector<int> clusterLastLine;

			char *lastLineStart;
			bool codebookSeen = false;

			while (0 != (bytesRead = fread(buffer + nextReadLoc, sizeof(char), availSpace, savedCodebook))) {
				// Plug a trailing zero into the end of the buffer, just in case.
				buffer[BufferSize - 1] = 0;

				size_t max = bytesRead + nextReadLoc;
				size_t pos = 0;
				lines.clear();
				lastLineStart = 0;
				pendingClusters.clear();
				clusterFirstLine.clear();
				clusterLastLine.clear();

				if (currentCluster) {
					pendingClusters.push_back(currentCluster);
					clusterFirstLine.push_back(0);
					// clusterLastLine lags one behind, do not add anything at this point.
				}

				while (pos < max) {
					while (pos < max && (buffer[pos] == '\n' || buffer[pos] == 0)) {
						pos++;
					}

					lastLineStart = buffer + pos;

					while (pos < max && !(buffer[pos] == '\n' || buffer[pos] == 0)) {
						pos++;
					}

					if (pos < max) {
						buffer[pos++] = 0;
						lines.push_back(lastLineStart);
						lastLineStart = buffer + pos;
					}
				}

				for (size_t i = 0; i < lines.size(); i++) {
					char *l = lines[i];

					if ((!codebookSeen) && match(l, "Codebook")) {
						codebookSeen = true;
						processCodebook(l);
					}
					else if ((!codebookSeen) && match(l, "Clusters")) {
						codebookSeen = true;
						processClusters(l);
					}
					else if (match(l, "Cluster")) {
						// Insert last line of previous cluster.
						if (currentCluster)
							clusterLastLine.push_back(i);

						// Create new cluster
						processCluster(l, currentCluster);

						// remember the new cluster
						pendingClusters.push_back(currentCluster);
						clusterFirstLine.push_back(i + 1);
					}
				}

				// The trailing cluster will be incomplete, but it extends to the end of the buffer.
				clusterLastLine.push_back((int)lines.size());

				for (size_t c = 0; c < pendingClusters.size(); c++) {
					pCluster currentCluster = pendingClusters[c];

					for (int i = clusterFirstLine[c]; i < clusterLastLine[c]; i++) {
						if (!ignoreInstances) processLine(lines[i], currentCluster);
					}
				}

				availSpace = lastLineStart - buffer;
				nextReadLoc = max - availSpace;
				memcpy(buffer, lastLineStart, nextReadLoc);
			}

			size_t kmerCount = GetKmerCount();

			if (splitClusterThreshold > 0) {
				SplitLargeClusters();
#if defined(PARANOID_CODEBOOK)
				size_t splitKmerCount = GetKmerCount();
				assert_equal(kmerCount, splitKmerCount);
#endif
			}

			CopyKmerDataNonRecursive((int)Alphabet::WordsRequiredToPack(kmerLength, charsPerWord));

			(cerr << codebook.size() << " clusters parsed, indexing " << kmerCount << " kmers.\n").flush();
		}

		virtual ~KmerCodebook() {
			for (uint i = 0; i < codebook.size(); i++) {
				auto cluster = codebook[i];

				if (cluster)
					delete cluster;
			}

			for (uint i = 0; i < literalSequences.size(); i++) {
				auto seq = literalSequences[i];

				if (seq)
					delete seq;
			}
		}

		void AllocateKmersToThreads(size_t numThreads) {
			for (pCluster cluster : codebook) {
				cluster->AllocateKmersToThreads(numThreads);
			}
		}

		size_t GetKmerCount() {
			size_t kmerCount = 0;

			for (auto cluster : *this) {
				auto x = cluster->Size();
				kmerCount += x;
			}

			return kmerCount;
		}

		size_t Size() {
			return codebook.size();
		}

	private:
		void SelectPrototypesUniform(
			//------------------------------------
			UniformRealRandom &rand
			//------------------------------------
		) {
			size_t kmerCount = kmerIndex.size();

			Selector s(rand, (uint)codebook_size, (uint)kmerCount);

			fprintf(stderr, "codebook_size = %zu\n", codebook_size);
			fprintf(stderr, "kmerCount = %zu\n", kmerCount);
			fflush(stderr);

			size_t total = 0;

			for (auto &pair : kmerIndex) {
				if (s.SelectThis()) {
					auto newCluster = new Cluster(pair.second, 0, distanceFunction);
					newCluster->index = codebook.size();
					codebook.push_back(newCluster);
				}
			}

			if (total != kmerCount) {
				fprintf(stderr, "total = %zu\n", total);
				fprintf(stderr, "This does not match the required number. Exiting now.\n");
				exit(1);
			}

			if (codebook.size() != codebook_size) {
				fprintf(stderr, "codebook.size() = %zu\n", codebook.size());
				fprintf(stderr, "This does not match the required number. Exiting now.\n");
				exit(1);
			}
		}

		void SelectPrototypesLinear(
			//------------------------------------
			UniformRealRandom &rand
			//------------------------------------
		) {
			throw NotImplementedException(FileAndLine);
		}

		void SelectPrototypesLog(
			//------------------------------------
			UniformRealRandom &rand
			//------------------------------------
		) {
			throw NotImplementedException(FileAndLine);
		}

		void SplitLargeClusters() {
			int C = int(codebook.size());

			for (int i = 0; i < C; i++) {
				pCluster cluster = (pCluster)codebook[i];
				auto &kmers = cluster->kmers;
				size_t cSize = cluster->Size();

				if (cSize >= splitClusterThreshold) {
					int slices = int(ceil(cSize / splitClusterThreshold));
					int itemsPerSlice = int(cSize / slices);

					for (int j = 0; j < slices - 1; j++) {
						auto p(cluster->prototype);

						pCluster newCluster = new Cluster(p, itemsPerSlice, distanceFunction);

						auto &newItems = newCluster->kmers;

						for (int k = 0; k < itemsPerSlice; k++) {
							auto item = kmers.back();
							kmers.pop_back();
							newItems.push_back(item);
						}

						newCluster->index = codebook.size();
						codebook.push_back(newCluster);
					}
				}
			}
		}

		using pCluster = Cluster * ;

		static bool match(const char *l, const char *s) {
			return (strncmp(l, s, strlen(s)) == 0);
		};

		// Advances the current character until a comma or semicolon is reached.
		// Upon return, l[i] is a comma. If no comma or semicolon is reached before
		// end of string, exception is thrown.
		static void gotoSymbol(const char *l, size_t &i, char symbol) {
			while (l[i] && l[i] != symbol)
				i++;

			if (!l[i]) {
				ostringstream str;
				str << "Expected '" << symbol << "' in string '" << l << "'";
				throw Exception(str.str(), FileAndLine);
			}
		};

		// Skips the current character, which is assumed to be a punctuation mark of some sort,
		// and parses an integer from the following positions. On return,
		// l[i] is the first non-digit encountered.
		static int parseInt(const char *l, size_t &i) {
			int y = 0;
			bool isNeg = l[i + 1] == '-';

			if (isNeg)
				i++;

			for (i++; isdigit(l[i]); i++) {
				y = y * 10 + (l[i] - '0');
			}

			return isNeg ? -y : y;
		};

		void processCodebook(const char *l) {
			size_t i = 0;
			gotoSymbol(l, i, ',');
			uint k = parseInt(l, i);

			if (k != kmerLength) {
				throw Exception("Invalid kmerLength", FileAndLine);
			}
		}

		void processClusters(const char *l) {
			// We don't really care how many clusters are present.
			// The "Clusters" line serves only to tell us that we are reading a codebook.
		}

		void processLine(
			char *l,
			pCluster &currentCluster
		) {
			size_t i = 0;

			char * seqId = l + i;
			gotoSymbol(l, i, ':');
			l[i] = 0;

			Distance d = numeric_limits<Distance>::max();

			if (strcmp("distance", seqId) == 0) {
				d = parseInt(l, i);

				if (l[i] == ';') {
					i++;
					seqId = l + i;
					gotoSymbol(l, i, ':');
					l[i] = 0;
				}
			}

			auto seqItem = dbIndex.find(seqId);

			if (seqItem == dbIndex.end()) {
				throw Exception(string("sequence with Id ") + seqId + " cannot be found in database.", FileAndLine);
			}

			EncodedFastaSequence *seq = seqItem->second;
			uint offset = parseInt(l, i);

			Substring kmerText(seq->Sequence().data(), offset, kmerLength);

			try {
				auto kmer = kmerIndex.at(kmerText);

				currentCluster->Add(kmer, d);

				// if (l[i] == ';') i++;
				// The kmer should be identified by the text content of the first instance, so it
				// is only necessary to add the kmer once. All other instances should be connected
				// to the pointer obtained from the lookup table.
			}
			catch (out_of_range & ex) {
				ostringstream s;
				s << "kmer " << kmerText << " not found in index!";
				throw Exception(s.str(), FileAndLine);
			}
		}

		void processCluster(
			char *l,
			pCluster &currentCluster
		) {
			size_t len = strlen(l);
			size_t i = 0;
			gotoSymbol(l, i, ',');

			size_t expectedSize = parseInt(l, i);

			assert_equal(',', l[i]);
			i++;

			char *seqId = l + i;

			// cerr << "seqId = " << seqId << "\n";

			gotoSymbol(l, i, ':');

			l[i] = 0;

			auto item{ protoIndex.find(seqId) };

			KmerClusterPrototype *seq = 0;

			if (item != protoIndex.end()) {
				seq = (KmerClusterPrototype *)item->second;
			}
			else {
				ostringstream s;
				s << "Unable to find prototype " << seqId << " in the prototype index.";
				throw Exception(s.str(), FileAndLine);
			}

			uint offset = parseInt(l, i);

			while (l[i] != 0 && l[i] != ';') i++;

			if (false) {
				(cerr << "l: " << l
					<< "\noffset: " << offset
					<< "\nexpectedSize: " << expectedSize
					<< "\n")
					.flush();
			}

			K *prototype = new K(*seq, 0, kmerLength);
			currentCluster = new Cluster(prototype, expectedSize, distanceFunction);
			currentCluster->index = codebook.size();
			codebook.push_back(currentCluster);

			if (false) {
				(cerr << "currentCluster->prototype: " << currentCluster->prototype
					<< "\ncurrentCluster->expectedSize: " << currentCluster->expectedSize
					<< "\n")
					.flush();
			}

			while (i < len - 1) {
				i++;

				if (i >= len) break;

				char * key = l + i;
				gotoSymbol(l, i, ':');
				l[i++] = 0;

				if (i >= len) break;

				char * val = l + i;
				gotoSymbol(l, i, ';');
				l[i] = 0;

				currentCluster->AddMetadata(key, val);
			}
		}

		void CopyKmerData(int wordsPerKmer) {
			int N = (int)codebook.size();
			kmerData.resize(N, wordsPerKmer);

#pragma omp parallel for
			for (int i = 0; i < N; i++) {
				pCluster cluster = (pCluster)codebook[i];
				memcpy(kmerData.row(i), cluster->PackedEncoding(), wordsPerKmer * sizeof(KmerWord));
				cluster->CopyKmerData(wordsPerKmer);
			}
		}

		void CopyKmerDataNonRecursive(int wordsPerKmer) {
			int64_t N = (int64_t)codebook.size();
			kmerData.resize(N, wordsPerKmer);

#pragma omp parallel for
			for (int64_t i = 0; i < N; i++) {
				auto cluster = codebook[i];
				memcpy((KmerWord *)kmerData.row(i), cluster->prototype->PackedEncoding(), wordsPerKmer * sizeof(KmerWord));
			}
		}

		void AddWordsToClusters2(Distance threshold = 0) {
			size_t n = kmerIndex.size();
			vector<pKmer> allKmers;

			for (auto &pair : kmerIndex) {
				allKmers.push_back(pair.second);
			}

#if USE_OMP
			vector<omp_lock_t> clusterLock(codebook_size);

#pragma omp parallel for
			for (uint c = 0; c < codebook_size; c++) {
				omp_init_lock(&clusterLock[c]);
				this->codebook[c]->kmers.reserve(n / codebook_size);
			}
#else
			for (uint c = 0; c < codebook_size; c++) {
				this->codebook[c]->kmers.reserve(n / codebook_size);
			}
#endif

			auto N = (int64_t)n;

#if USE_OMP
#pragma omp parallel for schedule(dynamic, 10)
			for (int64_t i = 0; i < N; i++) {
				int threadId = omp_get_thread_num();
				pKmer kmer = allKmers[i];
				Distance distance = numeric_limits<Distance>::max();
				Cluster *closest = FindNearestCluster(*kmer, distance);

				if (threshold == 0 || distance <= threshold) {
					omp_set_lock(&clusterLock[closest->index]);
					closest->expectedSize++;
					closest->Add(*kmer, distance);
					omp_unset_lock(&clusterLock[closest->index]);
				}
			}
#else
			for (int64_t i = 0; i < N; i++) {
				int threadId = 0;
				pKmer kmer = allKmers[i];
				Distance dist = numeric_limits<Distance>::max();
				KmerCluster<D, K> *closest = FindNearestCluster(*kmer, dist);

				if (threshold == 0 || dist <= threshold) {
					// omp_set_lock(&clusterLock[closest->index]);
					closest->expectedSize++;
					closest->Add(*kmer, dist);
					// omp_unset_lock(&clusterLock[closest->index]);
				}
			}
#endif

#if USE_OMP
#pragma omp parallel for
			for (uint c = 0; c < codebook_size; c++) {
				omp_destroy_lock(&clusterLock[c]);
			}
#endif
		}

	public:
		vector<Cluster *> &Codebook() { return codebook; }

		Cluster *FindNearestCluster(K &kmer, Distance &dist) {
			Cluster *nearestCluster = (pCluster)codebook[0];
			size_t nearestIdx = 0;
			KmerWord *seqStr = kmer.PackedEncoding();

			dist = distanceFunction(seqStr, kmerData.row(0), kmerLength);

			for (size_t i = 1; i < codebook_size; i++) {
				Distance d = distanceFunction(seqStr, kmerData.row(i), kmerLength);

				if (d < dist) {
					dist = d;
					nearestCluster = (pCluster)codebook[i];
					nearestIdx = i;
				}
			}

			if (nearestIdx == 0 && dist < 100) {
#pragma omp critical
				{
					cerr << "nearest distance:\t" << dist << "\n";
					cerr << "prototype sequence:\t" << codebook[nearestIdx]->prototype.Substr() << "\n";
					cerr << "kmer sequence:\t" << kmer.Substr() << "\n";

					cerr << "prototype encoding:";
					for (int t = 0; t < (int)Alphabet::WordsRequiredToPack(kmerLength, charsPerWord); t++) {
						KmerWord x = kmerData.row(nearestIdx)[t];
						cerr << "\t" << x;
					}
					cerr << "\n";

					cerr << "kmer encoding:";
					for (int t = 0; t < (int)Alphabet::WordsRequiredToPack(kmerLength, charsPerWord); t++) {
						KmerWord y = seqStr[t];
						cerr << "\t" << y;
					}
					cerr << "\n";

					cerr.flush();
				}
			}

			return nearestCluster;
		}

		ostream &Write(ostream &out) const {
			out << "Codebook," << kmerLength << "," << codebook_size << endl;
			for (auto cluster : codebook) {
				out << cluster;
			}
			return out;
		}

		typename vector<Cluster *>::iterator begin() {
			return codebook.begin();
		}

		typename vector<Cluster *>::iterator end() {
			return codebook.end();
		}

		typename vector<Cluster *>::reverse_iterator rbegin() {
			return codebook.rbegin();
		}

		typename vector<Cluster *>::reverse_iterator rend() {
			return codebook.rend();
		}

	public:
		Cluster *ClusterAt(size_t pos) {
			return codebook[pos];
		}

		uint KmerLength() {
			return kmerLength;
		}

		// Assign all single-class sequences belonging to each cluster to the
		// cluster by inserting them into mapping seqs.

		void GetSequencesPerCluster(
			vector<SequenceAssignment<DistanceFunction, KmerType>> &seqs
			//
		) {
			const uint N = codebook.size();

#if USE_OMP
#pragma omp parallel for
#endif
			for (uint i = 0; i < N; i++) {
				pCluster cluster{ codebook[i] };
				seqs.emplace_back(cluster);
				SequenceSet &seqList{ seqs.back() };

				for (auto &kmer : cluster->kmers) {
					if (kmer.Sequence()->classNumbers.size() == 1) {
						seqList.insert(kmer.Sequence());
					}
				}
			}
		}
	};
} // namespace QutBio
