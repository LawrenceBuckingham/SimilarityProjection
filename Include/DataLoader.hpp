#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <mutex>
#include <vector>
#include <map>
#include <iostream>

#include <FastaSequence.hpp>
#include <EncodedFastaSequence.hpp>
#include <KmerClusterPrototype.hpp>
#include <KmerCodebook.hpp>
#include <Domain.hpp>
#include <KmerCodebook.hpp>
#include <Homologs.hpp>
#include "KmerIndex.hpp"

using namespace std;

namespace QutBio {
	/**
	 *	Convenience functions to load the most commonly used data
	 *	collections.
	 */
	struct Load {
		// TODO: Convert this to uniq_ptr
		static vector<FastaSequence *> Fasta(
			const string & fileName,
			size_t idIndex,
			Alphabet * alphabet
		) {
			return FastaSequence::Read(fileName, idIndex, alphabet);
		}

		static vector<FastaSequence *> Fasta(
			istream & reader,
			size_t idIndex,
			Alphabet * alphabet
		)
		{
			return FastaSequence::Read(reader, idIndex, alphabet);
		}

		static void Fasta(
			istream & reader,
			size_t idIndex,
			Alphabet * alphabet,
			vector<FastaSequence *> &result
		)
		{
			FastaSequence::Read(reader, idIndex, alphabet, result);
		}

		static vector<EncodedFastaSequence *> Encoded(
			vector<FastaSequence *> db,
			int classIndex,
			Alphabet * alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			Symbol defaultSymbol
		) {
			vector<EncodedFastaSequence *> encoded;

			for (auto seq : db) {
				encoded.push_back(new EncodedFastaSequence(seq, classIndex, alphabet, kmerLength, charsPerWord, defaultSymbol));
			}

			return encoded;
		}

		static vector<KmerClusterPrototype *> Prototypes(
			vector<FastaSequence *> & db,
			Alphabet * alphabet,
			size_t kmerLength,
			size_t charsPerWord
		) {
			vector<KmerClusterPrototype *> encoded;

			for (auto seq : db) {
				auto p = new KmerClusterPrototype(seq, -1, alphabet, kmerLength, charsPerWord, alphabet->DefaultSymbol());
				encoded.push_back(p);
			}

			return encoded;
		}

		template<typename T>
		static Index<T> GetSequenceIndex(vector<T *> collection) {
			Index<T> idx(collection);
			return idx;
		}

		template<typename T>
		static KmerIndex GetKmerIndex(vector<T *> collection, size_t k) {
			KmerIndex idx(collection, k);
			return idx;
		}

		static map<string, Domain> Domains(const string & domFileName) {
			map<string, Domain> domains;
			ifstream domFile(domFileName);
			Domain::Load(domFile, domains);
			domFile.close();
			return domains;
		}

		template <typename DistanceFunction>
		static KmerCodebook<DistanceFunction, Kmer> * Codebook(
			const string & inFile,
			Alphabet * alphabet,
			DistanceFunction & distanceFunction,
			size_t kmerLength,
			size_t charsPerWord,
			Index<EncodedFastaSequence> & dbIndex,
			Index<KmerClusterPrototype> & protoIndex,
			KmerIndex &kmerIndex
		) {
			KmerCodebook<DistanceFunction, Kmer> *codebook = 0;

			FILE * f = fopen(inFile.c_str(), "r");

			if (f) {
				codebook = new KmerCodebook<DistanceFunction, Kmer>(
					alphabet, distanceFunction, charsPerWord, kmerLength, dbIndex,
					protoIndex, kmerIndex, f
					);

				fclose(f);
			}
			else {
				ostringstream err;
				err << "Unable to open codebook file '" << inFile << "'\n";
				throw Exception(err.str(), FileAndLine);
			}

			return codebook;
		}

		static map<string, set<string>> Homologs(
			const string & fileName,
			char separator = ' '
		) {
			ifstream inFile(fileName);
			return Homologs::Parse(inFile, separator);
		}
	};
}
