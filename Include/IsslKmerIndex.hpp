#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Kmer.hpp"
#include "AllocatedKmer.hpp"
#include "KmerIndex.hpp"
#include "PackedArray.hpp"

namespace QutBio {
	/**
	*	ISSL index of kmers in a dataset. This is not to be confused with a
	*	Sequence Signature index, but it does have some resemblance. Instead of
	*	signatures, the sequences are stored in packed arrays. DNA takes 2 bits per
	*	symbol; AA takes 5 bits per symbol. The structure strives to support fast lookup
	*	based on kmer identity.
	*
	*	Version 1 will use "naive" object oriented techniques to establish proof of concept.
	*	Version 2 will optimise (e.g. replace list of packed arrays with a matrix representation).
	*	Version 3 will support serialisation and de-serialisation of the index structure, and
	*		possibly move the bit-packing into the Substr object rather than doing it as an afterthought.
	*/

	class IsslKmerIndex : public KmerIndex {
		size_t slices;
		vector<PackedArray<char>> allSignatures;
		vector<vector<vector<Kmer *>>> issl;

	public:
		/**
		*	Maintain the class invariant:
		*		allKmers.length = allSignatures.length
		*		for all i: dom(allKmers) -> allSignatures[i] = PackedArray(allKmers[i],bitsPerChar);
		*
		*	@param db			A list of sequences which are to be indexed.
		*	@param kmerLength	The size of kmers to be extracted.
		*	@param bitsPerChar	The number of bits used to pack characters.
		*	@param inverse	    Mapping from char (e.g. "acgt") to 0-origin byte values.
		*	@param slices		The number of slices to use for
		*/
		IsslKmerIndex(
			vector<pEncodedFastaSequence> db,
			size_t kmerLength,
			size_t bitsPerChar,
			const uint8_t inverse[],
			size_t slices
		) :
			KmerIndex(db, kmerLength),
			slices(slices)
			//
		{
			// Tell superclass to instantiate the list of Kmers. 
			GetKmers();

			cerr << this->allKmers.size() << " kmers indexed.\n";

			Pack(bitsPerChar, inverse);
			BuildIndex();

			if (false) {
				uint slice = 0;

				for (auto & sliceVec : issl) {
					uint sliceVal = 0;

					for (auto & sliceValSet : sliceVec) {
						cerr << slice << "," << sliceVal << "," << sliceValSet.size() << "\n";
						sliceVal++;
					}

					slice++;
				}
			}
		}

	private:
		void Pack(
			size_t bitsPerChar,
			const uint8_t inverse[]
		) {
			// cerr << "bitsPerChar:" << bitsPerChar << "\n";

			for (auto & kmer : *this) {
				// PackedArray<char> sig(kmer.Substr(), bitsPerChar, inverse);
				// allSignatures.push_back(sig);
				allSignatures.emplace_back(kmer.second->Substr(), bitsPerChar, inverse);
			}
		}

		void BuildIndex() {
			const auto & firstSig = allSignatures.front();
			const size_t sigBits = firstSig.BitsPerItem() * firstSig.Length();

			issl.resize(slices);

			for (uint slice = 0; slice < slices; slice++) {
				size_t bitsPerSlice = (slice + 1) * sigBits / slices - slice * sigBits / slices;
				size_t sliceVals = 1 << bitsPerSlice;
				issl[slice].resize(sliceVals);
			}

			for (uint i = 0; i < allSignatures.size(); i++) {
				auto & sig = allSignatures[i];

				for (uint slice = 0; slice < slices; slice++) {
					size_t sliceValue = sig.GetSlice(slice, slices);
					issl[slice][sliceValue].push_back(this->allKmers[i]);
				}
			}
		}

	public:
		vector<PackedArray<char>> & AllSignatures() {
			return allSignatures;
		}

		using BaseType = KmerIndex;

		size_t size() {
			return BaseType::size();
		}

		vector<Kmer *> & GetKmers() {
			return BaseType::GetKmers();
		}

		PackedArray<char> & Signature(size_t i) {
			return allSignatures[i];
		}

		size_t Slices() {
			return slices;
		}

		vector<vector<vector<Kmer *>>> & Issl() {
			return issl;
		}
	};

}
