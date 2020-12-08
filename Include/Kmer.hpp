#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <ios>
#include <unordered_map>
#include <unordered_set>

#include "EncodedFastaSequence.hpp"
#include "EncodedKmer.hpp"
#include "Substring.hpp"
#include "DistanceType.hpp"

namespace QutBio {

	using pKmer = class Kmer *;

	class Kmer {
	public:
		/*
		** <summary>The location of an instance of the Kmer.</summary>
		*/
		typedef struct Instance {
			// The Id of the containing sequence.
			EncodedFastaSequence & sequence;

			// The offset of the first character of the kmer from the start of the sequence.
			size_t kmerPosition;

			Instance(EncodedFastaSequence & sequence, size_t kmerPosition) : sequence(sequence), kmerPosition(kmerPosition) {}

			friend ostream & operator<<(ostream & str, const Instance & instance) {
				str << instance.sequence.IdStr() << ":" << instance.kmerPosition;
				return str;
			}

			EncodedKmer PackedEncoding() const {
				return sequence.GetEncodedKmer(kmerPosition);
			}

			EncodedKmer UnpackedEncoding() const {
				return sequence.GetEncodedKmer1(kmerPosition);
			}

			static Instance & Zero() {
				static Instance zero(EncodedFastaSequence::Zero(), 0);
				return zero;
			}

			friend bool operator==(const Instance & lhs, const Instance & rhs) {
				return &(lhs.sequence) == &(rhs.sequence) && lhs.kmerPosition == rhs.kmerPosition;
			}

			friend bool operator!=(const Instance & lhs, const Instance & rhs) {
				return &(lhs.sequence) != &(rhs.sequence) || lhs.kmerPosition != rhs.kmerPosition;
			}
		} Instance;

	private:
		// The string representation of the kmer.
		Substring substring;

		// The locations of all instances of the current kmer.
		vector<Instance> instances;

		size_t serialNumber = 0;

	public:
		/**
		**	Summary:
		**		Construct a kmer belonging to a sequence.
		**
		**	Parameters:
		**		seq:          a sequence containing a prototypical instance of the kmer;
		**		kmerPosition: the offset of the kmer from the start of the sequence;
		**		kmerLength:   the length of the kmer;
		**		dist:         a distance value associated with the kmer, used when clustering
		**			          Arguably this should not be present in the base class.
		*/
		Kmer(
			EncodedFastaSequence & seq,
			size_t kmerPosition,
			size_t kmerLength
		) :
			substring(seq.Sequence().data(), kmerPosition, kmerLength)
			//
		{
			Add(seq, kmerPosition);
		}

		virtual ~Kmer() {

		}

		/**
		**	Summary:
		**		Adds a new instance to the current kmer.
		**
		**	Parameters:
		**		seq:          a sequence containing a prototypical instance of the kmer;
		**		kmerPosition: the offset of the kmer from the start of the sequence;
		**		kmerLength:   the length of the kmer;
		**		dist:         a distance value associated with the kmer, used when clustering
		**			          Arguably this should not be present in the base class.
		*/
		void Add(EncodedFastaSequence & seq, size_t kmerPosition) {
			instances.emplace_back(seq, kmerPosition);
		}

		void Add(vector<Instance> & other) {
			for (auto & i : other) {
				Add(i.sequence, i.kmerPosition);
			}
		}

		const Substring & Substr() const {
			return substring;
		}

		/// <summary>
		/// Returns a vector of byte containing a copy of the kmer.
		/// </summary>

		vector<Symbol> Word() const {
			vector<Symbol> ret(substring.Chars(), substring.Chars() + substring.Length());
			return ret;
		}

		const vector<Instance> & Instances() const {
			return instances;
		}

		// Gets the address of the first word in the packed numerically encoded kmer
		// array.
		//	*	This is currently a unit16_t array containing 1, 2, or 3 symbols
		//		packed into a single number.
		//	*	If I switch to a 64-bit dual embedding this will become a uint64_t
		//		array.
		//	*	I may have to keep both formats alive side-by-side for a while.

		EncodedKmer PackedEncoding() const {
			return instances.size() == 0 ? 0 : instances[0].PackedEncoding();
		}

		// Gets the address of the first word in the unpacked numerically encoded kmer
		// array.
		//	*	This is a unit16_t array containing 1 symbol in each element of the
		//		EncodedKmer.

		EncodedKmer UnpackedEncoding() {
			return instances.size() == 0 ? 0 : instances[0].UnpackedEncoding();
		}

		friend ostream & operator << (ostream & str, const Kmer & kmer) {
			for (auto & instance : kmer.instances) {
				str << instance << ";";
			}

			return str;
		}

		Kmer & operator=(const Kmer & other) = delete;

		friend bool operator==(const Kmer & lhs, const Kmer & rhs) {
			return lhs.substring == rhs.substring;
		}

		friend bool operator!=(const Kmer & lhs, const Kmer & rhs) {
			return lhs.substring != rhs.substring;
		}

		friend bool operator<(const Kmer & lhs, const Kmer & rhs) {
			return lhs.substring < rhs.substring;
		}

		ostream & write(ostream & stream) const {
			return stream << substring;
		}

		/**
		 *	<summary>
		 *		Returns the number of kmers required to tile the longest sequence in a dataset.
		 *	</summary>
		 *	<param name="db">A list of sequences.</param>
		 *	<param name="kmerLength">The kmer length.</param>
		 */
		static size_t GetMaxKmerCount(vector<EncodedFastaSequence *> & db, size_t kmerLength) {
			size_t maxKmerCount = 0;

			for (auto seq : db) {
				size_t K = seq->KmerCount(kmerLength);

				if (K > maxKmerCount) {
					maxKmerCount = K;
				}
			}

			return maxKmerCount;
		}

		EncodedFastaSequence & Sequence() {
			return instances.size() > 0 ? instances[0].sequence : EncodedFastaSequence::Zero();
		}

		size_t KmerPosition() { return instances.size() > 0 ? instances[0].kmerPosition : 0; }

		Instance & FirstInstance() { return instances.size() > 0 ? instances[0] : Instance::Zero(); }

		size_t Length() const { return substring.Length(); }

		size_t SerialNumber() const { return serialNumber; }

		Kmer & SetSerialNumber(size_t value) { serialNumber = value;  return *this; }
	};
}

