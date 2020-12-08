#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cmath>
#include <functional>

#include "Alphabet.hpp"
#include "Array.hpp"
#include "EncodedKmer.hpp"
#include "Exception.hpp"
#include "Types.hpp"

using namespace QutBio;

namespace QutBio {

	/// <summary> A call-back function used to process document fragments.
	/// </summary>
	/// <param name="fragmentIndex">
	///		The zero-origin index number of the current fragment within the seq.
	/// </param>
	/// <param name="fragmentCount">
	///		The number of fragments in the document.
	/// </param>
	/// <param name="subSequence">
	///		The text content of the current fragment.
	/// </param>

	using FragmentProcessor = function<void(
		size_t fragmentIndex,
		size_t fragmentCount,
		const char * subSequence,
		size_t len
		)>;

	/// <summary> A call-back function used to process sub-sequences.
	/// </summary>
	/// <param name="offset">
	///		The zero-origin index of the current fragment within the seq.
	/// </param>
	/// <param name="fragmentIndex">
	///		The zero origin index number of the subsequence as a fragment of 
	///		the containing sequence.
	/// </param>
	/// <param name="length">
	///		The length of the subsequence.
	/// </param>
	/// <param name="fragmentCount">
	///		The number of fragments in the document.
	/// </param>
	/// <param name="sequence">
	///		The sequence that contains the current fragment.
	/// </param>

	template <typename T>
	using SubsequenceProcessor = function<void(
		vector<T> & sequence,
		int offset,
		int length,
		int fragmentIndex,
		int fragmentCount
		)>;

	/// <summary> A set of functions used to facilitate document fragmentation.
	/// </summary>

	class Fragment {
	public:
		/// <summary>
		/// 
		/// </summary>
		/// <param name="seq">
		///		The sequence to be processed.
		/// </param>
		/// <param name="fragmentLength">
		///		The length of each fragment.
		/// </param>
		/// <param name="stepSize">
		///		The number of characters between fragment start points.
		/// </param>
		/// <param name="minLength">
		///		The minimum permitted fragment size.
		/// </param>
		/// <param name="process">
		///		A procedure to be carried out once for 
		/// </param>

		static void DissectString(
			const char * seq,
			size_t seqLen,
			uint fragmentLength,
			uint stepSize,
			uint minLength,
			FragmentProcessor process
		) {
			if (fragmentLength < minLength) {
				throw new Exception("Argument Exception: fragmentLength", __FILE__, __LINE__);
			}

			if (stepSize <= 0) {
				throw new Exception("Argument Exception: stepSize", __FILE__, __LINE__);
			}

			if (minLength <= 0) {
				throw new Exception("Argument Exception: minLength", __FILE__, __LINE__);
			}

			size_t fragmentCount = GetCount(seqLen, fragmentLength);

			for (uint fragmentIndex = 0; fragmentIndex < fragmentCount; fragmentIndex++) {
				size_t start = fragmentIndex * stepSize;
				process(fragmentIndex, fragmentCount, seq + start, fragmentLength);
			}
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="seq">
		///		The sequence to be processed.
		/// </param>
		/// <param name="fragmentLength">
		///		The length of each fragment.
		/// </param>
		/// <param name="stepSize">
		///		The number of characters between fragment start points.
		/// </param>
		/// <param name="minLength">
		///		The minimum permitted fragment size.
		/// </param>
		/// <param name="process">
		///		A procedure to be carried out once for 
		/// </param>
		/// <typeparam name="T">
		///		The element type of the sequence.
		/// </typeparam>

		template<typename T>
		static void DissectSequence(
			vector<T> & seq,
			int fragmentLength,
			int stepSize,
			int minLength,
			SubsequenceProcessor<T> process
		) {
			if (fragmentLength < minLength) {
				throw new Exception("Argument Exception: fragmentLength", __FILE__, __LINE__);
			}

			if (stepSize <= 0) {
				throw new Exception("Argument Exception: stepSize", __FILE__, __LINE__);
			}

			if (minLength <= 0) {
				throw new Exception("Argument Exception: minLength", __FILE__, __LINE__);
			}

			int fragmentCount = GetCount(seq.Length, fragmentLength, stepSize, minLength);

			for (int fragmentIndex = 0; fragmentIndex < fragmentCount; fragmentIndex++) {
				int start = fragmentIndex * stepSize;
				process(seq, start, fragmentLength, fragmentIndex, fragmentLength);
			}
		}

		/**
		 *	Gets the minimum number of fragments required to fully cover a sequence with {seqLength}
		 *	elements by sub-sequences containing {fragmentLength} contiguous elements.
		 *
		 *	Notes:
		 *	*	You can use more than this if you want to have overlapping fragments.
		 *	*	This scheme ensures that no part of the sequence will be completely omitted, that all
		 *		fragments have the same number of elements. However, fragments will inevitably overlap,
		 *		especially in short sequences, so that some elements (those in the overlapping regions)
		 *		will have greater weight than others. I find this to be an acceptable compromise, as
		 *		it tends to smooth the data rather than make it more jagged.
		 *	Returns:
		 *		ceil({seqLength}�{fragLength}) : {fragLength} < {seqLength}
		 *		1                              : otherwise.
		 */

		static size_t GetCount(size_t seqLength, size_t fragmentLength) {
			uint fragmentCount = 0;

			if (fragmentLength < seqLength) {
				fragmentCount = (uint)ceil((double)seqLength / fragmentLength);
			}
			else {
				fragmentCount = 1;
			}

			return fragmentCount;
		}

		/**
		 *	Gets the real-valued step size required to fully cover a sequence with {seqLength}
		 *	elements by {fragCount} sub-sequences containing {fragmentLength} contiguous elements.
		 *
		 *	Notes:
		 *	*	This scheme ensures that no part of the sequence will be completely omitted, that all
		 *		fragments have the same number of elements. However, fragments will inevitably overlap,
		 *		especially in short sequences, so that some elements (those in the overlapping regions)
		 *		will have greater weight than others. I find this to be an acceptable compromise, as
		 *		it tends to smooth the data rather than make it more jagged.
		 *
		 *	Returns:
		 *		1                           : {fragCount} = 1
		 *		{seqLength} � {fragCount}     : otherwise
		 */

		static double GetRealStepSize(size_t seqLength, size_t fragLength, size_t fragCount) {
			return fragCount == 1 ? fragLength : (double)seqLength / fragCount;
		}

		/**
		 *	Gets the integer-valued starting point for step {idx}, of nominal real-value length 
		 *	{stepSize}, in a sequence of {kmerCount} items.
		 *
		 *	Notes:
		 *	*	This scheme ensures that no part of the sequence will be completely omitted, that all
		 *		fragments have the same number of elements. However, fragments will inevitably overlap,
		 *		especially in short sequences, so that some elements (those in the overlapping regions)
		 *		will have greater weight than others. I find this to be an acceptable compromise, as
		 *		it tends to smooth the data rather than make it more jagged.
		 *
		 *	Returns:
		 *		min(round({idx} � {stepSize}), {kmerCount}).
		 */

		static uint GetFragmentStart(uint idx, double stepSize, size_t kmerCount) {
			return min((uint)round(idx * stepSize), (uint)kmerCount);
		}

		/**
		 *	Gets a list containing the successive start locations of all steps for a fragment-covering
		 *	of sequence with {kmerCount} items, using fragments of length {fragmentLength}.
		 *
		 *	Returns a vector v such that:
		 *		v[i] = GetFragmentStart(i, stepSize, {kmerCount}), where
		 *		stepSize = GetRealStepSize({kmerCount}, {fragmnetLength}, fragCount), with
		 *		fragCount = GetCount({kmerCount}, {fragmentLength}).
		 */

		static vector<uint> GetFragmentStartList( uint kmerCount, uint fragmentLength ) {
			auto fragCount = GetCount(kmerCount, fragmentLength);
			auto stepSize = GetRealStepSize(kmerCount, fragmentLength, fragCount );
			uint start = 0;
			vector<uint> fragStarts;

			for (uint i = 0; i < fragCount; i++ ) {
				fragStarts.push_back(start);
				start = Fragment::GetFragmentStart(i + 1, stepSize, kmerCount);
			}

			return fragStarts;
		}

		/**
		 *	Traverses a grid of paired indices corresponding to the starting points of fragment
		 *	coverings of a pair of sequences, invoking a callback function at each point in the grid,
		 *	and another at the end of each row.
		 */

		static void PartitionSequencePair(
			uint fragmentLength,
			uint queryKmerCount,
			uint queryFragCount,
			uint subjectKmerCount,
			uint subjectFragCount,
			function<void(uint qFragIdx, uint qStart, uint qEnd, uint sFragIdx, uint sStart, uint sEnd)> processCell,
			function<void(uint qFragIdx)> processEndOfRow
		) {
			double qStepSize = Fragment::GetRealStepSize(queryKmerCount, fragmentLength, queryFragCount);
			double sStepSize = Fragment::GetRealStepSize(subjectKmerCount, fragmentLength, subjectFragCount);

			uint qStart = 0;

			for (uint i = 0; i < queryFragCount; i++) {
				uint qEnd = Fragment::GetFragmentStart(i + 1, qStepSize, queryKmerCount);
				uint sStart = 0;

				for (uint j = 0; j < subjectFragCount; j++) {
					uint sEnd = Fragment::GetFragmentStart(j + 1, sStepSize, subjectKmerCount);
					processCell( i, qStart, qEnd, j, sStart, sEnd);
					sStart = sEnd;
				}

				qStart = qEnd;
				processEndOfRow(i);
			}
		}


#ifdef WANT_ENCODED_KMER
		/// <summary> Partitions a sequence, and tiles the resulting sub-sequences into kmers.
		/// </summary>
		/// <param name="sequence"></param>
		/// <param name="fragLength"></param>
		/// <param name="fragInterval"></param>
		/// <param name="kmerLength"></param>
		/// <param name="charsPerWord">The number of characters packed into each word of the encoded kmer.</param>
		/// <param name="fragList">A preallocated array that will contain the fragments.</param>
		/// <param name="tileList">A preallocated two-dimensional array that will contain the collection of kmer tilings corresponding to the fragments.</param>
		/// <param name="alphabet">The alphabet used to encode kmers.</param>

		static void SplitAndTile(
			const char * sequence,
			size_t seqLen,
			Alphabet * alphabet,
			int fragLength,
			int fragInterval,
			int kmerLength,
			int charsPerWord,
			Array<const char *> & fragList,
			Array<Array<EncodedKmer>> & tileList

		) {

#pragma warning( disable:4100 )
			DissectString(sequence, seqLen, fragLength, fragInterval, kmerLength,
				[&fragList](int fragmentIndex, int fragmentCount, const char * subSequence, size_t len) {
				fragList[fragmentIndex] = subSequence;
			}
			);

			for (uint i = 0; i < fragList.Length(); i++) {
				Array<EncodedKmer> & tiles = tileList[i];
				DissectString(fragList[i], fragLength, kmerLength, 1, 1,
					[&tiles, charsPerWord, alphabet](int fragmentIndex, int fragmentCount, const char * subSequence, size_t len) {
					alphabet->Encode(subSequence, len, charsPerWord, tiles[fragmentIndex]);
				}
				);
			}
		}

#endif // WANT_ENCODED_KMER

	};

}
