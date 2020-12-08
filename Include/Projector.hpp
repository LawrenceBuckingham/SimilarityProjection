#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Alphabet.hpp"
#include "SimilarityMatrix.hpp"
#include "FastaSequence.hpp"
#include "FreeList.hpp"
#include "Util.hpp"
#include "SequenceDistanceFunction.hpp"
#include "DiagonalGenerator.hpp"

#include <string>
#include <ctime>
#include <cfloat>

namespace QutBio {

	class Projector : public SequenceDistanceFunction {
	protected:
		const FragmentAggregationMode * fragMode;

	public:
		Projector(
			SimilarityMatrix * matrix,
			uint kmerLength,
			const FragmentAggregationMode * fragMode
		) :
			SequenceDistanceFunction(matrix, kmerLength),
			fragMode(fragMode) 
		{
			
		}

		virtual ~Projector() {}

		virtual double ComputeDistance(
			FastaSequence & querySeq,
			FastaSequence & subjectSeq
		) {
			// TODO: padding is no good if alphabet is not Amino Acid.
			auto padding = Alphabets::AA()->Encode('x');
			querySeq.Pad(kmerLength, padding);
			size_t queryKmerCount = querySeq.KmerCount(kmerLength);

			subjectSeq.Pad(kmerLength, padding);
			size_t subjectKmerCount = subjectSeq.KmerCount(kmerLength);

			vector<Distance> rowMinima(queryKmerCount);
			vector<Distance> colMinima(subjectKmerCount);

			ComputeDistanceMatrix(
				// querySeq, subjectSeq,
				querySeq.Sequence().data(),
				subjectSeq.Sequence().data(),
				rowMinima,
				colMinima,
				queryKmerCount,
				subjectKmerCount
			);

			double distance = GetSequenceDistance(rowMinima.data(), colMinima.data(), queryKmerCount, subjectKmerCount);

			return distance;
		}

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>

	protected:
		typedef struct {
			Distance * rowMinima;
			Distance * colMinima;
		} DataPacket;

		static void ProcessDistance(
			void * processingObject,
			size_t queryPos,
			size_t subjectPos,
			Distance distance
		) {
			DataPacket * data = (DataPacket *)processingObject;
			Distance * rowMinima = data->rowMinima;
			Distance * colMinima = data->colMinima;
			if (distance < rowMinima[queryPos]) rowMinima[queryPos] = distance;
			if (distance < colMinima[subjectPos]) colMinima[subjectPos] = distance;
		}

	public: virtual void ComputeDistanceMatrix(
		// EncodedFastaSequence &querySeq,
		// EncodedFastaSequence &subjectSeq,
		const Symbol* queryChars,
		const Symbol* subjectChars,
		vector<Distance> & rowMinima,
		vector<Distance> & colMinima,
		size_t queryKmerCount,
		size_t subjectKmerCount
	) {
		DataPacket data{ rowMinima.data(), colMinima.data() };

		Util::Fill(rowMinima, numeric_limits<Distance>::max());
		Util::Fill(colMinima, numeric_limits<Distance>::max());

		DiagonalGenerator generator;

		generator.GenerateDistances(
			// &querySeq, &subjectSeq, 
			queryChars,
			subjectChars,
			kmerLength,
			queryKmerCount,
			subjectKmerCount,
			distanceLookup,
			&data,
			ProcessDistance
		);
	}

	protected: virtual double GetSequenceDistance(
		Distance * rowMinima,
		Distance * colMinima,
		size_t queryKmerCount,
		size_t subjectKmerCount
	) {
		int rowTotal = 0;

		for (uint i = 0; i < queryKmerCount; i++) {
			rowTotal += rowMinima[i];
		}

		int colTotal = 0;

		for (uint i = 0; i < subjectKmerCount; i++) {
			colTotal += colMinima[i];
		}

		if (fragMode == FragmentAggregationMode::HausdorffAverageAverage()) {
			return ((double)rowTotal / queryKmerCount + (double)colTotal / subjectKmerCount) / 2;
		}
		else if (fragMode == FragmentAggregationMode::HausdorffAverage()) {
			static bool deja = false;

			if (!deja) {
				deja = true;
				cerr << "FragmentAggregationMode::HausdorffAverage()" << endl;
			}

			return max((double)rowTotal / queryKmerCount, (double)colTotal / subjectKmerCount);
		}
		else if (fragMode == FragmentAggregationMode::Hausdorff()) {
			double rowMax = 0;
			double colMax = 0;

			for (uint i = 0; i < queryKmerCount; i++) {
				if (rowMinima[i] > rowMax) {
					rowMax = rowMinima[i];
				}
			}

			for (uint i = 0; i < subjectKmerCount; i++) {
				if (colMinima[i] > colMax) {
					colMax = colMinima[i];
				}
			}

			return max(rowMax, colMax);
		}
		else if (fragMode == FragmentAggregationMode::BestOfBest()) {
			double min = HUGE_VAL;

			for (uint i = 0; i < queryKmerCount; i++) {
				if (rowMinima[i] < min) {
					min = rowMinima[i];
				}
			}

			for (uint i = 0; i < subjectKmerCount; i++) {
				if (colMinima[i] < min) {
					min = colMinima[i];
				}
			}

			return min;
		}
		else {
			throw NotImplementedException(FileAndLine);
		}
	}
	};
}
