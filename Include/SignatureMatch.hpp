#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream>
#include <string>
#include <unordered_map>

#include "Exception.hpp"
#include "Signature.hpp"


namespace QutBio {
	class SignatureMatch {
		private: string subjectId;
		private: string subjectClass;
		private: int queryFragment;
		private: int subjectFragment;
		private: double similarity;

		public: static int FIELD_COUNT() { return 5; }

		public: SignatureMatch() : similarity(-1) {}

		public: SignatureMatch( vector<string> & record, int start ) {
			int i = start;
			this->subjectId = record[i++];
			this->subjectClass = record[i++];
			this->queryFragment = Int::Parse( record[i++] );
			this->subjectFragment = Int::Parse( record[i++] );
			this->similarity = Double::Parse( record[i++] );
		}

		public: SignatureMatch(
			string subjectId,
			string subjectClass,
			int queryFragment,
			int subjectFragment,
			double similarity
			) {
			this->subjectId = subjectId;
			this->subjectClass = subjectClass;
			this->queryFragment = queryFragment;
			this->subjectFragment = subjectFragment;
			this->similarity = similarity;
		}

		public: vector<string> ToStringArray() {
			return vector<string>( {
				subjectId,
				subjectClass,
				Int::ToString( queryFragment ),
				Int::ToString( subjectFragment ),
				Double::ToString( similarity ),
			} );
		}

				/// <summary> Returns the similarity between two lists of signatures.
				/// <para>
				///		For lists, similarity is the maximum pairwise similarity between 
				///		members of the lists.
				/// </para>
				/// </summary>
				/// <param name="other"></param>
				/// <returns></returns>

		public: static SignatureMatch BestMatch(
			vector<Signature> & querySignatures,
			vector<Signature> & subjectSignatures
			) {
			double bestSimilarity = -1e100;
			int bestQueryFragment = -1;
			int bestSubjectFragment = -1;

			for (auto querySignature_ = querySignatures.begin(); querySignature_ != querySignatures.end(); querySignature_++) {
				Signature & querySignature( *querySignature_ );

				for (auto subjectSignature_ = subjectSignatures.begin(); subjectSignature_ != subjectSignatures.end(); subjectSignature_++) {
					Signature & subjectSignature( *subjectSignature_ );

					double similarity = querySignature.Similarity( subjectSignature );

					// cerr << querySignature.FragmentIndex() << ", " << subjectSignature.FragmentIndex() << " --> " << similarity << endl;

					if (similarity > bestSimilarity) {
						bestSimilarity = similarity;
						bestQueryFragment = querySignature.FragmentIndex();
						bestSubjectFragment = subjectSignature.FragmentIndex();
					}
				}
			}

			return SignatureMatch(
				subjectSignatures[0].GI(),
				subjectSignatures[0].ClassLabel(),
				bestQueryFragment,
				bestSubjectFragment,
				bestSimilarity
			);
		}

		public: string & SubjectId() { return subjectId; }
		public: string & SubjectClass() { return subjectClass; }
		public: int QueryFragment() { return queryFragment; }
		public: int SubjectFragment() { return subjectFragment; }
		public: double Similarity() { return similarity; }
	};
}
