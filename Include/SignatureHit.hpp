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
#include "SignatureMatch.hpp"


namespace QutBio {
	class SignatureHit {
		private: string queryId;
		private: string queryClass;
		private: int fragmentLength;
		private: int interval;
		private: vector<SignatureMatch> matches;

		public: string & QueryId() { return queryId; }

		private: void  QueryId( string & value ) { queryId = value; }

		public: string & QueryClass() { return queryClass; }

		private: void  QueryClass( string & value ) { queryClass = value; }

		public: int FragmentLength() { return fragmentLength; }

		private: void FragmentLength( int value ) {
			if ( value <= 0 ) throw Exception( "value must be strictly positive.", __FILE__, __LINE__ );

			fragmentLength = value;
		}

		public: int Interval() { return interval; }

		private: void  Interval( int value ) {
			if ( value <= 0 ) throw Exception( "value must be strictly positive.", __FILE__, __LINE__ );
			
			interval = value;
		}

		public: vector<SignatureMatch> & Matches() { return matches; }

		private: void Matches( vector<SignatureMatch> & value ) { matches = value; }

		public: SignatureHit( vector<string> & record ) {
			int j = 0;
			this->queryId = record[j++];
			this->queryClass = record[j++];
			this->fragmentLength = Int::Parse( record[j++] );
			this->interval = Int::Parse( record[j++] );

			int matchCount = Int::Parse( record[j++] );
			matches.resize( matchCount );

			int fieldCount = SignatureMatch::FIELD_COUNT();

			for (int i = 0; i < matchCount; i++) {
				int start = 5 + i * fieldCount;
				SignatureMatch match( record, start );
				matches.push_back( match );
			}
		}

		public: SignatureHit(
			string & queryId,
			string & queryClass,
			int fragmentLength,
			int interval
			) {
			this->queryId = queryId;
			this->queryClass = queryClass;
			this->fragmentLength = fragmentLength;
			this->interval = interval;
		}

		public: SignatureHit( Signature & query, vector<Signature> & hits ) :
			queryId( query.GI() ),
			queryClass( query.ClassLabel() ),
			fragmentLength( query.FragmentLength() ),
			interval( query.Interval() )

		{
			for (auto hitIter = hits.begin(); hitIter != hits.end(); hitIter++) {
				Signature & subject( *hitIter );
				// TODO: this is wrong, I don't know what was going on here. Cross reference against the C# version before using this code.
				throw Exception( "TODO: this is wrong, I don't know what was going on here. Cross reference again the C# version before using this code.", __FILE__, __LINE__ );
				SignatureMatch match( subject.GI(), subject.ClassLabel(), 0, 0, 0 );
				matches.push_back( match );
			}
		}

		public: static vector<string> Headings() {
			string fields[] = {
				"QueryId",
				"QueryClass",
				"FragmentLength",
				"Interval",
				"HitCount",
			};
			vector<string> result( fields, fields + 5 );
			return result;
		}

		public: vector<string> ToStringArray() {
			vector<string> result;
			result.push_back( queryId );
			result.push_back( queryClass );
			result.push_back( Int::ToString( fragmentLength ) );
			result.push_back( Int::ToString( interval ) );
			result.push_back( Int::ToString( (int) matches.size() ) );

			for (auto matchIter = matches.begin(); matchIter != matches.end(); matchIter++) {
				SignatureMatch &match( *matchIter );
				vector<string> record = match.ToStringArray();
				result.insert( result.end(), record.begin(), record.end() );
			}

			return result;
		}


				/// <summary> Incrementally populates the matches array by adding the 
				/// <para>
				///		Inserts a match representing the similarity between two lists of 
				///		signatures which respectively tile a single query and single subject sequence.
				/// </para>
				/// </summary>
				/// <param name="k">
				///		The number of elements desired in the array.
				/// </param>
				/// <param name="neighbours">
				///		The collection of nearest neighbours.
				/// </param>

		public: void InsertKnnMatch(
			vector<Signature> & queries,
			vector<Signature> & subjects,
			size_t k
			) {
			if (matches.size() != k) matches.resize( k );

			SignatureMatch match = SignatureMatch::BestMatch( queries, subjects );
			double similarity = match.Similarity();

			for ( size_t i = k; i > 0; i--) {
				SignatureMatch t = matches[i - 1];

				if (t.Similarity() < similarity) {
					matches[i - 1] = match;

					if (i - 1 < k - 1) {
						matches[i] = t;
					}
				}
			}
		}

	};
}
