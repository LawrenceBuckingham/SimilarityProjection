#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include <string>
#include "Exception.hpp"
#include "Types.hpp"
#include "CsvIO.hpp"
#include "SparseSet.hpp"

using namespace std;

namespace QutBio {
	class SparseSignature :
		public virtual SparseSet,
		public virtual ICsvReader,
		public virtual ICsvWriter {

	private:
		const FastaSequence *sequence;

	public:
		SparseSignature( const FastaSequence *sequence = 0 ) : sequence( sequence ) {}

		virtual ~SparseSignature() {}

		const FastaSequence * Sequence() const {
			return sequence;
		}

		void SetSequence( const FastaSequence * value ) {
			sequence = value;
		}

		static vector<SparseSignature *> Read(
			string &sigFile
		) {
			ifstream sigStream( sigFile );

			if ( sigStream.fail() ) {
				cerr << "File " << sigFile << " did not open properly\n";
				throw Exception( "Error reading file " + sigFile, FileAndLine );
			}

			return Read( sigStream );
		}

		/// Dependency injection: function that gets sequences from their ID.
		static function<const FastaSequence *(const string & seqId)> Lookup;

		/// Parse a collection of signatures. 
		static vector<SparseSignature *> Read( ifstream &sigStream ) {
			vector<SparseSignature *> signatures;
			string seqId;

			while ( !sigStream.eof() ) {
				sigStream >> seqId;

				if ( seqId.length() == 0 ) {
					break;
				}

				const FastaSequence * seq = Lookup( seqId );
				SparseSignature * sig = new SparseSignature( seq );
				sigStream >> *sig;
			}

			return signatures;
		}

		static void CreatePostingList(
			const vector<SparseSignature *> &dbSigs,
			vector<vector<uint>> & index
		) {
			const uint D = dbSigs.size();

			uint max_index = 0;

			for ( uint d = 0; d < D; d++ ) {
				for ( auto i : dbSigs[d]->features ) {
					if ( i > max_index ) max_index = i;
				}
			}

			index.resize( max_index + 1 );

			for ( uint d = 0; d < D; d++ ) {
				for ( auto i : dbSigs[d]->features ) {
					index[i].push_back( d );
				}
			}
		}

		/**
		 *	Derive a posting list from sparse signature.
		 *	@pre (dbSigs.size() == M)
		 *		 and ((i in dom(selectedSignatures)) => (0 <= selectedSignatures[i] < M))
		 *	@post ((i in dom(index) and (j in ran(index[i]))) => (i in ran(dbSigs[j].indices)))
		 */

		static void CreatePostingList(
			const vector<SparseSignature> &dbSigs,
			const vector<size_t> &selectedSignatures,
			vector<vector<size_t>> & index
		) {
			const uint D = selectedSignatures.size();

			uint max_index = 0;

			for ( auto d : selectedSignatures ) {
				for ( auto i : dbSigs[d].features ) {
					if ( i > max_index ) max_index = i;
				}
			}

			index.clear();
			index.resize( max_index + 1 );

			for ( auto d : selectedSignatures ) {
				for ( auto i : dbSigs[d].features ) {
					index[i].push_back( d );
				}
			}
		}

		void Write( CsvWriter &w ) const override {
			w << sequence->IdStr();
			SparseSet::Write(w);
		}

		void Read( CsvReader &r ) override {
			string idStr;
			r >> idStr;

			sequence = Lookup( idStr );
			SparseSet::Read(r);
		}

	};
}
