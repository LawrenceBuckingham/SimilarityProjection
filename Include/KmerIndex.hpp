#pragma once

#include <unordered_map>
#include <vector>

#include "Kmer.hpp"
#include <cstdio>
#include <chrono>
#include <thread>

using namespace std;

namespace QutBio {

	class KmerIndex : public unordered_map<Substring, Kmer *, Substring::Hash> {
		// TODO: this is not a valid use of inheritance. Convert the map to a member variable.

	protected:
		vector<Kmer *> allKmers;

	public:
		using BaseType = unordered_map<Substring, Kmer *, Substring::Hash>;

		/*
			**	Summary:
			**		Constructs a KmerIndex.
			**
			**	Parameters:
			**		dataset:    a list of sequences;
			**		kmerLength: the word length for tiling;
			*/

		KmerIndex(
			vector<EncodedFastaSequence *> & dataset,
			size_t kmerLength
		) {
			size_t length = dataset.size();

			for ( size_t i = 0; i < length; i++ ) {
				auto & seq = *dataset[i];
				const Symbol*residues{ seq.Sequence().data() };
				const size_t seqLen{ seq.Sequence().size() };

				if ( seqLen < kmerLength ) continue;

				uint kmerCount = seq.KmerCount( kmerLength );

				for ( size_t kmerPos = 0; kmerPos < kmerCount; kmerPos++ ) {
					AddKmer( residues, kmerPos, kmerLength, seq );
				}
			}
		}

		KmerIndex( const KmerIndex && other ) {
			allKmers.clear();
			allKmers = move( other.allKmers );
		}

		KmerIndex & operator= ( const KmerIndex && other ) {
			allKmers.clear();
			allKmers = move( other.allKmers );
			return *this;
		}

		void AddKmer(
			const Symbol* residues,
			const size_t &kmerPos,
			const size_t &kmerLength,
			QutBio::EncodedFastaSequence & seq
		) {
			Substring key( residues, kmerPos, kmerLength );
			auto item = this->find( key );

			if ( item == this->end() ) {
				auto value = new Kmer( seq, kmerPos, kmerLength );
				std::pair<Substring, Kmer *> p( key, value );
				BaseType::insert( p );
			}
			else {
				item->second->Add( seq, kmerPos );
			}
		}

		KmerIndex(
			vector<Subsequence> & substrings,
			size_t kmerLength
		) {
			for ( auto & substring : substrings ) {
				auto & seq = *substring.source;
				const Symbol*residues{ seq.Sequence().data() };
				uint kmerCount = seq.KmerCount( kmerLength );

				for ( size_t kmerPos = substring.start;
					kmerPos < kmerCount && kmerPos + kmerLength <= substring.start + substring.length;
					kmerPos++
					) {
					AddKmer( residues, kmerPos, kmerLength, seq );
				}
			}
		}

		KmerIndex(
			const Subsequence * substrings,
			size_t length,
			size_t kmerLength
			//
		) {
			for ( size_t i = 0; i < length; i++ ) {
				auto & substring = substrings[i];
				auto & seq = *substring.source;
				auto residues{ seq.Sequence().data() };
				auto kmerCount = seq.KmerCount( kmerLength );

				for ( auto kmerPos = substring.start; kmerPos < kmerCount && kmerPos + kmerLength <= substring.start + substring.length; kmerPos++ ) {
					AddKmer( residues, kmerPos, kmerLength, seq );
				}
			}
		}

		/**
			**	Summary:
			**		Destructor.
			*/
		virtual ~KmerIndex() {
			for ( auto &p : *this ) {
				delete p.second;
			}
		}

		/**
			**	Summary:
			**		Updates an internal cache of kmer pointers and then returns it by reference.
			*/
		vector<Kmer *> &GetKmers() {
			size_t id = 0;

			if ( allKmers.size() != this->size() ) {
				allKmers.clear();

				for ( auto &pair : *this ) {
					auto kmer = pair.second;
					kmer->SetSerialNumber( id++ );
					allKmers.push_back( kmer );
				}
			}

			return allKmers;
		}

		size_t size() {
			return BaseType::size();
		}
	};

	/**
	 *	Alternate version for use when the kmer encoding implemented in Alphabet
	 *	yields a perfect hash. That is,
	 *		|Sigma|^k-1 <= numeric_limits<size_t>::max();
	 */
	class KmerHashIndex : public unordered_map<size_t, Kmer *> {
		// TODO: this is not a valid use of inheritance. Convert the map to a member variable.

	protected:
		vector<Kmer *> allKmers;
		size_t kmerLength;
		const Alphabet * alphabet;

	public:
		using BaseType = unordered_map<size_t, Kmer *>;

		/*
			**	Summary:
			**		Constructs a KmerIndex.
			**
			**	Parameters:
			**		dataset:    a list of sequences;
			**		kmerLength: the word length for tiling;
			**		alphabet:   reference to the alphabet used to generate
			**		            hash codes.
			*/

		KmerHashIndex(
			vector<EncodedFastaSequence *> & dataset,
			size_t kmerLength,
			const Alphabet * alphabet
		) :
			kmerLength( kmerLength ),
			alphabet( alphabet )
			//
		{
			size_t length = dataset.size();

			for ( size_t i = 0; i < length; i++ ) {
				auto & seq = *dataset[i];
				const Symbol*residues{ seq.Sequence().data() };
				const size_t seqLen{ seq.Sequence().size() };

				if ( seqLen < kmerLength ) continue;

				uint kmerCount = seq.KmerCount( kmerLength );

				for ( size_t kmerPos = 0; kmerPos < kmerCount; kmerPos++ ) {
					AddKmer( residues, kmerPos, seq );
				}
			}
		}

		KmerHashIndex( const KmerHashIndex && other ) :
			kmerLength( other.kmerLength ),
			alphabet( other.alphabet )
		{
			allKmers.clear();
			allKmers = move( other.allKmers );
		}

		KmerHashIndex & operator= ( const KmerHashIndex && other ) {
			kmerLength = other.kmerLength;
			alphabet = other.alphabet;
			allKmers.clear();
			allKmers = move( other.allKmers );
			return *this;
		}

		void AddKmer(
			const Symbol* residues,
			const size_t kmerPos,
			EncodedFastaSequence & seq
		) {
			Substring substr( residues, kmerPos, kmerLength, alphabet );
			size_t key = substr.HashCode();
			auto item = this->find( key );

			if ( item == this->end() ) {
				auto value = new Kmer( seq, kmerPos, kmerLength );
				std::pair<size_t, Kmer *> p( key, value );
				BaseType::insert( p );
			}
			else {
				item->second->Add( seq, kmerPos );
			}
		}

		KmerHashIndex(
			vector<Subsequence> & substrings,
			size_t kmerLength,
			const Alphabet * alphabet
		) :
			kmerLength( kmerLength ),
			alphabet( alphabet )
			//
		{
			for ( auto & substring : substrings ) {
				auto & seq = *substring.source;
				const Symbol*residues{ seq.Sequence().data() };
				uint kmerCount = seq.KmerCount( kmerLength );

				for ( size_t kmerPos = substring.start;
					kmerPos < kmerCount && kmerPos + kmerLength <= substring.start + substring.length;
					kmerPos++
					) {
					AddKmer( residues, kmerPos, seq );
				}
			}
		}

		KmerHashIndex(
			const Subsequence * substrings,
			size_t length,
			size_t kmerLength,
			const Alphabet * alphabet
		) :
			kmerLength( kmerLength ),
			alphabet( alphabet )
			//
		{
			for ( size_t i = 0; i < length; i++ ) {
				auto & substring = substrings[i];
				auto & seq = *substring.source;
				auto residues{ seq.Sequence().data() };
				auto kmerCount = seq.KmerCount( kmerLength );

				for ( auto kmerPos = substring.start; kmerPos < kmerCount && kmerPos + kmerLength <= substring.start + substring.length; kmerPos++ ) {
					AddKmer( residues, kmerPos, seq );
				}
			}
		}

		/**
			**	Summary:
			**		Destructor.
			*/
		virtual ~KmerHashIndex() {
			for ( auto &p : *this ) {
				delete p.second;
			}
		}

		/**
			**	Summary:
			**		Updates an internal cache of kmer pointers and then returns it by reference.
			*/
		vector<Kmer *> &GetKmers() {
			size_t id = 0;

			if ( allKmers.size() != this->size() ) {
				allKmers.clear();

				for ( auto &pair : *this ) {
					auto kmer = pair.second;
					kmer->SetSerialNumber( id++ );
					allKmers.push_back( kmer );
				}
			}

			return allKmers;
		}

		size_t size() {
			return BaseType::size();
		}
	};
} // namespace QutBio
