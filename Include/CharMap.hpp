#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cstdint>
#include "Alphabet.hpp"

namespace QutBio {

	typedef struct BitRep {
		uint64_t lo, hi;

		void Clear() {
			lo = 0;
			hi = 0;
		}
	} BitRep;

#define BITS_PER_WORD (sizeof(uint64_t)<<3)

	class CharMap : public ByteIdxArray<BitRep> {
	public:

		/// <summary>
		/// Gets a 64-bit half-embedding, 27 bits set for each residue, for BLOSUM62.
		/// </summary>
		/// <returns>A 64-bit half-embedding, 27 bits set for each residue, for BLOSUM62.</returns>
		static const CharMap& Blosum62QueryEncoding() {
			static CharMap* encoding = 0;

			if ( !encoding ) {
				encoding = new CharMap();
				auto a = Alphabets::AA();
				( *encoding )[a->Encode( 'a' )].lo = 4196281838917878893ull;
				( *encoding )[a->Encode( 'r' )].lo = 7650073181085339229ull;
				( *encoding )[a->Encode( 'n' )].lo = 16820669735176575068ull;
				( *encoding )[a->Encode( 'd' )].lo = 13974388523644329108ull;
				( *encoding )[a->Encode( 'c' )].lo = 3754921625820924652ull;
				( *encoding )[a->Encode( 'q' )].lo = 2966207013620391484ull;
				( *encoding )[a->Encode( 'e' )].lo = 3615058547148921981ull;
				( *encoding )[a->Encode( 'g' )].lo = 7081679552086086861ull;
				( *encoding )[a->Encode( 'h' )].lo = 11541817753105046620ull;
				( *encoding )[a->Encode( 'i' )].lo = 1989466549711871335ull;
				( *encoding )[a->Encode( 'l' )].lo = 2034503094600777063ull;
				( *encoding )[a->Encode( 'k' )].lo = 8693214589846654589ull;
				( *encoding )[a->Encode( 'm' )].lo = 8078335720694865167ull;
				( *encoding )[a->Encode( 'f' )].lo = 1334197304103321434ull;
				( *encoding )[a->Encode( 'p' )].lo = 2323611251589552409ull;
				( *encoding )[a->Encode( 's' )].lo = 3613470385805040733ull;
				( *encoding )[a->Encode( 't' )].lo = 3560253876602510204ull;
				( *encoding )[a->Encode( 'w' )].lo = 134217727ull;
				( *encoding )[a->Encode( 'y' )].lo = 2454539073132404596ull;
				( *encoding )[a->Encode( 'v' )].lo = 10528300240591231349ull;
				( *encoding )[a->Encode( 'b' )].lo = 16242056208945323541ull;
				( *encoding )[a->Encode( 'z' )].lo = 2417730536003701791ull;
				( *encoding )[a->Encode( 'x' )].lo = 4840157387973009236ull;
				( *encoding )[a->Encode( 'A' )].lo = 4196281838917878893ull;
				( *encoding )[a->Encode( 'R' )].lo = 7650073181085339229ull;
				( *encoding )[a->Encode( 'N' )].lo = 16820669735176575068ull;
				( *encoding )[a->Encode( 'D' )].lo = 13974388523644329108ull;
				( *encoding )[a->Encode( 'C' )].lo = 3754921625820924652ull;
				( *encoding )[a->Encode( 'Q' )].lo = 2966207013620391484ull;
				( *encoding )[a->Encode( 'E' )].lo = 3615058547148921981ull;
				( *encoding )[a->Encode( 'G' )].lo = 7081679552086086861ull;
				( *encoding )[a->Encode( 'H' )].lo = 11541817753105046620ull;
				( *encoding )[a->Encode( 'I' )].lo = 1989466549711871335ull;
				( *encoding )[a->Encode( 'L' )].lo = 2034503094600777063ull;
				( *encoding )[a->Encode( 'K' )].lo = 8693214589846654589ull;
				( *encoding )[a->Encode( 'M' )].lo = 8078335720694865167ull;
				( *encoding )[a->Encode( 'F' )].lo = 1334197304103321434ull;
				( *encoding )[a->Encode( 'P' )].lo = 2323611251589552409ull;
				( *encoding )[a->Encode( 'S' )].lo = 3613470385805040733ull;
				( *encoding )[a->Encode( 'T' )].lo = 3560253876602510204ull;
				( *encoding )[a->Encode( 'W' )].lo = 134217727ull;
				( *encoding )[a->Encode( 'Y' )].lo = 2454539073132404596ull;
				( *encoding )[a->Encode( 'V' )].lo = 10528300240591231349ull;
				( *encoding )[a->Encode( 'B' )].lo = 16242056208945323541ull;
				( *encoding )[a->Encode( 'Z' )].lo = 2417730536003701791ull;
				( *encoding )[a->Encode( 'X' )].lo = 4840157387973009236ull;
			}

			return ( *encoding );
		}

		/// <summary>
		/// Gets a 64-bit complementary embedding, 27 bits set for each residue, for BLOSUM62.
		/// </summary>
		/// <returns>A 64-bit complementary embedding, 27 bits set for each residue, for BLOSUM62.</returns>
		static const CharMap& Blosum62SubjectEncoding() {
			static CharMap* encoding = 0;

			if ( !encoding ) {
				encoding = new CharMap();
				auto a = Alphabets::AA();
				( *encoding )[a->Encode( 'a' )].lo = 2863761771407970925ull;
				( *encoding )[a->Encode( 'r' )].lo = 7651199062198035261ull;
				( *encoding )[a->Encode( 'n' )].lo = 14505852547472661084ull;
				( *encoding )[a->Encode( 'd' )].lo = 3595913551146720277ull;
				( *encoding )[a->Encode( 'c' )].lo = 3755053567216261860ull;
				( *encoding )[a->Encode( 'q' )].lo = 3006730097971289629ull;
				( *encoding )[a->Encode( 'e' )].lo = 12874384598663773244ull;
				( *encoding )[a->Encode( 'g' )].lo = 7658265648044020940ull;
				( *encoding )[a->Encode( 'h' )].lo = 11541819024448920664ull;
				( *encoding )[a->Encode( 'i' )].lo = 269102453885837161ull;
				( *encoding )[a->Encode( 'l' )].lo = 584915626282040166ull;
				( *encoding )[a->Encode( 'k' )].lo = 6558930587529087837ull;
				( *encoding )[a->Encode( 'm' )].lo = 8073852185476959501ull;
				( *encoding )[a->Encode( 'f' )].lo = 1334828286049501018ull;
				( *encoding )[a->Encode( 'p' )].lo = 7007073522020817209ull;
				( *encoding )[a->Encode( 's' )].lo = 4262410801802746462ull;
				( *encoding )[a->Encode( 't' )].lo = 8316072681063168622ull;
				( *encoding )[a->Encode( 'w' )].lo = 134217727ull;
				( *encoding )[a->Encode( 'y' )].lo = 2455735375069421426ull;
				( *encoding )[a->Encode( 'v' )].lo = 17516751889262022129ull;
				( *encoding )[a->Encode( 'b' )].lo = 7054334882014501973ull;
				( *encoding )[a->Encode( 'z' )].lo = 2390568716419798137ull;
				( *encoding )[a->Encode( 'x' )].lo = 5930836213530205298ull;
				( *encoding )[a->Encode( 'A' )].lo = 2863761771407970925ull;
				( *encoding )[a->Encode( 'R' )].lo = 7651199062198035261ull;
				( *encoding )[a->Encode( 'N' )].lo = 14505852547472661084ull;
				( *encoding )[a->Encode( 'D' )].lo = 3595913551146720277ull;
				( *encoding )[a->Encode( 'C' )].lo = 3755053567216261860ull;
				( *encoding )[a->Encode( 'Q' )].lo = 3006730097971289629ull;
				( *encoding )[a->Encode( 'E' )].lo = 12874384598663773244ull;
				( *encoding )[a->Encode( 'G' )].lo = 7658265648044020940ull;
				( *encoding )[a->Encode( 'H' )].lo = 11541819024448920664ull;
				( *encoding )[a->Encode( 'I' )].lo = 269102453885837161ull;
				( *encoding )[a->Encode( 'L' )].lo = 584915626282040166ull;
				( *encoding )[a->Encode( 'K' )].lo = 6558930587529087837ull;
				( *encoding )[a->Encode( 'M' )].lo = 8073852185476959501ull;
				( *encoding )[a->Encode( 'F' )].lo = 1334828286049501018ull;
				( *encoding )[a->Encode( 'P' )].lo = 7007073522020817209ull;
				( *encoding )[a->Encode( 'S' )].lo = 4262410801802746462ull;
				( *encoding )[a->Encode( 'T' )].lo = 8316072681063168622ull;
				( *encoding )[a->Encode( 'W' )].lo = 134217727ull;
				( *encoding )[a->Encode( 'Y' )].lo = 2455735375069421426ull;
				( *encoding )[a->Encode( 'V' )].lo = 17516751889262022129ull;
				( *encoding )[a->Encode( 'B' )].lo = 7054334882014501973ull;
				( *encoding )[a->Encode( 'Z' )].lo = 2390568716419798137ull;
				( *encoding )[a->Encode( 'X' )].lo = 5930836213530205298ull;
			}

			return ( *encoding );
		}

		static std::pair<const CharMap&, const CharMap&> Get27BitBlosum62DualEmbedding() {
			std::pair<const CharMap&, const CharMap&> res{
				Blosum62QueryEncoding(), Blosum62SubjectEncoding()
			};

			return res;
		}
	};
}
