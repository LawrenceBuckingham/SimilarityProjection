#pragma once
#include <string>
#include <vector>
#include <cctype>

#include "Exception.hpp"
#include "Types.hpp"
#include "LineReader.hpp"
#include "Registry.hpp"
#include "String.hpp"
#include "SubstitutionMatrix.hpp"
#include "Util.hpp"

using namespace std;

namespace QutBio {
	class Sequence {
		using Register = Registry<string, std::hash<string>>;

	protected:
		size_t id;
		vector<string> metadata;
		vector<size_t> classes;
		string chars;
		vector<Symbol> seq;
		vector<Digram> digrams;

		void UpdateSeq( const QutBio::SubstitutionMatrix& matrix ) {
			matrix.Encode( chars.begin(), chars.end(), seq );
			Interleave( seq, matrix.Size(), 2, digrams );
		}

	public:
		Sequence() {}
		Sequence( const Sequence& other ) = delete;
		Sequence( Sequence&& other ) = delete;
		virtual ~Sequence() {}

		size_t Id() const { return id; }

		void Id( size_t value ) {
			if ( value >= Sequence::IdReg().Size() ) {
				throw ArgumentException( NAMEOF_VARIABLE( value ), FileAndLine );
			}

			id = value;
		}

		const string & IdString() const {
			auto& reg = Sequence::IdReg();
			auto& s = reg.At( id );
			return s;
		}

		void IdString( const string& value ) {
			id = Sequence::GetId( value );
		}

		const vector<size_t>& Classes() const { return classes; }

		void Classes( const vector<size_t>& value ) {
			classes = value;
			std::sort( classes.begin(), classes.end() );
		}

		void Classes( const vector<string>& value ) {
			classes.resize( value.size() );
			size_t i = 0;

			for ( auto& s : value ) {
				classes[i++] = Sequence::GetClass( s );
			}

			std::sort( classes.begin(), classes.end() );
		}

		vector<string> ClassStrings() {
			vector<string> res{ classes.size() };
			size_t i = 0;

			for ( auto c : classes ) {
				res[i++] = Sequence::ClassReg().At( c );
			}

			return res;
		}

		const string& Chars() const { return chars; }

		string DefLine() const {
			return String::Join(metadata, "|");
		}

		void Chars( const string& value, const SubstitutionMatrix& matrix ) {
			chars = value;
			UpdateSeq( matrix );
		}

		const vector<Symbol>& Seq() const { return seq; }

		const vector<Digram>& Digrams() const { return digrams; }

		bool Parse( const string& defLine, const string& charData, const SubstitutionMatrix& matrix, int idIndex = -1, int classIndex = -1 ) {
			if ( defLine[0] != '>' ) return false;

			metadata.clear();

			if ( idIndex >= 0 || classIndex >= 0 ) {
				metadata = String::Split( defLine, ">|" );
			}

			IdString( idIndex < 0 ? defLine : metadata.at( idIndex ) );

			if ( classIndex >= 0 ) {
				auto classNames = String::Split( metadata.at( classIndex ), ';' );
				Classes( classNames );
			}

			chars.clear();

			for ( auto c : charData ) {
				if ( matrix.IsDefined(c) ) chars.push_back( c );
			}

			UpdateSeq( matrix );

			return true;
		}

		bool ParseFasta( LineReader& in, const SubstitutionMatrix& matrix, int idIndex = -1, int classIndex = -1 ) {
			if ( !in.Ok() ) return false;

			if ( in.CurrentLine()[0] != '>' ) return false;

			string defLine = in.CurrentLine();
			string charData;

			while ( in.Advance() ) {
				if ( in.CurrentLine()[0] == '>' ) break;
				charData += in.CurrentLine();
			}

			return Parse(defLine, charData, matrix, idIndex, classIndex);
		}

		static vector<Sequence*> ParseAllFasta( LineReader& in, const SubstitutionMatrix& matrix, int idIndex = -1, int classIndex = -1  ) {
			vector<Sequence*> res;

			while ( in.Ok() ) {
				auto p = new Sequence();

				if ( !p->ParseFasta( in, matrix, idIndex, classIndex ) ) {
					delete p;
					break;
				}

				res.push_back( p );
			}

			return res;
		}

		static Register& IdReg() {
			static Register reg;
			return reg;
		}

		static Register& ClassReg() {
			static Register reg;
			return reg;
		}

		static size_t GetId( const string& s ) {
			auto& reg = IdReg();
			return reg( s );
		}

		static size_t GetClass( const string& s ) {
			auto& reg = ClassReg();
			return reg( s );
		}

		template<typename T>
		static void Interleave( const vector<Symbol>& seq, int radix, int stepSize, vector<T>& res ) {
			const size_t n = seq.size();
			res.resize( n + 1 - stepSize );

			for ( size_t i = 0; i < res.size(); i++ ) {
				T t{seq[i].value};
				for ( int j = 1; j < stepSize; j++ ) {
					t += radix * seq[i + j].value;
				}
				res[i] = t;
			}
		}
	};
}

