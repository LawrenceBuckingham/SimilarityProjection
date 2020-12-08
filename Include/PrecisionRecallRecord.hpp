#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include <ios>

#include "Exception.hpp"
#include "String.hpp"
#include "Util.hpp"

using namespace std;

namespace QutBio {
	/// <summary>
	///	Record type which represents a Precision-Recall pair.
	/// </summary>
	class PrecisionRecall {
	private: double precision;
	private: double recall;

	public: PrecisionRecall() {}

	public: PrecisionRecall( size_t relevantItemsRetrieved, size_t itemsRetrieved, size_t relevantDocumentCount ) {
		SetPrecision( (double) relevantItemsRetrieved / itemsRetrieved );
		SetRecall( (double) relevantItemsRetrieved / relevantDocumentCount );
	}

	public: double Precision() {
		return precision;
	}

	public: void SetPrecision( double value ) {
		if ( value < 0 || value > 1 ) throw Exception( "ArgumentException: value must be between 0 and 1.", FileAndLine );

		precision = value;
	}

	public: double Recall() { return recall; }

	public: void SetRecall( double value ) {
		if ( value < 0 || value > 1 ) throw Exception( "ArgumentException: value must be between 0 and 1.", FileAndLine );

		recall = value;
	}

	public: string ToString() {
		stringstream s;
		s << "(" << precision << "," << recall << ")";
		return s.str();
	}

	public: friend ostream & operator <<( ostream & out, const PrecisionRecall & p ) {
		return ( out << "(" << p.precision << "," << p.recall << ")" );
	}

	public: static PrecisionRecall Parse( string s ) {
		s = String::Trim( s );

		if ( s[0] != '(' ) throw Exception( "FormatException: PrecisionRecall must start with '('.", FileAndLine );
		if ( s[s.size() - 1] != ')' ) throw Exception( "FormatException: PrecisionRecall must start with ')'.", FileAndLine );

		s = s.substr( 1, s.size() - 2 );
		vector<string> parts = String::Split( s, ',' );

		if ( parts.size() != 2 ) throw Exception( "FormatException: PrecisionRecall must consist of two comma-separated fields.", FileAndLine );

		PrecisionRecall pr;
		pr.precision = Double::Parse( String::Trim( parts[0] ) );
		pr.recall = Double::Parse( String::Trim( parts[1] ) );

		return pr;
	}
	};

	class PrecisionRecallRecord {
	public:
		string queryId;
		string queryClass;
		size_t relevantDocumentCount;
		vector<PrecisionRecall> kmers;

		PrecisionRecallRecord( const string & queryId, const string & queryClass, size_t relevantDocumentCount )
			: queryId( queryId ), queryClass( queryClass ), relevantDocumentCount( relevantDocumentCount ) {}

		friend ostream & operator << ( ostream & out, const PrecisionRecallRecord & p ) {
			out << p.queryId << ( p.queryClass.size() > 0 ? "|" + p.queryClass : "" );
			out << "," << p.relevantDocumentCount;

			for ( size_t i = 0; i < p.kmers.size(); i++ ) {
				out << ",\"" << p.kmers[i] << "\"";
			}

			out << endl;

			return out;
		}

	public: static PrecisionRecallRecord Parse( vector<string> & csvRecord ) {
		size_t fieldNum = 0;
		auto id_class = String::Split( csvRecord[fieldNum++], "|" );

		PrecisionRecallRecord prRecord( id_class[0], id_class.size() > 1 ? id_class[1] : "", Int::Parse( csvRecord[fieldNum++] ) );

		while ( fieldNum < csvRecord.size() ) {
			prRecord.kmers.push_back( PrecisionRecall::Parse( csvRecord[fieldNum] ) );
			fieldNum++;
		}

		if ( prRecord.kmers.size() > prRecord.relevantDocumentCount ) {
			prRecord.relevantDocumentCount = prRecord.kmers.size();
		}

		return prRecord;
	}
	};

}
