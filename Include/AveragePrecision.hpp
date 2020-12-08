#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif


#include <string>
#include <vector>

#include "Delegates.hpp"
#include "Exception.hpp"
#include "String.hpp"
#include "Util.hpp"

namespace QutBio {
	class AveragePrecision {
	private:
		string queryId;
		string queryClass;
		double averagePrecision;
		int numberDetected;
		int relevantDocumentCount;
	public:
		AveragePrecision(
			const string & queryId,
			const string & queryClass,
			double averagePrecision,
			int numberDetected,
			int relevantDocumentCount
			) :
			queryId( queryId ),
			queryClass( queryClass ),
			averagePrecision( averagePrecision ),
			numberDetected( numberDetected ),
			relevantDocumentCount( relevantDocumentCount ) {}

		const string & QueryId() { return queryId; }
		void SetQueryId( const string & value ) { queryId = value; }

		const string &  QueryClass() { return queryClass; }
		void SetQueryClass( const string & value ) { queryClass = value; }

		int RelevantDocumentCount() { return relevantDocumentCount; }
		void SetRelevantDocumentCount( int value ) { relevantDocumentCount = value; }

		double AvgPrecision() { return averagePrecision; }
		void SetAvgPrecision( double value ) { averagePrecision = value; }

		int NumberDetected() { return numberDetected; }
		void SetNumberDetected( int value ) { numberDetected = value; }

		void ToStringArray( vector<string> & result ) {
			result.clear();
			result.push_back( queryId );
			result.push_back( queryClass );
			result.push_back( Double::ToString( averagePrecision ) );
			result.push_back( Int::ToString( numberDetected ) );
			result.push_back( Int::ToString( relevantDocumentCount ) );
		}

		string ToString() {
			vector<string> record;
			ToStringArray( record );
			return String::Join( record );
		}

		template<typename Collection>
		static AveragePrecision * Parse( const Collection & record ) {
			return new AveragePrecision(
				record[0],
				record[1],
				Double::Parse( record[2] ),
				Int::Parse( record[3] ),
				Int::Parse( record[4] )
				);
		}


#if 0
		/// <summary> Processes a table of AveragePrecision records, one record at a time.
		/// </summary>
		/// <param name="file">The file name (or in HTML, a File object).</param>
		/// <param name="processRecord">A callback function that will be invoked to process each record.</param>
		/// <param name="loadComplete">A callback that is invoked once all records have been processed.</param>
		/// <param name="onError">Error handler.</param>

		static void StreamFile(
			string file,
			Action1<AveragePrecision &> processRecord,
			Action loadComplete,
			Action1<Exception> onError
			) {

		}

#endif // 0

	};

}
