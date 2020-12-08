#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <algorithm>
#include <ios>
#include <string>

#include "IArrayParser.hpp"
#include "EncodedFastaSequence.hpp"

namespace QutBio {

	/// <summary> A ranked document result.
	/// </summary>

	struct Ranking : public IArrayParser {
		/// <summary> The document ID of the ranked document. 
		/// </summary>

	public:
		const string * queryId;
		const string * subjectId;
		double distance;
		size_t rank;
		int hits;

	public:
		Ranking(
			const string & queryId,
			const string & subjectId,
			double distance,
			size_t rank,
			int hits
		) :
			queryId(&queryId),
			subjectId(&subjectId),
			distance(distance),
			rank(rank),
			hits(hits) {}

		// For serialisation only.
		Ranking() {}

		/// <summary> Orders two ranking objects in ascending order of distance.
		/// </summary>
		/// <param name="x"></param>
		/// <param name="y"></param>
		/// <returns></returns>

		static bool AscendingDistance(const Ranking & x, const Ranking & y) {
			return x.distance < y.distance;
		}

		friend bool operator< (const Ranking & x, const Ranking & y) {
			return x.distance < y.distance;
		}

		string ToString() {
			ostringstream writer;
			CsvWriter csv(writer);
			vector<string> fields;
			ToStringArray(fields);
			csv.Write(fields);
			return writer.str();
		}

		/// <summary> Returns an array of string which contains the fields, in order. 
		/// <para>Any missing values are represented by a single '?' character.</para>
		/// </summary>
		/// <returns></returns>

		void ToStringArray(vector<string> & record) {
			record.clear();
			record.push_back(*queryId);
			record.push_back(*subjectId);
			record.push_back(std::to_string(distance));
			record.push_back(std::to_string(rank));
			record.push_back(std::to_string(hits));
		}

		/**
		 ** Serialise a ranking object to text format suitable for
		 **	trec_eval. Distance is negated to give a "similarity".
		 */

		friend ostream & operator<<(ostream & out, const Ranking & ranking) {
			out << *(ranking.queryId)
				<< " 0 "
				<< *(ranking.subjectId)
				<< " 0 "
				<< (-ranking.distance)
				<< " "
				<< ranking.hits;

			return out;
		}

		/// <summary>
		/// Returns a new Ranking record extracted from an array of strings which represent
		/// the fields of a Ranking record.
		/// </summary>
		/// <param name="record"></param>
		/// <returns></returns>

		void Parse(vector<string> & record) {
			throw NotImplementedException(FileAndLine);
		}

		static int CompareByDistance(const void * r1, const void * r2) {
			Ranking **x = (Ranking **)r1;
			Ranking **y = (Ranking **)r2;
			double xd = (**x).distance;
			double yd = (**y).distance;

			return xd < yd ? -1 : xd > yd ? +1 : 0;
		}

		template<typename Collection, typename Node>
		static void SerialiseCompact( Collection & collection, function<const Ranking *(const Node & n)> project, ostream & out ) {
			const string * previousQuery = 0;

			for ( Node & node: collection ) {
				const Ranking * ranking = project(node);

				// (cerr << "Node = (" << node.first << "," << node.second << ")\n").flush();

				if ( ranking->queryId != previousQuery ) {
					if ( previousQuery ) out << "\n";

					out << *(ranking->queryId);
					previousQuery = ranking->queryId; 
				}

				out << " " << *(ranking->subjectId) << " " << (-ranking->distance);
			}

			out << "\n";
		}

	};

	typedef vector<vector<Ranking>> Rankings;
}
