#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <ios>
#include <set>
#include <unordered_map>

#include "CsvIO.hpp"
#include "FastaSequence.hpp"
#include "Types.hpp"

namespace QutBio {
	class Homologs {
	public:
		/**
		*	<summary>
		*	Reads a lookup table which tells us which of the genes present in
		*	the current database under consideration are homologs of each query gene.
		*	</summary>
		*	<param name="homologStream">
		*		An input stream that contains the CSV-formatted homolog lookup table.
		*		The lookup table contains of one row for each sequence in the database; the
		*		first field in each record is the numeric Id of a sequence (this is padded
		*		on the left with zeros and has the prefix letter "g" attached, in line with
		*		the naming scheme used by Wayne when he created these datasets.
		*	</param>
		*	<param name="dbIndex">
		*		A lookup table that tells us which genes are present in the database.
		*		This might be obtained by calling FastaSequence::GetIndex on a list of
		*		sequences.
		*	</param>
		*/

		template <typename T1, typename T2>
		static void Parse(
			istream & homologStream,
			Index<T1> &queryIndex,
			Index<T2> & dbIndex
		) {
			CsvReader reader(homologStream);

			auto processRecord = [&](vector<string> & record) {
				// string MakeGeneId( string & s );
				string queryId(record[0]);

				auto found = queryIndex.find(queryId);

				if (found == queryIndex.end()) return true;

				auto seq = found->second;

				for (size_t i = 1; i < record.size(); i++) {
					string homologId(record[i]);
					auto homolog = dbIndex[homologId];
					seq->homologs.push_back(homolog);
				}

				return true;
			};

			auto loadComplete = []() {};

			reader.StreamRecords(processRecord, loadComplete);
		}

		/**
		 *	Gets a lookup table which can be used to conveniently
		 *	determine if two sequences are related.
		 */
		static map<string, set<string>> Parse(
			istream & homologStream,
			char separator = ' '
		) {
			CsvReader reader(homologStream,separator);
			map<string, set<string>> lookupTable;

			auto processRecord = [&](vector<string> & record) {
				lookupTable[record[0]].insert(record.begin(),record.end());
				return true;
			};

			auto loadComplete = []() {};

			reader.StreamRecords(processRecord, loadComplete);

			return lookupTable;
		}

		/**
		**	<summary>
		**	Reads a lookup table which tells us which of the genes present in
		**	the current database under consideration are homologs of each query gene.
		**	</summary>
		**	<param name="qrelsStream">
		**		An input stream that contains the CSV-formatted homolog lookup table.
		**		The lookup table contains of one row for each sequence in the database; the
		**		first field in each record is the numeric Id of a sequence (this is padded
		**		on the left with zeros and has the prefix letter "g" attached, in line with
		**		the naming scheme used by Wayne when he created these datasets.
		**	</param>
		**	<param name="dbIndex">
		**		A lookup table that tells us which genes are present in the database.
		**		This might be obtained by calling FastaSequence::GetIndex on a list of
		**		sequences.
		**	</param>
		**	<typeparam name="T">
		**		A type which has member function (string Id()) and field (vector<T *> homologs).
		**	</typeparam>
		**/

		template<typename T>
		static void ParseQrels(
			FILE * qrelsStream,
			Index<T> &queryIndex,
			Index<T> & dbIndex
		) {
			const int BufferSize = 1000000;
			vector<char> _buffer(BufferSize);
			char * buffer = _buffer.data();
			size_t nextReadLoc = 0;
			size_t availSpace = BufferSize - 1;
			size_t bytesRead;
			char * lastLineStart;

			int linesProcessed = 0;

			auto processLine = [&](char * line) {
				linesProcessed++;

				size_t i = 0;

				// Skip leading spaces
				while (line[i] && isspace(line[i])) i++;
				if (!line[i]) return;

				// Get query Id (first word)
				char * queryId = line + i;
				while (line[i] && !isspace(line[i])) i++;
				line[i++] = 0;

				// Skip one number and intervening spaces
				while (line[i] && isspace(line[i])) i++;
				while (line[i] && !isspace(line[i])) i++;
				while (line[i] && isspace(line[i])) i++;

				if (!line[i]) {
					cerr << "Malformed qrels record starting with '" << queryId << "'" << endl;
					throw Exception("Malformed qrels record", FileAndLine);
				}

				// Get subjectId
				char * subjectId = line + i;
				while (line[i] && !isspace(line[i])) i++;
				line[i++] = 0;

				auto queryLocation = queryIndex.find(queryId);
				if (queryLocation == queryIndex.end()) return;

				auto subjectLocation = dbIndex.find(subjectId);
				if (subjectLocation == dbIndex.end()) return;

				auto querySeq = queryLocation->second;
				auto subjectSeq = subjectLocation->second;

				querySeq->homologs.push_back(subjectSeq);
			};

			while (0 != (bytesRead = fread(buffer + nextReadLoc, sizeof(char), availSpace, qrelsStream))) {
				// Plug a trailing zero into the end of the buffer, just in case.
				buffer[BufferSize - 1] = 0;

				size_t max = bytesRead + nextReadLoc;
				size_t pos = 0;
				lastLineStart = 0;

				while (pos < max) {
					while (pos < max && (buffer[pos] == '\n' || buffer[pos] == 0)) {
						pos++;
					}

					lastLineStart = buffer + pos;

					while (pos < max && !(buffer[pos] == '\n' || buffer[pos] == 0)) {
						pos++;
					}

					if (pos < max) {
						buffer[pos++] = 0;
						processLine(lastLineStart);
						lastLineStart = buffer + pos;
					}
				}

				availSpace = lastLineStart - buffer;
				nextReadLoc = max - availSpace;
				memcpy(buffer, lastLineStart, nextReadLoc);
			}

			(cerr << linesProcessed << " qrels lines processed." << endl).flush();
		}

	private:

		/// <summary>
		///		Inserts a prefix starting with the letter 'g' followed by 
		///		sufficient zeroes to bring the length of the string up to
		///		zeroPaddedLength + 1 characters or more.
		///	</summary>

		static void MakeGeneId(string & s, size_t zeroPaddedLength) {
			while (s.size() < zeroPaddedLength) {
				s.insert(0, 1, '0');
			}
			s.insert(0, 1, 'g');
		}

	};
}
