#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <mutex>
#include <vector>
#include <map>
#include <iostream>

#include <EnumBase.hpp>
#include <EncodedFastaSequence.hpp>
#include <Index.hpp>

using namespace std;

namespace QutBio {
	struct Domain {
		struct Extent {
			uint begin;
			uint end;

			Extent(uint begin, uint end)
				: begin(begin), end(end) {
				if (begin > end) throw Exception("begin > end!", FileAndLine);
			}
		};

		struct Entry {
			string seqId;
			uint seqLen;
			vector<Extent> extents;
		};

		string pfamId;
		string pfamDesc;
		map<string, Entry> entries;

		static bool Parse(istream & str, map<string, Domain> &domains) {
			string s;

			do {
				std::getline(str, s);
				String::TrimInPlace(s);
			}
			while (!str.eof() && !str.fail() && s.length() == 0);

			if (s.length() == 0) return false;

			vector<string> fields = String::Split(s, " \t");

			if (fields[0].length() < 1 || fields[0][0] != '>') {
				throw Exception("First letter of sequence domain list should be '>'.", FileAndLine);
			}

			string seqId = fields[0].substr(1);
			uint seqLen = atoi(fields[fields.size() - 2].c_str());

			// cerr << "\n" << "seqId = " << seqId << "\nseqLen = " << seqLen << "\n";

			while (true) {
				std::getline(str, s);
				String::TrimInPlace(s);

				if (str.eof() || str.fail() || s.length() == 0) break;

				fields = String::Split(s, " \t");

				// cerr << "fields.size() = " << fields.size() << "\n";

				uint pfamIdx = 0;
				string pfamId;

				for (size_t i = 1; i < fields.size(); i++) {
					string s = fields[i];

					if (s.find("PF") == 0) {
						pfamIdx = i;
						pfamId = s.substr(0, s.find("."));
						break;
					}
				}

				if (pfamIdx == 0) {
					throw Exception("Expected 'PF' somewhere in line: s = '" + s + "'", FileAndLine);
				}

				uint numExtents = atoi(fields[1].c_str());

				string pfamDesc = fields[pfamIdx + 1];

				for (uint i = pfamIdx + 2; i + numExtents < fields.size(); i++) {
					pfamDesc += " ";
					pfamDesc += fields[i];
				}

				//cerr << "pfamId = " << pfamId << "\n";
				//cerr << "pfamDesc = " << pfamDesc << "\n";

				Domain & domain = domains[pfamId];
				domain.pfamId = pfamId;
				domain.pfamDesc = pfamDesc;

				Entry & entry = domain.entries[seqId];
				entry.seqId = seqId;
				entry.seqLen = seqLen;
				entry.extents.clear();

				for (uint i = fields.size() - numExtents; i < fields.size(); i++) {
					istringstream str(fields[i]);
					uint begin, end;
					char gap;
					str >> begin >> gap >> end;

					//cerr << "begin = " << begin << "\n";
					//cerr << "gap = " << gap << "\n";
					//cerr << "end = " << end << "\n";

					// Note that Swissprot uses 1-origin for positions, but I use 0-origin.
					entry.extents.emplace_back(begin - 1, end - 1);
				}
			}

			return true;
		}

		friend ostream & operator<<(ostream & out, const Domain & domain) {
			string desc = domain.pfamDesc;

			for (uint i = 0; i < desc.length(); i++) {
				if (desc[i] == ' ') desc[i] = '~';
			}

			out << domain.pfamId << ' ' << desc << ' ' << domain.entries.size() << '\n';

			for (auto & entry : domain.entries) {
				out << entry.second.seqId << ' ' << entry.second.seqLen << ' ' << entry.second.extents.size();

				for (auto & extent : entry.second.extents) {
					out << ' ' << extent.begin << ' ' << extent.end;
				}

				out << '\n';
			}

			return out;
		}

		friend istream & operator>>(istream & in, Domain & domain) {
			uint numEntries;

			in >> domain.pfamId >> domain.pfamDesc >> numEntries;

			for (uint i = 0; i < domain.pfamDesc.length(); i++) {
				if (domain.pfamDesc[i] == '~') domain.pfamDesc[i] = ' ';
			}

			for (uint i = 0; i < numEntries; i++) {
				string seqId;
				uint seqLen, numExtents;

				in >> seqId >> seqLen >> numExtents;

				Entry & entry = domain.entries[seqId];
				entry.seqId = seqId;
				entry.seqLen = seqLen;
				entry.extents.clear();

				for (uint j = 0; j < numExtents; j++) {
					uint begin, end;
					in >> begin >> end;
					entry.extents.emplace_back(begin, end);
				}
			}

			return in;
		}

		static void Load(istream & in, map<string, Domain> &domains) {
			uint numDomains;
			string prefix;

			in >> prefix >> numDomains;

			for (uint i = 0; i < numDomains && !(in.eof() || in.fail()); i++) {
				Domain d;
				in >> d;
				domains[d.pfamId] = d;
			}
		}

		void GetInstances(
			vector<Subsequence> & domainInstances,
			const Index<EncodedFastaSequence> &dbIdx
		) const {
			//cerr << "Processing domain " << this->pfamId << "\n";

			for (auto & entryRec : entries) {
				const string & seqId = entryRec.first;

				auto seqPos = dbIdx.find(seqId);

				if (seqPos == dbIdx.end()) continue;

				for (auto & extent : entryRec.second.extents) {
					// Remember: biologists count from 1, not zero.
					Subsequence sub{ seqPos->second, extent.begin - 1, extent.end - extent.begin + 1 };
					domainInstances.push_back(sub);
				}
			}

			//cerr << "domainInstances.size() = " << domainInstances.size() << "\n";
		}
	};
}
