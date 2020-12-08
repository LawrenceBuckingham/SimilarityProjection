#pragma once

#include "Kmer.hpp"
#include <ios>

namespace QutBio {

	// Represents an ungapped pairwise alignment between two sequences.
	struct Edge {
		string id1, id2;
		int start1, stop1, start2, stop2;

		Edge() {}

		Edge(const Kmer::Instance & k1, const Kmer::Instance & k2, size_t length) {
			id1 = k1.sequence.IdStr();
			id2 = k2.sequence.IdStr();
			start1 = k1.kmerPosition + 1;
			start2 = k2.kmerPosition + 1;
			stop1 = start1 + length - 1;
			stop2 = start2 + length - 1;
		}

		bool operator<(Edge const & other) const {
			if (id1 < other.id1) {
				return true;
			}
			else if (id1 == other.id1) {
				if (id2 < other.id2) {
					return true;
				}
				else if (id2 == other.id2) {
					if (start1 < other.start1) {
						return true;
					}
					else if (start1 == other.start1) {
						return start2 < other.start2;
					}
					else {
						return false;
					}
				}
				else {
					return false;
				}
			}
			else {
				return false;
			}
		}

		friend ostream & operator<<(ostream & out, Edge const & edge) {
			out << edge.id1 << "\t"
				<< edge.start1 << "\t"
				<< edge.stop1 << "\t"
				<< edge.id2 << "\t"
				<< edge.start2 << "\t"
				<< edge.stop2;
			return out;
		}
	};

}
