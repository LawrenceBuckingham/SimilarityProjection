#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <functional>
#include <numeric>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <csignal>
#include <cstdio>

#include "Kmer.hpp"

#if defined(SHOW_PROGRESS)
#define PROGRESS(x) x
#else
#define PROGRESS(x)
#endif

#if USE_OMP
#include <omp.h>
#endif

namespace QutBio {
	/**
	 **	<summary>
	 **		Template class representing cluster of Kmers with a "central" prototype.
	 **	</summary>
	 */
	template<typename DistanceFunction>
	struct CovTree {
		// using DistanceFunction = function<Distance(const Kmer *, const Kmer *)>;

		CovTree * parent;
		const Kmer * prototype;
		const Distance * threshold;
		const DistanceFunction &distance;
		vector<CovTree *> subTrees;
		vector<Kmer *> storedKmers;

		/**
		 *	Initialise a CovTreeNode.
		 *
		 *	Parameters:
		 *		parent:    The node that covers this branch and its siblings.
		 *		prototype: Reference to the central attractor for a branch of the tree.
		 *		threshold: Zero-terminated monotonic decreasing list of threshold values for this node and its subTrees.
		 *		distance:	The distance function used for comparison of kmers.
		 */

		CovTree(
			CovTree *parent,
			const Kmer * prototype,
			const Distance * threshold,
			const DistanceFunction & distanceFunction
			//
		) :
			parent(parent),
			prototype(prototype),
			threshold(threshold),
			distance(distanceFunction) {
			if (!prototype) {
				throw Exception("prototype may not be null!", FileAndLine);
			}
		}

		/**
		 *	Release memory in child nodes.
		 */

		virtual ~CovTree() {
			for (auto child : subTrees) {
				delete child;
			}
		}

		/**
		 *	Returns true if and only if the tree is a leaf, i.e. a terminal branch which
		 *	stores kmers other than the prototype.
		 */

		bool IsLeaf() {
			return threshold[1] == 0;
		}

		/**
		 *	Traverses the tree and finds the deepest node that covers the
		 *	designated Kmer. Does not "find" the kmer.
		 *
		 *	Parameters:
		 *		kmer: The address of a kmer.
		 *
		 *	Returns:
		 *		If the kmer is not covered by the current sub-tree, returns NULL.
		 *		Otherwise, returns the address of the first sub-tree detected which
		 *		covers the point.
		 */
		CovTree * GetDeepestCover(Kmer * queryPoint) {
			if (!Covers(queryPoint)) return 0;

			CovTree * nearbyChild = 0;

			for (auto child : subTrees) {
				nearbyChild = child->GetDeepestCover(queryPoint);

				if (nearbyChild) break;
			}

			return nearbyChild ? nearbyChild : this;
		}

		/**
		 *	Predicate returns true if and only if the current sub-tree covers
		 *	the kmer.
		 */

		bool Covers(Kmer * kmer) {
			auto d = distance(kmer, prototype);
			return d <= threshold[0];
		}

		/**
		 *	Inserts a kmer into the tree, where necessary extending the
		 *	list of cover trees at each threshold level to ensure there is
		 *	a path leading to a storage node where the kmer can be found.
		 */

		void Insert(Kmer * kmer) {
			CovTree * nearestCover = GetDeepestCover(kmer);

			if (nearestCover) {
				nearestCover->InsertDirectly(kmer);
			}
			else {
				parent->InsertDirectly(kmer);
			}
		}

		/**
		 *	Inserts a kmer into the tree, adding new cover trees at each threshold level
		 *	as needed to ensure there is a path leading to a storage node where the kmer
		 *	can be found. The newly added kmer will be the prototype of all new sub-trees
		 *	created in this process.
		 */

		void InsertDirectly(Kmer * point) {
			if (IsLeaf()) {
				storedKmers.push_back(point);
			}
			else {
				CovTree * subtree = new CovTree(this, point, threshold + 1, distance);
				subTrees.push_back(subtree);
				subtree->InsertDirectly(point);
			}
		}

		void Print(ostream & out) __attribute__((used)) {
			PrintAux(out, 0);
		}

		void PrintAux(ostream & out, int nestingLevel) {
			Indent(out, nestingLevel);
			out << "<Tree threshold=\"" << threshold[0] << "\" prototype=\"" << prototype->Substr() << "\">\n";
			for (auto child : subTrees) {
				child->PrintAux(out, nestingLevel + 1);
			}
			for (auto kmer : storedKmers) {
				Indent(out, nestingLevel + 1);
				out << kmer->Substr() << "\n";
			}
			Indent(out, nestingLevel);
			out << "</Tree>\n";
		}

		void Indent(ostream & out, int nestingLevel) {
			for (int i = 0; i < nestingLevel; i++) {
				out << "\t";
			}
		}

		void Serialize(ostream & out) {
			out << "CoverTree;";

			if (parent == 0) {
				out << "thresholds:";

				for (auto t = threshold; *t; t++) {
					out << *t << ";";
				}

				out << "0;";
			}

			out << subTrees.size() << ";" << storedKmers.size() << "\n";

			for (auto subTree : subTrees) {
				subTree->Serialize(out);
			}

			for (auto kmer : storedKmers) {
				out << (*kmer) << "\n";
			}
		}

		static Kmer & KmerZero() {
			static EncodedFastaSequence seqZero("", "", "", "", Alphabet::AA(), 1, 1, 'a');
			static Kmer zero(&seqZero, 0, 1, 1);
			return zero;
		}

		// template<typename KmerType>
		CovTree * Parse(
			istream & in,
			const EncodedFastaSequence::Index & seqIdx,
			const KmerIndex & kmerIdx
			// 
		) {
			Parse(in, seqIdx, kmerIdx, 0);
		}

		//template<typename KmerType>
		//CovTree * Parse(
		//	istream & in,
		//	const FastaSequence::Index & seqIdx,
		//	const KmerIndex<Kmer> & kmerIdx,
		//	CovTree * parent
		//	// 
		//)
		//{
		//	string s;
		//	getline( in, s );

		//	if ( s.find( "Cover;" ) == 0 )
		//	{
		//		Distance threshold;
		//		uint branchCount;
		//		uint leafCount;
		//		vector<string> parts = String::Split( s, ";" );
		//		threshold = atoi( parts[0].c_str() );
		//		branchCount = atoi( parts[1].c_str() );
		//		leafCount = atoi( parts[2].c_str() );

		//		Parse( in, seqIdx, kmerIdx, tree );
		//	}
		//	else {
		//		 
		//	}
		//}
	};
} // name space QutBio
