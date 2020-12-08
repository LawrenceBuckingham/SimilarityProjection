#pragma once

#include "FastaSequence.hpp"
#include "SequenceWrapper.hpp"
#include "DistanceType.hpp"

namespace QutBio {
	using EncodingMatrix = vector<vector<KmerWord>>;

	typedef class EncodedFastaSequence *pEncodedFastaSequence;

	class EncodedFastaSequence: public SequenceWrapper {
	protected:
		string classLabel;

		vector<uint64_t> embedding;
		pAlphabet alphabet;
		size_t charsPerWord = 0;
		size_t kmerLength;
		char defaultSymbol;

	public:
		// I know, this is messy as anything, but pending a total rewrite, this is what we have for now.
		// Position, rowMinima, colMinima, fragDistAuxData and hit are used by database sequences to store
		// information used in the collection-wide kmer sweep.
		size_t position;
		vector<Distance> rowMinima;
		vector<Distance> colMinima;
		void *fragDistAuxData = 0;
		function<void(void *p)> deleteFragDistAuxData;
		vector<pEncodedFastaSequence> homologs;
		vector<int> classNumbers;

		// When I'm using the cached 2-mer and/or 3-mer score tables, this is the packed image of the kmers in the collection.
		// Update (hacky!) When switching to the centroid-based codebook, I need to access one-char and two-char encodings.
		EncodingMatrix encoding1, encoding2;

		EncodedFastaSequence(const EncodedFastaSequence & other) = delete;

		EncodedFastaSequence & operator= (const EncodedFastaSequence & other) = delete;

		EncodedFastaSequence & operator= (const EncodedFastaSequence && other) {
			classLabel = std::move(other.classLabel);
			embedding = std::move(other.embedding);
			alphabet = other.alphabet;
			charsPerWord = other.charsPerWord;
			kmerLength = other.kmerLength;
			defaultSymbol = other.defaultSymbol;
			base = other.base;
			position = other.position;
			rowMinima = std::move(other.rowMinima);
			colMinima = std::move(other.colMinima);
			fragDistAuxData = other.fragDistAuxData;
			deleteFragDistAuxData = other.deleteFragDistAuxData;
			homologs = std::move(other.homologs);
			classNumbers = std::move(classNumbers);
			encoding1 = std::move(other.encoding1);
			encoding2 = std::move(other.encoding2);
			return *this;
		}

		EncodedFastaSequence(const EncodedFastaSequence && other) :
			SequenceWrapper(base),
			classLabel( std::move(other.classLabel)),
			embedding(std::move(other.embedding)),
			alphabet(other.alphabet),
			charsPerWord(other.charsPerWord),
			kmerLength(other.kmerLength),
			defaultSymbol(other.defaultSymbol),
			position(other.position),
			rowMinima(std::move(other.rowMinima)),
			colMinima(std::move(other.colMinima)),
			fragDistAuxData(other.fragDistAuxData),
			deleteFragDistAuxData(other.deleteFragDistAuxData),
			homologs(std::move(other.homologs)),
			classNumbers(std::move(classNumbers)),
			encoding1(std::move(other.encoding1)),
			encoding2(std::move(other.encoding2))
		{}

		static EncodedFastaSequence & Zero() {
			static EncodedFastaSequence zero(FastaSequence::Zero(), -1, nullptr, 0, 0, Symbol::From(0));
			return zero;
		}

		const string &ClassLabel() const { return classLabel; }
		const vector<uint64_t> Embedding() { return embedding; }

		static unordered_map<string, int> &ClassNumberRegistry() {
			static unordered_map<string, int> classNumberRegistry;
			return classNumberRegistry;
		}

		static vector<string> &ClassNameRegistry() {
			static vector<string> classNameRegistry;
			return classNameRegistry;
		}

		static int GetClassId(const string &classLabel) {
			auto &&classNumberRegistry = ClassNumberRegistry();
			int classNumber = -1;

			auto existingRecord = classNumberRegistry.find(classLabel);

			if (existingRecord == classNumberRegistry.end()) {
				classNumber = (int)classNumberRegistry.size();
				classNumberRegistry.emplace(classLabel, classNumber);
				auto &&classNameRegistry = ClassNameRegistry();
				classNameRegistry.push_back(classLabel);
			}
			else {
				classNumber = existingRecord->second;
			}

			return classNumber;
		}

		EncodedFastaSequence(
			FastaSequence * charData,
			int classIndex,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			Symbol defaultSymbol
		) : SequenceWrapper(charData) {
			Init(classIndex, alphabet, kmerLength, charsPerWord, defaultSymbol);
		}

		void Init(
			int classIndex,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			Symbol defaultSymbol
		) {
			this->classLabel = classIndex >= 0 ? base->Metadata(classIndex) : "";

			// cerr << "classLabel = '" << classLabel << "'\n";

			if (classLabel.size() > 0) {
				auto splitClassLabels = String::Split(classLabel, ';');

				for (string &classLabel : splitClassLabels) {
					int classNumber = GetClassId(classLabel);
					classNumbers.push_back(classNumber);
				}
			}

			Encode(alphabet, kmerLength, charsPerWord, defaultSymbol);
		}

		virtual ~EncodedFastaSequence() {
			if (fragDistAuxData) {
				deleteFragDistAuxData(fragDistAuxData);
			}
		}

		/**
		*	<summary>
		*		Returns true iff a designated sequence appears in the
		*		list of homologs of this sequence.
		*		(Needless to say, this is not to be used in search/classification.)
		*	</summary>
		*	<param name="other">A candidate homolog.</param>
		*/
		bool IsHomolog(const EncodedFastaSequence &other) {
			if (homologs.size() > 0) {
				return find(homologs.begin(), homologs.end(), &other) != homologs.end();
			}
			else {
				for (auto i : classNumbers) {
					for (auto j : other.classNumbers) {
						if (i == j)
							return true;
					}
				}
			}

			return false;
		}

		bool IsHomologAny(const vector<const EncodedFastaSequence *> &others) {
			for (auto other : others) {
				if (IsHomolog(*other)) {
					return true;
				}
			}

			return false;
		}

		void SetEmbedding(const CharMap &charMap) {
			auto len = base->Length();
			auto sequence = base->Sequence();
			embedding.resize(len);

			for (size_t i = 0; i < len; i++) {
				embedding[i] = charMap.at(sequence[i]).lo;
			}
		}

		/***
		*	<summary>
		*		Uses the supplied selector to choose zero or more kmers from this sequence and pass them back for processing.
		*	</summary>
		*	<param name="kmerLength">The size of the desired kmers.</param>
		*	<param name="selector">A reference to a Selector that determines which, if any, kmers are processed.</param>
		*	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		***/

		void SelectKmers(
			size_t kmerLength,
			Selector &selector,
			function<void(EncodedFastaSequence *seq, size_t pos, size_t length)> process
		) {
			auto length = base->Length();

			if (length < kmerLength)
				return;

			size_t n = length - kmerLength + 1;

			for (size_t i = 0; i < n; i++) {
				if (selector.SelectThis()) {
					process(this, i, kmerLength);
				}
			}
		}

		/***
		*	<summary>
		*		Iterates over the kmers in this sequence and passes them back for processing.
		*	</summary>
		*	<param name="kmerLength">The size of the desired kmers.</param>
		*	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		***/

		void SelectKmers(size_t kmerLength, function<void(EncodedFastaSequence *seq, size_t pos, size_t length)> process) {
			size_t n = base->KmerCount(kmerLength);

			for (size_t i = 0; i < n; i++) {
				process(this, i, kmerLength);
			}
		}

		/**
		*	<summary>
		*		Packs each kmer in the present sequence into a list of KmerWord arrays with the
		*		designated length and density. Each kmer is then represented as one row in the
		*		encoding member of this sequence.
		*
		*		As of January 2017, I need to have simultaneous access to both 1- and 2-letter
		*		numeric codes (1 for the KmerCLusterAL centroid distance calc, and 2 for the
		*		Kmer-vs-Kmer distance calc. So I create two encodings.
		*
		*	</summary>
		*/

		void Encode(Alphabet *alphabet, size_t kmerLength, size_t charsPerWord, Symbol defaultSymbol ) {
			base->EnsureLengthAtLeast( kmerLength, alphabet->Decode(defaultSymbol) );

			auto sequence = base->Sequence().data();
			auto len = base->Length();

			this->charsPerWord = charsPerWord;
			this->kmerLength = kmerLength;
			base->Pad(kmerLength, defaultSymbol);
			alphabet->Encode(sequence, len, kmerLength, 1, encoding1);

			if (charsPerWord > 1) {
				alphabet->Encode(sequence, len, kmerLength, charsPerWord, encoding2);
			}
		}

		EncodedKmer GetEncodedKmer(size_t pos) {
			return charsPerWord == 0 ? GetEncodedKmerError(pos) : charsPerWord == 1 ? GetEncodedKmer1(pos) : charsPerWord == 2 ? GetEncodedKmer2(pos) : charsPerWord == 3 ? GetEncodedKmer3(pos) : GetEncodedKmerGeneral(pos);
		}

		EncodedKmer GetEncodedKmerGeneral(size_t pos) {
			return kmerLength <= charsPerWord
				? &encoding2[0][pos]
				: &encoding2[pos % charsPerWord][pos / charsPerWord];
		}

		EncodedKmer GetEncodedKmer1(size_t pos) {
			return &encoding1[0][pos];
		}

		EncodedKmer GetEncodedKmer2(size_t pos) {
			return kmerLength <= charsPerWord ? GetEncodedKmer1(pos) : &encoding2[pos % 2][pos / 2];
		}

		EncodedKmer GetEncodedKmer3(size_t pos) {
			return kmerLength <= charsPerWord ? GetEncodedKmer1(pos) : &encoding2[pos % 3][pos / 3];
		}

		EncodedKmer GetEncodedKmerError(size_t pos) {
			(cerr << "GetEncodedKmerError -- this function should never be called.\n").flush();
			return 0;
		}
	};

	struct pFastaHash {
		/**
		 *	Gets a hash code from the Id of the base sequence.
		 */
		size_t operator()(const pEncodedFastaSequence &__val) const noexcept {
			auto &s = __val->base->IdStr();
			return _Hash_impl::hash((void *)s.c_str(), s.length());
		}
	};

	struct Subsequence {
		EncodedFastaSequence * source;
		size_t start;
		size_t length;
	};
}
