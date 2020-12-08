#pragma once

#define WANT_FRAG 1

#include <FL/Enumerations.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl.H>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <thread>
#include <algorithm>

#include <GridLayout.hpp>
#include <FlowLayout.hpp>
#include <Arg.hpp>
#include <Graph.hpp>
#include <Random.hpp>
#include <Histogram.hpp>
#include <IEventHandler.hpp>
#include <IPropertyChangedEventHandler.hpp>
#include <MarginLayout.hpp>
#include <Rectangle.hpp>
#include <SimpleKmer.hpp>
#include <Button.hpp>
#include <SparseFeatureVector.hpp>
#include <BitSet.hpp>
#include "Helpers.hpp"
#include "RankingBrowserBase.hpp"
#include "RankingBrowser.hpp"
#include "RankingDetailBrowser.hpp"
#include "PrecisionRecallBrowser.hpp"
#include "PrecisionDistanceBrowser.hpp"
#include "PrecisionRecallSummary.hpp"
#include "Page.hpp"
#include "StringChooser.hpp"

using namespace std;
using namespace QutBio;

/**
<summary>
	DbLoader gathers the parameters for an interactive version of DnaClust,
	with some additional analytic displays to render the software a little more
	usable.
</summary>
 */
class RankSimProjIndexed :
	public virtual Page,
	public virtual LBFL::IPropertyChangedEventHandler
	//
{
#if WANT_FRAG
	const int rows = 12;
#else
	const int rows = 11;
#endif

	int rowHeight = 25;

#pragma region MaxMatches Property
private:
	Arg<int, Fl_Int_Input> maxMatches;

public:
	int MaxMatches() {
		return maxMatches.Value();
	}

	void SetMaxMatches( const int val ) {
		if ( val == MaxMatches() ) return;

		this->maxMatches.SetValue( val );
	}
#pragma endregion

#pragma region KmerLength Property
private:
	Arg<uint, Fl_Int_Input> kmerLength;

public:
	uint KmerLength() {
		return kmerLength.Value();
	}

	void SetKmerLength( const uint val ) {
		if ( val == KmerLength() ) return;

		this->kmerLength.SetValue( val );
	}
#pragma endregion

#if WANT_FRAG
#pragma region FragLength Property
private:
	Arg<uint, Fl_Int_Input> fragLength;

public:
	uint FragLength() {
		return fragLength.Value();
	}

	void SetFragLength( const uint val ) {
		if ( val == FragLength() ) return;

		this->fragLength.SetValue( val );
	}
#pragma endregion
#endif

#pragma region ReduceMode Property
private:
	TypedArg<string, StringChooser> reduceMode;

public:
	string ReduceMode() {
		return reduceMode.Value();
	}

	void SetReduceMode( const string & val ) {
		if ( val == ReduceMode() ) return;

		this->reduceMode.SetValue( val );
	}
#pragma endregion

#pragma region CombineMode Property
private:
	TypedArg<string, StringChooser> combineMode;

public:
	string CombineMode() {
		return combineMode.Value();
	}

	void SetCombineMode( const string & val ) {
		if ( val == CombineMode() ) return;

		this->combineMode.SetValue( val );
	}
#pragma endregion

#pragma region RankingFileName Property
private:
	TypedArg<string, OutFileChooser> rankingFileName;

public:
	string RankingFileName() {
		return rankingFileName.Value();
	}

	void SetOutFileName( const string & val ) {
		if ( val == RankingFileName() ) return;

		this->rankingFileName.SetValue( val );
	}
#pragma endregion 

#pragma region private data members
private:
	GridLayout stack;
	Fl_Box titleBox;
	FlowLayout buttons;

	TextDisplay disp;
	Button dispButton;

	mutex mute;
#pragma endregion

#pragma region Data flow connections
	const vector<FastaSequence *> *dbSeqs = 0;
	const LookupTable_<size_t, const FastaSequence> * dbIndex = 0;
	const vector<size_t> * testSetIndex = 0;
	const vector<size_t> * trainingSetIndex = 0;
	const SimilarityMatrix * matrix = 0;
	const vector<SparseSignature> * signatures = 0;
	const vector<vector<size_t>> *featurePostingList = 0;
#pragma endregion

#pragma region Results
	vector<::Ranking> rankings;
#pragma endregion

public:
	RankSimProjIndexed(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Rank SimProj" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, w(), rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),

		rankingFileName(
			"rankingFileName", "", "Ranking File", Requirement::Required,
			"The path to a file which will be overwritten with rankings, in CSV \n"
			"format.",
			labelWidth, 0, width, rowHeight ),

		maxMatches(
			"maxMatches", 100, "Max Matches", Requirement::Required,
			"The maximum number of database sequences to include in ranking report \n"
			"for each query.",
			labelWidth, 0, width, rowHeight ),

		kmerLength(
			"kmerLength", 16, "K-mer Length", Requirement::Required,
			"The word size used for k-mer tiling.",
			labelWidth, 0, width, rowHeight ),
#if WANT_FRAG
		fragLength(
			"fragLength", 1, "Fragment Length", Requirement::Required,
			"The nominal fragment size for fragmentation.",
			labelWidth, 0, width, rowHeight ),
#endif
		reduceMode(
			"reduceMode", "", "Reduce Mode", Requirement::Required,
			"The reduction mode used to collapse row and column minima to produce the respective summary value.",
			labelWidth, 0, width, rowHeight ),

		combineMode(
			"combineMode", "", "Combine Mode", Requirement::Required,
			"The reduction mode used to combine the reduced row and column summaries.",
			labelWidth, 0, width, rowHeight ),

		useDigrams(
			"useDigrams", 0, "Use Digrams?", Requirement::Required,
			"Check this box to use digrams for kmer distance calculations.",
			labelWidth, 0, width, rowHeight ),

		saveMatrix(
			"saveMatrix", 0, "Save Diagnostics?", Requirement::Required,
			"Check this box to save a HUGE AMOUNT OF DATA in file 'SimProj.diagnostics.csv'.",
			labelWidth, 0, width, rowHeight ),

		numThreads(
			"numThreads", 0, "OMP # Threads", Requirement::Required,
			"Number of OpenMP threads for parallel regions.",
			labelWidth, 0, width, rowHeight ),

		autoLoad(
			"autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load the partition file if it exists.",
			labelWidth, 0, width, rowHeight ),

		buttons( 0, 0, w(), rowHeight ),

		disp( 0, 0, w(), 100 ),
		dispButton( 0, 0, 150, rowHeight, "Log", [this]( Button *b ) { ShowCentrePanel( &disp ); } )
		//
	{
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &kmerLength );
#if WANT_FRAG
		stack.add( &fragLength );
#endif
		stack.add( &useDigrams );
		stack.add( &saveMatrix );
		stack.add( &reduceMode );
		stack.add( &combineMode );
		stack.add( &numThreads );
		stack.add( &maxMatches );
		stack.add( &rankingFileName );
		stack.add( &autoLoad );

		vector<string> reduceModes{ "min", "max", "mean" };

		reduceMode.inputField.SetValues( reduceModes );
		combineMode.inputField.SetValues( reduceModes );

		rankingFileName.AddPropertyChangedEventHandler( this );
		kmerLength.AddPropertyChangedEventHandler( this );
#if WANT_FRAG
		fragLength.AddPropertyChangedEventHandler( this );
#endif
		useDigrams.AddPropertyChangedEventHandler( this );
		reduceMode.AddPropertyChangedEventHandler( this );
		combineMode.AddPropertyChangedEventHandler( this );
		maxMatches.AddPropertyChangedEventHandler( this );
		numThreads.AddPropertyChangedEventHandler( this );
		saveMatrix.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &dispButton );

		add( &disp );

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		ShowCentrePanel( &disp );
	}

	virtual ~RankSimProjIndexed() {}

#pragma region Data flow connections to antecedents
	function<decltype(dbSeqs) ()>             getDbSeqs;
	function<decltype(dbIndex) ()>            getDbIndex;
	function<decltype(testSetIndex) ()>       getTestSetIndex;
	function<decltype(trainingSetIndex) ()>   getTrainingSetIndex;
	function<decltype(matrix)()>              getMatrix;
	function<const vector<SparseSignature> * ()>	getSignatures;
	function<const vector<vector<size_t>> *()>      getPostingList;
#pragma endregion

	void GainFocus() override {
		Page::GainFocus();
		dbSeqs = getDbSeqs();
		dbIndex = getDbIndex();
		testSetIndex = getTestSetIndex();
		trainingSetIndex = getTrainingSetIndex();
		matrix = getMatrix();
		signatures = getSignatures();
		featurePostingList = getPostingList();


		UpdateBrowsers();
		UpdateReadiness();
	}

	YesNoCancel CheckLoad() {
		auto fileName = String::Trim( RankingFileName() );

		if ( fileName.length() > 0 && File::Exists( fileName ) ) {
			if ( AutoLoad() ) {
				return YesNoCancel::Yes;
			}

			int answer = fl_choice(
				"The selected file already exists.\n"
				"Would you like to load those rankings?",
				"Cancel operation.",
				"Yes, load the rankings.",
				"No, I want to redo the rankings."
			);

			return (YesNoCancel) answer;
		}
		else {
			return YesNoCancel::No;
		}
	}

	void Run() {
		SetReady( false );
		disp.Clear();
		auto choice = CheckLoad();

		if ( choice == YesNoCancel::Yes ) {
			ResetBrowsers();
			thread t( LoadImpl_, this );
			t.detach();
		}
		else if ( choice == YesNoCancel::No ) {
			ResetBrowsers();
			thread t( RunImpl_, this );
			t.detach();
		}
		else {
			UpdateReadiness();
			disp( "Action cancelled.\n" );
		}
	}

	static void RunImpl_( RankSimProjIndexed * self ) {
		self->RunImpl();
	}

	void RunImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		uint numThreads = std::max( std::min( NumThreads(), thread::hardware_concurrency() ), (uint) 1 );

#if USE_OMP
		omp_set_num_threads( numThreads );
#endif

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

#if INTERLEAVE_SIG
		try {
			Page::runTime = -omp_get_wtime();
			ofstream outFile( RankingFileName() );
			Rank( &outFile );
			Page::runTime += omp_get_wtime();
			Page::saveTime = 0;;
		}
		catch ( Exception & ex ) {
			Page::saveTime += omp_get_wtime();
			ShowProgress( "Error in interleaved ranking: %s\n", ex.what() );
		}
#else
		Page::runTime = -omp_get_wtime();
		Rank( 0 );
		Page::runTime += omp_get_wtime();

		try {
			Page::saveTime = -omp_get_wtime();
			Save();
			Page::saveTime += omp_get_wtime();
		}
		catch ( Exception & ex ) {
			Page::saveTime += omp_get_wtime();
			ShowProgress( "Unable to save rankings: %s\n", ex.what() );
		}
#endif

		UpdateGui();

		mute.unlock();
	}

	void UpdateGui() {
		UPDATE_GUI(
			UpdateReadiness(); \
			NotifyRunComplete();
		);
	}

	void ResetBrowsers() {}

	void UpdateBrowsers() {}

	/**
	 *	Initialise a ranking set ready for RankSignatures.
	 *
	 *	@param[in] database A list of sparse binary signatures, some (or possibly all) of
	 *		which will be selected as query points.
	 *
	 *	@param[in] queryIndices
	 *
	 *	@param[in] maxMatches The number of results required for each query.
	 *
	 *	@param[out] rankings A list of ::Ranking objects which will be prepared
	 *		for RankSignatures.
	 *
	 *	@post rankings.size = queryIndices.size and
	 *		for all i in dom rankings,
	 *			rankings[i].sequence = database[i].sequence and
	 *			rankings[i].knn.capacity = maxMatches and
	 *			rankings[i].precision.size = maxMatches and
	 *			rankings[i].recall.size = maxMatches
	 */
	static void Setup(
		const decltype(RankSimProjIndexed::dbSeqs) database,
		const decltype(RankSimProjIndexed::testSetIndex) queryIndices,
		const decltype(RankSimProjIndexed::trainingSetIndex) referenceIndices,
		int maxMatches,
		vector<::Ranking> & rankings,
		size_t & maxQueryLength,
		size_t & maxReferenceLength
	) {
		rankings.clear();

		KnnVector<const FastaSequence *, double> knn( maxMatches, -1 );
		vector<double> precision;
		vector<double> recall;

		for ( auto i : *queryIndices ) {
			auto & sig = (*database)[i];
			rankings.emplace_back( sig, knn, precision, recall );
			auto & ranking = rankings.back();

			auto &precision = ranking.precision;
			precision.clear();
			precision.reserve( maxMatches );

			auto &recall = ranking.recall;
			recall.clear();
			recall.reserve( maxMatches );
		}

		GetMaxLength( maxQueryLength, queryIndices, database );
		GetMaxLength( maxReferenceLength, referenceIndices, database );
	}

	static void GetMaxLength(
		size_t & maxQueryLength,
		const decltype(RankSimProjIndexed::testSetIndex) queryIndices,
		const decltype(RankSimProjIndexed::dbSeqs) database
	) {
		maxQueryLength = 0;
		for ( auto i : *queryIndices ) {
			auto len = (*database)[i]->Length();

			if ( len > maxQueryLength ) {
				maxQueryLength = len;
			}
		}
	}

	/**
	 *	Identify the k-nearest signatures in the training set (the members of
	 *	which are indexed in a posting list) and rank them in ascending order
	 *	of distance from the designated query. This is intended to be called
	 *	from the body of a parallel loop, so all data structures are required
	 *	to be initialised with necessary capacity before entry.
	 *
	 *	The function will also compute interpolated precision and recall values
	 *	at each recovered document if class labels are available.
	 *
	 *	@param[in] querySig The query signature.
	 *
	 *	@param[in] signatures A list of sparse binary feature vectors.
	 *		Let (s == signatures[i]) be a signature.
	 *		Then (j in s.elements) iff (sequence i contains feature j).
	 *
	 *	@param[in] sigPostingList An inverted index mapping feature numbers to the
	 *		sequences that contain them.
	 *		Let f be a feature.
	 *		Then (i in sigPostingList[f]) iff (f in signatures[i].elements).
	 *
	 *	@param[in] classPostingList An inverted index mapping class labels to the sequences
	 *		that are annotated with those labels.
	 *		Let c be a class.
	 *		Then (i in classPostingList[c]) iff (c in signatures[i].sequence->classes).
	 *
	 *	@param[in/out] processed A bit-set which has been previously initialised
	 *		with capacity to record all values in [0,signatures.size()). Used for
	 *		working storage.
	 *
	 *	@param[in/out] ranking A ::Ranking object with pre-allocated members sufficient
	 *		to store the required number of ranked results, together with any precision
	 *		and recall values that may be computed. The number of results returned
	 *		is equal to the capacity of ranking.knn.
	 *
	 *	@typeparam SigPostingList Object implementing `vector<size_t> & at(size_t)`.
	 */
	template <typename Reduce, typename Combine>
	static void RankQuery(
		const FastaSequence * query,
		const SparseSignature & querySig,
		const decltype(RankSimProjIndexed::dbSeqs) dbSeqs,
		const decltype(RankSimProjIndexed::trainingSetIndex) trainingSetIndex,
		const decltype(RankSimProjIndexed::featurePostingList) featurePostingList,
		uint kmerLength,
		const SimilarityMatrix * matrix,
#if WANT_FRAG
		uint fragLength,
#endif
		Reduce && reduce,
		Combine && combine,
		bool useDigrams,
		ofstream * saveMatrixFile,
		::Ranking & ranking,
		vector<Distance> & rowMin,
		vector<Distance> & colMin,
		BitSet &processed
	) {
		ranking.knn.clear();

#if WANT_FRAG
		auto queryKmerCount = query->KmerCount( kmerLength );
		auto queryFragCount = Fragment::GetCount( queryKmerCount, fragLength );
		auto queryStepSize = Fragment::GetRealStepSize( queryKmerCount, fragLength, queryFragCount );

		rowMin.resize( queryFragCount );
		processed.Clear();
#endif
		for ( auto feature : querySig ) {
			for ( auto i : (*featurePostingList)[feature] ) {
				if ( processed.Contains( i ) ) continue;

				processed.Insert( i );

				auto subject = (*dbSeqs)[i];
#if WANT_FRAG
				auto subjectKmerCount = subject->KmerCount( kmerLength );
				auto subjectFragCount = Fragment::GetCount( subjectKmerCount, fragLength );
				auto subjectStepSize = Fragment::GetRealStepSize( subjectKmerCount, fragLength, subjectFragCount );

				colMin.resize( subjectFragCount );
#else
				vector<Distance> rowMin( query->KmerCount( kmerLength ) );
				vector<Distance> colMin( subject->KmerCount( kmerLength ) );
#endif
				std::fill( rowMin.begin(), rowMin.end(), numeric_limits<Distance>::max() );
				std::fill( colMin.begin(), colMin.end(), numeric_limits<Distance>::max() );

				if ( useDigrams ) {
					if ( saveMatrixFile ) {
						FlatMatrix<Distance> distanceMatrix;
						PopulateDigrams( query->Digrams(), subject->Digrams(), kmerLength, matrix, rowMin, colMin, distanceMatrix );
						SaveDiagnosticData( saveMatrixFile, query, subject, distanceMatrix, rowMin, colMin );
					}
					else {
						PopulateDigrams( query->Digrams(), subject->Digrams(), kmerLength, matrix, rowMin, colMin );
					}
				}
				else {
					if ( saveMatrixFile ) {
						FlatMatrix<Distance> distanceMatrix;
						PopulateBytes( query->Sequence(), subject->Sequence(), kmerLength, matrix, rowMin, colMin, distanceMatrix );
						SaveDiagnosticData( saveMatrixFile, query, subject, distanceMatrix, rowMin, colMin );
					}
					else {
#if WANT_FRAG
						PopulateBytes( query->Sequence(), subject->Sequence(), kmerLength, matrix,
							[&rowMin, &colMin, queryStepSize, subjectStepSize]( long long r, long long c, Distance distance ) {
							auto queryFragIdx = (size_t) (r / queryStepSize);
							auto subjectFragIdx = (size_t) (c / subjectStepSize);
							if ( distance < rowMin[queryFragIdx] ) rowMin[queryFragIdx] = distance;
							if ( distance < colMin[subjectFragIdx] ) colMin[subjectFragIdx] = distance;
						} );
#else
						PopulateBytes( query->Sequence(), subject->Sequence(), kmerLength, matrix, rowMin, colMin );
#endif
					}
				}

				double distance = combine( reduce( rowMin ), reduce( colMin ) );

				if ( ranking.knn.canPush( distance ) ) {
					ranking.knn.push( subject, distance );
				}
			}
		}

		ranking.knn.sort();
	}

	static void SaveDiagnosticData(
		ofstream * saveMatrixFile,
		const FastaSequence * query,
		const FastaSequence * subject,
		const FlatMatrix<Distance> & distanceMatrix,
		const vector<Distance> & rowMin,
		const vector<Distance> & colMin
	) {
#if USE_OMP
#pragma omp critical
#endif
		{
			(*saveMatrixFile) << query->IdStr() << "," << subject->IdStr() << "," << distanceMatrix.rows() << "," << distanceMatrix.cols() << "\n" << distanceMatrix << "\n";
			(*saveMatrixFile) << query->IdStr() << "," << subject->IdStr() << "," << "RowMinima";
			for ( auto x : rowMin ) {
				(*saveMatrixFile) << "," << x;
			}
			(*saveMatrixFile) << "\n" << query->IdStr() << "," << subject->IdStr() << "," << "ColMinima";
			for ( auto x : colMin ) {
				(*saveMatrixFile) << "," << x;
			}
			(*saveMatrixFile) << "\n\n";
		}
	}

	static void PopulateDigrams(
		const vector<Digram> & query,
		const vector<Digram> & subject,
		uint kmerLength,
		const SimilarityMatrix * matrix,
		vector<Distance> & rowMin,
		vector<Distance> & colMin
	) {
		std::fill( rowMin.begin(), rowMin.end(), numeric_limits<Distance>::max() );
		std::fill( colMin.begin(), colMin.end(), numeric_limits<Distance>::max() );

		// Special treatment necessary for digrams: note +2 instead of +1
		const long long m = query.size() - kmerLength + 2;
		const long long n = subject.size() - kmerLength + 2;

		if ( m < 1 || n < 1 ) {
			return;
		}

		for ( long long r = 0; r < m; r++ ) {
			// if ( r == m - 1 ) {
			// 	cerr << "here\n";
			// }

			const int c_upper = r == 0 ? (int) n : 1;

			//	A diagonal run starts at (r,c).
			//	r traverses the query sequence.
			//	When r == 0, c traverses the subject sequence.
			//	When r > 0, we just start a run from the beginning of the subject series
			for ( int c = 0; c < c_upper; c++ ) {
				long long diagLength = std::min( m - r, n - c );

				if ( diagLength < 0 ) break;

				for ( int interleave = 0; interleave <= 1 && interleave < diagLength; interleave++ ) {

					Distance buffer[1000];
					const Digram * a = subject.data() + c + interleave;
					const Digram * b = query.data() + r + interleave;
					Distance distance = 0;

					// Prime the circular buffer with the first kmer in the query
					for ( size_t t = 0; t < kmerLength / 2; t++, a += 2, b += 2 ) {
						Distance currentTerm = matrix->DigramDifference( *a, *b );
						distance += currentTerm;
						buffer[t] = currentTerm;
					}

					// process( processingObject, r, c, distance );
					if ( distance < rowMin[r + interleave] ) rowMin[r + interleave] = distance;
					if ( distance < colMin[c + interleave] ) colMin[c + interleave] = distance;

					for ( long long offset = 2, buffptr = 0;
						offset + interleave < diagLength;
						a += 2, b += 2, offset += 2, buffptr++
						) {
						if ( buffptr >= kmerLength / 2 ) {
							buffptr = 0;
						}

						distance -= buffer[buffptr];
						Distance currentTerm = matrix->DigramDifference( *a, *b );
						buffer[buffptr] = currentTerm;
						distance += currentTerm;

						// process( processingObject, r + offset, c + offset, distance );
						if ( distance < rowMin[r + offset + interleave] ) rowMin[r + offset + interleave] = distance;
						if ( distance < colMin[c + offset + interleave] ) colMin[c + offset + interleave] = distance;
					}

				}

			}
		}
	}

	static void PopulateDigrams(
		const vector<Digram> & query,
		const vector<Digram> & subject,
		uint kmerLength,
		const SimilarityMatrix * matrix,
		vector<Distance> & rowMin,
		vector<Distance> & colMin,
		FlatMatrix<Distance> & distanceMatrix
	) {
		const auto outrageous = numeric_limits<Distance>::max();
		distanceMatrix.resize( query.size(), subject.size() );
		distanceMatrix.fill( outrageous );
		std::fill( rowMin.begin(), rowMin.end(), outrageous );
		std::fill( colMin.begin(), colMin.end(), outrageous );

		// Special treatment necessary for digrams: note +2 instead of +1
		const long long m = query.size() - kmerLength + 2;
		const long long n = subject.size() - kmerLength + 2;

		if ( m < 1 || n < 1 ) {
			return;
		}

		for ( long long r = 0; r < m; r++ ) {
			// if ( r == m - 1 ) {
			// 	cerr << "here\n";
			// }

			const int c_upper = r == 0 ? (int) n : 1;

			//	A diagonal run starts at (r,c).
			//	r traverses the query sequence.
			//	When r == 0, c traverses the subject sequence.
			//	When r > 0, we just start a run from the beginning of the subject series
			for ( int c = 0; c < c_upper; c++ ) {
				long long diagLength = std::min( m - r, n - c );

				if ( diagLength < 0 ) break;

				for ( int interleave = 0; interleave <= 1 && interleave < diagLength; interleave++ ) {

					Distance buffer[1000];
					const Digram * a = subject.data() + c + interleave;
					const Digram * b = query.data() + r + interleave;
					Distance distance = 0;

					// Prime the circular buffer with the first kmer in the query
					for ( size_t t = 0; t < kmerLength / 2; t++, a += 2, b += 2 ) {
						Distance currentTerm = matrix->DigramDifference( *a, *b );
						distance += currentTerm;
						buffer[t] = currentTerm;
					}

					// process( processingObject, r, c, distance );
					distanceMatrix( r + interleave, c + interleave ) = distance;
					if ( distance < rowMin[r + interleave] ) rowMin[r + interleave] = distance;
					if ( distance < colMin[c + interleave] ) colMin[c + interleave] = distance;

					for ( long long offset = 2, buffptr = 0;
						offset + interleave < diagLength;
						a += 2, b += 2, offset += 2, buffptr++
						) {
						if ( buffptr >= kmerLength / 2 ) {
							buffptr = 0;
						}

						distance -= buffer[buffptr];
						Distance currentTerm = matrix->DigramDifference( *a, *b );
						buffer[buffptr] = currentTerm;
						distance += currentTerm;

						// process( processingObject, r + offset, c + offset, distance );
						distanceMatrix( r + offset + interleave, c + offset + interleave ) = distance;
						if ( distance < rowMin[r + offset + interleave] ) rowMin[r + offset + interleave] = distance;
						if ( distance < colMin[c + offset + interleave] ) colMin[c + offset + interleave] = distance;
					}
				}
			}
		}
	}

	static void PopulateBytes(
		const vector<Symbol> & query,
		const vector<Symbol> & subject,
		uint kmerLength,
		const SimilarityMatrix * matrix,
#if WANT_FRAG
		function<void( long long r, long long c, Distance distance )> process
#else
		vector<Distance> & rowMin,
		vector<Distance> & colMin
#endif
	) {
#if ! WANT_FRAG
		std::fill( rowMin.begin(), rowMin.end(), numeric_limits<Distance>::max() );
		std::fill( colMin.begin(), colMin.end(), numeric_limits<Distance>::max() );
#endif

		const long long m = query.size() - kmerLength + 1;
		const long long n = subject.size() - kmerLength + 1;

		if ( m < 1 || n < 1 ) {
			return;
		}

		for ( long long r = 0; r < m; r++ ) {
			const int c_upper = r == 0 ? (int) n : 1;

			//	A diagonal run starts at (r,c).
			//	r traverses the query sequence.
			//	When r == 0, c traverses the subject sequence.
			//	When r > 0, we just start a run from the beginning of the subject series
			for ( int c = 0; c < c_upper; c++ ) {
				Distance buffer[1000];
				const Symbol * a = subject.data() + c;
				const Symbol * b = query.data() + r;
				Distance distance = 0;

				size_t diagLength = std::min( m - r, n - c );

				// Prime the circular buffer with the first kmer in the query
				for ( size_t t = 0; t < kmerLength; t++, a++, b++ ) {
					Distance currentTerm = matrix->Difference( *a, *b );
					distance += currentTerm;
					buffer[t] = currentTerm;
				}

#if WANT_FRAG
				process( r, c, distance );
#else
				if ( distance < rowMin[r] ) rowMin[r] = distance;
				if ( distance < colMin[c] ) colMin[c] = distance;
#endif

				for ( size_t offset = 1, buffptr = 0;
					offset < diagLength;
					a++, b++, offset++, buffptr++
					) {
					if ( buffptr >= kmerLength ) {
						buffptr = 0;
					}

					distance -= buffer[buffptr];
					Distance currentTerm = matrix->Difference( *a, *b );
					buffer[buffptr] = currentTerm;
					distance += currentTerm;

#if WANT_FRAG
					process( r + offset, c + offset, distance );
#else
					if ( distance < rowMin[r + offset] ) rowMin[r + offset] = distance;
					if ( distance < colMin[c + offset] ) colMin[c + offset] = distance;
#endif
				}
			}

		}
	}

	static void PopulateBytes(
		const vector<Symbol> & query,
		const vector<Symbol> & subject,
		uint kmerLength,
		const SimilarityMatrix * matrix,
		vector<Distance> & rowMin,
		vector<Distance> & colMin,
		FlatMatrix<Distance> & distanceMatrix
	) {
		const auto outrageous = numeric_limits<Distance>::max();
		distanceMatrix.resize( query.size(), subject.size() );
		distanceMatrix.fill( outrageous );
		std::fill( rowMin.begin(), rowMin.end(), outrageous );
		std::fill( colMin.begin(), colMin.end(), outrageous );

		const long long m = query.size() - kmerLength + 1;
		const long long n = subject.size() - kmerLength + 1;

		if ( m < 1 || n < 1 ) {
			return;
		}

		for ( long long r = 0; r < m; r++ ) {
			const int c_upper = r == 0 ? (int) n : 1;

			//	A diagonal run starts at (r,c).
			//	r traverses the query sequence.
			//	When r == 0, c traverses the subject sequence.
			//	When r > 0, we just start a run from the beginning of the subject series
			for ( int c = 0; c < c_upper; c++ ) {
				Distance buffer[1000];
				const Symbol * a = subject.data() + c;
				const Symbol * b = query.data() + r;
				Distance distance = 0;

				auto diagLength = std::min( m - r, n - c );

				// Prime the circular buffer with the first kmer in the query
				for ( uint t = 0; t < kmerLength; t++, a++, b++ ) {
					Distance currentTerm = matrix->Difference( *a, *b );
					distance += currentTerm;
					buffer[t] = currentTerm;
				}

				// process( processingObject, r, c, distance );
				distanceMatrix( r, c ) = distance;
				if ( distance < rowMin[r] ) rowMin[r] = distance;
				if ( distance < colMin[c] ) colMin[c] = distance;

				for ( long long offset = 1, buffptr = 0;
					offset < diagLength;
					a++, b++, offset++, buffptr++
					) {
					if ( buffptr >= kmerLength ) {
						buffptr = 0;
					}

					distance -= buffer[buffptr];
					Distance currentTerm = matrix->Difference( *a, *b );
					buffer[buffptr] = currentTerm;
					distance += currentTerm;

					// process( processingObject, r + offset, c + offset, distance );
					distanceMatrix( r + offset, c + offset ) = distance;
					if ( distance < rowMin[r + offset] ) rowMin[r + offset] = distance;
					if ( distance < colMin[c + offset] ) colMin[c + offset] = distance;
				}
			}

		}
	}

	static double ReduceMin( const vector<Distance> & minima ) {
		if ( minima.size() == 0 ) return numeric_limits<Distance>::max();

		auto begin = minima.begin();
		auto end = minima.end();
		double ret = *begin;

		for ( auto p = begin + 1; p != end; p++ ) {
			double t = *p;
			if ( t < ret ) ret = t;
		}

		return ret;
	}

	static double ReduceMax( const vector<Distance> & minima ) {
		if ( minima.size() == 0 ) return numeric_limits<Distance>::max();

		auto begin = minima.begin();
		auto end = minima.end();
		double ret = *begin;

		for ( auto p = begin + 1; p != end; p++ ) {
			double t = *p;
			if ( t > ret ) ret = t;
		}

		return ret;
	}

	static double ReduceMean( const vector<Distance> & minima ) {
		if ( minima.size() == 0 ) return numeric_limits<Distance>::max();

		auto begin = minima.begin();
		auto end = minima.end();
		double ret = *begin;

		for ( auto p = begin + 1; p != end; p++ ) {
			double t = *p;
			ret += t;
		}

		return ret / minima.size();
	}

	static double CombineMin( double x, double y ) {
		return std::min( x, y );
	}

	static double CombineMax( double x, double y ) {
		return std::max( x, y );
	}

	static double CombineMean( double x, double y ) {
		return (x + y) / 2;
	}

	using ProgressFunc = function<void( int & )>;

	template<typename R, typename C>
	void Rank( ProgressFunc && progress, R && reduce, C && combine, CsvWriter * rankingWriter ) {
		int done = 0;
		bool useDigrams = UseDigrams();
#if WANT_FRAG
		auto fragLength = FragLength();
#endif

#if USE_OMP
		auto maxThreads = std::thread::hardware_concurrency();
		auto numThreads = NumThreads();

		if ( numThreads <= 0 || numThreads > maxThreads ) {
			numThreads = maxThreads;
		}

		omp_set_num_threads( numThreads );
		disp( "numThreads = %u\n", numThreads );

#pragma omp parallel
#endif
		{
			vector<Distance> rowMin( maxQuerySize );
			vector<Distance> colMin( maxReferenceSize );
			BitSet processed( dbSeqs->size() );

#if USE_OMP
#pragma omp for schedule(dynamic, 1)
#endif
			for ( size_t q = 0; q < testSetIndex->size(); q++ ) {
				auto queryIdx = (*testSetIndex)[q];
				auto query = (*dbSeqs)[queryIdx];
				auto & querySig = (*signatures)[queryIdx];

				RankQuery(
					query, querySig, dbSeqs, trainingSetIndex, featurePostingList, kmerLength.Value(), matrix,
#if WANT_FRAG
					fragLength,
#endif
					reduce, combine, useDigrams, saveMatrixFile, rankings[q], rowMin, colMin, processed
				);

				if ( rankingWriter ) {
#pragma omp critical
					SigRank::SaveRanking( *rankingWriter, rankings[q] );
				}

				progress( done );
			}
		}
	}

	ofstream * saveMatrixFile = 0;

	size_t maxQuerySize;
	size_t maxReferenceSize;

	void Rank( ofstream * outFile = 0 ) {
		CsvWriter * rankingWriter = outFile ? new CsvWriter( *outFile ) : 0;

		const int maxMatches = MaxMatches();
		const size_t Q = (*testSetIndex).size();

		if ( SaveMatrix() ) {
			saveMatrixFile = new ofstream( "SimProj.diagnostics.csv" );
		}

		Setup( dbSeqs, testSetIndex, trainingSetIndex, maxMatches, rankings, maxQuerySize, maxReferenceSize );

		auto progress = [this, Q]( int & done ) {
			done++;
#if WANT_FEEDBACK
			int thisPercent = done * 100 / Q, lastPercent = (done - 1) * 100 / Q;

			if ( thisPercent != lastPercent ) {
				UpdateGuiProgress( thisPercent );
			}
#endif
		};

		if ( ReduceMode() == "min" && CombineMode() == "min" )   Rank( progress, ReduceMin, CombineMin, rankingWriter );
		else if ( ReduceMode() == "min" && CombineMode() == "max" )   Rank( progress, ReduceMin, CombineMax, rankingWriter );
		else if ( ReduceMode() == "min" && CombineMode() == "mean" )  Rank( progress, ReduceMin, CombineMean, rankingWriter );
		else if ( ReduceMode() == "max" && CombineMode() == "min" )   Rank( progress, ReduceMax, CombineMin, rankingWriter );
		else if ( ReduceMode() == "max" && CombineMode() == "max" )   Rank( progress, ReduceMax, CombineMax, rankingWriter );
		else if ( ReduceMode() == "max" && CombineMode() == "mean" )  Rank( progress, ReduceMax, CombineMean, rankingWriter );
		else if ( ReduceMode() == "mean" && CombineMode() == "min" )  Rank( progress, ReduceMean, CombineMin, rankingWriter );
		else if ( ReduceMode() == "mean" && CombineMode() == "max" )  Rank( progress, ReduceMean, CombineMax, rankingWriter );
		else if ( ReduceMode() == "mean" && CombineMode() == "mean" ) Rank( progress, ReduceMean, CombineMean, rankingWriter );

		UpdateBrowsers();

		if ( saveMatrixFile ) {
			saveMatrixFile->close();
			delete saveMatrixFile;
			saveMatrixFile = 0;
		}
	}

	void UpdateGuiProgress( int thisPercent ) {
		UPDATE_GUI( disp( "Ranking training set: %d%%\n", thisPercent ); )
	}

	void Save() {
		auto fileName = String::Trim( RankingFileName() );
		SigRank::Save( rankings, fileName );
	}

	static void LoadImpl_( RankSimProjIndexed * self ) {
		self->LoadImpl();
	}

	void LoadImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

		try {
			Page::loadTime = -omp_get_wtime();
			Load();
			Page::loadTime += omp_get_wtime();
		}
		catch ( Exception & ex ) {
			ShowProgress( "Error loading rankings: %s\n", ex.what() );
		}

		UpdateGui();
		mute.unlock();
	}

	void Load() {
		auto fileName = String::Trim( RankingFileName() );
		auto errorMsg = [this]( const string & tag ) {
			UPDATE_GUI(
				disp( "Expected 'rankings' tag in dataset, but found %s.\nLoad abandoned.\n", tag.c_str() ); \
				UpdateReadiness();
			);
		};
		auto setMaxMatches = [this]( size_t maxMatches ) {
			UPDATE_GUI( SetMaxMatches( maxMatches ); disp( "%zu rankings loaded.\n", rankings.size() ); );
		};

		ShowProgress( "Attempting to load rankings from %s:\n", fileName.c_str() );
		SigRank::Load( fileName, *dbIndex, 0, errorMsg, setMaxMatches, rankings );

		UpdateBrowsers();
	}

	void ShowProgress( const char * format, ... ) {
		Fl::lock();
		va_list args;
		va_start( args, format );
		disp( format, args );
		va_end( args );
		Fl::unlock();
		Fl::awake();
	}

	void UpdateReadiness() {
		auto trainingSetSize = trainingSetIndex == 0 ? 0 : trainingSetIndex->size();
		auto testSetSize = testSetIndex == 0 ? 0 : testSetIndex->size();
		auto dbSeqSize = dbSeqs == 0 ? 0 : dbSeqs->size();
		auto fileName = String::Trim( RankingFileName() );

		SetReady( trainingSetSize > 0 && testSetSize > 0 && fileName.length() > 0 && matrix != 0 );
		SetOk( Ready() && rankings.size() == testSetSize );
	}

	void PropertyChanged( void *sender, const string &propertyName ) override {
		// Handle property changed events...
		UpdateReadiness();
		redraw();
	}

	/// Get parameters so they can be serialised when program finishes.
	/// @param parms A Param list to which the component should append its current parameters. 
	void GetParams( set<Param> &parms ) override {
#define addPar(x) parms.emplace( Name(), STRING(x), (x).Value() )
#define addVar(x,fmt) parms.emplace( Name(), STRING(x), String::Format(fmt,x) )
		addPar( maxMatches );
		addPar( rankingFileName );
		addPar( autoLoad );
		addPar( useDigrams );
		addPar( saveMatrix );
		addPar( kmerLength );
#if WANT_FRAG
		addPar( fragLength );
#endif
		addPar( reduceMode );
		addPar( combineMode );
		addPar( numThreads );
		addVar( runTime, "%0.17g" );
		addVar( loadTime, "%0.17g" );
		addVar( saveTime, "%0.17g" );
#undef add_parm
#undef addVar
	}

	/// Set parameters which have been de-serialised when program starts.
	/// @param parms A Param list containing de-serialised parameters.
	void SetParams( const set<Param> & parms ) override {
#define setPar(x,t) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value<t>() )
#define setVar(x,t) if ( p.ParamName() == STRING(x) ) (x) = p.Value<t>()
		for ( auto & p : parms ) {
			if ( p.ComponentName() != Name() ) continue;
			setPar( maxMatches, int );
			setPar( rankingFileName, string );
			setPar( autoLoad, int );
			setPar( useDigrams, int );
			setPar( saveMatrix, int );
			setPar( kmerLength, uint );
#if WANT_FRAG
			setPar( fragLength, uint );
#endif
			setPar( numThreads, uint );
			setVar( runTime, double );
			setVar( loadTime, double );
			setVar( saveTime, double );
			setPar( reduceMode, string );
			setPar( combineMode, string );
		}

		UpdateReadiness();
#undef setPar
#undef setVar
	}

public:
	const decltype(rankings) * Rankings() {
		return &rankings;
	}

private:
	TypedArg<int, Fl_Check_Button> autoLoad;

public:
	int AutoLoad() {
		return autoLoad.Value();
	}

	void SetAutoLoad( const int val ) {
		if ( val == AutoLoad() ) return;

		this->autoLoad.SetValue( val );
		NotifyPropertyChanged( NAMEOF( AutoLoad ) );
	}


#pragma region NumThreads Property
private:
	Arg<uint, Fl_Int_Input> numThreads;

public:
	uint NumThreads() {
		return numThreads.Value();
	}

	void SetNumThreads( const uint val ) {
		if ( val == NumThreads() ) return;

		this->numThreads.SetValue( val );
	}
#pragma endregion


#pragma region UseDigrams Property
private:
	TypedArg<int, Fl_Check_Button> useDigrams;

public:
	int UseDigrams() {
		return useDigrams.Value();
	}

	void SetUseDigrams( const int val ) {
		if ( val == UseDigrams() ) return;

		this->useDigrams.SetValue( val );
	}
#pragma endregion


#pragma region SaveMatrix Property
private:
	TypedArg<int, Fl_Check_Button> saveMatrix;

public:
	int SaveMatrix() {
		return saveMatrix.Value();
	}

	void SetSaveMatrix( const int val ) {
		if ( val == SaveMatrix() ) return;

		this->saveMatrix.SetValue( val );
	}
#pragma endregion


public:
	void Reset() {
		ResetBrowsers();

		dbSeqs = 0;
		dbIndex = 0;
		testSetIndex = 0;
		trainingSetIndex = 0;
		matrix = 0;
	}
};
