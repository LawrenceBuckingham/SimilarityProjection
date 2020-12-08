#pragma once

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
#include "Helpers.hpp"
#include "RankingBrowserBase.hpp"
#include "RankingBrowser.hpp"
#include "RankingDetailBrowser.hpp"
#include "PrecisionRecallBrowser.hpp"
#include "PrecisionDistanceBrowser.hpp"
#include "PrecisionRecallSummary.hpp"
#include "Page.hpp"

using namespace std;
using namespace QutBio;

/**
<summary>
	DbLoader gathers the parameters for an interactive version of DnaClust,
	with some additional analytic displays to render the software a little more
	usable.
</summary>
 */
class RankBagOfWords :
	public virtual Page,
	public virtual LBFL::IPropertyChangedEventHandler
	//
{
	const int rows = 6;
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

	PrecisionRecallSummary summary;
	Button summaryButton;

	RankingBrowser browser;
	Button browserButton;

	PrecisionRecallBrowser precisionRecallBrowser;
	Button precisionRecallBrowserButton;

	PrecisionDistanceBrowser precisionDistanceBrowser;
	Button precisionDistanceBrowserButton;

	mutex mute;
#pragma endregion

#pragma region Data flow connections
	const vector<SparseFeatureVector> * featureVectors = 0;
	const vector<FastaSequence *> *dbSeqs = 0;
	const LookupTable_<size_t, const FastaSequence> * dbIndex = 0;
	const vector<size_t> * testSetIndex = 0;
	const vector<size_t> * trainingSetIndex = 0;
	const vector<vector<size_t>> *featurePostingList = 0;
	const unordered_map<size_t, vector<size_t>> *classPostingList = 0;
#pragma endregion

#pragma region Results
	vector<::Ranking> rankings;
	PrecisionRecallStats stats;
	double meanAveragePrecision;
#pragma endregion

public:
	RankBagOfWords(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Rank D2" ),
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

		autoLoad(
			"autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load the partition file if it exists.",
			labelWidth, 0, width, rowHeight ),

		precisionSteps(
			"precisionSteps", 100, "Precision Steps", Requirement::Required,
			"The number of steps for precision-recall summary.",
			labelWidth, 0, width, rowHeight ),

		buttons( 0, 0, w(), rowHeight ),

		disp( 0, 0, w(), 100 ),
		dispButton( 0, 0, 150, rowHeight, "Log", [this]( Button *b ) { ShowCentrePanel( &disp ); } ),

		browser(),
		browserButton( 0, 0, 150, rowHeight, "Overview", [this]( Button *b ) { ShowCentrePanel( &browser ); } ),

		precisionRecallBrowser(),
		precisionRecallBrowserButton( 0, 0, 150, rowHeight, "Prec vs Recall", [this]( Button *b ) { ShowCentrePanel( &precisionRecallBrowser ); } ),

		precisionDistanceBrowser(),
		precisionDistanceBrowserButton( 0, 0, 150, rowHeight, "Prec vs Distance", [this]( Button *b ) { ShowCentrePanel( &precisionDistanceBrowser ); } ),

		summary( stats ),
		summaryButton( 0, 0, 150, rowHeight, "Summary", [this]( Button *b ) { ShowCentrePanel( &summary ); } )
		//
	{
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &maxMatches );
		stack.add( &rankingFileName );
		stack.add( &autoLoad );

		rankingFileName.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &dispButton );
		buttons.add( &summaryButton );
		buttons.add( &browserButton );
		buttons.add( &precisionRecallBrowserButton );
		buttons.add( &precisionDistanceBrowserButton );

		add( &disp );
		add( &summary );
		add( &browser );
		add( &precisionRecallBrowser );
		add( &precisionDistanceBrowser );

		browser.AddPropertyChangedEventHandler( this );
		precisionRecallBrowser.AddPropertyChangedEventHandler( this );
		precisionDistanceBrowser.AddPropertyChangedEventHandler( this );

		auto getRankings = [this]() { return &rankings; };
		browser.getListItems = getRankings;
		precisionRecallBrowser.getListItems = getRankings;
		precisionDistanceBrowser.getListItems = getRankings;

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		ShowCentrePanel( &disp );
	}

	virtual ~RankBagOfWords() {}

#pragma region Data flow connections to antecedents
	function<decltype(featureVectors) ()>     getFeatureVectors;
	function<decltype(dbSeqs) ()>             getDbSeqs;
	function<decltype(dbIndex) ()>            getDbIndex;
	function<decltype(testSetIndex) ()>       getTestSetIndex;
	function<decltype(trainingSetIndex) ()>   getTrainingSetIndex;
	function<decltype(featurePostingList) ()> getPostingList;
	function<decltype(classPostingList) ()>   getClassPostingList;
#pragma endregion

	void GainFocus() override {
		Page::GainFocus();
		featureVectors = getFeatureVectors();
		dbSeqs = getDbSeqs();
		dbIndex = getDbIndex();
		testSetIndex = getTestSetIndex();
		trainingSetIndex = getTrainingSetIndex();
		featurePostingList = getPostingList();
		classPostingList = getClassPostingList();

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

	static void RunImpl_( RankBagOfWords * self ) {
		self->RunImpl();
	}

	void RunImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

		Page::runTime = -omp_get_wtime();
		Rank();
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

		stats.Update( rankings, 100 );
		meanAveragePrecision = stats.meanAveragePrecision;
		UpdateGui();

		mute.unlock();
	}

	void UpdateGui() {
		UPDATE_GUI(
			summary.Refresh(); \
			browser.Refresh(); \
			precisionRecallBrowser.Refresh(); \
			precisionDistanceBrowser.Refresh();  \
			UpdateReadiness(); \
			NotifyRunComplete();
		);
	}

	void ResetBrowsers() {
		summary.Clear();
		browser.Clear();
		precisionRecallBrowser.Clear();
		precisionDistanceBrowser.Clear();
	}

	void UpdateBrowsers() {
		summary.GainFocus();
		browser.GainFocus();
		precisionRecallBrowser.GainFocus();
		precisionDistanceBrowser.GainFocus();
	}

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
		const vector<SparseFeatureVector> & database,
		const vector<size_t> & queryIndices,
		int maxMatches,
		vector<::Ranking> & rankings
	) {
		rankings.clear();

		KnnVector<const FastaSequence *, double> knn( maxMatches, -1 );
		vector<double> precision;
		vector<double> recall;

		for ( auto i : queryIndices ) {
			auto & sig = database[i];
			rankings.emplace_back( sig.Sequence(), knn, precision, recall );
			auto & ranking = rankings.back();

			auto &precision = ranking.precision;
			precision.clear();
			precision.reserve( maxMatches );

			auto &recall = ranking.recall;
			recall.clear();
			recall.reserve( maxMatches );
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
	static void RankFeatureVectors(
		const SparseFeatureVector & querySig,
		const vector<SparseFeatureVector> & signatures,
		const vector<vector<size_t>> & sigPostingList,
		const unordered_map<size_t, vector<size_t>> & classPostingList,
		BitSet & processed,
		::Ranking & ranking
	) {
		processed.Clear();
		ranking.knn.clear();

		auto eol = sigPostingList.end();

		for ( auto c : querySig ) {
			for ( size_t d : sigPostingList[c.key] ) {
				if ( !processed.Contains( d ) ) {
					processed.Insert( d );
					auto &dbSig = signatures[d];
					double distance = -querySig.Dot( dbSig );

					if ( ranking.knn.canPush( distance ) ) {
						ranking.knn.push( dbSig.Sequence(), distance );
					}
				}
			}
		}

		ranking.knn.sort();

		SigRank::ComputePrecisionRecall( classPostingList, processed, ranking );
	}

	void Rank() {
		const int maxMatches = MaxMatches();
		const size_t Q = (*testSetIndex).size();

		Setup( *featureVectors, *testSetIndex, maxMatches, rankings );

		int done = 0;

#if USE_OMP
#pragma omp parallel
#endif
		{
			BitSet processed( featureVectors->size() );

#if USE_OMP
#pragma omp for
#endif
			for ( size_t q = 0; q < Q; q++ ) {
				auto &querySig = (*featureVectors)[(*testSetIndex)[q]];
				RankFeatureVectors( querySig, *featureVectors, *featurePostingList, *classPostingList, processed, rankings[q] );
#if USE_OMP
#pragma omp critical
#endif
				{
					done++;
					int thisPercent = done * 100 / Q, lastPercent = (done - 1) * 100 / Q;

					if ( thisPercent != lastPercent ) {
						UpdateGuiProgress( thisPercent );
					}
				}
			}

		}

		UpdateBrowsers();
	}

	void UpdateGuiProgress( int thisPercent ) {
		UPDATE_GUI( browser.Refresh(); precisionRecallBrowser.Refresh(); disp( "Ranking training set: %d%%\n", thisPercent ); )
	}

	void Save() {
		auto fileName = String::Trim( RankingFileName() );
		SigRank::Save( rankings, fileName );
	}

	static void LoadImpl_( RankBagOfWords * self ) {
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

		stats.Update( rankings, PrecisionSteps() );
		meanAveragePrecision = stats.meanAveragePrecision;
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
		SigRank::Load( fileName, *dbIndex, classPostingList, errorMsg, setMaxMatches, rankings );

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

		SetReady( trainingSetSize > 0 && testSetSize > 0 && dbSeqSize <= (trainingSetSize + testSetSize) && fileName.length() > 0 );
		SetOk( Ready() && rankings.size() == testSetSize );
	}

	void PropertyChanged( void *sender, const string &propertyName ) override {
		// Handle property changed events...
		long Focus = -1;

		if ( propertyName == NAMEOF( Focus ) ) {
			auto source = (RankingBrowserBase *) (sender);
			Focus = source->Focus();

			if ( Focus >= 0 ) {
				browser.SetFocus( Focus );
				precisionRecallBrowser.SetFocus( Focus );
				precisionDistanceBrowser.SetFocus( Focus );
			}
		}

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
		addPar( precisionSteps );
		addVar( runTime, "%0.17g" );
		addVar( loadTime, "%0.17g" );
		addVar( saveTime, "%0.17g" );
		addVar( meanAveragePrecision, "%0.17g" );
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
			setPar( precisionSteps, double );
			setVar( runTime, double );
			setVar( loadTime, double );
			setVar( saveTime, double );
			setVar( meanAveragePrecision, double );
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

private:
	Arg<int, Fl_Int_Input> precisionSteps;

public:
	int PrecisionSteps() {
		return precisionSteps.Value();
	}

	void SetPrecisionSteps( int value ) {
		if ( value == PrecisionSteps() ) return;

		this->precisionSteps.SetValue( value );
		NotifyPropertyChanged( NAMEOF( PrecisionSteps ) );
	}

public:
	void Reset() {
		stats.Clear();
		ResetBrowsers();

		featureVectors = 0;
		dbSeqs = 0;
		dbIndex = 0;
		testSetIndex = 0;
		trainingSetIndex = 0;
		featurePostingList = 0;
		classPostingList = 0;
	}
};
