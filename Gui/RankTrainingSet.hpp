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
#include <SimilarityMatrix.hpp>
#include <SimpleKmer.hpp>
#include <Centroid.hpp>
#include <Button.hpp>
#include "Helpers.hpp"
#include "RankingBrowserBase.hpp"
#include "RankingBrowser.hpp"
#include "RankingDetailBrowser.hpp"
#include "PrecisionRecallBrowser.hpp"
#include "PrecisionDistanceBrowser.hpp"
#include "SigRank.hpp"
#include "PrecisionRecallSummary.hpp"

using namespace std;
using namespace QutBio;

/**
<summary>
	DbLoader gathers the parameters for an interactive version of DnaClust,
	with some additional analytic displays to render the software a little more
	usable.
</summary>
 */
class RankTrainingSet :
	public virtual Page,
	public virtual LBFL::IPropertyChangedEventHandler
	//
{
	const int rows = 8;
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

	RankingDetailBrowser rankingDetailBrowser;
	Button rankingDetailBrowserButton;

	mutex mute;
#pragma endregion

#pragma region Data flow connections
	const vector<SparseSignature> * signatures;
	const vector<FastaSequence *> *dbSeqs;
	const LookupTable_<size_t, const FastaSequence> * dbIndex;
	const vector<size_t> * testSetIndex;
	const vector<size_t> * trainingSetIndex;
	const vector<vector<size_t>> *featurePostingList;
	const unordered_map<size_t, vector<size_t>> *classPostingList;

	AlignInfo alignInfo;
#pragma endregion

#pragma region Results
	vector<::Ranking> rankings;
	PrecisionRecallStats stats;
	double meanAveragePrecision;
#pragma endregion

public:
	RankTrainingSet(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Rank Training Set" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, w(), rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),

		rankingFileName(
			"rankingFileName", "*.rankings.csv", "Ranking File", Requirement::Required,
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

		interleave(
			"interleave", 0, "Interleave?", Requirement::Required,
			"Check this box to interleave output to disk with result calculation. \n"
			"This allows calculations to proceed while one thread at a time is \n"
			"waiting on output, which improves overall throughput.",
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

		rankingDetailBrowser(),
		rankingDetailBrowserButton( 0, 0, 150, rowHeight, "Details", [this]( Button *b ) { ShowCentrePanel( &rankingDetailBrowser ); } ),

		summary( stats ),
		summaryButton( 0, 0, 150, rowHeight, "Summary", [this]( Button *b ) { ShowCentrePanel( &summary ); } )
		//
	{
		rankingDetailBrowser.Init( &alignInfo );

		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &maxMatches );
		stack.add( &rankingFileName );
		stack.add( &autoLoad );
		stack.add( &interleave );
		stack.add( &precisionSteps );

		rankingFileName.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &dispButton );
		buttons.add( &summaryButton );
		buttons.add( &browserButton );
		buttons.add( &precisionRecallBrowserButton );
		buttons.add( &precisionDistanceBrowserButton );
		buttons.add( &rankingDetailBrowserButton );

		add( &disp );
		add( &summary );
		add( &browser );
		add( &precisionRecallBrowser );
		add( &precisionDistanceBrowser );
		add( &rankingDetailBrowser );

		browser.AddPropertyChangedEventHandler( this );
		precisionRecallBrowser.AddPropertyChangedEventHandler( this );
		precisionDistanceBrowser.AddPropertyChangedEventHandler( this );
		rankingDetailBrowser.AddPropertyChangedEventHandler( this );

		auto getRankings = [this]() { return &rankings; };
		browser.getListItems = getRankings;
		precisionRecallBrowser.getListItems = getRankings;
		precisionDistanceBrowser.getListItems = getRankings;
		rankingDetailBrowser.getListItems = getRankings;

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		ShowCentrePanel( &disp );
	}

	virtual ~RankTrainingSet() {}

#pragma region Data flow connections to antecedents
	function<const vector<SparseSignature> * ()>	getSignatures;
	function<const vector<FastaSequence *> *()>		getDbSeqs;
	function<const LookupTable_<size_t, const FastaSequence> * ()> getDbIndex;
	function<const vector<size_t> * ()>				getTestSetIndex;
	function<const vector<size_t> * ()>				getTrainingSetIndex;
	function<const vector<vector<size_t>> *()>      getPostingList;
	function<const unordered_map<size_t, vector<size_t>> *()> getClassPostingList;

	function<decltype(alignInfo.kmerLength) ()> getKmerLength;
	function<decltype(alignInfo.matrix) ()> getMatrix;
	function<decltype(alignInfo.sigIndex) ()> getSigIndex;
	function<decltype(alignInfo.threshold) ()> getThreshold;
	function<decltype(alignInfo.vocab) ()> getVocab;
#pragma endregion

	void GainFocus() override {
		Page::GainFocus();
		signatures = getSignatures();
		dbSeqs = getDbSeqs();
		dbIndex = getDbIndex();
		testSetIndex = getTestSetIndex();
		trainingSetIndex = getTrainingSetIndex();
		featurePostingList = getPostingList();
		classPostingList = getClassPostingList();

		alignInfo.sigIndex = getSigIndex();
		alignInfo.kmerLength = getKmerLength();
		alignInfo.vocab = getVocab();
		alignInfo.matrix = getMatrix();
		alignInfo.sigIndex = getSigIndex();
		alignInfo.threshold = getThreshold();

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

	static void RunImpl_( RankTrainingSet * self ) {
		self->RunImpl();
	}

	void RunImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

		if ( Interleave() ) {
			Page::runTime = -omp_get_wtime();

			try {
				Rank( RankingFileName() );
			}
			catch ( Exception & ex ) {
				ShowProgress( "Error during interleaved ranking: %s\n", ex.what() );
			}

			Page::runTime += omp_get_wtime();
			Page::saveTime = 0;
		}
		else {
			Page::runTime = -omp_get_wtime();
			Rank();
			Page::runTime += omp_get_wtime();

			Page::saveTime = -omp_get_wtime();
			try {
				Save();
			}
			catch ( Exception & ex ) {
				ShowProgress( "Unable to save rankings: %s\n", ex.what() );
			}
			Page::saveTime += omp_get_wtime();
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
			rankingDetailBrowser.Refresh();  \
			UpdateReadiness(); \
			NotifyRunComplete();
		);
	}

	void ResetBrowsers() {
		summary.Clear();
		browser.Clear();
		precisionRecallBrowser.Clear();
		precisionDistanceBrowser.Clear();
		rankingDetailBrowser.Clear();
	}

	void UpdateBrowsers() {
		summary.GainFocus();
		browser.GainFocus();
		precisionRecallBrowser.GainFocus();
		precisionDistanceBrowser.GainFocus();
		rankingDetailBrowser.GainFocus();
	}

	void Rank() {
		const int maxMatches = MaxMatches();
		const size_t Q = (*testSetIndex).size();

		SigRank::Setup( *signatures, *testSetIndex, maxMatches, rankings );

#if WANT_FEEDBACK
		int done = 0;
#endif

#if USE_OMP
#pragma omp parallel
#endif
		{
			BitSet processed( signatures->size() );

#if USE_OMP
#pragma omp for schedule(dynamic,1)
#endif
			for ( size_t q = 0; q < Q; q++ ) {
				auto &querySig = (*signatures)[(*testSetIndex)[q]];
				SigRank::RankSignatures( querySig, *signatures, *featurePostingList, *classPostingList, processed, rankings[q] );

				if ( querySig.Sequence()->ClassIndex() >= 0 ) {
					SigRank::ComputePrecisionRecall( *classPostingList, processed, rankings[q] );
				}

#if false || WANT_FEEDBACK
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
#endif
			}

		}

		UpdateBrowsers();
	}

	/**
	 *	Rank training set with Interleave. This DOES NOT DO PRECISION/RECALL.
	 *	You'll have to use a stand-alone post-processing step to derive stats
	 *	from the saved rankings if you use this step.
	 */

	void Rank( const string & fileName ) {
		ofstream rankingFile(fileName);
		CsvWriter rankingWriter(rankingFile);

		const int maxMatches = MaxMatches();
		const size_t Q = (*testSetIndex).size();

		rankingWriter << "rankings" << Q << '\n';

		SigRank::Setup( *signatures, *testSetIndex, maxMatches, rankings );

#if WANT_FEEDBACK
		int done = 0;
#endif

#if USE_OMP
#pragma omp parallel
#endif
		{
			BitSet processed( signatures->size() );

#if USE_OMP
#pragma omp for schedule(dynamic,1)
#endif
			for ( size_t q = 0; q < Q; q++ ) {
				auto &querySig = (*signatures)[(*testSetIndex)[q]];
				SigRank::RankSignatures( querySig, *signatures, *featurePostingList, *classPostingList, processed, rankings[q] );

				#pragma omp critical
				{
					SigRank::SaveRanking( rankingWriter, rankings[q] );
				}

#if false || WANT_FEEDBACK
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
#endif
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

	static void LoadImpl_( RankTrainingSet * self ) {
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
			);
		};

		auto setMaxMatches = [this]( size_t maxMatches ) {
			UPDATE_GUI( SetMaxMatches( maxMatches ); disp( "%zu rankings loaded.\n", rankings.size() ); );
		};

		ShowProgress( "Attempting to load rankings from %s:\n", fileName.c_str() );
		SigRank::Load( fileName, *dbIndex, classPostingList, errorMsg, setMaxMatches, rankings );

		UpdateReadiness();
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
				rankingDetailBrowser.SetFocus( Focus );
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
		addPar( interleave );
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
			setPar( interleave, int );
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
	TypedArg<int, Fl_Check_Button> interleave;

public:
	int Interleave() {
		return interleave.Value();
	}

	void SetInterleave( const int val ) {
		if ( val == Interleave() ) return;

		this->interleave.SetValue( val );
		NotifyPropertyChanged( NAMEOF( Interleave ) );
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

		signatures = 0;
		dbSeqs = 0;
		dbIndex = 0;
		testSetIndex = 0;
		trainingSetIndex = 0;
		featurePostingList = 0;
		classPostingList = 0;

		alignInfo.sigIndex = 0;
		alignInfo.kmerLength = 0;
		alignInfo.vocab = 0;
		alignInfo.matrix = 0;
		alignInfo.sigIndex = 0;
		alignInfo.threshold = 0;
	}
};
