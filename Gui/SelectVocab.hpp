#pragma once

#include <FL/Enumerations.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl.H>
#include <FL/Fl_Chart.H>
#include <FL/Fl_Tabs.H>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <thread>

#include <FlowLayout.hpp>
#include <Graph.hpp>
#include <Random.hpp>
#include <Histogram.hpp>
#include <IEventHandler.hpp>
#include <IPropertyChangedEventHandler.hpp>
#include <MarginLayout.hpp>
#include <Rectangle.hpp>
#include <SimilarityMatrix.hpp>
#include <SimpleKmer.hpp>
#include "IncludeAll.hpp"
#include <Centroid.hpp>
#include <Button.hpp>
#include <kNearestNeighbours.hpp>
#include "Helpers.hpp"

using namespace LBFL;
using namespace QutBio;
using namespace LBGraph;
using namespace std;

/**
<summary>
	DbLoader gathers the parameters for an interactive version of DnaClust,
	with some additional analytic displays to render the software a little more
	usable.
</summary>
 */
class SelectVocab : public virtual Page,
	public virtual IPropertyChangedEventHandler
	//
{
	const int rows = 7;
	int rowHeight = 0;

#pragma region Seed
private:
	Arg<int, Fl_Int_Input> seed;

public:
	int Seed() {
		return seed.Value();
	}

	void SetSeed( const int val ) {
		if ( val == Seed() ) return;

		this->seed.SetValue( val );
	}
#pragma endregion Seed

#pragma region VocabSize
private:
	Arg<uint, Fl_Int_Input> vocabSize;

public:
	uint VocabSize() {
		return vocabSize.Value();
	}

	void SetVocabSize( const uint val ) {
		if ( val == VocabSize() ) return;

		this->vocabSize.SetValue( val );
	}
#pragma endregion VocabSize

#pragma region SelectionMode
private:
	TypedArg<int, Fl_Choice> selectionMode;

public:
	int SelectionMode() {
		return selectionMode.Value();
	}

	void SetSelectionMode( const int val ) {
		if ( val == SelectionMode() ) return;

		this->selectionMode.SetValue( val );
	}
#pragma endregion SelectionMode

#pragma region VocabFileName
private:
	TypedArg<string, OutFileChooser> vocabFileName;

public:
	string VocabFileName() {
		return vocabFileName.Value();
	}

	void SetOutFileName( const string & val ) {
		if ( val == VocabFileName() ) return;

		this->vocabFileName.SetValue( val );
	}
#pragma endregion VocabFileName

#pragma region
private:
	GridLayout stack;
	Fl_Box titleBox;
	FlowLayout buttons;

	TextDisplay note;
	Button noteButton;

	TextDisplay disp;
	Button dispButton;

	UniformIntRandom<size_t> rand;

	mutex mute;

	/// Input: list of prototypes from which the vocabulary will be selected.
	vector<Centroid> * prototypes;

	/// Input: word length for kmer parsing.
	size_t kmerLength;

	/// Input: Sequence index.
	const LookupTable_<size_t, const FastaSequence> * seqIndex;

	/// Output: Selected k-mers (initially obtained from prototypes, but after 
	/// load may be disconnected. 
	vector<Centroid> vocab;

#pragma endregion  // private data members

public:
	SelectVocab(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
#pragma region Ctor init
		Page( left, top, width, height, "Select Vocabulary" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, w(), rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),
		seed(
			"seed", time( nullptr ), "RNG Seed", Requirement::Required,
			"Cardinal number; required. Seed for random number generator used\n"
			"to sample k-mer population to compute empirical distances, and \n"
			"available for downstream components if necessary.",
			labelWidth, 0, width, rowHeight ),

		vocabSize(
			"vocabSize", 10000, "Vocabulary Size",
			Requirement::Required,
			"Cardinal number; required. Number of kmers that will be added to the \n"
			"vocabulary.",
			labelWidth, 0, width, rowHeight ),

		vocabFileName(
			"vocabFileName", "*.vocab.csv", "Vocabulary File", Requirement::Required,
			"The path to a file which will be overwritten with k-mers representing \n"
			"vocabulary terms, in CSV format.",
			labelWidth, 0, width, rowHeight ),

		autoLoad( "autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load the partition file if it exists.",
			labelWidth, 0, width, rowHeight ),

		buttons( 0, 0, w(), rowHeight ),
		note( 0, 0, w(), 100 ),
		disp( 0, 0, w(), 100 ),
		rand( (unsigned long) time( 0 ) ),

		selectionMode( "selectionMode", 0, "Selection mode", Requirement::Required,
			"Vocabulary selection mode:\n"
			"First N -- select the N first prototypes identified by greedy incremental clustering;\n"
			"MinHash -- select the N items with lowest hash codes;\n"
			"Uniform -- select N items uniformly at random.",
			labelWidth, 0, width, rowHeight
		),

		noteButton( 0, 0, 100, rowHeight, "Notes", [this]( Button *b ) { ShowCentrePanel( &note ); } ),
		dispButton( 0, 0, 100, rowHeight, "Log", [this]( Button *b ) { ShowCentrePanel( &disp ); } )
		//
#pragma endregion
	{
#pragma region Ctor code
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		selectionMode.inputField.add( "First N" );
		selectionMode.inputField.add( "MinHash" );
		selectionMode.inputField.add( "Uniform" );
		SetSelectionMode( 0 );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &seed );
		stack.add( &vocabSize );
		stack.add( &selectionMode );
		stack.add( &vocabFileName );
		stack.add( &autoLoad );

		seed.AddPropertyChangedEventHandler( this );
		vocabSize.AddPropertyChangedEventHandler( this );
		vocabFileName.AddPropertyChangedEventHandler( this );
		selectionMode.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &noteButton );
		buttons.add( &dispButton );

		add( &note );
		add( &disp );

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		note.textfont( FL_COURIER );
		note.textsize( 16 );
		note( "Select prototypes to form a compressed vocabulary.\n\n"
			"<... more will be added here ...>\n\n"
			"Meantime, switch to 'Log' and use 'Update' to run this operation."
		);

		ShowCentrePanel( &disp );
#pragma endregion
	}

	virtual ~SelectVocab() {}

#pragma region Data flow connections to antecedents
	function<vector<Centroid> *()> getPrototypes;
	function<size_t()> getKmerLength;
	function<const LookupTable_<size_t, const FastaSequence> *()> getSeqIndex;
#pragma endregion Data flow connections to antecedents

	void GainFocus() override {
		Page::GainFocus();
		prototypes = getPrototypes();
		kmerLength = getKmerLength();
		seqIndex = getSeqIndex();
		UpdateReadiness();
	}

	YesNoCancel CheckLoad() {
		auto fileName = String::Trim( VocabFileName() );

		if ( fileName.length() > 0 && File::Exists( fileName ) ) {
			if ( AutoLoad() ) {
				return YesNoCancel::Yes;
			}

			int answer = fl_choice(
				"The selected file already exists.\n"
				"Would you like to load that vocabulary?",
				"Cancel operation.",
				"Yes, load the vocabulary.",
				"No, I want to redo the vocabulary."
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
			thread t( LoadImpl_, this );
			t.detach();
		}
		else if ( choice == YesNoCancel::No ) {
			thread t( RunImpl_, this );
			t.detach();
		}
		else {
			UpdateReadiness();
			disp( "Action cancelled.\n" );
		}
	}

	static void RunImpl_( SelectVocab * self ) {
		self->RunImpl();
	}

	void RunImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

		Page::runTime = -omp_get_wtime();
		Select();
		Page::runTime += omp_get_wtime();

		try {
			Page::saveTime = -omp_get_wtime();
			Save();
			Page::saveTime += omp_get_wtime();
		}
		catch ( Exception & ex ) {
			ShowProgress( "Unable to save vocabulary: %s\n", ex.what() );
		}

		UPDATE_GUI( UpdateReadiness();  NotifyRunComplete(); );
		mute.unlock();
	}

	void Select() {
		auto seed = Seed();
		auto vocabSize = VocabSize();
		auto selectionMode = SelectionMode();
		auto wanted = std::min( prototypes->size(), (size_t) vocabSize );

		vocab.clear();

		if ( selectionMode == 0 ) {
			for ( size_t i = 0; i < wanted; i++ ) {
				vocab.push_back( prototypes->at( i ) );
			}
		}
		else if ( selectionMode == 2 ) {
			UniformRealRandom rand( seed );
			Selector sel( rand, wanted, prototypes->size() );

			for ( auto & centroid : *prototypes ) {
				if ( sel.SelectThis() ) {
					vocab.push_back( centroid );
				}
			}
		}
		else /* Minhash */ {
			KnnVector<const Centroid *, size_t> knn( wanted, numeric_limits<size_t>::min() );
			Substring::Hash hash;

			for ( auto &centroid : *prototypes ) {
				auto code = _Hash_impl::hash( (void *) centroid.centroid->Bytes(), kmerLength );

				if ( knn.canPush( code ) ) {
					knn.push( &centroid, code );
				}
			}

			for ( auto &p : knn ) {
				vocab.push_back( *(p.second) );
			}
		}

		ShowProgress( "%zu terms added to vocabulary.\n", vocab.size() );
	}

	void Save() {
		auto fileName = String::Trim( VocabFileName() );

		ofstream f( fileName );
		CsvWriter w( f );

		w << "vocabulary" << vocab.size() << '\n';
		for ( auto & term : vocab ) {
			w << term << '\n';
		}
	}

	static void LoadImpl_( SelectVocab * self ) {
		self->LoadImpl();
	}

	void LoadImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		Helpers::UpdateGui( [this]() {
			SetOk( false );
			disp.Clear();
		} );

		try {
			Page::loadTime = -omp_get_wtime();
			Load();
			Page::loadTime += omp_get_wtime();
		}
		catch ( Exception & ex ) {
			ShowProgress( "Error loading vocabulary: %s\n", ex.what() );
		}

		UPDATE_GUI( UpdateReadiness();  NotifyRunComplete(); );
		mute.unlock();
	}

	void Load() {
		auto fileName = String::Trim( VocabFileName() );
		vocab.clear();

		ShowProgress( "Attempting to load vocabulary terms from %s:\n", fileName.c_str() );
		ifstream f( fileName );
		CsvReader r( f );
		string s;
		string kmerDef;
		size_t n;

		r >> s >> n;
		vocab.clear();
		vocab.resize( n );

		for ( size_t i = 0; i < n && !r.IsEOF(); i++ ) {
			Helpers::GetCentroid( r, kmerDef, kmerLength, *seqIndex, vocab[i] );
		}

		Helpers::UpdateGui( [this]() {
			SetVocabSize( vocab.size() );
			disp( "%zu vocabulary terms loaded.\n", vocab.size() );
		} );
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
		auto vocabSize = VocabSize();
		auto fileName = String::Trim( VocabFileName() );

		SetReady( vocabSize > 0 && prototypes != 0 && fileName.length() > 0 );
		SetOk( Ready() && vocab.size() == std::min( prototypes->size(), (size_t) VocabSize() ) );
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
		addPar( seed );
		addPar( vocabSize );
		addPar( selectionMode );
		addPar( vocabFileName );
		addPar( autoLoad );
		addVar( runTime, "%0.17g" );
		addVar( loadTime, "%0.17g" );
		addVar( saveTime, "%0.17g" );
#undef addPar
#undef addVar
	}

	/// Set parameters which have been de-serialised when program starts.
	/// @param parms A Param list containing de-serialised parameters.
	void SetParams( const set<Param> & parms ) override {
#define setPar(x,t) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value<t>() )
#define setVar(x,t) if ( p.ParamName() == STRING(x) ) (x) = p.Value<t>()
		for ( auto & p : parms ) {
			if ( p.ComponentName() != Name() ) continue;
			setPar( seed, int );
			setPar( vocabSize, uint );
			setPar( selectionMode, int );
			setPar( vocabFileName, string );
			setPar( autoLoad, int );
			setVar( runTime, double );
			setVar( loadTime, double );
			setVar( saveTime, double );
		}

		UpdateReadiness();
#undef setPar
#undef setVar
	}

#pragma region Vocab
public:
	const vector<Centroid> * Vocab() {
		return &vocab;
	}
#pragma endregion

#pragma region AutoLoad Property
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
#pragma endregion

	void Reset() override {
		vocab.clear();
		prototypes = 0;
		kmerLength = 0;
		seqIndex = 0;
	}
};
