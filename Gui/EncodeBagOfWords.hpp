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
#include "Helpers.hpp"
#include <SparseFeatureVector.hpp>
#include "FeatureVectorBrowser.hpp"

#include <ComboBox.hpp>

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
class EncodeBagOfWords : public virtual Page,
	public virtual IPropertyChangedEventHandler
	//
{
	const int rows = 6;
	int rowHeight = 0;

	GridLayout stack;
	Fl_Box titleBox;
	FlowLayout buttons;

	TextDisplay disp;
	Button dispButton;

	FeatureVectorBrowser featureVectorBrowser;
	Button sigButton;

	unordered_map<uint, Registry<Substring, Substring::Hash>> kmerRegistry;

	UniformIntRandom<size_t> rand;

	mutex mute;

	/// Input: Sequence index.
	const LookupTable_<size_t, const FastaSequence> * seqIndex = 0;

	/// Input: list of sequences to be encoded.
	const vector<pFastaSequence> *  dbSeqs = 0;

	/// Input: Locations of the training set in the sequence dataset.
	const vector<size_t> * trainingSetIndex = 0;

	/// Output: list of sequence signatures.
	vector<SparseFeatureVector> featureVectors;

	/// Map feature number to list of sequences containing feature. 
	vector<vector<size_t>> featurePostingList;

	/// Internal: feature selection Ok?
	bool featuresOk = false;

	size_t vocabSize;

#pragma endregion 

public:
	EncodeBagOfWords(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Encode Bag Of Words" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, w(), rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),

		kmerLength(
			"kmerLength", 3, "K-mer Length", Requirement::Required,
			"The word size (k) for k-mer tiling.",
			labelWidth, 0, width, rowHeight ),

		outFileName(
			"outFileName", "", "Feature Vector File", Requirement::Required,
			"The path to a file which will be overwritten with a list of sparse \n"
			"feature vectors. A feature vector is a list of (feature,weight) pairs.\n"
			"Feature is a natural number obtained from the sequence database, so \n"
			"if the database changes you will need to rebuild the database feature \n"
			"vectors.",
			labelWidth, 0, width, rowHeight ),

		autoLoad(
			"autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load from the feature vector file, \n"
			"if it exists.",
			labelWidth, 0, width, rowHeight ),

		buttons( 0, 0, w(), rowHeight ),
		disp( 0, 0, w(), 100 ),
		rand( (unsigned long) time( 0 ) ),
		dispButton( 0, 0, 100, rowHeight, "Log", [this]( Button *b ) { ShowCentrePanel( &disp ); } ),
		sigButton( 0, 0, 100, rowHeight, "Features", [this]( Button *b ) { ShowCentrePanel( &featureVectorBrowser ); } ),
		featureVectorBrowser( featureVectors, vocabSize )
		//
	{
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &kmerLength );
		stack.add( &outFileName );
		stack.add( &autoLoad );

		outFileName.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &dispButton );
		buttons.add( &sigButton );

		add( &disp );
		add( &featureVectorBrowser );

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		ShowCentrePanel( &disp );
	}

	virtual ~EncodeBagOfWords() {}

#pragma region Data flow connections to antecedents
	function<const vector<FastaSequence *> *()> getDatabase;
	function<const LookupTable_<size_t, const FastaSequence> *()> getSeqIndex;
	function<const vector<size_t> * ()> getTrainingSetIndex;
#pragma endregion

	void GainFocus() override {
		// cerr << "GainFocus\n";

		Page::GainFocus();

		dbSeqs = getDatabase();
		seqIndex = getSeqIndex();
		trainingSetIndex = getTrainingSetIndex();

		featureVectorBrowser.SetDbSeqs( dbSeqs );
		featureVectorBrowser.GainFocus();

		UpdateReadiness();
	}

	YesNoCancel CheckLoad() {
		auto fileName = String::Trim( OutFileName() );

		if ( fileName.length() > 0 && File::Exists( fileName ) ) {
			if ( AutoLoad() ) {
				return YesNoCancel::Yes;
			}

			int answer = fl_choice(
				"The selected file already exists.\n"
				"Would you like to load those signatures?",
				"Cancel operation.",
				"Yes, load the signatures.",
				"No, I want to redo the signatures."
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
			ostringstream str;
			if ( KmerLength() == 0 ) str << "K-mer length must be greater than zero.\n";
			if ( !dbSeqs ) str << "Sequence dataset not initialised!\n";
			auto s = str.str();

			if ( s.length() > 0 ) {
				fl_message( "%s\nRun abandoned.\n", s.c_str() );
				UpdateReadiness();
				disp( "Run abandoned.\n" );
			}
			else {
				thread t( RunImpl_, this );
				t.detach();
			}
		}
		else {
			UpdateReadiness();
			disp( "Action cancelled.\n" );
		}
	}

	void IndexDatabase() {
		if ( !featuresOk ) return;

		SparseFeatureVector::CreatePostingList( featureVectors, *trainingSetIndex, featurePostingList );
	}

	static void RunImpl_( EncodeBagOfWords * self ) {
		self->RunImpl();
	}

	void RunImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		auto fileName = String::Trim( OutFileName() );

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

		Page::runTime = -omp_get_wtime();
		GenerateSignatures();
		IndexDatabase();
		Page::runTime += omp_get_wtime();

		try {
			Page::saveTime = -omp_get_wtime();
			Save( fileName );
			Page::saveTime += omp_get_wtime();
		}
		catch ( Exception & ex ) {
			ShowProgress( "Unable to save vocabulary: %s\n", ex.what() );
		}

		UPDATE_GUI( UpdateReadiness(); NotifyRunComplete(); );

		mute.unlock();
	}

	void GenerateFeatureVector(
		const FastaSequence & seq,
		uint kmerLength,
		SparseFeatureVector & featureVector
	) {
		auto &reg = kmerRegistry[kmerLength];

		featureVector.SetSequence( &seq );

		size_t K = seq.KmerCount( kmerLength );
		auto bytes = seq.Sequence().data();

		for ( size_t i = 0; i < K; i++ ) {
			Substring kmer( bytes, i, kmerLength );
			size_t kmerId = 0;
			kmerId = reg( kmer );
			featureVector.Add( kmerId );
		}

		featureVector.Sort();
	}

	void GenerateSignatures() {
		featuresOk = false;
		const size_t N = dbSeqs->size();
		size_t done = 0;

		vocabSize = 0;

		featureVectors.clear();
		featureVectors.resize( N );

		for ( size_t i = 0; i < N; i++ ) {
			pFastaSequence seq = dbSeqs->at( i );
			GenerateFeatureVector( *seq, KmerLength(), featureVectors[i] );

			auto maxKey = featureVectors[i].MaxKey();

			if ( maxKey + 1 > vocabSize ) {
				vocabSize = maxKey + 1;
			}

			// Periodically show progress
			{
				done++;
				int thisPercent = done * 100 / N, lastPercent = (done - 1) * 100 / N;

				if ( thisPercent != lastPercent ) {
					UPDATE_GUI( disp( "Encoding sequences: %d%%\n", thisPercent ); )
				}
			}
		}

		featuresOk = true;
		featureVectorBrowser.GainFocus();
	}

	const string FileTag = "signatures";

	void Save( const string & fileName ) {
		ofstream f( fileName );
		CsvWriter w( f );

		w << FileTag << featureVectors.size() << '\n';

		for ( auto &c : featureVectors ) {
			w << c.Sequence()->IdStr() << c;
		}
	}

	static void LoadImpl_( EncodeBagOfWords * self ) {
		self->LoadImpl();
	}

	void LoadImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		auto fileName = String::Trim( OutFileName() );

		{ UPDATE_GUI( SetOk( false ); disp.Clear();) }

		Page::loadTime = -omp_get_wtime();
		Load( fileName );
		IndexDatabase();
		Page::loadTime += omp_get_wtime();

		{ UPDATE_GUI( UpdateReadiness(); NotifyRunComplete(); ) }
		mute.unlock();
	}

	void Load( const string & fileName ) {
		featuresOk = false;

		ifstream f( fileName );
		CsvReader r( f );

		string scratch;
		size_t N;
		r >> scratch >> N;

		if ( scratch != FileTag ) {
			throw Exception( "Expected tag " + FileTag + " not found in file " + fileName, FileAndLine );
		}

		vocabSize = 0;

		featureVectors.clear();
		featureVectors.resize( N );

		for ( size_t i = 0; (!r.IsEOF()) && (i < N); i++ ) {
			auto &sig = featureVectors[i];
			r >> scratch >> sig;
			sig.SetSequence( seqIndex->at( FastaSequence::Register( scratch ) ) );

			auto maxKey = featureVectors[i].MaxKey();

			if ( maxKey + 1 > vocabSize ) {
				vocabSize = maxKey + 1;
			}
		}

		featuresOk = true;
		featureVectorBrowser.GainFocus();
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
		SetReady( String::Trim( OutFileName() ).length() > 0 );
		SetOk( Ready() && featuresOk );
	}

	void PropertyChanged( void *sender, const string &propertyName ) override {
		// Handle property changed events...

		//if ( (void *) &kmerLength == sender ) {
		//	for (auto & kReg : kmerRegistry ) {
		//		kReg.Clear();
		//	}
		//}

		UpdateReadiness();
		redraw();
	}

	/// Get parameters so they can be serialised when program finishes.
	/// @param parms A Param list to which the component should append its current parameters. 
	void GetParams( set<Param> &parms ) override {
#define addPar(x) parms.emplace( Name(), STRING(x), (x).Value() )
#define addVar(x,fmt) parms.emplace( Name(), STRING(x), String::Format(fmt,x) )
		addPar( kmerLength );
		addPar( outFileName );
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
			setPar( outFileName, string );
			setPar( kmerLength, uint );
			setPar( autoLoad, int );
			setVar( runTime, double );
			setVar( loadTime, double );
			setVar( saveTime, double );
		}

		UpdateReadiness();
#undef setPar
#undef setVar
	}


#pragma region Miscellaneous properties
	const decltype(featureVectors) * FeatureVectors() const { return &featureVectors; }

	const decltype(featurePostingList) *PostingList() const { return &featurePostingList; }
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

#pragma region AutoLoad Property
private:
	TypedArg<string, OutFileChooser> outFileName;

public:
	string OutFileName() {
		return outFileName.Value();
	}

	void SetOutFileName( const string & val ) {
		if ( val == OutFileName() ) return;

		this->outFileName.SetValue( val );
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
		// cerr << "Reset\n";

		featureVectorBrowser.Clear();
		featureVectors.clear();
		dbSeqs = 0;
		seqIndex = 0;
		trainingSetIndex = 0;
	}
};
