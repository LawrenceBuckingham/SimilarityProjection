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
#include <SparseSignature.hpp>
#include "SignatureBrowser.hpp"
#include "SigRank.hpp"

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
class EncodePartition : public virtual Page,
	public virtual IPropertyChangedEventHandler
	//
{
	const int rows = 5;
	int rowHeight = 0;

#pragma region DbSigFileName
private:
	TypedArg<string, OutFileChooser> dbSigFileName;

public:
	string DbSigFileName() {
		return dbSigFileName.Value();
	}

	void SetOutFileName( const string & val ) {
		if ( val == DbSigFileName() ) return;

		this->dbSigFileName.SetValue( val );
	}
#pragma endregion DbSigFileName

#pragma region Private data members
private:
	GridLayout stack;
	Fl_Box titleBox;
	FlowLayout buttons;

	TextDisplay disp;
	Button dispButton;

	SignatureBrowser sigBrowser;
	Button sigButton;

	UniformIntRandom<size_t> rand;

	mutex mute;

	/// Input: list of prototypes from which the vocabulary will be selected.
	const vector<Centroid> * vocab;

	/// Input: word length for kmer parsing.
	size_t kmerLength;

	/// Input: Sequence index.
	const LookupTable_<size_t, const FastaSequence> * seqIndex;

	pSimilarityMatrix matrix;

	/// Input: list of sequences to be encoded.
	const vector<pFastaSequence> *  dbSeqs;

	/// Input: Locations of the trainig set in the sequence dataset.
	const vector<size_t> * trainingSetIndex;

	/// Input: cluster radius for feature labelling.
	Distance threshold;

	/// Output: list of sequence signatures.
	vector<SparseSignature> signatures;

	/// Map feature number to list of sequences containing feature. 
	vector<vector<size_t>> featurePostingList;

	/// Map sequence unique id to signature.
	unordered_map<size_t, SparseSignature *> sigIndex;

	/// Internal: feature selection Ok?
	bool featuresOk = false;

	int classIndex;
	int nameIndex;

#pragma endregion 

public:
	EncodePartition(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Encode Partition" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, w(), rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),

		dbSigFileName(
			"dbSigFileName", "*.signatures.csv", "Signature File", Requirement::Required,
			"The path to a file which will be overwritten with sparse binary \n"
			"signatures derived from the vocabulary.",
			labelWidth, 0, width, rowHeight ),

		autoLoad(
			"autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load the partition file if it exists.",
			labelWidth, 0, width, rowHeight ),

		buttons( 0, 0, w(), rowHeight ),
		disp( 0, 0, w(), 100 ),
		rand( (unsigned long) time( 0 ) ),
		dispButton( 0, 0, 100, rowHeight, "Log", [this]( Button *b ) { ShowCentrePanel( &disp ); } ),
		sigButton( 0, 0, 100, rowHeight, "Signatures", [this]( Button *b ) { ShowCentrePanel( &sigBrowser ); } )
		//
	{
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &dbSigFileName );
		stack.add( &autoLoad );

		dbSigFileName.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &dispButton );
		buttons.add( &sigButton );

		add( &disp );
		add( &sigBrowser );

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		ShowCentrePanel( &disp );
	}

	virtual ~EncodePartition() {}

#pragma region Data flow connections to antecedents
	function<pSimilarityMatrix()> getMatrix;
	function<const vector<Centroid> *()> getVocab;
	function<const vector<FastaSequence *> *()> getDatabase;
	function<size_t()> getKmerLength;
	function<const LookupTable_<size_t, const FastaSequence> *()> getSeqIndex;
	function<Distance()> getThreshold;
	function<int()> getClassIndex;
	function<int()> getNameIndex;
	function<const vector<size_t> * ()> getTrainingSetIndex;
#pragma endregion

	void GainFocus() override {
		// cerr << "GainFocus\n";

		Page::GainFocus();
		matrix = getMatrix();
		vocab = getVocab();
		dbSeqs = getDatabase();
		kmerLength = getKmerLength();
		seqIndex = getSeqIndex();
		threshold = getThreshold();
		classIndex = getClassIndex();
		nameIndex = getNameIndex();
		trainingSetIndex = getTrainingSetIndex();
		UpdateReadiness();

		sigBrowser.SetClassIndex( classIndex );
		sigBrowser.SetNameIndex( nameIndex );
		sigBrowser.SetDbSeqs( dbSeqs, &signatures );
		sigBrowser.SetVocabSize( vocab->size() );
	}

	YesNoCancel CheckLoad() {
		auto fileName = String::Trim( DbSigFileName() );

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
			if ( threshold == 0 ) str << "Threshold must be greater than zero.\n";
			if ( kmerLength == 0 ) str << "K-mer length must be greater than zero.\n";
			if ( !dbSeqs ) str << "Sequence dataset not initialised!\n";
			if ( vocab->size() == 0 ) str << "Vocabulary contains no terms.\n";
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

		SparseSignature::CreatePostingList( signatures, *trainingSetIndex, featurePostingList );

		sigIndex.clear();

		for ( size_t i = 0; i < signatures.size(); i++ ) {
			auto & sig = signatures[i];
			sigIndex[sig.Sequence()->Id()] = &sig;
		}
	}

	static void RunImpl_( EncodePartition * self ) {
		self->RunImpl();
	}

	void RunImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		auto fileName = String::Trim( DbSigFileName() );

		UPDATE_GUI( SetOk( false ); disp.Clear(); );

#if INTERLEAVE_SIG
		try {
			Page::runTime = -omp_get_wtime();
			GenerateSignatures( fileName );
			IndexDatabase();
			Page::runTime += omp_get_wtime();

			Page::saveTime = 0;
		}
		catch ( Exception & ex ) {
			ShowProgress( "Unable to encode and save signatures: %s\n", ex.what() );
		}
#else
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
			ShowProgress( "Unable to save signatures: %s\n", ex.what() );
		}
#endif

		UPDATE_GUI( UpdateReadiness(); NotifyRunComplete(); );

		mute.unlock();
	}

#if INTERLEAVE_SIG
	void GenerateSignatures( const string & fileName ) {
#else
	void GenerateSignatures() {
#endif
		featuresOk = false;
		const size_t N = dbSeqs->size();

#if INTERLEAVE_SIG
		ofstream f( fileName );
		CsvWriter w( f );

		w << FileTag << N << '\n';
#endif

#if WANT_FEEDBACK
		size_t done = 0;
#endif

		signatures.clear();
		signatures.resize( N );

#if USE_OMP
#pragma omp parallel for schedule(dynamic,1)
#endif
		for ( size_t i = 0; i < N; i++ ) {
			// SigRank::GenerateSignature( *((*dbSeqs)[i]), kmerLength, *vocab, threshold, *matrix, signatures[i] );
			SigRank::GenerateSignatureDigrams( *((*dbSeqs)[i]), kmerLength, *vocab, threshold, *matrix, signatures[i] );

#if INTERLEAVE_SIG
#if USE_OMP
#pragma omp critical
			{
				w << signatures[i];
			}
#endif
#endif

#if false || WANT_FEEDBACK
#if USE_OMP
#pragma omp critical
#endif
			{
				done++;
				int thisPercent = done * 100 / N, lastPercent = (done - 1) * 100 / N;

				if ( thisPercent != lastPercent ) {
					UPDATE_GUI( disp( "Encoding sequences: %d%%\n", thisPercent ); )
				}
			}
#endif
		}

		featuresOk = true;
#undef symbolDistance
	}

	const string FileTag = "signatures";

	void Save( const string & fileName ) {
		ofstream f( fileName );
		CsvWriter w( f );

		w << FileTag << signatures.size() << '\n';

		for ( auto &c : signatures ) {
			w << c;
		}
	}

	static void LoadImpl_( EncodePartition * self ) {
		self->LoadImpl();
	}

	void LoadImpl() {
		if ( !mute.try_lock() ) {
			return;
		}

		auto fileName = String::Trim( DbSigFileName() );

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

		SparseSignature::Lookup = [this]( const string & id ) {
			auto idNumber = FastaSequence::Register( id );
			auto i = seqIndex->find( idNumber );
			return i == seqIndex->end() ? nullptr : i->second;
		};

		signatures.clear();
		signatures.resize( N );

		for ( size_t i = 0; (!r.IsEOF()) && (i < N); i++ ) {
			auto &sig = signatures[i];
			r >> sig;

			auto vocabSize = vocab->size();
			auto sigMax = sig.Size() == 0 ? 0 : sig.Max();

			if ( sigMax >= vocabSize ) {
				UPDATE_GUI( disp( "Feature out of range for vocabulary: load abandoned.\n" ) );

				for ( size_t j = 0; j <= i; j++ ) {
					signatures[j].Clear();
				}

				return;
			}
		}

		featuresOk = true;
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
		SetReady( String::Trim( DbSigFileName() ).length() > 0 );
		SetOk( Ready() && featuresOk );
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
		addPar( dbSigFileName );
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
			setPar( dbSigFileName, string );
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
	const vector<SparseSignature> * Signatures() const { return &signatures; }

	const vector<vector<size_t>> *PostingList() const { return &featurePostingList; }

	const unordered_map<size_t, SparseSignature *> * SigIndex() const { return &sigIndex; }
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

		sigBrowser.Clear();
		signatures.clear();
		matrix = 0;
		vocab = 0;
		dbSeqs = 0;
		kmerLength = 0;
		seqIndex = 0;
		threshold = 0;
		classIndex = 0;
		nameIndex = 0;
		trainingSetIndex = 0;
	}
};
