#pragma once

#include <FL/Enumerations.H>
#include <FL/Fl.H>
#include <FL/Fl_Chart.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Int_Input.H>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>

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
#include <Index.hpp>

/**
<summary>
	DbLoader gathers the parameters for an interactive version of DnaClust,
	with some additional analytic displays to render the software a little more
	usable.
</summary>
 */
class DbLoader :
	public Page,
	public virtual IPropertyChangedEventHandler
	//
{

private:
	const uint rows = 10;
	int rowHeight;
	Fl_Box titleBox;

#pragma region DbFileName
private:
	Arg<string, InFileChooser> dbFileName;

public:
	string DbFileName() {
		return dbFileName.Value();
	}

	void SetDbFileName( const string & val ) {
		if ( val == DbFileName() ) return;

		this->dbFileName.SetValue( val );
	}
#pragma endregion DbFileName

#pragma region Alphabet
private:
	TypedArg<pAlphabet, AlphabetChooser> alphabet;

public:
	pAlphabet Alphabet() {
		return alphabet.Value();
	}

	void SetAlphabet( const pAlphabet val ) {
		if ( val == Alphabet() ) return;

		this->alphabet.SetValue( val );
	}
#pragma endregion Alphabet

#pragma region MatrixFileName Property
private:
	Arg<string, InFileChooser> matrixFileName;

public:
	string MatrixFileName() {
		return matrixFileName.Value();
	}

	void SetMatrixFileName( const string & val ) {
		if ( val == MatrixFileName() ) return;

		this->matrixFileName.SetValue( val );
	}
#pragma endregion

#pragma region IdIndex
private:
	Arg<uint, Fl_Int_Input> idIndex;

public:
	uint IdIndex() {
		return idIndex.Value();
	}

	void SetIdIndex( const uint val ) {
		if ( val == IdIndex() ) return;

		this->idIndex.SetValue( val );
	}
#pragma endregion IdIndex

#pragma region ClassIndex
private:
	Arg<int, Fl_Int_Input> classIndex;

public:
	int ClassIndex() {
		return classIndex.Value();
	}

	void SetClassIndex( const int val ) {
		if ( val == ClassIndex() ) return;

		this->classIndex.SetValue( val );
	}
#pragma endregion ClassIndex

#pragma region NameIndex Property
private:
	Arg<int, Fl_Int_Input> nameIndex;

public:
	int NameIndex() {
		return nameIndex.Value();
	}

	void SetNameIndex( const int val ) {
		if ( val == NameIndex() ) return;

		this->nameIndex.SetValue( val );
	}
#pragma endregion

#pragma region ReverseComplement
private:
	TypedArg<int, Fl_Check_Button> reverseComplement;

public:
	int ReverseComplement() {
		return reverseComplement.Value();
	}

	void SetReverseComplement( const int val ) {
		if ( val == ReverseComplement() ) return;

		this->reverseComplement.SetValue( val );
	}
#pragma endregion ReverseComplement

#pragma region PreviewSize
private:
	Arg<uint, Fl_Int_Input> previewSize;

public:
	uint PreviewSize() {
		return previewSize.Value();
	}

	void SetPreviewSize( const uint val ) {
		if ( val == PreviewSize() ) return;

		this->previewSize.SetValue( val );
	}
#pragma endregion PreviewSize

#pragma region PreviewSeqLength
private:
	Arg<uint, Fl_Int_Input> previewSeqLength;

public:
	uint PreviewSeqLength() {
		return previewSeqLength.Value();
	}

	void SetPreviewSeqLength( const uint val ) {
		if ( val == PreviewSeqLength() ) return;

		this->previewSeqLength.SetValue( val );
	}
#pragma endregion PreviewSeqLength

#pragma region DbSeqs
private:
	vector<pFastaSequence> dbSeqs;

public:
	vector<pFastaSequence> *DbSeqs() {
		return &dbSeqs;
	}
#pragma endregion DbSeqs

#pragma region SeqIndex
private:
	LookupTable_<size_t, const FastaSequence> seqIndex;

public:
	LookupTable_<size_t, const FastaSequence> * SeqIndex() {
		return &seqIndex;
	}
#pragma endregion SeqIndex

	SimilarityMatrix matrix;

	GridLayout stack;
	TextDisplay disp;

	vector<map<string, size_t>> metadataCounters;

public:
	DbLoader( int left, int top, int width, int height, int labelWidth = 150, int rowHeight = 25, int vgap = 5
	) :
		Page( left, top, width, height, "Load Seqs" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, width, rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),

		dbFileName( "dbFileName", "", "[in] DB File", Requirement::Optional,
			"The path to a list of DNA/RNA sequences represented in FASTA format.\n"
			"The definition line is assumed to contain a list of fields, \n"
			"including the sequence Id and possibly a list of class names.\n"
			"The fields are traditionally separated by '|' (pipe) symbols. \n"
			"These are numbered from 0. Particular fields of note are idIndex\n"
			"and classIndex.",
			labelWidth, 0, width, rowHeight ),

		alphabet( "alphabet", Alphabets::DNA(), "alphabet", Requirement::Required,
			"The alphabet from which symbols in the input sequences are drawn.",
			labelWidth, 0, width, rowHeight ),

		matrixFileName( "matrixFileName", "", "[in] Matrix File", Requirement::Optional,
			"The path to a plain text file containing a tabulated symbol \n"
			"similarity matrix.",
			labelWidth, 0, width, rowHeight ),

		idIndex( "idIndex", 0, "Sequence Id Index", Requirement::Optional,
			"The zero-origin field number of the sequence ID in the pipe-separated multi-field definition line.",
			labelWidth, 0, width, rowHeight ),

		classIndex( "classIndex", -1, "Sequence Class Index", Requirement::Optional,
			"If there is a field containing class-names, this classIndex tells\n"
			"us which one it is (classIndex >= 0). Class names are strings,\n"
			"separated by semicolons. If there are no class names, use a\n"
			"negative value for classIndex.",
			labelWidth, 0, width, rowHeight ),

		nameIndex( "nameIndex", -1, "Sequence Name Index", Requirement::Optional,
			"If there is a field containing sequence names, this field tells\n"
			"us which one it is (nameIndex >= 0). Sequence names are strings.\n"
			"If there is no sequence name, use a negative value for nameIndex.",
			labelWidth, 0, width, rowHeight ),

		reverseComplement( "reverseComplemnt", 0, "Complement?", Requirement::Required,
			"Should we add the reverse-complement of each sequence to the database?",
			labelWidth, 0, width, rowHeight ),

		previewSize( "previewSize", 10, "Preview Size", Requirement::Required,
			"The number of sequences to preview in the status display.",
			labelWidth, 0, width, rowHeight ),

		previewSeqLength(
			"previewSeqLength", 100, "Preview Seq Len", Requirement::Required,
			"The number of leading charac.",
			labelWidth, 0, width, rowHeight ),

		disp( 0, 0, 0, 0 )
		//
	{
		add( &stack, Location::North );

		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );
		reverseComplement.SetValue( 1 );

		stack.add( &titleBox );
		stack.add( &dbFileName );
		stack.add( &alphabet );
		stack.add( &matrixFileName );
		stack.add( &idIndex );
		stack.add( &classIndex );
		stack.add( &nameIndex );
		stack.add( &reverseComplement );
		stack.add( &previewSize );
		stack.add( &previewSeqLength );

		dbFileName.AddPropertyChangedEventHandler( this );
		alphabet.AddPropertyChangedEventHandler( this );
		matrixFileName.AddPropertyChangedEventHandler( this );
		idIndex.AddPropertyChangedEventHandler( this );
		classIndex.AddPropertyChangedEventHandler( this );
		nameIndex.AddPropertyChangedEventHandler( this );
		reverseComplement.AddPropertyChangedEventHandler( this );
		previewSize.AddPropertyChangedEventHandler( this );
		previewSeqLength.AddPropertyChangedEventHandler( this );

		add( &disp, Location::Centre );
		disp.textfont( FL_COURIER );
		disp.textsize( 16 );
	}

	virtual ~DbLoader() {
		Util::Free( dbSeqs );
	}

	void PropertyChanged( void* sender, const string& propertyName ) override {
		if ( sender == (void *) &classIndex ) {
			UpdateClasses();
		}

		if ( sender == (void *) &nameIndex ) {
			UpdateNames();
		}

		UpdateReadiness();
	}

	/// Update the class sets for database sequences.
	void UpdateClasses() {
		for ( auto seq : dbSeqs ) {
			seq->SetClassIndex( ClassIndex() );
		}
	}

	/// Update the name index field of database sequences.
	void UpdateNames() {
		for ( auto seq : dbSeqs ) {
			seq->SetNameIndex( NameIndex() );
		}
	}

	/**
	 *	Returns true if and only if the input fields validate and the output fields contain useble results.
	 */
	void UpdateReadiness() {
		SetReady( (DbFileName().length() > 0) && (Alphabet() != 0) && (MatrixFileName().length() > 0) );

		SetOk( Ready() &&
			(dbSeqs.size() > 0) &&
			(metadataCounters.size() > IdIndex()) &&
			(metadataCounters[IdIndex()].size() == dbSeqs.size())
		);
	}

	/**
	 *	Attempts to load the designated dataset and pack it into EncodedSequence form.
	 */
	void Run() {
		double runStart = omp_get_wtime();
		double loadStart = omp_get_wtime();

		SetReady( false );
		SetOk( false );
		disp.Clear();

		auto dbFileName = DbFileName();
		auto alphabet = Alphabet();
		auto matrixFileName = MatrixFileName();
		auto idIndex = IdIndex();
		auto classIndex = ClassIndex();
		auto nameIndex = NameIndex();
		auto reverseComplement = ReverseComplement();
		auto previewSize = PreviewSize();
		auto previewSeqLength = PreviewSeqLength();

		try {
			ifstream matrixFile( matrixFileName );
			matrix.SetAlphabet(alphabet).Parse(matrixFile);
		}
		catch ( Exception & ex ) {
			disp( "Error while loading similarity matrix:\n%s\n", ex.what() );
			return;
		}

		try {
			ifstream fasta( dbFileName );
			Util::Free( dbSeqs );
			dbSeqs.clear();
			Load::Fasta( fasta, idIndex, alphabet, dbSeqs );
			seqIndex.clear();

			for ( auto seq : dbSeqs ) {
				seqIndex[seq->Id()] = seq;
			}
			Page::saveTime = omp_get_wtime() - loadStart;
		}
		catch ( Exception &ex ) {
			disp( "Error while loading sequence data:\n%s\n", ex.what() );
			return;
		}

		int N = dbSeqs.size();

		if ( N == 0 ) {
			disp( "Sequence file '%s' contains no sequences.\n", dbFileName.c_str() );
		}
		else {
			disp( "%d sequences loaded from %s.\n", N, dbFileName.c_str() );
		}

		if ( reverseComplement ) {
			for ( int i = 0; i < N; i++ ) {
				auto seq = dbSeqs[i];
				vector<string> metadata;

				for ( size_t i = 0; i < seq->MetaCount(); i++ ) {
					if ( i == idIndex ) {
						metadata.push_back( seq->Metadata( i ) + '_' );
					}
					else {
						metadata.push_back( seq->Metadata( i ) );
					}
				}

				auto complementDefLine = String::Join( metadata, "|" );
				auto s = alphabet->ReverseComplement( seq->CharData() );
				auto comp = new FastaSequence( complementDefLine, s, idIndex, alphabet );
				dbSeqs.push_back( comp );
				seqIndex[comp->Id()] = comp;
			}

			disp( "%d reverse-complement sequences added to database.\n", N );
			N *= 2;
		}

		UpdateClasses();
		UpdateNames();

		size_t totalChars = 0;

		for ( auto seq : dbSeqs ) {
			totalChars += seq->Length();
		}

		disp( "Total database size: %d sequences, %zu characters.\n", N, totalChars );

		FastaSequence::GetMetadataCounts( dbSeqs, metadataCounters );

		disp( "Metadata Column\tDistinct Values\n" );

		for ( size_t i = 0; i < metadataCounters.size(); i++ ) {
			disp( "%d\t%zu\n", i, metadataCounters[i].size() );
		}

		for ( uint i = 0; i < previewSize && i < dbSeqs.size(); i++ ) {
			FastaSequence *seq = dbSeqs[i];
			disp( "\n>%s\n", seq->DefLine().substr( 0, previewSeqLength ).c_str() );
			string s = seq->CharData();
			disp( "%s\n", s.substr( 0, std::min( (size_t) previewSeqLength, s.length() ) ).c_str() );
			disp( "Id = %zu\nIdStr = %s\nName = %s\nClass = %s\n",
				seq->Id(),
				seq->IdStr().c_str(),
				(nameIndex >= 0 && nameIndex < (int) seq->MetaCount() ? seq->Metadata( nameIndex ).c_str() : "n/a"),
				(classIndex >= 0 && classIndex < (int) seq->MetaCount() ? seq->Metadata( classIndex ).c_str() : "n/a")
			);
			if ( classIndex >= 0 && classIndex < (int) seq->MetaCount() ) {
				for ( auto classe : seq->Classes() ) {
					disp( "Member of class %u\n", classe );
				}
			}
		}

		Page::runTime = omp_get_wtime() - loadStart;

		UpdateReadiness();
		NotifyRunComplete();
	}

	void GetParams( set<Page::Param> & parms ) override {
#define add_parm(x) parms.emplace( Name(), STRING(x), (x).Value() )
#define add_parm_name(x) parms.emplace( Name(), STRING(x), (x).Value()->Name() )
#define add_val(x) parms.emplace( Name(), STRING(x), Util::ToString(x) )

		add_parm_name( alphabet );
		add_parm( dbFileName );
		add_parm( matrixFileName );
		add_parm( idIndex );
		add_parm( classIndex );
		add_parm( nameIndex );
		add_parm( reverseComplement );
		add_parm( previewSize );
		add_parm( previewSeqLength );
		add_val( runTime );
		add_val( loadTime );
		add_val( saveTime );

#undef add_parm
#undef add_parm_name
#undef add_val
	}

	void SetParams( const set<Param> & parms ) {
#define setParm(x,t) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value<t>() )
#define setString(x) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value() )
#define setVal(x,v) if ( p.ParamName() == STRING(x) ) (x).SetValue( v )
#define get_val(x, T) if ( p.ParamName() == STRING(x) ) (x) = p.Value<T>()

		for ( auto & p : parms ) {
			if ( p.ComponentName() != Name() ) continue;

			setVal( alphabet, Alphabets::ByName( p.Value() ) );
			setString( dbFileName );
			setString( matrixFileName );
			setParm( idIndex, int );
			setParm( classIndex, int );
			setParm( nameIndex, int );
			setParm( reverseComplement, int );
			setParm( previewSize, uint );
			setParm( previewSeqLength, uint );
			get_val( runTime, double );
			get_val( loadTime, double );
			get_val( saveTime, double );
		}

		UpdateReadiness();
#undef setParm
#undef setString
#undef setVal
#undef get_val
	}

	void Reset() override {
		dbSeqs.clear();
		UpdateReadiness();
	}

	SimilarityMatrix * Matrix() { return &matrix; }
};
