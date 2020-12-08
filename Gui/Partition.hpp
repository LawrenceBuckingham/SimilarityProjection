#pragma once

#include <FL/Enumerations.H>
#include <FL/Fl.H>
#include <FL/Fl_Chart.H>
#include <FL/Fl_Tabs.H>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

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
#include <TextDisplay.hpp>

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
class Partition :
	public virtual Page,
	public virtual IPropertyChangedEventHandler {
public:
	/// User inputs and computed outputs from this component.
	struct State {
		/// User-supplied: seed for RNG.
		ulong seed;

		/// User-supplied: sample size.
		size_t sampleSize;

		/// User-supplied: trainingSetFilename
		string partitionFileName;

		/// Address of input dataset.
		const vector<FastaSequence *> *dbSeqs;

		/// Computed: subset of
		const vector<size_t> * trainingSetIndices;

		/// Computed:
		const vector<size_t> * testSetIndices;

		/**
		 *	Constructor for GetKmerDistances::State. If the supplied gui is not
		 *	null, then its contents are loaded into the new State.
		 *
		 *	@param[in] gui	Address of a GetKmerDistances control.
		 */
		State( const Partition *gui ) {
			if ( gui ) Snapshot( gui );
		}

		/**
		 *	Compare two objects for equality.
		 *
		 *	@param	other	Reference to another state object for comparison.
		 *	@returns		Returns true if and only if the user-supplied fields
		 *					are all equal.
		 */
		bool Equals( const State &other ) {
			return seed == other.seed && sampleSize == other.sampleSize;
		}

		/**
		 *	Loads the instance with values obtained from the host control.
		 */
		void Snapshot( const Partition *gui ) {
			// Get user-supplied values.
			this->seed = gui->seed->Value();
			this->sampleSize = gui->sampleSize->Value();
			this->partitionFileName = gui->partitionFileName->Value();

			// Get computed values.
			this->dbSeqs = gui->dbSeqs;
			this->trainingSetIndices = &gui->trainingSetIndices;
			this->testSetIndices = &gui->testSetIndices;
		}
	};

#pragma region
private:
	const int rows = 5;
	int rowHeight;
	Arg<ulong, Fl_Input> *seed = 0;
	Arg<size_t, Fl_Input> *sampleSize = 0;
	TypedArg<string, OutFileChooser> *partitionFileName = 0;
	GridLayout *stack = 0;
	Fl_Box *titleBox = 0;
	Fl_Text_Buffer *buff = 0;
	TextDisplay *disp = 0;

	State state;

	const vector<FastaSequence *> *dbSeqs = 0;
	UniformRealRandom *rand = 0;
	vector<size_t> trainingSetIndices;
	vector<size_t> testSetIndices;
#pragma endregion  // private data members

public:
	Partition( int left, int top, int width, int height,
		int labelWidth = 150, int rowHeight = 25, int vgap = 5 )
		: Page( left, top, width, height, "Partition" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		state( nullptr ),
		rand( nullptr ),
		autoLoad( "autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load the partition file if it exists.",
			labelWidth, 0, width, rowHeight )
		//
	{
		add( stack = new GridLayout( 0, 0, w(), rows * rowHeight, rows, 1 ),
			Location::North );

		stack->add( titleBox = new Fl_Box( 0, 0, width, rowHeight, name.c_str() ) );
		titleBox->align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		stack->add( seed = new Arg<ulong, Fl_Input>(
			"seed", time( nullptr ), "RNG Seed", Requirement::Required,
			"Cardinal number; required. Seed for random number generator used\n"
			"to sample k-mer population to compute empirical distances, and \n"
			"available for downstream components if necessary.",
			labelWidth, 0, width, rowHeight ) );

		stack->add( sampleSize = new Arg<size_t, Fl_Input>(
			"sampleSize", 1000, "Sample Size",
			Requirement::Required,
			"Cardinal number; required. Number of ",
			labelWidth, 0, width, rowHeight ) );

		stack->add( partitionFileName = new TypedArg<string, OutFileChooser>(
			"partitionFileName", "", "[out] Partition File", Requirement::Required,
			"The path to a file which will be overwritten with the partitioned dataset \n"
			"tabulated in CSV format. The training set will be written first, followed \n"
			"by the test set.",
			labelWidth, 0, width, rowHeight ) );

		stack->add( &autoLoad );

		seed->AddPropertyChangedEventHandler( this );
		sampleSize->AddPropertyChangedEventHandler( this );
		partitionFileName->AddPropertyChangedEventHandler( this );
		partitionFileName->inputField.SetPattern("*.partition.csv");

		state.Snapshot( this );

		add( disp = new TextDisplay( 0, 0, w(), 100 ), Location::Centre );
		disp->textfont( FL_COURIER );
		disp->textsize( 16 );

		(*disp)("Partition Data Set: select the number of test points to hold out, then click \"Update\".\n");
	}

	virtual ~Partition() {
		delete this->rand;
		delete this->seed;
		delete this->sampleSize;
		delete this->buff;
		delete this->disp;
		delete stack;
	}

	function<const vector<FastaSequence *> *()> GetDatabase;

	void GainFocus() override {
		dbSeqs = GetDatabase();
		UpdateReadiness();
	}

	void UpdateReadiness() {
		SetReady( dbSeqs
			&& dbSeqs->size() > 0
			&& seed->HasValue()
			&& sampleSize->HasValue()
			&& String::Trim( partitionFileName->Value() ).length() > 0 );

		SetOk( Ready() && trainingSetIndices.size() > 0 && testSetIndices.size() > 0 );
	}

	State & GetState() {
		state.Snapshot( this );
		return state;
	}

	/**
	 *	Generates theoretical and empirical distance distributions. Allows user
	 *	to specify a p-value for the clustering step.
	 */
	void Run() {
		SetOk( false );

		if ( partitionFileName->Value() == "" ) {
			fl_message( "Output file name for partition must be specified." );
			return;
		}

		if ( !seed->HasValue() || !sampleSize->HasValue() ) {
			fl_message( "Seed and Sample size must have numeric values!" );
			return;
		}

		auto shouldLoad = CheckLoad();

		if ( shouldLoad == YesNoCancel::Yes ) {
			Load();
			return;
		}
		else if ( shouldLoad == YesNoCancel::Cancel ) {
			NotifyRunComplete();
			return;
		}
		else {
			RunImpl();
		}
	}

	void RunImpl() {
		double runStart = omp_get_wtime();

		disp->Clear();

		State currentState( this );

		delete rand;
		rand = new UniformRealRandom( currentState.seed );

		Selector sel( *rand, currentState.sampleSize, dbSeqs->size() );

		trainingSetIndices.clear();
		testSetIndices.clear();

		for ( size_t i = 0; i < dbSeqs->size(); i++ ) {
			if ( sel.SelectThis() ) {
				testSetIndices.push_back( i );
			}
			else {
				trainingSetIndices.push_back( i );
			}
		}

		(*disp)(
			"Database partitioned\n"
			"Training set sequence count:\t%zu\n"
			"Test set sequence count:\t%zu\n",
			trainingSetIndices.size(),
			testSetIndices.size()
			);

		state = currentState;
		UpdateClassPostingList();
		
		Page::runTime = omp_get_wtime() - runStart;

		SaveFiles();
		UpdateReadiness();
		NotifyRunComplete();
	}

	void SaveFiles() {
		double start = omp_get_wtime();
		ofstream f( partitionFileName->Value() );
		CsvWriter writer( f );
		(writer << Name()).Ln();
		(writer << "Seed" << seed->Value()).Ln();
		(writer << "SampleSize" << sampleSize->Value()).Ln();
		Save( writer, trainingSetIndices, "Training set" );
		Save( writer, testSetIndices, "Test set" );
		Page::saveTime = omp_get_wtime() - start;
	}

	void Save( CsvWriter & writer, const vector<size_t> & db, const string & tableName ) {
		writer << tableName;

		for ( auto seqIdx : db ) {
			writer << seqIdx;
		}

		writer.Ln();
	}

	YesNoCancel CheckLoad() {
		if ( String::Trim( partitionFileName->Value() ).length() > 0 && File::Exists( partitionFileName->Value() ) ) {
			if ( AutoLoad() ) {
				return YesNoCancel::Yes;
			}

			int answer = fl_choice(
				"The select file already exists\n"
				"Would you like to load that partition?",
				"Cancel",
				"Yes, load the partition.",
				"No, I want to overwrite the partition."
			);

			return (YesNoCancel) answer;
		}
		else {
			return YesNoCancel::No;
		}
	}

	virtual void PropertyChanged( void *sender, const string &propertyName ) override {
		//cerr << "Partition: Received PropertyChanged(" << propertyName << ")\n";

		//if (sender == static_cast<void *>(partitionFileName)) cerr << "sender = partitionFileName\n";
		//else if (sender == static_cast<void *>(sampleSize)) cerr << "sender = sampleSize\n";
		//else if (sender == static_cast<void *>(seed)) cerr << "sender = seed\n";
		//else cerr << "Dunno who sent this!\n";

		UpdateReadiness();
	}

	/**
	 *	Reads a record from the CSV save file, and ensures that the first word in the row matches a designated label.
	 *	Further check (optional) is made to ensure that the required number of field values have been obtained.
	 *	@param r Reference to a CsvReader from which records are obtained.
	 *	@param label The required contents of the first field in the record.
	 *	@param expectedLength The expected number of fields in the record, including the label. If this is zero, then an arbitrary number of fields is accepted.
	 *	@param[out] record A vector which will be populated with the fields obtained from the reader.
	 *	@throws Exception is thrown if record size is zero, first field is not equal to label, or expected length is greater than zero and record size is not equal to expected length.
	 */
	void Read(
		CsvReader & r,
		const string & label,
		size_t expectedLength,
		vector<string> & record
	) {
		record.clear();
		r.ReadRecord( record );
		auto tagNames = String::Split(label, '|');

		if ( record.size() == 0 || std::find(tagNames.begin(), tagNames.end(), record[0]) == tagNames.end() ) {
			throw Exception( String::Format( "Expected text '%s' not found. Load abandoned.\n", label.c_str() ), FileAndLine );
		}
		if ( expectedLength > 0 && record.size() != expectedLength ) {
			throw Exception( String::Format( "Expected %zu fields in record '%s', but found %zu.\n",
				expectedLength, label.c_str(), record.size() ), FileAndLine );
		}
	}

	void Load() {
		double loadStart = omp_get_wtime();

		auto fName = partitionFileName->Value();
		(*disp)("Attempting to load partition data from %s\n", fName.c_str());
		trainingSetIndices.clear();
		testSetIndices.clear();
		ifstream f( fName );
		CsvReader reader( f );

		try {
			vector<string> record;
			Read( reader, Name() + "|" + TagName(), 1, record );

			Read( reader, "Seed", 2, record );
			seed->SetValue( Ulong::Parse( record[1] ) );

			Read( reader, "SampleSize", 2, record );
			sampleSize->SetValue( Ulong::Parse( record[1] ) );

			Read( reader, "Training set", 0, record );

			for ( size_t i = 1; i < record.size(); i++ ) trainingSetIndices.push_back( Ulong::Parse( record[i] ) );

			Read( reader, "Test set", 0, record );

			for ( size_t i = 1; i < record.size(); i++ ) testSetIndices.push_back( Ulong::Parse( record[i] ) );

			(*disp)("%zu training set indices loaded.\n", trainingSetIndices.size());
			(*disp)("%zu test set indices loaded.\n", testSetIndices.size());
		}
		catch ( Exception & e ) {
			(*disp)("%s\nLoad abandoned.\n", e.what());
		}

		UpdateClassPostingList();

		Page::loadTime = omp_get_wtime() - loadStart;

		UpdateReadiness();
		NotifyRunComplete();
	}

	/// Get parameters so they can be serialised when program finishes.
	/// @param parms A Param list to which the component should append its current parameters. 
	void GetParams( set<Param> &parms ) override {
		parms.emplace( Name(), "seed", seed->Value() );
		parms.emplace( Name(), "sampleSize", sampleSize->Value() );
		parms.emplace( Name(), "partitionFileName", partitionFileName->Value() );
		parms.emplace( Name(), "autoLoad", autoLoad.Value() );

		parms.emplace( Name(), "runTime", Util::ToString(runTime) );
		parms.emplace( Name(), "loadTime",  Util::ToString(loadTime) );
		parms.emplace( Name(), "saveTime",  Util::ToString(saveTime) );
	}

	private:
		string tagName;

	public:
		const string & TagName() const {
			static string defaultVal = "Partition Database";
			return tagName.length() == 0 ? defaultVal : tagName;
		}

		void SetTagName( const string & value ) {
			tagName = value;
		}

	/// Set parameters which have been de-serialised when program starts.
	/// @param parms A Param list containing de-serialised parameters.
	void SetParams( const set<Param> & parms ) override {
		for ( auto & p : parms ) {
			if ( p.ComponentName() != Name() && p.ComponentName() != TagName() ) continue;

			if ( p.ParamName() == "seed" ) seed->SetValue( p.Value<ulong>() );
			if ( p.ParamName() == "sampleSize" ) sampleSize->SetValue( p.Value<size_t>() );
			if ( p.ParamName() == "partitionFileName" ) partitionFileName->SetValue( p.Value() );
			if ( p.ParamName() == "autoLoad" ) autoLoad.SetValue( p.Value<int>() );

			if ( p.ParamName() == "runTime" )  runTime = p.Value<double>();
			if ( p.ParamName() == "loadTime" ) loadTime = p.Value<double>();
			if ( p.ParamName() == "saveTime" ) saveTime = p.Value<double>();
		}

		UpdateReadiness();
	}

	const vector<size_t> * TrainingSetIndices() const { return &trainingSetIndices; }

	const vector<size_t> * TestSetIndices() const { return &testSetIndices; }

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

private:
	unordered_map<size_t, vector<size_t>> classPostingList;

	void UpdateClassPostingList() {
		classPostingList.clear();

		for ( auto idx: trainingSetIndices ) {
			auto seq = (*dbSeqs)[idx];

			for ( auto classId : seq->Classes() ) {
				classPostingList[classId].push_back(idx);
			}
		}
	}

public:
	/**
	 *	Gets an index structure enabling direct access to all sequences in the 
	 *	training set assigned to each class label.
	 *	@returns A map from classId to list of reference sequence indices. 
	 */
	const unordered_map<size_t, vector<size_t>> * ClassPostingList() { 
		return &classPostingList; 
	}

	void Reset() override {
		dbSeqs = 0;
		trainingSetIndices.clear();
		testSetIndices.clear();
		classPostingList.clear();
	}
};
