#include "IncludeAll.hpp"
#include <set>

#include <LBFL/DoubleWindow.hpp>

#include "DbLoader.hpp"
#include "GetKmerDistances.hpp"
#include "Partition.hpp"
#include "ClusterTrainingSet.hpp"
#include "SelectVocab.hpp"

#include "EncodePartition.hpp"
#include "RankTrainingSet.hpp"

#include <IPropertyChangedEventHandler.hpp>

#undef TRON
// #define TRON 1
#include <db.hpp>

struct DnaClustGui :
	virtual public DoubleWindow,
	virtual public IPropertyChangedEventHandler
	//
{
	const int rowHeight = 25;
	const int BW = 150;
	uint currentPage = 0;
	vector<Page *> pages;

	Button nextButton;
	Button prevButton;
	Button runButton;
	Button saveButton;
	Button loadButton;

	BorderLayout middle;
	FlowLayout bottom;

	DbLoader dbLoader;
	Partition partition;
	GetKmerDistances getKmerDistances;
	ClusterTrainingSet clusterTrainingSet;
	SelectVocab selectVocab;
	EncodePartition encodePartition;
	RankTrainingSet rankTrainingSet;

	bool autoRun = false;
	bool autoTerm = false;

	DnaClustGui( int W, int H, const char *l = 0 )
		: DoubleWindow( W, H, l ),
		nextButton( 0, 0, BW, rowHeight, "Next", [this]( Button *b ) { Forward(); } ),
		prevButton( 0, 0, BW, rowHeight, "Prev", [this]( Button *b ) { Back(); } ),
		runButton( 0, 0, BW, rowHeight, "Run", [this]( Button *b ) { Run(); } ),
		saveButton( 0, 0, BW, rowHeight, "Save Config", [this]( Button *b ) { SaveConfig(); } ),
		loadButton( 0, 0, BW, rowHeight, "Load Config", [this]( Button *b ) { LoadConfig(); } ),
		middle( 0, 0, W, H - rowHeight ),
		bottom( 0, H - rowHeight, W, rowHeight ),
		dbLoader( 0, 0, W, H - rowHeight ),
		getKmerDistances( 0, 0, W, H - rowHeight ),
		partition( 0, 0, W, H - rowHeight ),
		clusterTrainingSet( 0, 0, W, H - rowHeight ),
		selectVocab( 0, 0, W, H - rowHeight ),
		encodePartition( 0, 0, W, H - rowHeight ),
		rankTrainingSet( 0, 0, W, H - rowHeight )
		// 
	{
		add( &middle );
		resizable( &middle );
		add( &bottom );
		middle.add( &dbLoader );

		bottom.add( &loadButton );
		bottom.add( &saveButton );
		bottom.add( &prevButton );
		bottom.add( &runButton );
		bottom.add( &nextButton );

		pages.push_back( &dbLoader );
		pages.push_back( &partition );
		pages.push_back( &getKmerDistances );
		pages.push_back( &clusterTrainingSet );
		pages.push_back( &selectVocab );
		pages.push_back( &encodePartition );
		pages.push_back( &rankTrainingSet );

		function<void( Page * )> pageRunComplete = [this]( Page * page ) { PageRunComplete( page ); };

		for ( auto page : pages ) {
			page->RunComplete( pageRunComplete );
		}

		WireUp();

		Fl::add_handler( WindowClosing );
		LoadSettings();
		EnableButtons();
	}

	virtual ~DnaClustGui() {
		SaveSettings();
	}

	void PropertyChanged( void* sender, const string& propertyName ) override {
		// cerr << "DnaClustGui: Received PropertyChanged(," << propertyName << ")\n";

		bool Ok = 0;
		bool Ready = 0;

		if ( propertyName == NAMEOF( Ok ) || propertyName == NAMEOF( Ready ) ) {
			TRACE;
			EnableButtons();
		}
	}

	void WireUp() {
		dbLoader.AddPropertyChangedEventHandler( this );
		dbLoader.AddConsumer( &partition );

		auto getMatrix = [this]() { return dbLoader.Matrix(); };

		partition.AddPropertyChangedEventHandler( this );
		partition.GetDatabase = [this]() { return dbLoader.DbSeqs(); };
		partition.AddConsumer( &getKmerDistances );

		getKmerDistances.AddPropertyChangedEventHandler( this );
		getKmerDistances.getAlphabet = [this]() { return dbLoader.Alphabet(); };
		getKmerDistances.getMatrix = getMatrix;
		getKmerDistances.getDb = [this]() { return partition.GetState().dbSeqs; };
		getKmerDistances.getTrainingSetIndices = [this]() { return partition.GetState().trainingSetIndices; };
		getKmerDistances.AddConsumer( &clusterTrainingSet );

		clusterTrainingSet.AddPropertyChangedEventHandler( this );
		clusterTrainingSet.getMatrix = getMatrix;
		clusterTrainingSet.getKmerLength = [this]() { return getKmerDistances.KmerLength(); };
		clusterTrainingSet.getThreshold = [this]() { return getKmerDistances.Threshold(); };
		clusterTrainingSet.setThreshold = [this]( int t ) { getKmerDistances.SetThreshold( t ); };
		clusterTrainingSet.getClassIndex = [this]() { return dbLoader.ClassIndex(); };
		clusterTrainingSet.getIdIndex = [this]() { return dbLoader.IdIndex(); };
		clusterTrainingSet.getAlphabet = [this]() { return dbLoader.Alphabet(); };
		clusterTrainingSet.getSeqIndex = [this]() { return dbLoader.SeqIndex(); };
		clusterTrainingSet.getDb = [this]() { return partition.GetState().dbSeqs; };
		clusterTrainingSet.getTrainingSetIndices = [this]() { return partition.GetState().trainingSetIndices; };
		clusterTrainingSet.AddConsumer( &selectVocab );

		selectVocab.AddPropertyChangedEventHandler( this );
		selectVocab.getPrototypes = [this]() { return clusterTrainingSet.Prototypes(); };
		selectVocab.getKmerLength = [this]() { return getKmerDistances.KmerLength(); };
		selectVocab.getSeqIndex = [this]() { return dbLoader.SeqIndex(); };
		selectVocab.AddConsumer( &encodePartition );

		encodePartition.AddPropertyChangedEventHandler( this );
		encodePartition.getMatrix = getMatrix;
		encodePartition.getVocab = [this]() { return selectVocab.Vocab(); };
		encodePartition.getDatabase = [this]() { return dbLoader.DbSeqs(); };
		encodePartition.getSeqIndex = [this]() { return dbLoader.SeqIndex(); };
		encodePartition.getKmerLength = [this]() { return getKmerDistances.KmerLength(); };
		encodePartition.getThreshold = [this]() { return getKmerDistances.Threshold(); };
		encodePartition.getClassIndex = [this]() { return dbLoader.ClassIndex(); };
		encodePartition.getNameIndex = [this]() { return dbLoader.NameIndex(); };
		encodePartition.getTrainingSetIndex = [this]() { return partition.GetState().trainingSetIndices; };
		encodePartition.AddConsumer( &rankTrainingSet );

		rankTrainingSet.AddPropertyChangedEventHandler( this );
		rankTrainingSet.getDbSeqs = [this]() { return dbLoader.DbSeqs(); };
		rankTrainingSet.getDbIndex = [this]() { return dbLoader.SeqIndex(); };
		rankTrainingSet.getSignatures = [this]() { return encodePartition.Signatures(); };
		rankTrainingSet.getTestSetIndex = [this]() { return partition.TestSetIndices(); };
		rankTrainingSet.getTrainingSetIndex = [this]() { return partition.TrainingSetIndices(); };
		rankTrainingSet.getPostingList = [this]() { return encodePartition.PostingList(); };
		rankTrainingSet.getClassPostingList = [this]() { return partition.ClassPostingList(); };
		rankTrainingSet.getSigIndex = [this]() { return encodePartition.SigIndex(); };
		rankTrainingSet.getKmerLength = [this]() { return getKmerDistances.KmerLength(); };
		rankTrainingSet.getMatrix = getMatrix;
		rankTrainingSet.getThreshold = [this]() { return getKmerDistances.Threshold(); };
		rankTrainingSet.getVocab = [this]() { return selectVocab.Vocab(); };
		// rankTrainingSet.AddConsumer(&???);
	}

	Page *CurrentPage() {
		return (Page *) (this->pages[currentPage]);
	}

	/**
	 * Activate or deactivate the Next and Prev buttons.
	 * -	Next becomes inactive if current page is not ok, or if it is the
	 * 		last page in the workflow.
	 * -	Prev becomes inactive if the current page is the first page in the
	 * 		workflow.
	 */
	void EnableButtons() {
		if ( currentPage == this->pages.size() - 1 ) {
			nextButton.label( "Next" );
			nextButton.deactivate();
		}
		else if ( !CurrentPage()->Ok() ) {
			nextButton.label( this->pages[currentPage + 1]->Name().c_str() );
			nextButton.deactivate();
		}
		else {
			nextButton.label( this->pages[currentPage + 1]->Name().c_str() );
			nextButton.activate();
		}

		if ( currentPage == 0 ) {
			prevButton.label( "Prev" );
			prevButton.deactivate();
		}
		else {
			prevButton.label( pages[currentPage - 1]->Name().c_str() );
			prevButton.activate();
		}

		if ( pages[currentPage]->Ready() ) {
			runButton.activate();
		}
		else {
			runButton.deactivate();
		}

		nextButton.redraw();
		prevButton.redraw();
		runButton.redraw();
	}

	bool Forward() {
		bool ret = false;
		if ( CurrentPage()->Ok() && currentPage < this->pages.size() - 1 ) {
			middle.remove( CurrentPage() );
			CurrentPage()->LoseFocus();
			currentPage++;
			middle.add( CurrentPage(), BorderLayout::Location::Centre );
			CurrentPage()->GainFocus();
			middle.redraw();
			ret = true;
		}
		EnableButtons();
		return ret;
	}

	void Back() {
		if ( currentPage > 0 ) {
			middle.remove( CurrentPage() );
			CurrentPage()->LoseFocus();
			currentPage--;
			middle.add( CurrentPage(), BorderLayout::Location::Centre );
			CurrentPage()->GainFocus();
			middle.redraw();
		}
		EnableButtons();
	}

	void Run() {
		CurrentPage()->Run();
		EnableButtons();
	}

	void PageRunComplete( Page * page ) {
		if ( CurrentPage() == page ) {
			if ( Forward() && autoRun ) {
				Run();
			}
			else if ( autoTerm ) {
				exit( EXIT_SUCCESS );
			}
		}
	}

	static int WindowClosing( int event ) {
		// (cerr << "event = " << event << "\n").flush();

		if ( event == FL_SHORTCUT && Fl::event_key() == FL_Escape ) {
			return 1; // ignore Escape
		}
		else if ( event == Fl_Event::FL_CLOSE ) {
			if ( fl_choice( "Are you really absolutely positively certain you want to quit?", "Yes", "No", 0 ) ) {
				exit( 0 );
			}
		}

		return 0;
	}

	string settingName;
	string title;

	string SettingName() {
		return settingName;
	}

	void SetSettingName( const string & value ) {
		settingName = value;
		title = "DnaClust " + value;
		label( title.c_str() );
	}

	void SaveSettings() {
		set<Page::Param> newParms;

		for ( auto page : pages ) {
			page->GetParams( newParms );
		}

		for ( const auto & p: parms ) {
			newParms.insert( p );
		}

		ofstream f( SettingName() );
		CsvWriter w( f );

		for ( const ICsvWriter &parm : newParms ) {
			w.Write( parm ).Ln();
		}
	}

	set<Page::Param> parms;

	void LoadSettings() {
		ifstream f( SettingName() );
		CsvReader r( f );

		parms.clear();

		while ( !r.IsEOF() ) {
			Page::Param p;
			r >> p;
			parms.insert( p );
		}

		for ( auto page : pages ) {
			page->SetParams( parms );
		}
	}

	void LoadConfig() {
		Fl_File_Chooser flChooser( ".", "*.cfg", Fl_File_Chooser::SINGLE, "Open config file..." );
		flChooser.value( SettingName().c_str() );
		flChooser.show();

		while ( flChooser.shown() ) {
			Fl::wait();
		}

		if ( flChooser.value() ) {
			SetSettingName( flChooser.value( 1 ) );
			LoadSettings();
		}
	}

	void SaveConfig() {
		Fl_File_Chooser flChooser( ".", "*.cfg", Fl_File_Chooser::CREATE, "Save config to file..." );
		flChooser.value( SettingName().c_str() );
		flChooser.show();

		while ( flChooser.shown() ) {
			Fl::wait();
		}

		if ( flChooser.value() ) {
			SetSettingName( flChooser.value( 1 ) );
			SaveSettings();
		}
	}
};

using namespace QutBio;
using namespace std;

mutex FragmentAggregationMode::m;
mutex DistanceType::m;
function<const FastaSequence *(const string & seqId)> SparseSignature::Lookup;

DnaClustGui * dlg = 0;

/**
 *	Signal handler for ctrl-c to ensure screen is cleaned up properly.
 */
extern "C" void ctrl_c_handler( int signal_code ) {
	exit( 1 );
}

extern "C" void SaveSettings() {
	// This will save settings.
	dlg->SaveSettings();
}

static void AutoRun( void * userData ) {
	if ( dlg->autoRun && dlg->CurrentPage()->Ready() ) {
		dlg->Run();
	}
}

int main( int argc, char *argv[] ) {
	atexit( SaveSettings );
	signal( SIGINT, ctrl_c_handler );
	signal( SIGABRT, ctrl_c_handler );
	dlg = new DnaClustGui( 1024, 768 );
	dlg->SetSettingName( argc > 1 ? argv[1] : "DnaClust.cfg" );
	dlg->LoadSettings();
	dlg->autoRun = argc > 2 && strcasecmp( argv[2], "autoRun" ) == 0;
	dlg->autoTerm = argc > 3 && strcasecmp( argv[3], "autoTerm" ) == 0;
	Fl::add_timeout( 0.25, AutoRun, dlg );
	Fl::lock();
	dlg->show();
	Fl::run();
}
