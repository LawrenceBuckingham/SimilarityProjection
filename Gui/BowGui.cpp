#include "IncludeAll.hpp"

#include <LBFL/DoubleWindow.hpp>

#include "DbLoader.hpp"
#include "Partition.hpp"
#include "EncodeBagOfWords.hpp"
#include "RankBagOfWords.hpp"

#include <IPropertyChangedEventHandler.hpp>

#undef TRON
// #define TRON 1
#include <db.hpp>

struct BowGui :
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
	EncodeBagOfWords encodeBagOfWords;
	RankBagOfWords rankBagOfWords;

	bool autoRun = false;
	bool autoTerm = false;

	BowGui( int W, int H, const char *l = 0 )
		: DoubleWindow( W, H, l ),
		nextButton( 0, 0, BW, rowHeight, "Next", [this]( Button *b ) { Forward(); } ),
		prevButton( 0, 0, BW, rowHeight, "Prev", [this]( Button *b ) { Back(); } ),
		runButton( 0, 0, BW, rowHeight, "Run", [this]( Button *b ) { Run(); } ),
		saveButton( 0, 0, BW, rowHeight, "Save Config", [this]( Button *b ) { SaveConfig(); } ),
		loadButton( 0, 0, BW, rowHeight, "Load Config", [this]( Button *b ) { LoadConfig(); } ),
		middle( 0, 0, W, H - rowHeight ),
		bottom( 0, H - rowHeight, W, rowHeight ),
		dbLoader( 0, 0, W, H - rowHeight ),
		partition( 0, 0, W, H - rowHeight ),
		encodeBagOfWords( 0, 0, W, H - rowHeight ),
		rankBagOfWords( 0, 0, W, H - rowHeight )
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
		pages.push_back( &encodeBagOfWords );
		pages.push_back( &rankBagOfWords );

		function<void( Page * )> pageRunComplete = [this]( Page * page ) { PageRunComplete( page ); };

		for ( auto page : pages ) {
			page->RunComplete( pageRunComplete );
		}

		WireUp();

		Fl::add_handler( WindowClosing );
		LoadSettings();
		EnableButtons();
	}

	virtual ~BowGui() {
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

		auto getMatrix = [this]() { return dbLoader.Matrix(); };
		auto getDatabase = [this]() { return dbLoader.DbSeqs(); };
		auto getSeqIndex = [this]() { return dbLoader.SeqIndex(); };
		auto getTrainingSetIndex = [this]() { return partition.GetState().trainingSetIndices; };
		auto getFeatureVectors = [this]() { return encodeBagOfWords.FeatureVectors(); };
		auto getTestSetIndex = [this]() { return partition.TestSetIndices(); };
		auto getPostingList = [this]() { return encodeBagOfWords.PostingList(); };
		auto getClassPostingList = [this]() { return partition.ClassPostingList(); };

		dbLoader.AddPropertyChangedEventHandler( this );
		dbLoader.AddConsumer( &partition );

		partition.AddPropertyChangedEventHandler( this );
		partition.GetDatabase = getDatabase;
		partition.AddConsumer( &encodeBagOfWords );

		encodeBagOfWords.AddPropertyChangedEventHandler( this );
		encodeBagOfWords.getDatabase = getDatabase;
		encodeBagOfWords.getSeqIndex = getSeqIndex;
		encodeBagOfWords.getTrainingSetIndex = getTrainingSetIndex;
		encodeBagOfWords.AddConsumer( &rankBagOfWords );

		rankBagOfWords.AddPropertyChangedEventHandler( this );
		rankBagOfWords.getDbSeqs = getDatabase;
		rankBagOfWords.getDbIndex = getSeqIndex;
		rankBagOfWords.getFeatureVectors = getFeatureVectors;
		rankBagOfWords.getTestSetIndex = getTestSetIndex;
		rankBagOfWords.getTrainingSetIndex = getTrainingSetIndex;
		rankBagOfWords.getPostingList = getPostingList;
		rankBagOfWords.getClassPostingList = getClassPostingList;
			/// rankBagOfWords.AddConsumer(&???);

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
	set<Page::Param> parms;

	string SettingName() {
		return settingName;
	}

	void SetSettingName( const string & value ) {
		settingName = value;
		title = "DnaClust " + value;
		label( title.c_str() );
	}

	void SaveSettings() {
		set<Page::Param> parms;

		for ( auto page : pages ) {
			page->GetParams( parms );
		}

		parms.insert( this->parms.begin(), this->parms.end() );

		ofstream f( SettingName() );
		CsvWriter w( f );

		for ( const ICsvWriter &parm : parms ) {
			w.Write( parm ).Ln();
		}
	}

	void LoadSettings() {
		ifstream f( SettingName() );
		CsvReader r( f );

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

BowGui * dlg = 0;

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
	dlg = new BowGui( 1024, 768 );
	dlg->SetSettingName( argc > 1 ? argv[1] : "DnaClust.cfg" );
	dlg->LoadSettings();
	dlg->autoRun = argc > 2 && strcasecmp( argv[2], "autoRun" ) == 0;
	dlg->autoTerm = argc > 3 && strcasecmp( argv[3], "autoTerm" ) == 0;
	Fl::add_timeout( 0.25, AutoRun, dlg );
	Fl::lock();
	dlg->show();
	Fl::run();
}
