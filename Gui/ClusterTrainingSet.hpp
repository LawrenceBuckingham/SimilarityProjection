#pragma once

#include <FL/Enumerations.H>
#include <FL/Fl.H>
#include <FL/Fl_Chart.H>
#include <FL/Fl_Tabs.H>
#include <algorithm>
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
#include <Centroid.hpp>
#include <DiscreteDistribution.hpp>
#include "IncludeAll.hpp"
#include <Util.hpp>
#include "ClusterBrowser.hpp"

#include "Helpers.hpp"

#undef TRON
// #define TRON
#include <db.hpp>

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
class ClusterTrainingSet :
	public Page,
	// public virtual IEventHandler,
	public virtual IPropertyChangedEventHandler,
	public virtual ScatterPlot::MouseHandler
	//
{

	int rowHeight = 0;
	int rows;
	GridLayout stack;
	Fl_Box titleBox;
	TextDisplay disp;

	Arg<ulong, Fl_Input> seed;
	Arg<int, Fl_Input> numThreads;
	TypedArg<int, Fl_Check_Button> insertAll;
	TypedArg<string, OutFileChooser> clusterFileName;

	struct ButtonData {
		ClusterTrainingSet * context;
		Fl_Widget *targetPanel;

		ButtonData(
			ClusterTrainingSet * context,
			Fl_Widget *targetPanel
		) : context( context ), targetPanel( targetPanel ) {}

		virtual ~ButtonData() {}
	};

	FontSpec mono12;

	vector<Fl_Button *> buttonMutex;
	FlowLayout *buttons = 0;

	Fl_Button *showDisplay = 0;

	Fl_Button *showClusterSizes = 0;
	ScatterPlot *sizeChart = 0;
	Circle initialPopMarker;
	Series initialPop;
	Circle finalPopMarker;
	Series finalPop;

	Fl_Button *showPurity = 0;
	ScatterPlot *purityChart = 0;
	Square purityMarker;
	Series puritySeries;

	Fl_Button *showSizeVsPurity = 0;
	ScatterPlot *sizeVsPurityChart = 0;
	Series sizeVsPuritySeries;

	Fl_Button *showSizeDist = 0;
	ScatterPlot *sizeDistChart = 0;
	LineSpec sizeDistLine;
	Series sizeDistSeries;

	Fl_Button *showPurityDist = 0;
	ScatterPlot *purityDistChart = 0;
	LineSpec purityDistLine;
	Series purityDistSeries;

	/// Computed: prototypes
	vector<Centroid> prototypes;

	/// Computed: Assignment of allKmers to prototypes.
	vector<vector<shared_ptr<const SimpleKmer::Instance>>>  clusters;

	/// Computed: list of classes assigned to each prototype.
	vector<map<string, size_t>>	classCounters;

	ClusterBrowser clusterBrowser;
	Fl_Button showClusterBrowser;

	mutex mute;

	Alphabet * alphabet = 0;
	uint kmerLength = 0;
	SimilarityMatrix * matrix = 0;
	int threshold = 0;
	int classIndex = 0;
	int idIndex = 0;
	vector<pKmerInstance> allKmers;
	const LookupTable_<size_t, const FastaSequence> *seqIndex = 0;
	const vector<FastaSequence *> *dbSeqs = 0;
	const vector<size_t> *trainingSetIndices = 0;

public:
	/// Functional which will yield the k-mer length.
	function<size_t()> getKmerLength;

	/// Functional which will yield the threshold distance.
	function<int()> getThreshold;

	/// Functional which will allow threshold distance to be pushed back after load.
	function<void( int )> setThreshold;

	/** Functional which will yield a similarity matrix. */
	function<pSimilarityMatrix()> getMatrix;

	/** Functional which will yield the integer index of the class label in sequence metadata strings. */
	function<int()> getIdIndex;

	/** Functional which will yield the integer index of the class label in sequence metadata strings. */
	function<int()> getClassIndex;

	/// Get the alphabet;
	function<Alphabet *()> getAlphabet;

	/// Get the sequence index.
	function<const LookupTable_<size_t, const FastaSequence> *()> getSeqIndex;

	/** Functional which will yield a sequence dataset. */
	function<decltype(dbSeqs)()> getDb;

	/// Function which will yield a list of training set indices.
	function<decltype(trainingSetIndices)()> getTrainingSetIndices;

	ClusterTrainingSet(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Cluster Training Set" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		rows( 9 ),
		stack( 0, 0, 0, rows*rowHeight, rows, 1 ),
		titleBox( 0, 0, 0, rowHeight, name.c_str() ),
		disp( 0, 0, 0, 100 ),

		seed( "seed", time( nullptr ), "RNG Seed", Requirement::Required,
			"Cardinal number; required. Seed for random number generator used\n"
			"to sample k-mer population to compute empirical distances, and \n"
			"available for downstream components if necessary.",
			labelWidth, 0, width, rowHeight ),

		numThreads(
			"numThreads", 0, "OMP Threads", Requirement::Required,
			"Cardinal number; required. Number of threads for Open MP parallel sections.",
			labelWidth, 0, width, rowHeight ),

		clusterPercent(
			"clusterPercent", 95, "Cluster Percent", Requirement::Required,
			"Cardinal number; required. Percentage of training set to cover with greedy cluster\n"
			"algorithm. The last few percent take a long time, and yield little novelty, so it\n"
			"does little harm to skip them.",
			labelWidth, 0, width, rowHeight ),

		insertAll(
			"insertAll", 0, "Insert All Kmers", Requirement::Required,
			"Boolean. True if and only if you want to reinsert all allKmers to clusters\n"
			"(which is VERY SLOW). If this is not true, it also suppresses reloading \n"
			"of cluster members, which makes reuse of clusters practical even if they \n"
			"were fully populated.",
			labelWidth, 0, width, rowHeight ),

		clusterFileName(
			"clusterFileName", "*.protos.csv", "Cluster File", Requirement::Required,
			"The path to a file which will be overwritten with prototypes and k-mer cluster "
			"assignments, in CSV format.",
			labelWidth, 0, width, rowHeight ),

		autoLoad(
			"autoLoad", 0, "Auto Load?", Requirement::Required,
			"Check this box to automatically load the partition file if it exists.",
			labelWidth, 0, width, rowHeight ),

		mono12{ FL_SCREEN, 12 },

		initialPopMarker( nullptr, FL_BLUE, 5 ),
		initialPop( &initialPopMarker, nullptr, &mono12 ),

		finalPopMarker( nullptr, FL_GREEN, 3 ),
		finalPop( &finalPopMarker, nullptr, &mono12 ),

		purityMarker( nullptr, FL_GREEN, 3 ),
		puritySeries( &purityMarker, nullptr, &mono12 ),

		sizeVsPuritySeries( &purityMarker, nullptr, &mono12 ),

		sizeDistLine( FL_GREEN, 3, FL_SOLID ),
		sizeDistSeries( &purityMarker, &sizeDistLine, &mono12 ),

		purityDistLine( FL_GREEN, 3, FL_SOLID ),
		purityDistSeries( &purityMarker, &purityDistLine, &mono12 ),

		clusterBrowser( prototypes, clusters, [this]() { return kmerLength; } ),
		showClusterBrowser( 0, 0, 100, rowHeight, "Browse" )

		//
	{
		add( &stack, Location::North );

		stack.add( &titleBox );
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		stack.add( &seed );
		stack.add( &numThreads );
		stack.add( &clusterPercent );
		stack.add( &insertAll );
		stack.add( &clusterFileName );
		stack.add( &autoLoad );

		stack.add( buttons = new FlowLayout{ 0, 0, w(), rowHeight } );

		buttons->add( showDisplay = new Fl_Button{ 0, 0, 100, rowHeight, "Log" } );
		buttonMutex.push_back( showDisplay );
		add( &disp, Location::Centre );

		buttons->add( showClusterSizes = new Fl_Button{ 0, 0, 100, rowHeight, "Cluster sizes" } );
		buttonMutex.push_back( showClusterSizes );
		add( sizeChart = new ScatterPlot(), Location::Centre );
		InitChart( sizeChart, vector<Series *>{&initialPop, &finalPop} );

		buttons->add( showPurity = new Fl_Button{ 0, 0, 100, rowHeight, "Purity" } );
		buttonMutex.push_back( showPurity );
		add( purityChart = new ScatterPlot(), Location::Centre );
		InitChart( purityChart, vector<Series *>{&puritySeries} );

		buttons->add( showSizeVsPurity = new Fl_Button{ 0, 0, 100, rowHeight, "Purity v. size" } );
		buttonMutex.push_back( showSizeVsPurity );
		add( sizeVsPurityChart = new ScatterPlot(), Location::Centre );
		InitChart( sizeVsPurityChart, vector<Series *>{&sizeVsPuritySeries} );

		buttons->add( showSizeDist = new Fl_Button{ 0, 0, 100, rowHeight, "Size dist" } );
		buttonMutex.push_back( showSizeDist );
		add( sizeDistChart = new ScatterPlot(), Location::Centre );
		InitChart( sizeDistChart, vector<Series *>{&sizeDistSeries} );

		buttons->add( showPurityDist = new Fl_Button{ 0, 0, 100, rowHeight, "Purity dist" } );
		buttonMutex.push_back( showPurityDist );
		add( purityDistChart = new ScatterPlot(), Location::Centre );
		InitChart( purityDistChart, vector<Series *>{&purityDistSeries} );

		buttons->add( &showClusterBrowser );
		buttonMutex.push_back( &showClusterBrowser );
		add( &clusterBrowser, Location::Centre );

		vector<Fl_Widget *> panels{ &disp, sizeChart, purityChart, sizeVsPurityChart, sizeDistChart, purityDistChart, &clusterBrowser };

		for ( uint i = 0; i < buttonMutex.size(); i++ ) {
			buttonMutex[i]->callback( ShowOutputPanel, new ButtonData( this, panels[i] ) );
		}

		ShowOutputPanel( showDisplay );

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );

		seed.AddPropertyChangedEventHandler( this );
		numThreads.AddPropertyChangedEventHandler( this );
		clusterFileName.AddPropertyChangedEventHandler( this );
	}

	virtual ~ClusterTrainingSet() {
		for ( auto b : buttonMutex ) {
			delete (ButtonData *) b->user_data();
			if ( b != &showClusterBrowser ) delete b;
		}

		delete this->sizeChart;
		delete this->purityChart;
		delete this->sizeVsPurityChart;
		delete this->sizeDistChart;
		delete this->purityDistChart;
		delete this->buttons;
	}

	void InitChart( ScatterPlot * chart, vector<Series *> series ) {
		chart->SetXAxis( new LinearAxis{ 0, 2 * M_PI } );
		chart->SetYAxis( new LinearAxis{ -1, 1 } );
		chart->SetMargin( 75, 10, 30, 30 );
		chart->SetFillColour( FL_WHITE );
		chart->AddMouseHandler( this );

		for ( auto s : series ) {
			chart->Add( s );
		}
	}

	void ShowOutputPanel( Fl_Widget * control ) {
		Fl_Button * button = (Fl_Button *) control;

		for ( auto b : buttonMutex ) {
			auto bData = (ButtonData *) b->user_data();

			if ( b == button ) {
				b->deactivate();
				bData->targetPanel->show();
				bData->targetPanel->redraw();
			}
			else {
				b->activate();
				bData->targetPanel->hide();
			}
		}
	}

	static void ShowOutputPanel( Fl_Widget *sender, void *target_ ) {
		auto bData = (ButtonData *) target_;
		bData->context->ShowOutputPanel( sender );
	}

	void GainFocus() override {
		Page::GainFocus();
		kmerLength = getKmerLength();
		allKmers.clear();
		matrix = getMatrix();
		threshold = getThreshold();
		classIndex = getClassIndex();
		idIndex = getIdIndex();
		alphabet = getAlphabet();
		seqIndex = getSeqIndex();
		dbSeqs = getDb();
		trainingSetIndices = getTrainingSetIndices();
		UpdateReadiness();
	}

	void LoseFocus() override {
		Page::LoseFocus();
	}

	void UpdateReadiness() {
		bool isReady = seed.HasValue()
			&& numThreads.HasValue()
			&& String::Trim( clusterFileName.Value() ).size() > 0
			&& seqIndex != 0
			&& seqIndex->size() > 0
			;
		bool isOk = isReady
			&& prototypes.size() > 0;
		SetReady( isReady );
		SetOk( isOk );
	}

	/**
	 *	Generates theoretical and empirical distance distributions. Allows user
	 *	to specify a p-value for the clustering step.
	 */
	void Run() {
		double startTime = omp_get_wtime();

		auto confirm = CheckLoad();

		if ( confirm == YesNoCancel::Yes ) {
			SetReady( false );
			thread t( LoadImpl, this );
			t.detach();
		}
		else if ( confirm == YesNoCancel::No ) {
			if ( numThreads.Value() <= 0 || clusterFileName.Value() == "" ) {
				disp( "numThreads must be greater than zero; cluster file name must have value.\nAction abandoned.\n" );
				return;
			}
			SetReady( false );
			thread t( RunImpl, this );
			t.detach();
		}
		else {
			disp( "Action cancelled.\n" );
		}
	}

	static void RunImpl( ClusterTrainingSet * self ) {
		self->RunImpl_();
	}

	void RunImpl_() {
		if ( !mute.try_lock() ) {
			return;
		}

		UPDATE_GUI( SetOk( false ); disp( "Computing new clusters...\n" ); );

		Cluster();

		UPDATE_GUI( UpdateReadiness(); DrawChart(); NotifyRunComplete(); );
		mute.unlock();
	}

	static void LoadImpl( ClusterTrainingSet * self ) {
		self->LoadImpl_();
	}

	void LoadImpl_() {
		if ( !mute.try_lock() ) {
			return;
		}

		UPDATE_GUI( SetOk( false ); disp( "Computing new clusters...\n" ) );

		Page::loadTime = -omp_get_wtime();
		Load();
		Page::loadTime += omp_get_wtime();

		UPDATE_GUI( UpdateReadiness(); DrawChart(); NotifyRunComplete(); );

		mute.unlock();
	}

	YesNoCancel CheckLoad() {
		if ( String::Trim( clusterFileName.Value() ).length() > 0 && File::Exists( clusterFileName.Value() ) ) {
			if ( AutoLoad() ) {
				return YesNoCancel::Yes;
			}

			int answer = fl_choice(
				"The selected file already exists.\n"
				"Would you like to load those clusters?",
				"Cancel operation.",
				"Yes, load the clusters.",
				"No, I want to redo the clusters."
			);

			return (YesNoCancel) answer;
		}
		else {
			return YesNoCancel::No;
		}
	}

	void Cluster() {
		double startTime = omp_get_wtime();
		Page::runTime = -omp_get_wtime();
		// fprintf( stderr, "Run called at %f.\n", startTime );

		UPDATE_GUI( disp.clear(); disp( "Run called at %f.\n", startTime ); );

		prototypes.clear();

		size_t N = allKmers.size();

		if ( N == 0 ) {
			for ( auto i : *trainingSetIndices ) {
				N += (*dbSeqs)[i]->KmerCount( kmerLength );
			}

			allKmers.reserve( N );

			for ( auto i : *trainingSetIndices ) {
				auto seq = (*dbSeqs)[i];
				auto n = seq->KmerCount( kmerLength );

				for ( size_t j = 0; j < n; j++ ) {
					allKmers.push_back( pKmerInstance( new SimpleKmer::Instance( seq, j ) ) );
				}
			}
		}

		vector<pKmerInstance> unalloc = allKmers;

		int clusterPercent = ClusterPercent();

		// Shuffle the kmers
		UniformIntRandom<size_t> rand( seed.Value() );

		for ( size_t i = 0; i < N - 1; i++ ) {
			size_t idx = rand( i, N - 1 );
			std::swap( unalloc[i], unalloc[idx] );
		}

		struct IndexInterval {
			size_t begin = 0;
			size_t end = 0;
			size_t next = 0;
		};

		vector<IndexInterval> indexes;
		size_t previousCount = N;
		size_t terminateWhen = N - N * clusterPercent / 100;

#if USE_OMP
		int numThreads = this->numThreads.Value();
		int hardwareThreads = std::thread::hardware_concurrency();

		if ( numThreads < 1 || numThreads > hardwareThreads ) numThreads = hardwareThreads;

		if ( numThreads < 1 ) numThreads = 1;

		omp_set_num_threads( numThreads );
		indexes.resize( numThreads );

		for ( int i = 0; i < numThreads; i++ ) {
			indexes[i].begin = indexes[i].next = i * N / numThreads;
			indexes[i].end = (i + 1) * N / numThreads;
		}
#else
		int numThreads = 1;

		indexes.resize( 1 );
		indexes[0].begin = indexes[0].next = 0;
		indexes[0].end = N;
#endif

		for ( ;;) {
			pKmerInstance centre( 0 );
			size_t numRemaining = 0;

			for ( int i = 0; i < numThreads; i++ ) {
				numRemaining += indexes[i].end - indexes[i].next;
			}

			if ( numRemaining < terminateWhen ) break;

			if ( prototypes.size() > 0 ) {
				prototypes.back().initialClusterSize = previousCount - numRemaining;

#if WANT_FEEDBACK
				{
					int thisPercent = numRemaining * 100 / N, lastPercent = previousCount * 100 / N;

					if ( thisPercent != lastPercent ) {
						UPDATE_GUI( disp( "Cluster training set: %d%% remaining\n", thisPercent ); );
					}
				}
#endif

				previousCount = numRemaining;

				// fprintf( stderr, "Cluster %zu: %zu elements. %zu unallocated k-mers\n", prototypes.size(), prototypes.back().initialClusterSize, numRemaining );
				// fflush( stderr );
			}

			// Get a new prototype from the available pool.
			for ( int i = 0; i < numThreads; i++ ) {
				if ( indexes[i].next < indexes[i].end ) {
					centre = unalloc[indexes[i].next];
					prototypes.emplace_back( centre, 0, 0 );
					indexes[i].next++;
					break;
				}
			}

			if ( !centre ) break;

			// const unsigned char * centreBytes = centre->Bytes();
			const Digram * centreDigrams = centre->Digrams();

#if USE_OMP
#pragma omp parallel
#endif
			{
#if USE_OMP
				int threadId = omp_get_thread_num();
#else
				int threadId = 0;
#endif

				for ( size_t i = indexes[threadId].next; i < indexes[threadId].end; i++ ) {
					// Distance dist = matrix->Difference( centreBytes, unalloc[i]->Bytes(), kmerLength );
					Distance dist = matrix->DigramDifference( centreDigrams, unalloc[i]->Digrams(), kmerLength / 2 );

					if ( dist <= threshold ) {
						std::swap( unalloc[i], unalloc[indexes[threadId].next] );
						indexes[threadId].next++;
					}
				}
			}
		}

		if ( prototypes.size() > 0 ) {
			prototypes.back().initialClusterSize = previousCount;
		}

		std::sort( prototypes.begin(), prototypes.end(), Centroid::GreaterByInitSize );

		Page::runTime += omp_get_wtime();

		// fprintf( stderr, "\nPrototype selection finished in %fs.\n\n", (omp_get_wtime() - startTime) );
		ShowProgress( "Prototype selection finished in %fs.\n%zu centroids identified.\n", (omp_get_wtime() - startTime), prototypes.size() );

		{ UPDATE_GUI( DrawClusterSizeChart() ) }

		Page::runTime -= omp_get_wtime();

		clusters.resize( prototypes.size() );
		classCounters.resize( prototypes.size() );

		int classIndex = getClassIndex();

		if ( insertAll.Value() ) {
			size_t N = prototypes.size();

#if WANT_FEEDBACK
			size_t done = 0;
#endif

#if USE_OMP
#pragma omp parallel for
#endif
			for ( size_t i = 0; i < N; i++ ) {
				clusters[i].clear();
				classCounters[i].clear();

				Centroid & centre = prototypes[i];
				const Symbol * centreBytes = centre.centroid->Bytes();

				for ( auto kmer : allKmers ) {
					Distance dist = matrix->Difference( centreBytes, kmer->Bytes(), kmerLength );

					if ( dist <= threshold ) {
						clusters[i].push_back( kmer );
						centre.finalClusterSize++;
						centre.finalInstanceCount++;

						auto seq = kmer->sequence;
						const auto & classLabel = seq->Metadata( classIndex );
						classCounters[i][classLabel] ++;
					}
				}

				auto begin = classCounters[i].begin();
				auto end = classCounters[i].end();
				using tuple = pair<string, size_t>;

				auto pos = std::max_element( begin, end, []( tuple lhs, tuple rhs ) { return lhs.second < rhs.second; } );
				auto sum = std::accumulate( begin, end, size_t( 0 ), []( size_t lhs, tuple rhs ) { return lhs + rhs.second;  } );

				prototypes[i].purity = sum == 0 ? 0.0 : double( pos->second ) / sum;

#if WANT_FEEDBACK
#if USE_OMP
#pragma omp critical
#endif
				{
					done++;
					int thisPercent = done * 100 / N, lastPercent = (done - 1) * 100 / N;

					if ( thisPercent != lastPercent ) {
						UPDATE_GUI( disp( "Reinserting k-mers: %d%%\n", thisPercent ); DrawChart(); )
					}
				}
#endif

				std::sort( prototypes.begin(), prototypes.end(), Centroid::GreaterByFinalSize );
			}

			// fprintf( stderr, "\nReinsertion finished after %fs.\n\n", (omp_get_wtime() - startTime) );
			ShowProgress( "Reinsertion finished after %fs.\nSaving clusters...\n", (omp_get_wtime() - startTime) );
		}
		else {

#if USE_OMP
#pragma omp parallel for
#endif
			// Get rid of spurious stuff hanging around from previous run.
			for ( size_t i = 0; i < prototypes.size(); i++ ) {
				clusters[i].clear();
				classCounters[i].clear();
			}
		}

		Page::runTime += omp_get_wtime();

		Helpers::UpdateGui( [this]() {
			disp( "Saving clusters...\n" );
			SetOk( prototypes.size() > 0 && clusterFileName.Value() != "" );
		} );

		Page::saveTime = -omp_get_wtime();
		Save();
		Page::saveTime += omp_get_wtime();

		ShowProgress( "Clusters saved after %fs.\nDrawing charts...\n", (omp_get_wtime() - startTime) );
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

	void DrawChart() {
		DrawClusterSizeChart();
		DrawPurityChart();
		DrawSizeVsPurityChart();
		DrawSizeDistChart();
		DrawPurityDistChart();
		if ( clusterBrowser.Focus() == 0 ) clusterBrowser.Refresh(); else clusterBrowser.SetFocus( 0 );
	}

	void DrawClusterSizeChart() {
		initialPop.Clear();
		finalPop.Clear();

		for ( size_t i = 0; i < prototypes.size(); i++ ) {
			initialPop.Add( i, prototypes[i].initialClusterSize );
			finalPop.Add( i, prototypes[i].finalClusterSize );
		}

		SetAxes( sizeChart );
		sizeChart->redraw();
	}

	void DrawPurityChart() {
		puritySeries.Clear();

		for ( size_t i = 0; i < prototypes.size(); i++ ) {
			// cerr << "Added (" << i << "," << prototypes[i].purity << ") to purity series\n";
			puritySeries.Add( i, prototypes[i].purity );
		}

		SetAxes( purityChart );
		purityChart->redraw();
	}

	void DrawSizeVsPurityChart() {
		sizeVsPuritySeries.Clear();

		for ( size_t i = 0; i < prototypes.size(); i++ ) {
			// cerr << "Added (" << i << "," << prototypes[i].purity << ") to purity series\n";
			sizeVsPuritySeries.Add( prototypes[i].finalClusterSize, prototypes[i].purity );
		}

		SetAxes( sizeVsPurityChart );
		sizeVsPurityChart->redraw();
	}

	void DrawSizeDistChart() {
		Histogram<size_t> hist;

		for ( auto &c : prototypes ) {
			hist.Add( c.finalClusterSize );
		}

		IntegerDistribution dist( hist );
		sizeDistSeries.Clear();

		const int steps = 100;

		for ( size_t i = 0; i <= steps; i++ ) {
			double x = dist.Min() + i * (dist.Max() - dist.Min()) / steps;
			double y = dist.Cdf( x );
			sizeDistSeries.Add( x, y );
		}

		SetAxes( sizeDistChart );
		sizeDistChart->redraw();
	}

	void DrawPurityDistChart() {
		Histogram<double> hist;

		for ( auto &c : prototypes ) {
			hist.Add( c.purity );
		}

		// TODO: why is DiscreteDistribution called that? It seems to be a better name might be Empirical Distribution?
		DiscreteDistribution dist;
		dist.SetPmf( hist );
		purityDistSeries.Clear();

		auto keys = hist.GetKeys();
		double min, max;
		GetMin( keys, min );
		GetMax( keys, max );

		const int steps = 100;

		for ( size_t i = 0; i <= steps; i++ ) {
			double x = min + i * (max - min) / steps;
			double y = dist.Cdf( x );
			purityDistSeries.Add( x, y );
		}

		SetAxes( purityDistChart );
		purityDistChart->redraw();
	}

	static void SetAxes( ScatterPlot * chart ) {
		double xCrossesY = 0;

		Series *htics = chart->HTics();
		LineSpec tickLine;
		Plus tickMark( &tickLine, 5 );

		if ( htics->GetMarker() == 0 ) {
			htics->SetMarker( &tickMark );
		}
		htics->Clear();

		Series *vtics = chart->VTics();

		if ( vtics->GetMarker() == 0 ) {
			vtics->SetMarker( &tickMark );
		}

		vtics->Clear();

		double xMax = -1;
		double popMax = -1;

		for ( auto series : chart->DataSeries() ) {
			for ( auto p : series->Data() ) {
				if ( p->x > xMax ) xMax = p->x;
				if ( p->y > popMax ) popMax = p->y;
			}
		}

		double hMax = xMax; // (int)(pow(10, ceil(log10(xMax))));
		int numHTicks = 10; // std::min( 10, hMax );

		chart->XAxis()->SetMin( 0 );
		chart->XAxis()->SetMax( hMax );

		static FontSpec mono12{ FL_SCREEN, 12 };

		for ( int i = 0; i <= numHTicks; i++ ) {
			double x = i * hMax / numHTicks;
			Point *p = htics->Add( x, xCrossesY );

			p->SetLabel( Double::ToString( x ), 0, 7, 0.5, 0.0, &mono12 );
		}

		double yMax = popMax; // (int)(pow(10, ceil(log10(popMax))));
		int numYTicks = 10;

		chart->YAxis()->SetMax( 0 );
		chart->YAxis()->SetMin( yMax );

		for ( int i = 0; i <= numYTicks; i++ ) {
			double y = i * yMax / numYTicks;
			Point *p = vtics->Add( 0, y );
			p->SetLabel( Double::ToString( y ), -7, 0, 1.0, 0.5, &mono12 );
		}
	}

	bool Handle( const ScatterPlot::MouseEvent & me ) override {
		if ( me.eventCode == FL_RELEASE ) {
			auto datum = me.source->NearestTo( me.x, me.y );

			if ( me.source == sizeChart || me.source == sizeVsPurityChart || me.source == purityChart ) {
				size_t idx = datum.index;
#if SHOW_INTERMEDIATE_BROWSE_DIALOG
				auto & proto = (idx < prototypes.size()) ? prototypes[idx] : Centroid::Zero();
				ostringstream s;
				s << "Cluster around "
					<< proto.centroid->Sequence()->IdStr()
					<< "(" << proto.centroid->KmerPosition() << "):"
					<< "\nPattern: " << proto.centroid->Chars()
					<< "\nInitial cluster size: " << proto.initialClusterSize
					<< "\nFinal cluster size: " << proto.finalClusterSize
					<< "\nFinal instance count:" << proto.finalInstanceCount
					<< "\nClass purity: " << proto.purity
					<< "\nClass entropy: " << proto.entropy
					;
				int answer = fl_choice(
					s.str().c_str(),
					"Close",
					"Browse cluster",
					0
				);

				if ( answer == YesNoCancel::Yes ) {
					clusterBrowser.SetFocus( idx );
					ShowOutputPanel( &showClusterBrowser );
				}
#else
				clusterBrowser.SetFocus( idx );
				ShowOutputPanel( &showClusterBrowser );
				return true;
#endif
			}
			else {
				fl_message( "Point: %s:%d (%f, %f)", datum.series->Name().c_str(), datum.index, me.x, me.y );
				return true;
			}
		}

		return false;
	}

	void Save() {
		if ( !Ok() ) return;

		ofstream f( clusterFileName.Value() );
		CsvWriter w( f );

		w << Name() << '\n';
		w << "kmerLength" << kmerLength << '\n';
		w << "threshold" << threshold << '\n';
		w << "matrix" << *matrix << '\n';
		w << "classIndex" << classIndex << '\n';
		w << "prototypes" << prototypes.size() << '\n';

		for ( size_t i = 0; i < prototypes.size(); i++ ) {
			w << i << prototypes[i] << '\n';
		}

		if ( this->insertAll.Value() ) {
			w << "clusters" << '\n';

			for ( size_t i = 0; i < prototypes.size(); i++ ) {
				for ( auto & kmer : clusters[i] ) {
					w << i << *kmer << '\n';
				}
			}

			int sentinel = -1;
			w << sentinel << '\n'; // Sentinel for cluster list.
		}
	}

	void Expect( CsvReader & r, const string & val ) {
		string s;
		r >> s;
		if ( s != val ) {
			string message = "Expected text '" + val + "' not found.";
			throw Exception( message, FileAndLine );
		}
	}

	template<typename T>
	void Expect( CsvReader & r, const string & name, const T & expectedVal ) {
		Expect( r, name );
		T val;
		r >> val;

		if ( val != expectedVal ) {
			ostringstream ss;
			ss << "Value of " << name << " does not match current selection.";
			throw Exception( ss.str(), FileAndLine );
		}
	}

	template<typename T>
	void Expect( CsvReader & r, const string & name, const T & expectedVal, T & val ) {
		Expect( r, name );
		r >> val;

		if ( val != expectedVal ) {
			ostringstream ss;
			ss << "Value of " << name << " does not match current selection.";
			throw Exception( ss.str(), FileAndLine );
		}
	}

	bool Load() {
		auto fName = clusterFileName.Value();
		UPDATE_GUI( disp( "Attempting to load clusters from %s\n", fName.c_str() ) );
		prototypes.clear();
		clusters.clear();

		ifstream f( fName );
		CsvReader r( f );
		string buffer;
		TRACE;

		try {
			size_t prototypeCount;
			SimilarityMatrix matrix( alphabet );

			Expect( r, Name() );
			Expect( r, "kmerLength", kmerLength );
			Expect( r, "threshold" );
			r >> threshold;
			setThreshold( threshold );
			Expect( r, "matrix", *(this->matrix), matrix );

			r >> buffer;

			if ( buffer == "classIndex" ) {
				// Just skip the value; it is not required. 
				r >> buffer;
			}

			Expect( r, "prototypes" );
			r >> prototypeCount;
			TRACE;

			string kmerDef;

			while ( !r.IsEOF() && prototypes.size() < prototypeCount ) {
				size_t idx;
				r >> idx;

				if ( idx != prototypes.size() ) break;

				Centroid c;
				Helpers::GetCentroid( r, kmerDef, kmerLength, *seqIndex, c );
				prototypes.push_back( c );
			}

			size_t N = prototypes.size();
			clusters.resize( N );

			UPDATE_GUI( disp( "%zu prototypes loaded from %s.\n", prototypes.size(), clusterFileName.Value().c_str() ); DrawChart() );
			TRACE;

			if ( insertAll.Value() ) {
				// Expect( r, "clusters" );
				string s;

				while ( !r.IsEOF() ) {
					r >> s;
					if ( s == "clusters" ) break;
				}
				TRACE;

				if ( r.IsEOF() ) {
					UPDATE_GUI( disp( "End of file reached before clusters!", FileAndLine ); UpdateReadiness(); );
					return false;
				}
				TRACE;

				string kmerSeqId;
				size_t kmerOffset;
				int64_t prevPercent = 0;

				try {
					while ( !r.IsEOF() ) {
						int64_t idx;
						r >> idx;

						// PR( idx, %zd );

						if ( idx < 0 ) break;

						if ( idx >= (int64_t) clusters.size() ) {
							throw Exception( "Cluster index out of range.", FileAndLine );
						}

						r >> kmerSeqId >> kmerOffset;

						shared_ptr<const SimpleKmer::Instance> element = Helpers::GetKmer( kmerSeqId, kmerOffset, kmerLength, *seqIndex );
						clusters[idx].push_back( element );

						int64_t percent = (idx + 1) * 100 / N;

						if ( percent != prevPercent ) {
							prevPercent = percent;
							UPDATE_GUI( DrawChart(); disp( "Loading k-mer cluster assignments: %zu%%\n", percent ) )
						}
					}
				}
				catch ( Exception & ex ) {
					UPDATE_GUI( DrawChart(); UpdateReadiness(); disp( "Clusters loaded OK!\n" ) );
				}
			}

			std::sort( prototypes.begin(), prototypes.end(), Centroid::GreaterByFinalSize );

			UPDATE_GUI( DrawChart(); UpdateReadiness(); disp( "Clusters loaded OK!\n" ) );
			TRACE;
			return true;
		}
		catch ( Exception & e ) {
			UPDATE_GUI( (fl_message( "Configuration mismatch while loading clusters:\n%s", e.what() ),
				disp( "%s\nLoad abandoned.\n", e.what() )) );
			return false;
		}
	}

	void PropertyChanged( void *sender, const string &propertyName ) override {
		UpdateReadiness();
	}

	/// Get parameters so they can be serialised when program finishes.
	/// @param parms A Param list to which the component should append its current parameters. 
	void GetParams( set<Param> &parms ) override {
#define addPar(x) parms.emplace( Name(), STRING(x), (x).Value() )
#define addVar(x,fmt) parms.emplace( Name(), STRING(x), String::Format(fmt,x) )

		addPar( seed );
		addPar( numThreads );
		addPar( clusterFileName );
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
#define setStr(x) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value() )
#define setPar(x,t) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value<t>() )
#define setVar(x,t) if ( p.ParamName() == STRING(x) ) (x) = p.Value<t>()

		for ( auto & p : parms ) {
			if ( p.ComponentName() != Name() ) continue;
			setPar( seed, ulong );
			setPar( numThreads, int );
			setPar( autoLoad, int );
			setStr( clusterFileName );
			setVar( runTime, double );
			setVar( loadTime, double );
			setVar( saveTime, double );
		}

		UpdateReadiness();

#undef setPar
#undef setStr
#undef setVar
	}

#pragma region Prototypes
public:
	/// Get the address of the list of prototypes generated by this control.
	vector<Centroid> * Prototypes() {
		return &prototypes;
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


#pragma region ClusterPercent Property
private:
	Arg<uint, Fl_Int_Input> clusterPercent;

public:
	uint ClusterPercent() {
		uint ret = clusterPercent.Value();

		if ( ret > 100 ) ret = 95;

		return ret;
	}

	void SetClusterPercent( const uint val ) {
		auto val_ = std::min( val, 100u );

		if ( val_ == ClusterPercent() ) return;

		this->clusterPercent.SetValue( val_ );
	}
#pragma endregion

	void Reset() override {
		clusters.clear();
		prototypes.clear();
		allKmers.clear();
		kmerLength = 0;
		matrix = 0;
		threshold = 0;
		classIndex = 0;
		idIndex = 0;
		alphabet = 0;
		seqIndex = 0;
	}
};
