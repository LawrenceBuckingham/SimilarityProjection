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
#include "Helpers.hpp"


#undef TRON
// #define TRON
#include <db.hpp>

using namespace LBFL;
using namespace QutBio;
using namespace LBGraph;
using namespace std;

#include "DefaultChartHandler.hpp"

/**
<summary>
	DbLoader gathers the parameters for an interactive version of DnaClust,
	with some additional analytic displays to render the software a little more
	usable.
</summary>
 */
class GetKmerDistances : public virtual Page,
	public virtual IEventHandler,
	public virtual IPropertyChangedEventHandler,
	public virtual ScatterPlot::MouseHandler {
public:

#pragma region NumThreads Property
private:
	Arg<int, Fl_Int_Input> numThreads;

public:
	int NumThreads() {
		return numThreads.Value();
	}

	void SetNumThreads( const int val ) {
		if ( val == NumThreads() ) return;

		this->numThreads.SetValue( val );
	}
#pragma endregion

#pragma region
private:
	const int rows = 10;
	int rowHeight = 0;
	GridLayout stack;
	Fl_Box titleBox;
	Arg<uint, Fl_Input> kmerLength;
	Arg<ulong, Fl_Input> seed;
	Arg<size_t, Fl_Input> sampleSize;
	Arg<int, Fl_Input> threshold;
	Arg<double, Fl_Input> pValue;
	TypedArg<int, Fl_Check_Button> logScale;
	Arg<double, Fl_Input> logLower;
	TextDisplay disp;

	vector<Fl_Button *> buttonMutex;
	FlowLayout buttons;
	Fl_Button showDisplay;
	Fl_Button showHistograms;

	ScatterPlot chart;

	UniformRealRandom rand;
	Histogram<uint8_t> symbolDistribution;
	IntegerDistribution kmerDistanceTheoretical;
	IntegerDistribution kmerDistanceObserved;
	NormalDistribution normal;

	FontSpec mono12;

	LineSpec empiricalMarkerLine;
	LineSpec theoreticalMarkerLine;
	LineSpec normalMarkerLine;

	LineSpec empiricalCdfLine;
	Plus empiricalCdfMarker;
	Series empiricalCdf;

	LineSpec empiricalPdfLine;
	Plus empiricalPdfMarker;
	Series empiricalPdf;

	LineSpec theoreticalCdfLine;
	Cross theoreticalCdfMarker;
	Series theoreticalCdf;

	LineSpec theoreticalPdfLine;
	Cross theoreticalPdfMarker;
	Series theoreticalPdf;

	LineSpec normalCdfLine;
	Square normalCdfMarker;
	Series normalCdf;

	LineSpec normalPdfLine;
	Square normalPdfMarker;
	Series normalPdf;

	LineSpec crossHairLine;
	Marker crossHairMarker;
	Series crossHairSeriesH, crossHairSeriesV;

	mutex mute;

	pAlphabet alphabet = 0;
	pSimilarityMatrix matrix = 0;
	const vector<FastaSequence *> *dbSeqs = 0;
	const vector<size_t> *trainingSetIndices = 0;

	DefaultChartHandler handler;

#pragma endregion  // private data members

public:
	/** Functional which will yield an alphabet. */
	function<pAlphabet()> getAlphabet;

	/** Functional which will yield a similarity matrix. */
	function<pSimilarityMatrix()> getMatrix;

	/** Functional which will yield a sequence dataset. */
	function<const vector<FastaSequence *> *()> getDb;

	/// Function which will yield a list of training set indices.
	function<const vector<size_t> *()> getTrainingSetIndices;

	GetKmerDistances(
		int left,
		int top,
		int width,
		int height,
		int labelWidth = 125,
		int rowHeight = 25,
		int vgap = 5
	) :
		Page( left, top, width, height, "Kmer Distances" ),
		PropertyChangedEventSource( this ),
		rowHeight( rowHeight ),
		stack( 0, 0, w(), rows * rowHeight, rows, 1 ),
		titleBox( 0, 0, width, rowHeight, name.c_str() ),
		kmerLength(
			"kmerLength", 120, "k-mer Length", Requirement::Required,
			"Cardinal number; required. The number of symbols per k-letter word.",
			labelWidth, 0, width, rowHeight ),
		seed(
			"seed", time( nullptr ), "RNG Seed", Requirement::Required,
			"Cardinal number; required. Seed for random number generator used\n"
			"to sample k-mer population to compute empirical distances, and \n"
			"available for downstream components if necessary.",
			labelWidth, 0, width, rowHeight ),
		numThreads(
			"numThreads", 8, "Open MP Threads", Requirement::Required,
			"Number of Open MP threads to use.",
			labelWidth, 0, width, rowHeight ),
		sampleSize(
			"sampleSize", 10000, "Sample Size",
			Requirement::Required,
			"Cardinal number; required. Number of kmer pairs used to sample \n"
			"k-mer population to compute empirical distances, and available \n"
			"for downstream components if necessary.",
			labelWidth, 0, width, rowHeight ),
		threshold( "threshold", -1, "Cluster radius", Requirement::Optional,
			"Integer; required. Maximum distance from any k-mer in a cluster to\n"
			"the centroid of that cluster. Due to the peculiarity of biological\n"
			"sequence comparison, this may conceivably be negative, but usually\n"
			"it should be positive, and chosen so that clusters are tight enough\n"
			"to be unlikely to be the outcome of chance alone, yet not so tight\n"
			"that all clusters are singletons.",
			labelWidth, 0, width, rowHeight ),
		pValue( "seed", 0, "Threshold P-Value", Requirement::Required,
			"Probability that distance between random k-mers drawn from background\n"
			"distribution will be nearer than the threshold distance.",
			labelWidth, 0, width, rowHeight ),
		logScale( "logScale", 0, "Log Scale?", Requirement::Required,
			"Display probabilities on logarithmic scale.",
			labelWidth, 0, width, rowHeight ),
		logLower( "logLower", 1e-9, "Log low bound", Requirement::Required,
			"Lower limit for logarithmic scale.",
			labelWidth, 0, width, rowHeight ),
		buttons( 0, 0, w(), rowHeight ),
		showDisplay( 0, 0, 100, rowHeight, "Show Log" ),
		showHistograms( 0, 0, 100, rowHeight, "Show Charts" ),
		disp( 0, 0, w(), 100 ),
		rand( (unsigned long) time( 0 ) ),
		kmerDistanceTheoretical( 0, 0 ),
		kmerDistanceObserved( 0, 0 ),

		mono12( FL_SCREEN, 12 ),

		empiricalMarkerLine( FL_DARK_BLUE, 1, FL_SOLID ),
		theoreticalMarkerLine( FL_GREEN, 1, FL_SOLID ),
		normalMarkerLine( FL_RED, 1, FL_SOLID ),

		empiricalCdfLine( FL_DARK_BLUE, 3, FL_DASH ),
		empiricalCdfMarker( &empiricalMarkerLine, 5 ),
		empiricalCdf( &empiricalCdfMarker, &empiricalCdfLine, &mono12 ),

		empiricalPdfLine( FL_DARK_BLUE, 3, FL_DOT ),
		empiricalPdfMarker( &empiricalMarkerLine, 5 ),
		empiricalPdf( &empiricalPdfMarker, &empiricalPdfLine, &mono12 ),

		theoreticalCdfLine( FL_GREEN, 3, FL_DASH ),
		theoreticalCdfMarker( &theoreticalMarkerLine, 5 ),
		theoreticalCdf( &theoreticalCdfMarker, &theoreticalCdfLine, &mono12 ),

		theoreticalPdfLine( FL_GREEN, 3, FL_DOT ),
		theoreticalPdfMarker( &theoreticalMarkerLine, 5 ),
		theoreticalPdf( &theoreticalPdfMarker, &theoreticalPdfLine, &mono12 ),

		normalCdfLine( FL_RED, 1, FL_DASH ),
		normalCdfMarker( &normalMarkerLine, 3 ),
		normalCdf( &normalCdfMarker, &normalCdfLine, &mono12 ),

		normalPdfLine( FL_RED, 1, FL_DOT ),
		normalPdfMarker( &normalMarkerLine, 3 ),
		normalPdf( &normalPdfMarker, &normalPdfLine, &mono12 ),

		crossHairLine( FL_GRAY0, 3, FL_DOT ),
		crossHairMarker(),
		crossHairSeriesH( &crossHairMarker, &crossHairLine, &mono12 ),
		crossHairSeriesV( &crossHairMarker, &crossHairLine, &mono12 ),

		handler( &chart )

		//
	{
		titleBox.align( FL_ALIGN_CENTER | FL_ALIGN_INSIDE );

		add( &stack, Location::North );
		stack.add( &titleBox );
		stack.add( &kmerLength );
		stack.add( &seed );
		stack.add( &numThreads );
		stack.add( &sampleSize );
		stack.add( &threshold );
		stack.add( &pValue );
		stack.add( &logScale );
		stack.add( &logLower );

		kmerLength.AddPropertyChangedEventHandler( this );
		seed.AddPropertyChangedEventHandler( this );
		numThreads.AddPropertyChangedEventHandler( this );
		sampleSize.AddPropertyChangedEventHandler( this );
		threshold.AddPropertyChangedEventHandler( this );
		pValue.AddPropertyChangedEventHandler( this );
		logScale.AddPropertyChangedEventHandler( this );
		logLower.AddPropertyChangedEventHandler( this );

		stack.add( &buttons );
		buttons.add( &showDisplay );
		buttons.add( &showHistograms );

		showDisplay.callback( ShowDisplay, this );
		showHistograms.callback( ShowHistograms, this );

		add( &disp, Location::Centre );
		add( &chart, Location::Centre );
		chart.SetXAxis( new LinearAxis{ 0, 2 * M_PI } )
			.SetYAxis( new LinearAxis{ -1, 1 } )
			.Add( &empiricalCdf )
			.Add( &empiricalPdf )
			.Add( &theoreticalCdf )
			.Add( &theoreticalPdf )
			.Add( &normalCdf )
			.Add( &normalPdf )
			.Add( &crossHairSeriesH )
			.Add( &crossHairSeriesV )
			.SetMargin( 75, 10, 30, 30 )
			.SetFillColour( FL_WHITE )
			.AddMouseHandler( this )
			.AddMouseHandler( &handler );

		ShowDisplay( &showDisplay, this );

		disp.textfont( FL_COURIER );
		disp.textsize( 16 );
		disp( "Get Kmer Distance Distributions: Select k-mer length and sample size above, then click \"Update\".\n" );
	}

	virtual ~GetKmerDistances() {}

	static void ShowDisplay( Fl_Widget *sender, void *target_ ) {
		auto target = (GetKmerDistances *) target_;
		target->showDisplay.deactivate();
		target->showHistograms.activate();
		target->disp.show();
		target->chart.hide();
	}

	static void ShowHistograms( Fl_Widget *sender, void *target_ ) {
		auto target = (GetKmerDistances *) target_;
		target->showDisplay.activate();
		target->showHistograms.deactivate();
		target->disp.hide();
		target->chart.show();
		target->chart.redraw();
	}

	void SetYAxis() {
		delete chart.YAxis();

		double xCrossesY = 0;
		bool isLogY = logScale.Value();

		if ( isLogY ) {
			xCrossesY = logLower.Value();
			auto axis = new LogarithmicAxis( 1, xCrossesY );
			chart.SetYAxis( axis ).SetAxesCross( 0, axis->Max() );
		}
		else {
			auto axis = new LinearAxis( 1, 0 );
			chart.SetYAxis( axis ).SetAxesCross( 0, axis->Max() );
		}

		Series *htics = chart.HTics();
		LineSpec tickLine;
		Plus tickMark( &tickLine, 5 );

		if ( htics->GetMarker() == 0 ) {
			htics->SetMarker( &tickMark );
		}
		htics->Clear();

		Series *vtics = chart.VTics();

		if ( vtics->GetMarker() == 0 ) {
			vtics->SetMarker( &tickMark );
		}

		vtics->Clear();

		int max = chart.XAxis()->Max();
		int numTicks = std::min( 10, max );

		for ( int i = 0; i <= numTicks; i++ ) {
			int x = i * max / numTicks;
			Point *p = htics->Add( x, xCrossesY );
			p->SetLabel( Double::ToString( x ), 0, 7, 0.5, 0.0, &mono12 );
		}

		double yMax = std::max( chart.YAxis()->Max(), chart.YAxis()->Min() );
		double yMin = std::min( chart.YAxis()->Max(), chart.YAxis()->Min() );
		numTicks = yMin == 0 ? 10 : (int) round( std::log10( yMax / yMin ) );

		for ( int i = 0; i <= numTicks; i++ ) {
			double y = isLogY ? yMin * pow( 10, i ) : i * (yMax - yMin) / numTicks;
			Point *p = vtics->Add( 0, y );
			p->SetLabel( Double::ToString( y ), -7, 0, 1.0, 0.5, &mono12 );
		}

		UpdateCrossHair();
	}

	/**
	 *	Generates theoretical and empirical distance distributions. Allows user
	 *	to specify a p-value for the clustering step.
	 */
	void Run() {
		TRACE;
		SetReady( false );
		thread t( RunImpl_, this );
		t.detach();
	}

	static void RunImpl_( GetKmerDistances * self ) {
		TRACE;
		self->RunImpl();
	}

	void RunImpl() {
		TRACE;
		if ( !mute.try_lock() ) {
			return;
		}

		TRACE;
		UPDATE_GUI( SetOk( false ); disp.Clear(); );

		TRACE;
		double startTime = omp_get_wtime();
		// fprintf( stderr, "Run called at %f.\n", startTime );

		ShowProgress(
			"Symbol frequency distribution:\n"  //
			"------------------------------\n" );

		TRACE;
		Histogram<Symbol> symbolHist = FastaSequence::GetSymbolHistogram( *dbSeqs );
		pAlphabet alpha = getAlphabet();

		TRACE;
		UPDATE_GUI( {
			for ( auto key : symbolHist.GetKeys() ) {
				double v = symbolHist[key];
				disp( "%3c|", alpha->Decode( key ) );

				for ( int i = 0; i < v * 100; i++ ) {
					disp( "%c", '=' );
				}

				disp( " (p=%8.6f)\n", v );
			} }
		);

		uint kmerLength = KmerLength();

		function<Distance( Symbol x, Symbol y )> symbolDistance = [this]( Symbol x, Symbol y ) {
			return matrix->Difference( x, y );
		};

		kmerDistanceTheoretical = IntegerDistribution::GetKmerDistanceDistribution<Symbol, Distance>(
			symbolHist, symbolDistance, kmerLength );

		if ( 0 < PValue() && PValue() < 1 && Threshold() < 0 ) {
			UPDATE_GUI( SetThreshold( ThresholdForPval( PValue() ) ); );
		}

		uint sampleSize = this->sampleSize.Value();
		Distance max = (matrix->maxValue - matrix->minValue ) * kmerLength;

		double mu = kmerDistanceTheoretical.Mean();
		double sigma = kmerDistanceTheoretical.StdDev();
		normal.SetMean( mu );
		normal.SetStdDev( sigma );

		TRACE;
		UPDATE_GUI( {
			disp( "\n\nmu\t%0.17f\nsigma\t%0.17f\n", mu, sigma );
			disp( "\n\nTheoretical k-mer distance distribution.\n\nx\tp\tf\tnp\tnf\n" );

			for ( Distance x = 0; x <= max; x++ ) {
				double F = kmerDistanceTheoretical.Cdf( x );
				double p = kmerDistanceTheoretical.Pdf( x );
				double np = normal.Pdf( x );
				double nf = normal.Cdf( x );
				disp( "%d\t%0.10e\t%0.10e\t%0.10e\t%0.10e\n", x, p, F, np, nf );
			} }
		);

		TRACE;
		ShowProgress( "Theoretical distance calc done in %fs.\n", (omp_get_wtime() - startTime) );

		TRACE;
		vector<Distance> observedDistances;

		ShowProgress( "Computing empirical k-mer distance sample...\n" );

		rand.Reseed( Seed() );

		size_t totalKmerCount = 0;

		for ( auto i : *trainingSetIndices ) {
			totalKmerCount += (*dbSeqs)[i]->KmerCount( kmerLength );
		}

		Selector select( rand, 2 * sampleSize, totalKmerCount );
		vector<SimpleKmer::Instance> sample;

		for ( auto i : *trainingSetIndices ) {
			auto seq = (*dbSeqs)[i];
			for ( size_t j = 0; j < seq->KmerCount( kmerLength ); j++ ) {
				if ( select.SelectThis() ) {
					sample.emplace_back( seq, j );
				}
			}
		}

		observedDistances.clear();

		for ( uint i = 0; i < sampleSize; i++ ) {
			uint idx = i + uint( rand() * (sampleSize - i) );
			std::swap( sample[2 * i], sample[idx] );
			idx = i + 1 + uint( rand() * (sampleSize - i - 1) );
			std::swap( sample[2 * i + 1], sample[idx] );
			SimpleKmer::Instance k1 = sample[2 * i];
			SimpleKmer::Instance k2 = sample[2 * i + 1];
			Distance d = matrix->Difference( k1.Bytes(), k2.Bytes(), kmerLength );
			// Distance d = matrix->DigramDifference( k1.Digrams(), k2.Digrams(), kmerLength / 2 );
			observedDistances.push_back( d );

			//if ( d != d_ ) {
			//	cerr << "Distance calc error.\nk1 = " << k1.Chars(kmerLength) << "\nk2 = " << k2.Chars(kmerLength) << "\nd (bytes) = " << d << "\nd (digrams) = " << d_ << "\n";
			//}
		}

		ShowProgress( "Kmer distance sample done after %fs.\n", (omp_get_wtime() - startTime) );

		Histogram<Distance> empiricalHist;
		empiricalHist.Initialise( observedDistances );

		IntegerDistribution empiricalDist( empiricalHist );
		kmerDistanceObserved = empiricalDist;

		TRACE;
		UPDATE_GUI( {
			disp( "\n\nEmpirical k-mer distance distribution.\n\nx\tp\tf\n" );

			for ( Distance x = 0; x <= max; x++ ) {
				double F = empiricalDist.Cdf( x );
				double p = empiricalDist.Pdf( x );
				disp( "%d\t%0.10e\t%0.10e\n", x, p, F );
			} } );

		UPDATE_GUI( {
			disp( "Empirical distributions fitted in %fs.\nDrawing chart...", (omp_get_wtime() - startTime) );
			DrawChart();
			disp( "Task complete after %fs.\n", (omp_get_wtime() - startTime) );
			UpdateReadiness();
			} );
		TRACE;

		Page::runTime = (omp_get_wtime() - startTime);

		UPDATE_GUI( UpdateReadiness(); NotifyRunComplete(); );
		mute.unlock();
		TRACE;
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

	void GainFocus() override {
		Page::GainFocus();
		alphabet = getAlphabet();
		matrix = getMatrix();
		dbSeqs = getDb();
		trainingSetIndices = getTrainingSetIndices();
		UpdateReadiness();
	}

	void UpdateReadiness() {
		bool isReady = alphabet != 0 && matrix != 0 && dbSeqs != 0 && trainingSetIndices != 0;
		bool isOk = isReady
			&& (kmerDistanceTheoretical.Max() > kmerDistanceTheoretical.Min())
			&& threshold.Value() >= 0 /* == ThresholdForPval( pValue.Value() ) */
			;
		SetReady( isReady );
		SetOk( isOk );
	}

	void DrawChart() {
		double max = kmerDistanceTheoretical.Max();

		chart.XAxis()->SetMin( 0 );
		chart.XAxis()->SetMax( max );

		SetYAxis();

		theoreticalCdf.Clear();
		theoreticalPdf.Clear();

		for ( Distance x = 0; x <= max; x++ ) {
			theoreticalCdf.Add( x, kmerDistanceTheoretical.Cdf( x ) );
			theoreticalPdf.Add( x, kmerDistanceTheoretical.Pdf( x ) );
		}

		normalCdf.Clear();
		normalPdf.Clear();

		for ( Distance x = 0; x <= max; x++ ) {
			normalCdf.Add( x, normal.Cdf( x ) );
			normalPdf.Add( x, normal.Pdf( x ) );
		}

		empiricalCdf.Clear();
		empiricalPdf.Clear();

		for ( Distance x = 0; x <= max; x++ ) {
			empiricalCdf.Add( x, kmerDistanceObserved.Cdf( x ) );
			empiricalPdf.Add( x, kmerDistanceObserved.Pdf( x ) );
		}

		UpdateCrossHair();
		chart.redraw();
	}

	void UpdateCrossHair() {
		crossHairSeriesH.Clear();
		crossHairSeriesV.Clear();

		double p = pValue.Value();
		int t = threshold.Value();

		if ( p > 0 && p < 1 && t > 0 ) {
			crossHairSeriesH.Add( 0, p );
			crossHairSeriesH.Add( t, p );
			crossHairSeriesV.Add( t, p );
			crossHairSeriesV.Add( t, chart.YAxis()->Max() );
		}
	}

	bool HandleEvent( Fl_Widget *src, int eventCode ) override {
		// fprintf(stderr, "Received event %d\n", eventCode);

		if ( FL_RELEASE == eventCode ) {
			// fprintf(stderr, "Received release event.\n");

			int sx = Fl::event_x();
			int sy = Fl::event_y();
			Rectangle<int> plotArea = chart.PlotArea();

			double x = chart.XAxis()->ToWorld( sx, plotArea.left, plotArea.right );
			double y = chart.YAxis()->ToWorld( sy, plotArea.top, plotArea.bottom );

			this->threshold.SetValue( (int) (round( x )) );
			this->pValue.SetValue( y );

			UpdateCrossHair();
			UpdateReadiness();
			chart.redraw();
			return true;
		}

		return false;
	}

	bool Handle( const ScatterPlot::MouseEvent & event ) override {
		if ( event.eventCode == FL_RELEASE ) {
			threshold.SetValue( (int) round( event.x ) );
			pValue.SetValue( event.y );
			UpdateCrossHair();
			UpdateReadiness();
			chart.redraw();
			return true;
		}

		return false;
	}

	/// Compute the threshold distance corresponding to a designated probability.
	///	@param p The probability.
	///	@returns int(round(kmerDistanceTheoretical.InverseCdf(p)))
	int ThresholdForPval( double p ) { return int( round( kmerDistanceTheoretical.InverseCdf( p ) ) ); }

	/// Compute the p-value corresponding to a designated distance.
	/// @param threshold The distance.
	///	@returns  kmerDistanceTheoretical.Cdf(threshold)
	double PvalForThreshold( int threshold ) { return kmerDistanceTheoretical.Cdf( threshold ); }

	void PropertyChanged( void *sender, const string &propertyName ) override {
		if ( sender == &this->pValue ) {
			if ( pValue.HasValue() ) {
				double p = pValue.Value();

				if ( p < 0 || p > 1 ) return;

				threshold.SetValue( ThresholdForPval( p ) );
				pValue.SetValue( p );
			}
			else {
				pValue.SetValue( 0 );
				threshold.SetValue( 0 );
			}

			UpdateCrossHair();
			chart.redraw();
		}
		else if ( sender == &this->threshold ) {
			auto newThreshold = Threshold();

			if ( newThreshold < 0 ) return;

			auto oldInferredThreshold = ThresholdForPval( PValue() );

			if ( newThreshold != oldInferredThreshold ) {
				pValue.SetValue( PvalForThreshold( threshold.Value() ) );
				UpdateCrossHair();
				chart.redraw();
			}
		}
		else if ( sender == &this->logScale ) {
			SetYAxis();
			chart.redraw();
		}
		else if ( sender == &this->logLower ) {
			SetYAxis();
			chart.redraw();
		}
		// else {
		// 	fprintf( stderr, "some unknown source posted PropertyChanged(\"%s\")", propertyName.c_str() );
		// }

		UpdateReadiness();
		redraw();
	}

#pragma region TagName Property
private:
	string tagName;

public:
	const string & TagName() const {
		static string defaultVal = "Get Kmer Distances";
		return tagName.length() == 0 ? defaultVal : tagName;
	}

	void SetTagName( const string & value ) {
		tagName = value;
	}
#pragma endregion


	/// Get parameters so they can be serialised when program finishes.
	/// @param parms A Param list to which the component should append its current parameters. 
	void GetParams( set<Param> &parms ) override {
#define add_parm(x) parms.emplace( Name(), STRING(x), (x).Value() )
		add_parm( kmerLength );
		add_parm( seed );
		add_parm( numThreads );
		add_parm( sampleSize );
		add_parm( threshold );
		add_parm( pValue );
		add_parm( logScale );
		add_parm( logLower );

		parms.emplace( Name(), "runTime", Util::ToString( runTime ) );
		parms.emplace( Name(), "loadTime", Util::ToString( loadTime ) );
		parms.emplace( Name(), "saveTime", Util::ToString( saveTime ) );
#undef add_parm
	}

	/// Set parameters which have been de-serialised when program starts.
	/// @param parms A Param list containing de-serialised parameters.
	void SetParams( const set<Param> & parms ) override {
		ResetForm();

#define setParm(x,t) if ( p.ParamName() == STRING(x) ) (x).SetValue( p.Value<t>() )
		for ( auto & p : parms ) {
			if ( p.ComponentName() != Name() && p.ComponentName() != TagName() ) continue;
			setParm( kmerLength, uint );
			setParm( seed, ulong );
			setParm( numThreads, int );
			setParm( sampleSize, size_t );
			setParm( threshold, Distance );
			setParm( pValue, double );
			setParm( logScale, int );
			setParm( logLower, double );

			if ( p.ParamName() == "runTime" )  runTime = p.Value<double>();
			if ( p.ParamName() == "loadTime" ) loadTime = p.Value<double>();
			if ( p.ParamName() == "saveTime" ) saveTime = p.Value<double>();
		}

		UpdateReadiness();
#undef setParm
	}

	double PValue() const {
		return pValue.Value();
	}

	void SetPValue( double value ) {
		pValue.SetValue( value );
		PropertyChanged( &pValue, "Value" );
	}

	int Threshold() const {
		return threshold.Value();
	}

	void SetThreshold( int value ) {
		threshold.SetValue( value );
		PropertyChanged( &threshold, "Value" );
	}

	uint KmerLength() const {
		return kmerLength.Value();
	}

	void SetKmerLength( uint value ) {
		if ( KmerLength() == value ) return;
		kmerLength.SetValue( value );
		PropertyChanged( &kmerLength, "Value" );
	}

	ulong Seed() {
		return seed.Value();
	}

	void SetSeed( ulong val ) {
		if ( Seed() == val ) return;

		seed.SetValue( val );
		PropertyChanged( &seed, "Value" );
	}

	void Reset() override {
		alphabet = 0;
		matrix = 0;
		dbSeqs = 0;
		trainingSetIndices = 0;
	}

	void ResetForm() {
		kmerLength.Reset();
		seed.Reset();
		sampleSize.Reset();
		threshold.Reset();
		pValue.Reset();
		logScale.Reset();
		logLower.Reset();
	}

};
