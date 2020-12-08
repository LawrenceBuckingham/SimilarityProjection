#pragma once

#include <functional>
#include <set>

#include "BorderLayout.hpp"
#include "IRunnable.hpp"

using namespace std;
using namespace LBFL;

/**
 *	The common base class of pages in this App. I tried to do without this by
 *	separately sub-classing BorderLayou and IRunnable in the Page classes
 *	themselves, using IRunnable as an interface, but the pointers didn't work
 *	the way I expected them to.
 *
 */
class Page :
	public virtual BorderLayout,
	public virtual IRunnable {
protected:
	/** Gets the name of the page. */
	string name;

	/// List of event handlers for run complete.
	vector<function<void(Page *page)>> runComplete;

	/// Store the most recent elapsed run time here.
	double runTime = 0;

	/// Store the most recent load time here.
	double loadTime = 0;

	/// Store the most recent save time here. 
	double saveTime = 0;

public:
	/**
	 *	Constructs a new Page.
	 *
	 *	@param	left	The left edge of the control, nominally relative to the window.
	 *	@param	top		The top edge of the control, nominally relative to the window.
	 *	@param	width	The nominal width of the control.
	 *	@param	height	The nominal height of the control.
	 */
	Page(
		int left,
		int top,
		int width,
		int height,
		const string & name
	) : BorderLayout( left, top, width, height ), name( name ) {
		Fl_Group::current( 0 );
		Fl_Widget::box( FL_FLAT_BOX );
		RunComplete( [this](Page *p) { InvalidateConsumers(); } );
	}

	/// I don't want to clone pages.
	Page( const Page & other ) = delete;

	/// I don't want to clone pages.
	Page & operator=( const Page & other ) = delete;

	/**
	 * Gets the name of the control.
	 *
	 * @returns Returns a copy of the string containing the control name.
	 */
	const string & Name() const {
		return name;
	}

	/// 
	class Param :
		public virtual ICsvWriter,
		public virtual ICsvReader {
		string componentName, paramName, value;
	public:
		Param(
			const string & componentName = "",
			const string & paramName = "",
			const string & value = ""
		) :
			componentName( componentName ), paramName( paramName ), value( value ) {}

		template<typename T>
		Param(
			const string & componentName,
			const string & paramName,
			const T & value
		) :
			componentName( componentName ), paramName( paramName ), value( Util::ToString( value ) ) {}

		string ComponentName() const { return componentName; }

		string ParamName() const { return paramName; }

		string Value() const { return value; }

		template<typename T>
		T Value() const {
			T t;
			istringstream s( value );
			s >> t;
			return t;
		}

		void Write( CsvWriter & w ) const {
			w << componentName << paramName << value;
		}

		void Read( CsvReader &r ) {
			r >> componentName >> paramName >> value;
		}

		friend bool operator<( const Param & lhs, const Param & rhs ) {
			if ( lhs.componentName < rhs.componentName ) return true;
			else if ( rhs.componentName < lhs.componentName ) return false;
			else if ( lhs.paramName < rhs.paramName ) return true;
			else return false;
		}
	};

	/// Get parameters so they can be serialised when program finishes.
	/// @param parms A Param list to which the component should append its current parameters. 
	virtual void GetParams( set<Param> &parms ) = 0;

	/// Set parameters which have been de-serialised when program starts.
	/// @param parms A Param list containing de-serialised parameters.
	virtual void SetParams( const set<Param> & parms ) = 0;

	/// Hides all child controls in the Centre, other than the supplied argument.
	///	The child `panel` is (if found) made visible.
	/// @param panel The child control to make visible in the centre.
	void ShowCentrePanel( Fl_Widget * panel ) {
		for ( auto ix : this->Centre() ) {
			if ( ix == panel ) {
				ix->show();
			}
			else {
				ix->hide();
			}
		}
	}

	/**
	 *	Adds an action to the list executed by NotifyRunComplete.
	 */
	void RunComplete( function<void(Page * page)> action ) {
		runComplete.push_back(action);
	}

	/**
	 *	Calls the list of actions inserted via RunComplete. Preferably call this
	 *	from Run to tell the shell to advance the display. If called from a worker
	 *	thread you will want to use UPDATE_GUI to allow the UI to be updated.
	 */
	void NotifyRunComplete() {
		for ( auto & action : runComplete ) {
			action( this );
		}
	}
};
