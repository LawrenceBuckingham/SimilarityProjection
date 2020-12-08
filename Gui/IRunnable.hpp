#pragma once

#include <IPropertyChangedEventHandler.hpp>
#include <PropertyChangedEventSource.hpp>
#include <vector>

class IRunnable :
	public virtual PropertyChangedEventSource
{
private:
	/** Backing store for Ok property. */
	bool ok = false;
	bool ready = false;

protected:
	/** Observers for the invalidate event. */
	vector<IRunnable *> consumers;
public:
	virtual ~IRunnable() {}

	/// Called by containing application when the current page is brought to the foreground.
	///	Override to perform custom actions when focus gained. 
	///	Typical use case would be to fetch data from upstream as a precursor to Run. 
	virtual void GainFocus() {}

	/// Called by containing application when the update button is clicked with the current page in the foreground.
	/// Override to perform custom actions when the . 
	///	Typical use case would be to fetch data from upstream as a precursor to Run. 
	virtual void Run() = 0;

	virtual void LoseFocus() {}

	/**
	 *	Returns true if and only if this object is OK. That should be set AFTER 
	 *	it has successfully Run().
	 */
	bool Ok()
	{
		return ok;
	}

	/**
	 *	Sets a new value for Ok(), and transmits a PropertyChanged message to all
	 *	registered PropertyChanged event listeners. If the state does not change,
	 *	then the function does nothing, and returns without sending the message.
	 *	@param value The new value for Ok().
	 */
	void SetOk(bool value)
	{
		// cerr << "Received SetOk(" << value << ")\n";

		if (ok == value) return;

		ok = value;
		NotifyPropertyChanged(NAMEOF(Ok));
	}

	/// Return true if and only if the object is ready to Run. This is what determines 
	///	the state of the "Update" button in the application.
	bool Ready() {
		return ready;
	}

	/// Sets a new value for Ready(), and transmits a PropertyChanged message to all
	///	registered PropertyChanged event listeners. If the new value is the same as the 
	///	old value, no action is taken and no change message is broadcast.
	///	@param newValue The new state of the Ready() property.
	void SetReady(bool newValue)
	{
		if (Ready() == newValue) { return; }

		ready = newValue;
		NotifyPropertyChanged(NAMEOF(Ready));
	}

	/**
	 *	Adds a listener for the Invalidate event. These would normally be down-stream components which
	 *	depend on the current component.
	 *	@param consumer A runnable object which wants to know if the current component is invalidated.
	 */
	void AddConsumer(IRunnable * consumer)
	{
		consumers.push_back(consumer);
	}

	/**
	 *	Makes the current component "not ok" and propagates the Invalidate message to dependent
	 *	components to ensure they cannot execute without a valid antecedent (namely, this).
	 */
	void Invalidate()
	{
		SetOk(false);
		Reset();
		InvalidateConsumers();
	}

	void InvalidateConsumers() {
		for ( auto consumer : consumers ) {
			consumer->Invalidate();
		}
	}

	/// Reset this page to an effectively empty state.
	virtual void Reset () = 0;
};

