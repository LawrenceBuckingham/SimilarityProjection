#pragma once

#include "IPropertyChangedEventHandler.hpp"
#include <vector>

#define NAMEOF(x) #x

using namespace std;

namespace LBFL {
	/**
	 *	An object which broadcasts PropertyChanged messages to a list of registered listeners. 
	 */
	class PropertyChangedEventSource
	{
	private:
		vector<IPropertyChangedEventHandler *> listeners;
		void * eventSource;

	public:
		PropertyChangedEventSource( void * src ) : eventSource(src) {}

		virtual ~PropertyChangedEventSource() {}

		/**
		 *	Adds a PropertyChanged event listener.
		 */
		void AddPropertyChangedEventHandler(IPropertyChangedEventHandler * handler)
		{
			listeners.push_back(handler);
		}

		/**
		 *	Broadcasts a PropertyChanged message to all registered listeners.
		 */
		void NotifyPropertyChanged(const string & propertyName)
		{
			// cerr << "Received NotifyPropertyChanged(" << propertyName << ")\n";

			for (auto listener : listeners)
			{
				listener->PropertyChanged(eventSource, propertyName);
			}
		}
	};
}
