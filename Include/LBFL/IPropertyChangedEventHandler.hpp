#pragma once

namespace LBFL {
	class IPropertyChangedEventHandler {
	public:
		virtual ~IPropertyChangedEventHandler() {}
		virtual void PropertyChanged( void* sender, const string& propertyName ) = 0;
	};
}  // namespace LBFL
