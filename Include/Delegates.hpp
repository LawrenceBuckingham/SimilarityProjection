#pragma once


// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <functional>

namespace QutBio {
	using Action = std::function<void( void )>;

	template <typename T>
	using Action1 = std::function<void( T t )>;

	template <typename T1, typename T2>
	using Action2 = std::function<void( T1  t1, T2  t2 )>;

	template <typename T1, typename T2, typename T3>
	using Action3 = std::function<void( T1 t1, T2 t2, T3 t3 )>;

	template <typename T>
	using Func = std::function<T( void )>;

	template <typename T, typename Out>
	using Func1 = std::function<Out( T t )>;

	template <typename T1, typename T2, typename Out>
	using Func2 = std::function<Out( T1  t1, T2  t2 )>;

	template <typename T1, typename T2, typename T3, typename Out>
	using Func3 = std::function<Out( T1 t1, T2 t2, T3 t3 )>;
}
