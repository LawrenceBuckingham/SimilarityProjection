#include "TestFramework.h"

#include "TestArgs.hpp"
#include "TestCsv.hpp"
#include "TestEnums.hpp"
#include "TypedArray.hpp"

using namespace SystemQut;

int main() {
	TestFramework::RunAllTests( TestArgs::Tests() );
	TestFramework::RunAllTests( TestEnums::Tests() );
	TestFramework::RunAllTests( TestCsv::Tests() );
}