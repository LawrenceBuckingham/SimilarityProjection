#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>
#include <string>
#include <cstdlib>

using namespace std;

#include <Args.hpp>
#include <Domain.hpp>
#include <FileUtil.hpp>

using namespace QutBio;

namespace Utilities
{
	class ImportSwissprotDomains
	{
	public:

		/**
		Parses the swissprot domain text file swisspfam.
		*/
		static void Run(Args & args)
		{
			Parameters parms(args);

			ifstream in(parms.inFile);
			map<string, Domain> domains;

			for (bool ok = Domain::Parse(in, domains); ok; ok = Domain::Parse(in, domains)) {}

			in.close();

			ofstream out( parms.outFile );
			for ( auto & dom: domains ) {
				out << dom.second;
			}
		}

	private:
		struct Parameters
		{
			string inFile;
			string outFile;
			bool ok = true;

			Parameters(Args & args)
			{
				if (args.IsDefined("help"))
				{
					vector<string> text {
						"ImportSwissprotDomains: Parses the swisspfam text file and writes out a",
						"more accessible and compact list of pfam domain occurrences in all sequences.",
						"",
						"Results are ordered by PFam Id, with separate records for each sequence in the",
						"family ordered in ascending order of the Swissprot key.",
						"",
						"Arguments:",
						"--inFile   Required. The path of the swisspfam text file.",
						"",
						"--outFile  Required. The path of a file that will be overwritten with the list of",
						"           domains."
					};

					for (auto & s : text)
					{
						cerr << s << '\n';
					}
				}

				if (!args.Get("inFile", inFile))
				{
					cerr << "Argument 'inFile' not defined." << endl;
					ok = false;
				}

				if (!args.Get("outFile", outFile))
				{
					cerr << "Argument 'outFile' not defined." << endl;
					ok = false;
				}

				if (!ok)
				{
					throw Exception("Error in parameters.", FileAndLine);
				}
			}
		};

	};
};

using namespace Utilities;

int main(int argc, char ** argv)
{
	try
	{
		Args args(argc, argv);
		ImportSwissprotDomains::Run(args);
	}
	catch (Exception & ex)
	{
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
		return 1;
	}
	catch (runtime_error & err)
	{
		cerr << "Unhandled exception:" << endl << err.what() << endl;
		return 1;
	}
	return 0;
}
