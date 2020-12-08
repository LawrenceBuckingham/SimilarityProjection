#include <map>
#include <string>
#include <vector>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <FileUtil.hpp>
#include <DataLoader.hpp>

using namespace std;
using namespace QutBio;

namespace Utilities {
	struct GetDomainSubset {
		static void Run( int argc, char** argv ) {
			Args args( argc, argv );
			Params parms( args );

			map<string, Domain> allDomains;
			LoadDomains( parms.domainIn, allDomains );
			auto alphabet = Alphabets::AA();
			auto allFasta = Load::Fasta(parms.seqIn, parms.idIndex, alphabet );
			auto allSequences = Load::Encoded( allFasta, parms.classIndex, Alphabets::AA(), 1, 1, alphabet->DefaultSymbol());

			vector<const Domain *> selectedDomains;
			SelectDomains( parms.domainIds, allDomains, selectedDomains );

			map<string, EncodedFastaSequence *> selectedSequences;
			SelectSequences( allSequences, selectedDomains, selectedSequences );

			ofstream dOut( parms.domainOut );
			for ( auto dom : selectedDomains ) {
				dOut << (*dom);
			}

			ofstream sOut(parms.seqOut);
			for ( auto & p : selectedSequences ) {
				sOut << p.second;
			}
		}

		static void SelectDomains(
			const vector<string> & domainIds,
			const map<string, Domain> & domains,
			vector<const Domain *> & selectedDomains ) //
		{
			for ( auto & domainId : domainIds ) {
				auto dom = domains.find( domainId );

				if ( dom != domains.end() ) {
					selectedDomains.push_back( &dom->second );
				}
			}
		}

		static void SelectSequences(
			vector<EncodedFastaSequence *> & sequences,
			const vector<const Domain *> & domains,
			map<string, EncodedFastaSequence *> & selectedSequences //
		) {
			Index<EncodedFastaSequence> idx( sequences );

			for ( auto dom : domains ) {
				for ( auto & entry : dom->entries ) {
					const string & seqId = entry.second.seqId;
					auto seq = idx.find( seqId );

					if ( seq != idx.end() ) {
						selectedSequences[seqId] = seq->second;
					}
				}
			}
		}

		static void LoadDomains( const string & domFileName, map<string, Domain> & domains ) {
			OMP_TIMER_DECLARE( domLoad );
			OMP_TIMER_START( domLoad );

			ifstream domFile( domFileName );
			Domain::Load( domFile, domains );
			domFile.close();

			OMP_TIMER_END( domLoad );
			cerr << domains.size() << " domains loaded from " << domFileName << " in " << OMP_TIMER( domLoad ) << "s\n";
		}

		struct Params {
			bool ok;
			string domainIn, domainOut, seqIn, seqOut;
			vector<string> domainIds;
			int idIndex, classIndex;

			Params( Args & args ) {
				ok = true;

				if ( !args.Get( "domainIn", domainIn ) ) {
					cerr << "Argument 'domainIn' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "domainOut", domainOut ) ) {
					cerr << "Argument 'domainOut' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "seqIn", seqIn ) ) {
					cerr << "Argument 'seqIn' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "seqOut", seqOut ) ) {
					cerr << "Argument 'seqOut' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "domainIds", domainIds ) ) {
					cerr << "Argument 'domainIds' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "idIndex", idIndex ) ) {
					cerr << "Argument 'idIndex' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "classIndex", classIndex ) ) {
					cerr << "Argument 'classIndex' not supplied.\n";
					ok = false;
				}

				if ( !ok ) {
					cerr << "Example\n GetDomainSubset.exe --domainIn swissprot.domains --seqIn m:/swissprot/swissprot_n_ge_5/swissprot_n_ge_5.faa --domainOut PF01336.domains --seqOut PF01336.faa --domainIds PF01336 --idIndex 2 --classIndex 3" << '\n';
					throw Exception( "Invalid arguments.", FileAndLine );
				}
			}

		};

	};
}

using namespace Utilities;

int main( int argc, char** argv ) {
	try {
		GetDomainSubset::Run( argc, argv );
	}
	catch ( Exception & ex ) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
		return 1;
	}
	catch ( runtime_error & err ) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
		return 1;
	}
	return 0;
}

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
