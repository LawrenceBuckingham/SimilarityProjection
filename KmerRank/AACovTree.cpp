#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "KmerDistanceCache.hpp"
#include "FastaSequence.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "KmerIndex.hpp"
#include "OmpTimer.h"
#include "KmerCluster.hpp"

#define TRON 1
#include "db.hpp"

#include <cstdio>
#include <omp.h>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct AACovTree {
public:
    using DistanceFunction = KmerDistanceCache2;
    using Cluster = KmerCluster<DistanceFunction, Kmer>;
    using Seq = EncodedFastaSequence;

    static int Run() {
        OMP_TIMER_DECLARE( loadTime );
        OMP_TIMER_DECLARE( clusterTime );

        Params parms;

        if ( parms.ok == false || !parms.matrix ) {
            cerr << arguments->ProgName() << ": error - unable to construct similarity matrix. For help, run: AAClust --help\n";
            return 1;
        }

        Alphabet * alphabet = new Alphabet( parms.matrix );
        BlosumDifferenceFunction rawDistanceFunction( parms.matrix );
        DistanceFunction distanceFunction( alphabet, &rawDistanceFunction );
        UniformRealRandom rand( parms.seed );

        omp_set_num_threads( parms.numThreads );

        OMP_TIMER_START( loadTime );

        auto initSeq = [&]( Seq & instance ) {
            instance.Init( -1, alphabet, parms.wordLength,
                distanceFunction.CharsPerWord(), 'x' );
        };

        auto db = Seq::ReadSequences<Seq>(
            parms.fastaFile, parms.idIndex, initSeq
            );

        const int N = db.size();

        cerr << arguments->ProgName() << ": " << N << " sequences loaded.\n";

        KmerIndex kmerIndex( db, parms.wordLength );
        auto kmers = kmerIndex.GetKmers();
        auto K = kmers.size();

        OMP_TIMER_END( loadTime );
        OMP_TIMER_START( clusterTime );

        // Shuffle the kmers to avoid pathological selection.
        for ( size_t i = 0; i < K; i++ ) {
            auto newLoc = i + (size_t) ( rand() * ( K - i ) );
            std::swap( kmers[i], kmers[newLoc] );
        }

        vector<double> threshold = parms.threshold;
        sort( threshold.begin(), threshold.end(), []( double & lhs, double & rhs ) { return rhs < lhs; } );

        for ( auto t : threshold ) {
            cerr << arguments->ProgName() << ": threshold item: " << t << "\n";
        }

        omp_set_num_threads( parms.numThreads );

        OMP_TIMER_END( loadTime );

        Node * currentNode = new Node();
        size_t thresholdIdx = 0;

        ClusterKmersByThresholdRecursive( kmerIndex, alphabet, parms.wordLength, distanceFunction.CharsPerWord(), threshold, thresholdIdx, distanceFunction, parms.increment, parms.minClusterSize, rand, currentNode );


        OMP_TIMER_END( clusterTime );

        cerr << "Elapsed time loading: " << OMP_TIMER( loadTime ) << "\n";
        cerr << "Elapsed time clustering: " << OMP_TIMER( clusterTime ) << "\n";

        return 0;
    }

    struct Node {
        double threshold;
        vector<pEncodedFastaSequence> protos;
        vector<Cluster *> clusters;
        vector<Node *> children;

        ~Node() {
            for ( auto cluster : clusters ) {
                delete cluster;
            }
            for ( auto node : children ) {
                delete node;
            }
            for ( auto proto : protos ) {
                delete proto;
            }
        }
    };

    template<typename KmerPointerSource>
    static void ClusterKmersByThresholdRecursive(
        KmerPointerSource & kmerPointerSource,
        Alphabet * alphabet,
        size_t wordLength,
        size_t charsPerWord,
        vector<double> & threshold,
        size_t thresholdIndex,
        DistanceFunction & distanceFunction,
        size_t increment,
        size_t minClusterSize,
        UniformRealRandom & rand,
        Node * currentNode
    ) {
        vector<Kmer *> & allKmers = kmerPointerSource.GetKmers();

        ClusterKmersByThreshold( 
            allKmers, 
            alphabet, 
            wordLength, 
            distanceFunction.CharsPerWord(), 
            threshold[thresholdIndex], 
            distanceFunction, 
            increment, 
            rand, 
            currentNode->protos, 
            currentNode->clusters 
        );

        auto indent = [thresholdIndex]() { for ( uint i = 0; i < thresholdIndex; i++ ) cerr << "  "; };

        indent();
        cerr << "Kmers clustered at threshold " << threshold[thresholdIndex] << " under " << currentNode->protos.size() << " prototypes\n";

        if ( thresholdIndex < threshold.size() - 1 ) {
            for ( size_t i = 0; i < currentNode->clusters.size(); i++ ) {
                if ( currentNode->clusters[i]->kmers.size() < minClusterSize ) continue;

                // PR( currentNode->clusters[i]->Size(), %zu );
                indent();
                cerr << "Processing child cluster " << i << ": " << currentNode->clusters[i]->kmers.size() << " kmers\n";

                Node * newNode = new Node();
                currentNode->children.push_back( newNode );
                ClusterKmersByThresholdRecursive( *( currentNode->clusters[i] ), alphabet, wordLength, charsPerWord, threshold, thresholdIndex + 1, distanceFunction, increment, minClusterSize, rand, newNode );
            }
        }
        else {
            indent();
            cerr << "Leaf cluster would be saved to disk here." << "\n";
        }

        for ( auto c : currentNode->clusters ) {
            delete c;
        }

        currentNode->clusters.clear();
    }

    static void ClusterKmersByThreshold(
        vector<Kmer *> & allKmers,
        Alphabet * alphabet,
        size_t wordLength,
        size_t charsPerWord,
        double threshold,
        DistanceFunction & distanceFunction,
        size_t increment,
        UniformRealRandom & rand,
        vector<KmerClusterPrototype> & protos,
        vector<Cluster *> & clusters
    ) {
        auto createPrototype = [&]( Kmer *kmer ) {
            protos.emplace_back( kmer->Word(), alphabet, wordLength, charsPerWord);
            return protos.back();
        };

        Cluster::DoExhaustiveIncrementalClustering(
            allKmers,
            wordLength,
            threshold,
            alphabet->Size(),
            distanceFunction,
            rand,
            createPrototype,
            clusters
        );

        // Update prototype sizes.
        {
            for ( auto cluster : clusters ) {
                auto & proto = (KmerClusterPrototype &) cluster->prototype->Sequence();
                proto.Size( proto.Size() + cluster->InstanceCount() );
            }
        }
    }

    struct Params {
        bool ok = true;

        string fastaFile;
        int numThreads = 7;
        int wordLength = 30;
        vector<double> threshold;
        int  seed;
        int idIndex;
        string clusterOut;
        string protoOut;
        string treeOut;
        DistanceType * distanceType;
        SimilarityMatrix * matrix;
        bool isCaseSensitive = false;
        size_t increment = 1;
        size_t minClusterSize = 10;

        Params() {
            if ( arguments->IsDefined( "help" ) ) {
                vector<string> text{
    "AACovTree: Greedy k-tree-inspired tree structured clustering of Amino Acid kmers by substitution matrix.",
    "           Node = { KmerClusterPrototype prototype; uint level; }",
    "           StructuralNode = Node + { vector<Node *> children; }",
    "           StorageNodeNode = Node + { vector<Kmer> items; }",
    "           Tree = Node *;"
    "--help	Gets this text.",
    "--fastaFile    Required. A list of one or more file paths. Each file will be parsed as a FASTA file which contains DNA sequences to be clustered.",
    "--idIndex      Required. The 0-origin position of the sequence ID field in the pipe-separated definition line.",
    "--clusterOut   Required. The name the output file produced by the program.",
    "--protoOut     Required. The name the output file produced by the program.",
    "--threshold    Required. List of threshold distances used to assign kmers to clusters. Distance less than or equal to threshold corresponds to cluster membership. Kmers are added to clusters at the leaves which have the smallest threshold. Above leaf level, branches contain a prototype and a list of child branches which have a prototype within the ",
    "--seed         Required. The random number seed.",
    "--numThreads   Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
    "--wordLength   Optional; default value = 32. The word length used for kmer tiling.",
    "--matrixId     Optional, default = 62, but you need either --matrixId or --matrixFile. BLOSUM Matrix ID, one of { 35, 40, 45, 50, 62, 80, 100 }. This is ignored if a custom similarity matrix file is specified.",
    "--matrixFile   Optional, but you need either --matrixId or --matrixFile. File name for custom similarity matrix. Use this to specify some matrix other than BLOSUM, or if a custom alphabet is in use.",
    "--isCaseSensitive Optional, default = false. Should symbols be treated as case-sensitive.",
    "--increment    Required. Number of random kmers to recruit on each round. A figure of about 100 seems to provide a reasonable trade-off between speed and quality."
                };

                for ( auto s : text ) {
                    if ( s.find( "--" ) != string::npos ) {
                        cerr << "\n";
                    }

                    auto words = String::Split( s, ' ' );

                    int charsWritten = 0;

                    for ( auto & word : words ) {
                        if ( charsWritten + 1 + s.length() <= 80 ) {
                            if ( charsWritten > 0 ) cerr << ' ';
                            cerr << word;
                        }
                        else {
                            cerr << '\n' << word;
                            charsWritten = word.length();
                        }
                    }

                    cerr << "\n";
                }

                ok = false;
                return;
            }

            if ( !arguments->Get( "protoOut", protoOut ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--protoOut' not supplied.\n";
                ok = false;
            }

            if ( !arguments->Get( "fastaFile", fastaFile ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--fastaFile' not supplied.\n";
                ok = false;
            }

            if ( !arguments->Get( "idIndex", idIndex ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--idIndex' not supplied.\n";
                ok = false;
            }

            if ( !arguments->Get( "seed", seed ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--seed' not supplied.\n";
                ok = false;
            }

            if ( !arguments->Get( "increment", this->increment ) ) {
                cerr << arguments->ProgName() << ": note  - optional argument '--increment' not supplied. Using default value '" << increment << "'\n";
            }

            if ( !arguments->Get( "minClusterSize", this->minClusterSize ) ) {
                cerr << arguments->ProgName() << ": note  - optional argument '--minClusterSize' not supplied. Using default value '" << minClusterSize << "'\n";
            }

            if ( !arguments->Get( "threshold", threshold ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--threshold' not provided.\n";
                ok = false;
            }

            if ( !arguments->Get( "clusterOut", clusterOut ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--clusterOut' not provided.\n";
                ok = false;
            }

            if ( !arguments->Get( "treeOut", treeOut ) ) {
                cerr << arguments->ProgName() << ": error - required argument '--treeOut' not provided.\n";
                ok = false;
            }

            int matrixId;

            if ( arguments->IsDefined( "matrixId" ) ) {
                if ( !( arguments->Get( "matrixId", matrixId ) ) ) {
                    cerr << arguments->ProgName() << ": error - argument 'matrixId' not valid." << endl;
                    ok = false;
                }

                vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

                bool found = false;

                for ( auto x : matrices ) {
                    if ( x == matrixId ) { found = true; }
                }

                if ( !found ) {
                    cerr << arguments->ProgName() << ": error - matrix id not recognised." << endl;
                    ok = false;
                }
            }

            if ( !arguments->Get( "numThreads", numThreads ) ) {
                cerr << arguments->ProgName() << ": note  - optional argument '--numThreads' not set; running with default value " << numThreads << ".\n";
            }

            if ( !arguments->Get( "isCaseSensistive", isCaseSensitive ) ) {
                cerr << arguments->ProgName() << ": note  - optional argument '--isCaseSensistive' not set; running with default value " << isCaseSensitive << ".\n";
            }

            if ( !arguments->Get( "wordLength", wordLength ) ) {
                cerr << arguments->ProgName() << ": note  - optional argument '--wordLength' not set; running with default value " << wordLength << ".\n";
            }

            if ( !ok ) {
                cerr << "Invalid command line arguments supplied. For help, run: AAClust --help\n";
                return;
            }

            string matrixFile;
            distanceType = DistanceType::BlosumDistance();

            if ( arguments->IsDefined( "matrixFile" ) ) {
                arguments->Get( "matrixFile", matrixFile );
                distanceType = DistanceType::Custom();
                matrixId = -1;
            }

            bool isCaseSensitive = true;

            if ( arguments->IsDefined( "isCaseSensitive" ) ) {
                if ( !arguments->Get( "isCaseSensitive", isCaseSensitive ) ) {
                    cerr << arguments->ProgName() << ": error - Invalid data for argument 'isCaseSensitive'." << "\n";
                    ok = false;
                }
            }

            matrix = SimilarityMatrix::GetMatrix( distanceType, matrixId, matrixFile, isCaseSensitive );

            if ( !matrix ) {
                cerr << arguments->ProgName() << ": error - unable to construct similarity matrix.\n";
                cerr << arguments->ProgName() << ": error - you need to supply either matrixId or matrixFile arguments.\n";
                ok = false;
            }
        }
    };
};

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

int main( int argc, char *argv[] ) {
    try {
        Args args( argc, argv );

        arguments = &args;

        double start_time = omp_get_wtime();
        int retCode = AACovTree::Run();
        double end_time = omp_get_wtime();

        cout << "Elapsed time: " << ( end_time - start_time ) << "s" << endl;

        return retCode;
    }
    catch ( Exception ex ) {
        cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
        return 1;
    }
}


