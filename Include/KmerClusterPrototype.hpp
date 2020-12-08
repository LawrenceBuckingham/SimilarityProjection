#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <FastaSequence.hpp>
#include <Kmer.hpp>

namespace QutBio {
    /**
    ** <summary>
    **	An abuse of the FASTA format which allows a sequence
    **	which consists of a single kmer to be managed, along
    **	with additional properties, notably the size of a cluster.
    **	I might add other things later as needed.
    **
    **	<para>
    **	A prototype represents the centre of a kmer cluster that contains
    **	elements accumulated across sequences in multiple files. An
    **	individual cluster file will index kmers from multiple sequences
    **	in a single file. Prototypes unite clusters across files.
    **	</para>
    **
    ** </summary>
    */
    class KmerClusterPrototype : public EncodedFastaSequence {
    protected:
        size_t size = 0;
        size_t serialNumber = 0;
        Kmer thisKmer;

        const string & idPrefix() {
            static string s = "proto_";
            return s;
        }

        static size_t largestSerialNumber( size_t latest = 0 ) {
            static size_t value = 0;

#pragma omp critical
            {
                if ( latest > value ) {
                    value = latest;
                }
            }

            return value;
        }

    public:
        using Type = KmerClusterPrototype;
        using Pointer = Type * ;

        KmerClusterPrototype(
            FastaSequence * baseSequence,
            size_t classIndex,
            Alphabet * alphabet,
            size_t wordLength,
            size_t charsPerWord,
            Symbol defaultSymbol
        ) :
            EncodedFastaSequence( baseSequence, classIndex, alphabet, wordLength, charsPerWord, defaultSymbol ),
            size( 0 ),
            thisKmer( *this, 0, wordLength ) {

            istringstream idStream( IdStr().substr( idPrefix().length() ) );
            idStream >> serialNumber;

            largestSerialNumber( serialNumber );

			bool done = false;

            for ( uint i = 0; i < baseSequence->MetaCount(); i++ ) {
				const auto & meta = baseSequence->Metadata(i);

                if ( done ) break;

                if ( meta.find( "size=" ) != string::npos ) {
                    auto properties = String::Split( meta, "=" );
                    bool sizeFound = false;

                    for ( auto & p : properties ) {
                        if ( p == "size" ) {
                            sizeFound = true;
                        }
                        else if ( sizeFound ) {
                            istringstream( p ) >> size;
                            done = true;
                            break;
                        }
                    }
                }
            }
        }

        KmerClusterPrototype(
            size_t serialNumber,
            FastaSequence * kmerWord,
            Alphabet * alphabet,
            size_t wordLength,
            size_t charsPerWord,
            Symbol defaultSymbol
        ) :
            EncodedFastaSequence( kmerWord, -1, alphabet, wordLength, charsPerWord, defaultSymbol ),
            size( 0 ),
            serialNumber( serialNumber ),
            thisKmer( *this, 0, wordLength ) {

            UpdateDefLine();
            largestSerialNumber( serialNumber );
        }

        KmerClusterPrototype(
            FastaSequence * kmerWord,
            Alphabet * alphabet,
            size_t wordLength,
            size_t charsPerWord,
            Symbol defaultSymbol
        ) :
            EncodedFastaSequence( kmerWord, -1, alphabet, wordLength, charsPerWord, defaultSymbol ),
            size( 0 ),
            thisKmer( *this, 0, wordLength ) {
            
            largestSerialNumber( largestSerialNumber() + 1 );
            serialNumber = largestSerialNumber();
            UpdateDefLine();
        }

        /// <summary>Get the (total) size of the cluster (s) represented by this prototype.</summary>
        size_t Size() const {
            return size;
        }

        /// <summary>Set the size of the combined collection of indexed kmers.</summary>
        KmerClusterPrototype & Size( size_t size ) {
            this->size = size;
            UpdateDefLine();
            return *this;
        }

        /// <summary>Make the definition line consistent with the id and size.</summary>
        void UpdateDefLine() {
            ostringstream s;
            s << idPrefix() << serialNumber << "|size=" << size;
            base->SetDefLine( s.str() );
        }

        /**
        **	<summary>
        **	Get the singleton kmer represented by this prototype.
        **	</summary>
        */

        Kmer * SingletonKmer() {
            return &thisKmer;
        }

        EncodedKmer PackedEncoding() {
            return GetEncodedKmer( 0 );
        }

        /**
        **	<summary>
        **	Get an ID string for use when persisting to a FASTA file.
        **	</summary>
        */

        string GetId( size_t serialNumber ) {
            return idPrefix() + std::to_string( serialNumber );
        }
    };

    typedef KmerClusterPrototype * pKmerClusterPrototype;

}
