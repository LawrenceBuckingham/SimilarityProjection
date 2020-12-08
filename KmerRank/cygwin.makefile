BIN=bin
DEST=$(BIN)/cygwin

TARGETS= \
		$(DEST)/KmerRank \
		$(DEST)/GetRandomSubsetFasta \
		$(DEST)/AAClustGetRandomPrototypes \
		$(DEST)/AAD2 \
		$(DEST)/AAD2NoIndex \
		$(DEST)/AAD2DB \
		$(DEST)/AAClust \
		$(DEST)/AAClustDB \
		$(DEST)/AAClusterFirst \
		$(DEST)/AAClustRecluster \
		$(DEST)/AAClustSig \
		$(DEST)/AAClustSigEncode \
		$(DEST)/AAClustCountIndices \
		$(DEST)/AASP \
		$(DEST)/AASPDB \
		$(DEST)/SimProjDP

FLAGS=	-std=gnu++14 \
		-I $(SIG) \
		-I $(SIG)/HBFL \
		-I $(SIG)/HBFL/HBGraph \
		-I $(SIG)/Lmvq \
		-g \
		-D 'DEFAULT_THREADS=8' \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-D 'PARANOID_PACK=1' \
		-mpopcnt \
		-fopenmp \
		-W \
		-Werror \
		-Wno-unknown-pragmas \
		-fno-omit-frame-pointer \
		-lfltk \
		-D 'USE_OMP=1'

DBFLAGS=	-std=gnu++14 \
		-I $(SIG) \
		-I $(SIG)/HBFL \
		-I $(SIG)/HBFL/HBGraph \
		-I $(SIG)/Lmvq \
		-g \
		-D 'DEFAULT_THREADS=8' \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-D 'PARANOID_PACK=1' \
		-mpopcnt \
		-fopenmp \
		-W \
		-Werror \
		-Wno-unknown-pragmas \
		-fno-omit-frame-pointer \
		-lfltk \
		-D 'USE_OMP=0'

FLAGS_DB=	-std=c++14 \
		-I ../Include \
		-I ../Include/Lmvq \
		-g \
		-D 'DEFAULT_THREADS=8' \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-mpopcnt \
		-fopenmp \
		-Wall \
		-Werror \
		-Wno-unknown-pragmas \
		-fno-omit-frame-pointer

FLAGS_NO_OMP=	-std=c++14 \
		-I ../Include \
		-I ../Include/Lmvq \
		-g \
		-D 'DEFAULT_THREADS=8' \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-mpopcnt \
		-Wall \
		-Werror \
		-Wno-unknown-pragmas \
		-fno-omit-frame-pointer

SIG=../Include

FREQUENT=$(SIG)/Alphabet.hpp \
		$(SIG)/Args.hpp \
		$(SIG)/Array.hpp \
		$(SIG)/Assert.hpp \
		$(SIG)/AveragePrecision.hpp \
		$(SIG)/BestHitDetector.hpp \
		$(SIG)/CharMap.hpp \
		$(SIG)/Console.hpp \
		$(SIG)/Constants.h \
		$(SIG)/CsvIO.hpp \
		$(SIG)/db.hpp \
		$(SIG)/Delegates.hpp \
		$(SIG)/DiagonalGenerator.hpp \
		$(SIG)/DiscreteDistribution.hpp \
		$(SIG)/DistanceType.hpp \
		$(SIG)/Distribution.hpp \
		$(SIG)/EncodedKmer.hpp \
		$(SIG)/EnumBase.hpp \
		$(SIG)/Exception.hpp \
		$(SIG)/FastaSequence.hpp \
		$(SIG)/FileUtil.hpp \
		$(SIG)/FragDistProcessor.hpp \
		$(SIG)/Fragment.hpp \
		$(SIG)/FragmentAggregationMode.hpp \
		$(SIG)/FreeList.hpp \
		$(SIG)/HausdorffCalculator.hpp \
		$(SIG)/Random.hpp \
		$(SIG)/Homologs.hpp \
		$(SIG)/IArrayParser.hpp \
		$(SIG)/Kmer.hpp \
		$(SIG)/KmerCluster.hpp \
		$(SIG)/KmerClusterPrototype.hpp \
		$(SIG)/KmerCodebook.hpp \
		$(SIG)/KmerDistanceCache.hpp \
		$(SIG)/KmerDistanceCalculator.hpp \
		$(SIG)/KmerDistributions.hpp \
		$(SIG)/KmerIndex.hpp \
		$(SIG)/KmerSequenceRanker.hpp \
		$(SIG)/KmerSequenceRanker_Params.hpp \
		$(SIG)/kNearestNeighbours.hpp \
		$(SIG)/LookupTable.hpp \
		$(SIG)/NormalDistribution.hpp \
		$(SIG)/OmpTimer.h \
		$(SIG)/PackedArray.hpp \
		$(SIG)/Projector.hpp \
		$(SIG)/ProjectorBitEmbedding.hpp \
		$(SIG)/ProjectorSlice.hpp \
		$(SIG)/Ranking.hpp \
		$(SIG)/Selector.hpp \
		$(SIG)/SequenceDistanceFunction.hpp \
		$(SIG)/SimilarityMatrix.hpp \
		$(SIG)/String.hpp \
		$(SIG)/Substring.hpp \
		$(SIG)/Types.hpp \
		$(SIG)/Util.hpp \

all: $(TARGETS)

$(DEST)/AAClustGetRandomPrototypes: AAClustGetRandomPrototypes.cpp \
		$(SIG)/FastaSequence.hpp \
		$(SIG)/Kmer.hpp \
		$(SIG)/EncodedKmer.hpp \
		$(SIG)/Args.hpp \
		$(SIG)/DataLoader.hpp
	g++ AAClustGetRandomPrototypes.cpp \
		-o $@ \
		-D USE_OMP=1  \
		$(FLAGS) \
		-O3

$(DEST)/GetRandomSubsetFasta: GetRandomSubsetFasta.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h \
	$(SIG)/DataLoader.hpp \
	$(SIG)/kNearestNeighbours.hpp
	g++ GetRandomSubsetFasta.cpp \
		-o $@ \
		-D USE_OMP=1  \
		$(FLAGS) \
		-O3

$(DEST)/AAD2: AAD2.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h \
	$(SIG)/kNearestNeighbours.hpp \
	$(SIG)/Simproj.hpp
	g++ AAD2.cpp -o $@ \
		$(FLAGS)\
		-O3

$(DEST)/AAD2NoIndex: AAD2NoIndex.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h \
	$(SIG)/kNearestNeighbours.hpp
	g++ AAD2NoIndex.cpp -o $@ \
		$(FLAGS)\
		-O3

$(DEST)/AAD2DB: AAD2.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h \
	$(SIG)/kNearestNeighbours.hpp
	g++ AAD2.cpp -o $@ \
		$(DBFLAGS)\
		-O0

$(DEST)/AASP: AASP.cpp \
		$(SIG)/Args.hpp \
		$(SIG)/BitSet.hpp \
		$(SIG)/FastaSequence.hpp \
		$(SIG)/Alphabet.hpp \
		$(SIG)/Console.hpp \
		$(SIG)/Delegates.hpp \
		$(SIG)/Util.hpp \
		$(SIG)/OmpTimer.h \
		$(SIG)/kNearestNeighbours.hpp \
		$(SIG)/Sequence.hpp \
		$(SIG)/Simproj.hpp
	g++ AASP.cpp -o $@ \
		$(FLAGS)\
		-O3

$(DEST)/AASPDB: AASP.cpp \
		$(SIG)/Args.hpp \
		$(SIG)/BitSet.hpp \
		$(SIG)/Console.hpp \
		$(SIG)/Delegates.hpp \
		$(SIG)/Util.hpp \
		$(SIG)/OmpTimer.h \
		$(SIG)/kNearestNeighbours.hpp \
		$(SIG)/Sequence.hpp \
		$(SIG)/Simproj.hpp
	g++ AASP.cpp -o $@ \
		$(DBFLAGS)\
		-O0

$(DEST)/AAClustRecluster: AAClustRecluster.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerCluster.hpp \
	$(SIG)/KmerClusterPrototype.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustRecluster.cpp -o $@ \
		-D USE_OMP=1 \
		$(FLAGS)

$(DEST)/AAClustSigEncode: AAClustSigEncode.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustSigEncode.cpp -o $@ \
		-D USE_OMP=1 \
		$(FLAGS)

$(DEST)/AAClustSig: AAClustSig.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h \
	$(SIG)/SparseSignature.hpp \
	$(SIG)/kNearestNeighbours.hpp
	g++ AAClustSig.cpp \
		-o $@ \
		-D USE_OMP=1 \
		$(FLAGS) \
		-O3

$(DEST)/AAClustSigDB: AAClustSig.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h \
	$(SIG)/kNearestNeighbours.hpp
	g++ AAClustSig.cpp \
		-o $@ \
		$(FLAGS) \
		-D USE_OMP=0 \
		-O0

$(DEST)/AAClustCountIndices: AAClustCountIndices.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustCountIndices.cpp \
		-o $@ \
		-D USE_OMP=1 \
		$(FLAGS) \
		-O3

$(DEST)/KmerRank: KmerRank.cpp KmerRank.hpp $(FREQUENT)
	g++ KmerRank.cpp \
		-o $@ \
		-D PARANOID_CODEBOOK=1 \
		-D USE_OMP=1 \
		$(FLAGS) \
		-O3

$(DEST)/KmerRank_diag: KmerRank.cpp KmerRank.hpp  $(FREQUENT)
	g++ KmerRank.cpp \
		-o $@ \
		-O3 \
		-D USE_OMP=1 \
		-D WANT_DIAGNOSTIC_STREAM=1 \
		$(FLAGS)

$(DEST)/AAClust: AAClust.cpp  $(FREQUENT) $(SIG)/SubstitutionMatrix.hpp $(SIG)/DataLoader.hpp
	g++ AAClust.cpp \
		-o $@ \
		$(FLAGS) \
		-O3 

$(DEST)/AAClustDB: AAClust.cpp  $(FREQUENT) $(SIG)/DataLoader.hpp
	g++ AAClust.cpp \
		-o $@ \
		$(DBFLAGS)

$(DEST)/AAClusterFirst: AAClusterFirst.cpp  $(FREQUENT)
	g++ AAClusterFirst.cpp \
		-o $@ \
		-O3 \
		$(FLAGS)

clean:
	for f in $(TARGETS); do \
		if [ -f $${f}.exe ]; then rm $${f}.exe; fi; \
	done

rebuild: clean all

test: all
	gdb --args KmerRank \
		-dbFile m:/work-2015-10/wak4.faa \
		-queryFile m:/work-2015-10/wak4_test_951_1000_Query.faa \
		-statsFile q:/Temp/wak4_test_951_1000_Query.k3.f32.all.stats_gcc.csv \
		-rankingFile q:/Temp/wak4_test_951_1000_Query.k3.f32.all.rankings_gcc.csv \
		-kmerLength 3 \
		-kmerMode HausdorffAverage \
		-fragLength 32 \
		-interval 16 \
		-fragMode HausdorffAverage \
		-dist Blosum 
