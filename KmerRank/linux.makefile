BIN=bin
DEST=$(BIN)/linux

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
		$(DEST)/AASPDB

DEPRECATED = \
		AAClustEval \
		AAClusterSample \
		AAClusterSplit \
		AAClustInfo \
		AAClustJaccardMatrix \
		AAClustLearnDomains \
		AAClustSigDB \
		SimProjDP \
		SimProjDP_diag

all: $(TARGETS)

FLAGS=\
		-std=c++14 \
		-I $(SIG) \
		-I $(SIG)/Lmvq \
		-O3 \
		-static \
		-static-libgcc \
		-static-libstdc++ \
		-o $@ \
		-fopenmp \
		-lgomp \
		-lpthread \
		-lrt \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-D 'PARANOID_CODEBOOK=1' \
		-D 'USE_OMP=1' \
		-D 'DEFAULT_THREADS=72' \
		-mpopcnt

DB_FLAGS=\
		-g \
		-std=c++14 \
		-I $(SIG) \
		-I $(SIG)/Lmvq \
		-static \
		-o $@ \
		-fopenmp \
		-lgomp \
		-lpthread \
		-lrt \
		-static-libgcc \
		-static-libstdc++ \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-D 'PARANOID_CODEBOOK=1' \
		-D 'USE_OMP=0' \
		-D 'DEFAULT_THREADS=1' \
		-mpopcnt

SIG=../Include

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

$(DEST)/AAClustLearnDomains: AAClustLearnDomains.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustLearnDomains.cpp -o $@ \
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
	$(SIG)/kNearestNeighbours.hpp
	g++ AAD2.cpp -o $@ \
		-D USE_OMP=1  \
		$(FLAGS) \
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
		-D USE_OMP=0  \
		$(DB_FLAGS)\
		-O0

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
		-D USE_OMP=1  \
		$(FLAGS) \
		-O3

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
		-D USE_OMP=1  \
		$(FLAGS) \
		-O3

$(DEST)/AASPDB: AASP.cpp \
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
		$(DB_FLAGS) \
		-O0

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

$(DEST)/AAClustCountIndices: AAClustCountIndices.cpp \
	$(SIG)/Args.hpp \
	$(SIG)/BitSet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/EncodedFastaSequence.hpp \
	$(SIG)/DataLoader.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustCountIndices.cpp \
		-o $@ \
		-D USE_OMP=1 \
		$(FLAGS) \
		-O3

$(DEST)/AAClustSigEncode: AAClustSigEncode.cpp  \
	$(SIG)/Args.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerCluster.hpp \
	$(SIG)/KmerClusterPrototype.hpp \
	$(SIG)/KmerCodebook.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/EncodedFastaSequence.hpp \
	$(SIG)/DataLoader.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustSigEncode.cpp \
		$(FLAGS)

$(DEST)/AAClustJaccardMatrix: AAClustJaccardMatrix.cpp  \
	$(SIG)/Args.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustJaccardMatrix.cpp \
		$(FLAGS)

$(DEST)/AAClustRecluster: AAClustRecluster.cpp  \
		$(SIG)/Alphabet.hpp \
		$(SIG)/Args.hpp \
		$(SIG)/Console.hpp \
		$(SIG)/DataLoader.hpp \
		$(SIG)/Delegates.hpp \
		$(SIG)/EncodedFastaSequence.hpp \
		$(SIG)/FastaSequence.hpp \
		$(SIG)/KmerCluster.hpp \
		$(SIG)/KmerClusterPrototype.hpp \
		$(SIG)/KmerCodebook.hpp \
		$(SIG)/KmerDistanceCache.hpp \
		$(SIG)/OmpTimer.h \
		$(SIG)/SimilarityMatrix.hpp \
		$(SIG)/Util.hpp
	g++ AAClustRecluster.cpp \
		$(FLAGS)

$(DEST)/AAClustSig: AAClustSig.cpp  \
	$(SIG)/Args.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerCluster.hpp \
	$(SIG)/KmerClusterPrototype.hpp \
	$(SIG)/KmerCodebook.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/EncodedFastaSequence.hpp \
	$(SIG)/DataLoader.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/SparseSignature.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustSig.cpp \
		$(FLAGS) \
		-D USE_OMP=1 \
		-O3

$(DEST)/AAClustSigDB: AAClustSig.cpp  \
	$(SIG)/Args.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerCluster.hpp \
	$(SIG)/KmerClusterPrototype.hpp \
	$(SIG)/KmerCodebook.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/SparseSignature.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClustSig.cpp \
		$(FLAGS) \
		-D USE_OMP=0 \
		-O0

$(DEST)/AAClusterFirst: AAClusterFirst.cpp  \
	$(SIG)/Args.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerCluster.hpp \
	$(SIG)/KmerClusterPrototype.hpp \
	$(SIG)/KmerCodebook.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/OmpTimer.h
	g++ AAClusterFirst.cpp \
		$(FLAGS)

$(DEST)/SimProjDP: SimProjDP.cpp  \
	$(SIG)/Args.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/SimilarityMatrix.hpp \
	$(SIG)/KmerDistanceCache.hpp \
	$(SIG)/KmerDistanceCalculator.hpp \
	$(SIG)/KmerSequenceRanker.hpp \
	$(SIG)/Alphabet.hpp \
	$(SIG)/Console.hpp \
	$(SIG)/Delegates.hpp \
	$(SIG)/FastaSequence.hpp \
	$(SIG)/FragDistProcessor.hpp \
	$(SIG)/FreeList.hpp \
	$(SIG)/HausdorffCalculator.hpp \
	$(SIG)/BestHitDetector.hpp \
	$(SIG)/ProjectorBitEmbedding.hpp \
	$(SIG)/ProjectorSlice.hpp \
	$(SIG)/Util.hpp \
	$(SIG)/Ranking.hpp \
	$(SIG)/kNearestNeighbours.hpp \
	$(SIG)/OmpTimer.h
	g++ SimProjDP.cpp \
		$(FLAGS)

$(DEST)/KmerRank: KmerRank.cpp KmerRank.hpp $(SIG)/*.hpp
	g++ KmerRank.cpp \
		$(FLAGS)

$(DEST)/KmerRank_diag: KmerRank.cpp KmerRank.hpp $(SIG)/*.hpp
	g++ KmerRank.cpp \
		$(FLAGS) \
		-D 'WANT_DIAGNOSTIC_STREAM=1' \
		-D 'SUBVECTOR_BOUNDS=1'

$(DEST)/AAClust: AAClust.cpp $(SIG)/*.hpp
	g++ AAClust.cpp \
		$(FLAGS)

$(DEST)/AAClustDB: AAClust.cpp $(SIG)/*.hpp
	g++ AAClust.cpp \
		$(DB_FLAGS) \
		-O0 

clean:
	for f in $(TARGETS); do \
		if [ -f $${f} ]; then rm $${f}; fi; \
	done

rebuild: clean all

