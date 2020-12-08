DST=bin/cygwin

TARGETS= \
	$(DST)/BootstrapColumnAverage \
	$(DST)/ComputeGoodnessOfFit \
	$(DST)/GetLargestProtosByClass \
	$(DST)/AdHoc \
	$(DST)/GetDnaHammingDistribution \
	$(DST)/GetCdfInverse \
	$(DST)/GetKmerTheoreticalDistanceDistributions \
	$(DST)/GetKmerDistanceDistribution \
	$(DST)/ImportSwissprotDomains \
	$(DST)/GetDomainSubset \
	$(DST)/GetCodebookSubset \
	$(DST)/ExtractDomainSubsequences \
	$(DST)/KmerKMedoids \
	$(DST)/GetDomainCoverage \
	$(DST)/DomainKMedoids \
	$(DST)/trec_eval_tc \
	$(DST)/trec_eval_tc_compact \
	$(DST)/trec_eval_precision_at \
	$(DST)/SplitFastaHomologs \
	$(DST)/ComputeDistanceMatrix \
	$(DST)/SelectProtosMinHash

all: $(TARGETS)

clean:
	for f in $(TARGETS); do \
		if [ -f $${f}.exe ]; then rm $${f}.exe; fi; \
	done

rebuild: clean all

INCLUDE=../Include

FLAGS=	-I $(INCLUDE) \
		-I $(INCLUDE)/LBFL \
		-I $(INCLUDE)/LBFL/LBGraph \
		-I $(INCLUDE)/Lmvq \
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

$(DST)/BootstrapColumnAverage: BootstrapColumnAverage.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/CsvIO.hpp \
		$(INCLUDE)/Exception.hpp \
		$(INCLUDE)/Random.hpp
	g++ BootstrapColumnAverage.cpp \
		-o $@ \
		$(FLAGS) -O0

$(DST)/ComputeGoodnessOfFit: ComputeGoodnessOfFit.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/OmpTimer.h \
		$(INCLUDE)/Exception.hpp \
		$(INCLUDE)/IntegerDistribution.hpp
	g++ ComputeGoodnessOfFit.cpp \
		-o $@ \
		$(FLAGS) 

$(DST)/ComputeDistanceMatrix: ComputeDistanceMatrix.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/OmpTimer.h \
		$(INCLUDE)/Exception.hpp \
		$(INCLUDE)/FileUtil.hpp \
		$(INCLUDE)/HausdorffCalculator.hpp \
		$(INCLUDE)/EnumBase.hpp
	g++ ComputeDistanceMatrix.cpp \
		-o $@ \
		$(FLAGS) 

$(DST)/SplitFastaHomologs: SplitFastaHomologs.cpp
	g++ SplitFastaHomologs.cpp \
		-o $@ \
		$(FLAGS) 

$(DST)/trec_eval_tc: trec_eval_tc.cpp
	g++ trec_eval_tc.cpp \
		-o $@ \
		-O3 \
		-I ../Include
		-fopenmp -lgomp -lpthread -lrt -static-libgcc -static-libstdc++ \
		-Dalloca=__builtin_alloca \
		-DPOPCOUNT=__builtin_popcountll \
		-mpopcnt 

$(DST)/trec_eval_tc_compact: trec_eval_tc_compact.cpp
	g++ trec_eval_tc_compact.cpp \
		$(FLAGS) \
		-O3 \
		-o $@

$(DST)/trec_eval_precision_at: trec_eval_precision_at.cpp
	g++ trec_eval_precision_at.cpp \
		$(FLAGS) \
		-o $@ 

$(DST)/SelectProtosMinHash: SelectProtosMinHash.cpp \
	$(INCLUDE)/Args.hpp \
	$(INCLUDE)/BitSet.hpp \
	$(INCLUDE)/FastaSequence.hpp \
	$(INCLUDE)/SimilarityMatrix.hpp \
	$(INCLUDE)/KmerDistanceCache.hpp \
	$(INCLUDE)/KmerDistanceCalculator.hpp \
	$(INCLUDE)/KmerSequenceRanker.hpp \
	$(INCLUDE)/kNearestNeighbours.hpp \
	$(INCLUDE)/Alphabet.hpp \
	$(INCLUDE)/Console.hpp \
	$(INCLUDE)/Delegates.hpp \
	$(INCLUDE)/FastaSequence.hpp \
	$(INCLUDE)/Util.hpp \
	$(INCLUDE)/OmpTimer.h
	g++ SelectProtosMinHash.cpp -o $@ $(FLAGS) -D USE_OMP=1

$(DST)/GetLargestProtosByClass: GetLargestProtosByClass.cpp \
	$(INCLUDE)/Args.hpp \
	$(INCLUDE)/BitSet.hpp \
	$(INCLUDE)/FastaSequence.hpp \
	$(INCLUDE)/SimilarityMatrix.hpp \
	$(INCLUDE)/KmerDistanceCache.hpp \
	$(INCLUDE)/KmerDistanceCalculator.hpp \
	$(INCLUDE)/KmerSequenceRanker.hpp \
	$(INCLUDE)/Alphabet.hpp \
	$(INCLUDE)/Console.hpp \
	$(INCLUDE)/Delegates.hpp \
	$(INCLUDE)/FastaSequence.hpp \
	$(INCLUDE)/Util.hpp \
	$(INCLUDE)/OmpTimer.h
	g++ GetLargestProtosByClass.cpp -o $@ \
		-D USE_OMP=1 \
		$(FLAGS)

$(DST)/GetDomainCoverage: GetDomainCoverage.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/KmerIndex.hpp \
		$(INCLUDE)/Substring.hpp \
		$(INCLUDE)/Kmer.hpp \
		$(INCLUDE)/KmerClusterPrototype.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ GetDomainCoverage.cpp $(FLAGS) -O3 -o $@

$(DST)/DomainKMedoids: DomainKMedoids.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/KMedoids.hpp \
		$(INCLUDE)/Random.hpp
	g++ DomainKMedoids.cpp $(FLAGS) -O3 -o $@

$(DST)/KmerKMedoids: KmerKMedoids.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ KmerKMedoids.cpp $(FLAGS) -O3 -o $@

$(DST)/ExtractDomainSubsequences: ExtractDomainSubsequences.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp
	g++ ExtractDomainSubsequences.cpp $(FLAGS) -o $@

$(DST)/GetCodebookSubset: GetCodebookSubset.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerCodebook.hpp \
		$(INCLUDE)/KmerCluster.hpp \
		$(INCLUDE)/KmerClusterPrototype.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp
	g++ GetCodebookSubset.cpp $(FLAGS) -o $@

$(DST)/GetDomainSubset: GetDomainSubset.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp
	g++ GetDomainSubset.cpp $(FLAGS) -o $@

$(DST)/ImportSwissprotDomains: ImportSwissprotDomains.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp
	g++ ImportSwissprotDomains.cpp $(FLAGS) -o $@

$(DST)/AdHoc: Main.cpp *.hpp
	g++ Main.cpp $(FLAGS) -o $@

$(DST)/GetKmerDistanceDistribution: Main.cpp *.hpp 
	g++ Main.cpp -D Build_GetKmerDistanceDistribution=1 -o $@ $(FLAGS)

$(DST)/GetKmerTheoreticalDistanceDistributions: \
		GetKmerTheoreticalDistanceDistributions.cpp \
		$(INCLUDE)/IntegerDistribution.hpp \
		$(INCLUDE)/WeibullDistribution.hpp \
		$(INCLUDE)/NormalDistribution.hpp \
		$(INCLUDE)/Util.hpp \
		$(INCLUDE)/Distribution.hpp
	g++ GetKmerTheoreticalDistanceDistributions.cpp $(FLAGS) -o $@ 

$(DST)/GetDnaHammingDistribution: GetDnaHammingDistribution.cpp
	g++ GetDnaHammingDistribution.cpp  $(FLAGS) -o $@

$(DST)/GetCdfInverse: GetCdfInverse.cpp $(INCLUDE)/Histogram.hpp
	g++ GetCdfInverse.cpp  $(FLAGS) -o $@

