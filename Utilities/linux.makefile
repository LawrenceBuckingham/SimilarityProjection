# Makefile for Linux port. These programs are linked statically because I did not have access to an 
# up-to-date C++ compiler on the target server. I was successful compiling on UBUNTU and running
# the resulting statically linked executables on RHEL 7. If one has a better version of gcc on the
# target, use that instead so you can dynamically link.
#
# ~ Lawrence.

DST=bin/linux

TARGETS= \
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
	$(DST)/trec_eval_tc \
	$(DST)/trec_eval_tc_compact \
	$(DST)/trec_eval_precision_at \
	$(DST)/SplitFastaHomologs \
	$(DST)/DomainKMedoids \
	$(DST)/GetLargestProtosByClass \
	$(DST)/SelectProtosMinHash

all: $(TARGETS)

clean:
	for f in $(TARGETS); do \
		if [ -f $${f} ]; then rm $${f}; fi; \
	done

rebuild: clean all

FLAGS=	-g \
		-O3 \
		-I $(INCLUDE) \
		-I $(INCLUDE)/LBFL \
		-I $(INCLUDE)/LBFL/LBGraph \
		-I $(INCLUDE)/Lmvq \
		-D alloca=__builtin_alloca \
		-D POPCOUNT=__builtin_popcountll \
		-D USE_OMP=1 \
		-D DO_OMP_TIMER=1 \
		-fopenmp \
		-lgomp \
		-lpthread \
		-lrt \
		-lfltk \
		-static-libgcc \
		-static-libstdc++ \
		-mpopcnt

INCLUDE=../Include

$(DST)/trec_eval_tc: trec_eval_tc.cpp
	g++ trec_eval_tc.cpp -o $@ $(FLAGS)

$(DST)/trec_eval_tc_compact: trec_eval_tc_compact.cpp
	g++ trec_eval_tc_compact.cpp -o $@ $(FLAGS)

$(DST)/trec_eval_precision_at: trec_eval_precision_at.cpp
	g++ trec_eval_precision_at.cpp -o $@ $(FLAGS) 

$(DST)/SelectProtosMinHash: SelectProtosMinHash.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ SelectProtosMinHash.cpp $(FLAGS) -o $@

$(DST)/GetLargestProtosByClass: GetLargestProtosByClass.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ GetLargestProtosByClass.cpp $(FLAGS) -o $@

$(DST)/DomainKMedoids: DomainKMedoids.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ DomainKMedoids.cpp $(FLAGS) -o $@

$(DST)/SplitFastaHomologs: SplitFastaHomologs.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ SplitFastaHomologs.cpp $(FLAGS) -o $@

$(DST)/GetDomainCoverage: GetDomainCoverage.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ GetDomainCoverage.cpp $(FLAGS) -o $@

$(DST)/KmerKMedoids: KmerKMedoids.cpp \
		$(INCLUDE)/Args.hpp \
		$(INCLUDE)/Domain.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/KmerDistanceCache.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Random.hpp
	g++ KmerKMedoids.cpp $(FLAGS) -o $@

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

$(DST)/AdHoc: *.cpp *.hpp
	g++ Main.cpp -o $@ $(FLAGS)


$(DST)/GetKmerDistanceDistribution: Main.cpp *.hpp
	g++ Main.cpp -D Build_GetKmerDistanceDistribution=1 -o $@ $(FLAGS)

$(DST)/GetKmerTheoreticalDistanceDistributions: GetKmerTheoreticalDistanceDistributions.cpp
	g++ GetKmerTheoreticalDistanceDistributions.cpp $(FLAGS) -o $@

$(DST)/GetDnaHammingDistribution: GetDnaHammingDistribution.cpp
	g++ GetDnaHammingDistribution.cpp -o $@ $(FLAGS)

$(DST)/GetCdfInverse: GetCdfInverse.cpp
	g++ GetCdfInverse.cpp -o $@ $(FLAGS)

