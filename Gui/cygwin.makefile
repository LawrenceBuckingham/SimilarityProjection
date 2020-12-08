BIN=bin
DEST=$(BIN)/cygwin

FEEDBACK_TARGETS=\
	$(DEST)/BrowseRankingsFB \
	$(DEST)/BowGuiFB \
	$(DEST)/DnaClustGuiFB \
	$(DEST)/SimProjGuiFB \
	$(DEST)/SimProjIndexedGuiFB
feedback: $(FEEDBACK_TARGETS)

RELEASE_TARGETS=\
	$(DEST)/BrowseRankings \
	$(DEST)/BowGui \
	$(DEST)/DnaClustGui \
	$(DEST)/SimProjGui \
	$(DEST)/SimProjIndexedGui
release: $(RELEASE_TARGETS)

DEBUG_TARGETS=\
	$(DEST)/DnaClustGuiDB \
	$(DEST)/BrowseRankingsDB \
	$(DEST)/BowGuiDB \
	$(DEST)/SimProjGuiDB \
	$(DEST)/SimProjIndexedGuiDB
debug: $(DEBUG_TARGETS)

TARGETS =  $(FEEDBACK_TARGETS) $(RELEASE_TARGETS) $(DEBUG_TARGETS)

all: $(TARGETS)

INCLUDE=../Include

DEPS= \
		*.hpp \
		$(INCLUDE)/Assert.hpp \
		$(INCLUDE)/Alphabet.hpp \
		$(INCLUDE)/Centroid.hpp \
		$(INCLUDE)/CsvIO.hpp \
		$(INCLUDE)/FastaSequence.hpp \
		$(INCLUDE)/Index.hpp \
		$(INCLUDE)/KNearestNeighbours.hpp \
		$(INCLUDE)/Registry.hpp \
		$(INCLUDE)/SimilarityMatrix.hpp \
		$(INCLUDE)/SimpleKmer.hpp \
		$(INCLUDE)/SparseFeatureVector.hpp \
		$(INCLUDE)/SparseSignature.hpp \
		$(INCLUDE)/String.hpp \
		$(INCLUDE)/Util.hpp \
		$(INCLUDE)/LBFL/*.hpp \
		$(INCLUDE)/LBFL/LBGraph/*.hpp


FLAGS=	-std=gnu++14 \
		-I $(INCLUDE) \
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

DBFLAGS=	-std=gnu++14 \
		-I $(INCLUDE) \
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
		-D 'USE_OMP=0'

$(DEST)/DnaClustGui: DnaClustGui.cpp $(DEPS) 
	g++ DnaClustGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=0' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/DnaClustGuiFB: DnaClustGui.cpp $(DEPS) 
	g++ DnaClustGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

# profile flags. -pg \

$(DEST)/DnaClustGuiDB: DnaClustGui.cpp  $(DEPS) 
	g++ DnaClustGui.cpp -o $@ -g $(DBFLAGS) -Wa,-mbig-obj -O0 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/BowGui: BowGui.cpp $(DEPS) 
	g++ BowGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=0' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/BowGuiFB: BowGui.cpp $(DEPS) 
	g++ BowGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/BowGuiDB: BowGui.cpp  $(DEPS) 
	g++ BowGui.cpp  -o $@ -g $(DBFLAGS) -Wa,-mbig-obj -O0 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/SimProjGui: SimProjGui.cpp $(DEPS) 
	g++ SimProjGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=0' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/SimProjGuiFB: SimProjGui.cpp $(DEPS) 
	g++ SimProjGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/SimProjGuiDB: SimProjGui.cpp  $(DEPS) 
	g++ SimProjGui.cpp  -o $@ -g $(DBFLAGS) -Wa,-mbig-obj -O0 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/SimProjIndexedGui: SimProjIndexedGui.cpp $(DEPS) 
	g++ SimProjIndexedGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=0' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/SimProjIndexedGuiFB: SimProjIndexedGui.cpp $(DEPS) 
	g++ SimProjIndexedGui.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/SimProjIndexedGuiDB: SimProjIndexedGui.cpp  $(DEPS) 
	g++ SimProjIndexedGui.cpp  -o $@ -g $(DBFLAGS) -Wa,-mbig-obj -O0 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_INCLUDE=1'

$(DEST)/BrowseRankings: BrowseRankings.cpp $(DEPS) 
	g++ BrowseRankings.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=0'

$(DEST)/BrowseRankingsFB: BrowseRankings.cpp $(DEPS) 
	g++ BrowseRankings.cpp -o $@ -g $(FLAGS) -Wa,-mbig-obj -O3 -D 'WANT_FEEDBACK=1'

$(DEST)/BrowseRankingsDB: BrowseRankings.cpp  $(DEPS) 
	g++ BrowseRankings.cpp  -o $@ -g $(DBFLAGS) -Wa,-mbig-obj -O0 -D 'WANT_FEEDBACK=1'

clean:
	for f in $(TARGETS); do \
		if [ -f $${f}.exe ]; then rm $${f}.exe; fi; \
	done

rebuild: clean all

