BIN=bin
DEST=$(BIN)/linux

RELEASE = \
	$(DEST)/BrowseRankings \
	$(DEST)/BowGui \
	$(DEST)/DnaClustGui \
	$(DEST)/SimProjGui \
	$(DEST)/SimProjIndexedGui

FEEDBACK = \
	$(DEST)/BrowseRankingsFB \
	$(DEST)/BowGuiFB \
	$(DEST)/DnaClustGuiFB \
	$(DEST)/SimProjGuiFB

DEBUG = \
	$(DEST)/BrowseRankingsDB \
	$(DEST)/BowGuiDB \
	$(DEST)/DnaClustGuiDB \
	$(DEST)/SimProjGuiDB

TARGETS = $(RELEASE) $(FEEDBACK) $(DEBUG)

all: $(TARGETS)

release: $(RELEASE)

feedback: $(FEEDBACK)

debug: $(DEBUG)

SIG=../Include

FLAGS=	-std=gnu++14 \
		-I $(SIG) \
		-I $(SIG)/LBFL \
		-I $(SIG)/LBFL/LBGraph \
		-I $(SIG)/Lmvq \
		-g \
		-D 'DEFAULT_THREADS=8' \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-D 'PARANOID_PACK=1' \
		-mpopcnt \
		-W \
		-Werror \
		-Wno-unknown-pragmas \
		-fno-omit-frame-pointer \
		-lfltk \
		-D 'USE_OMP=1' \
		-fopenmp \
		-lgomp \
		-lpthread \
		-lrt \
		-static-libgcc \
		-static-libstdc++

DBFLAGS=	-std=gnu++14 \
		-I $(SIG) \
		-I $(SIG)/LBFL \
		-I $(SIG)/LBFL/LBGraph \
		-I $(SIG)/Lmvq \
		-g \
		-D 'DEFAULT_THREADS=8' \
		-D 'alloca=__builtin_alloca' \
		-D 'POPCOUNT=__builtin_popcountll' \
		-D 'DO_OMP_TIMER=1' \
		-D 'PARANOID_PACK=1' \
		-mpopcnt \
		-W \
		-Werror \
		-Wno-unknown-pragmas \
		-fno-omit-frame-pointer \
		-lfltk \
		-D 'USE_OMP=0' \
		-fopenmp \
		-lgomp \
		-lpthread \
		-lrt \
		-static-libgcc \
		-static-libstdc++

DEPS = *.hpp \
		$(SIG)/*.hpp \
		$(SIG)/LBFL/*.hpp \
		$(SIG)/LBFL/LBGraph/*.hpp

$(BIN):
	if [ ! -d $(BIN) ] ; then mkdir $(BIN); fi

$(DEST): $(BIN)
	if [ ! -d $(DEST) ] ; then mkdir $(DEST); fi

$(DEST)/DnaClustGui: DnaClustGui.cpp $(DEPS)	
	g++ DnaClustGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=0' -D 'INTERLEAVE_SIG=1'

$(DEST)/DnaClustGuiFB: DnaClustGui.cpp $(DEPS)	
	g++ DnaClustGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_SIG=1'

$(DEST)/DnaClustGuiDB: DnaClustGui.cpp $(DEPS)
	g++ DnaClustGui.cpp -o $@ -g $(DBFLAGS) -O0 -D 'WANT_FEEDBACK=1' -D 'INTERLEAVE_SIG=1'

$(DEST)/BowGui: BowGui.cpp $(DEPS)	
	g++ BowGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=0'
	
$(DEST)/BowGuiFB: BowGui.cpp $(DEPS)	
	g++ BowGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=1'

$(DEST)/BowGuiDB: BowGui.cpp $(DEPS)
	g++ BowGui.cpp -o $@ -g $(DBFLAGS) -O0 -D 'WANT_FEEDBACK=1'

$(DEST)/SimProjGui: SimProjGui.cpp $(DEPS)	
	g++ SimProjGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=0'

$(DEST)/SimProjGuiFB: SimProjGui.cpp $(DEPS)	
	g++ SimProjGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=1'

$(DEST)/SimProjGuiDB: SimProjGui.cpp $(DEPS)
	g++ SimProjGui.cpp -o $@ -g $(DBFLAGS) -O0 -D 'WANT_FEEDBACK=1'

$(DEST)/SimProjIndexedGui: SimProjIndexedGui.cpp $(DEPS)	
	g++ SimProjIndexedGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=0'

$(DEST)/SimProjIndexedGuiFB: SimProjIndexedGui.cpp $(DEPS)	
	g++ SimProjIndexedGui.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=1'

$(DEST)/SimProjIndexedGuiDB: SimProjIndexedGui.cpp $(DEPS)
	g++ SimProjIndexedGui.cpp -o $@ -g $(DBFLAGS) -O0 -D 'WANT_FEEDBACK=1'

$(DEST)/BrowseRankings: BrowseRankings.cpp $(DEPS)	
	g++ BrowseRankings.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=0'

$(DEST)/BrowseRankingsFB: BrowseRankings.cpp $(DEPS)	
	g++ BrowseRankings.cpp -o $@ -g $(FLAGS) -O3 -D 'WANT_FEEDBACK=1'

$(DEST)/BrowseRankingsDB: BrowseRankings.cpp $(DEPS)
	g++ BrowseRankings.cpp -o $@ -g $(DBFLAGS) -O0 -D 'WANT_FEEDBACK=1'

clean:
	for f in $(TARGETS); do \
		if [ -f $${f}     ]; then rm $${f}    ; fi; \
	done

rebuild: clean all

