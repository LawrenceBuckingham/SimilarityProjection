#pragma once 

#include "db.hpp"

#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "DnaDistance.hpp"
#include "FastaSequence.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "DataLoader.hpp"

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>

#include <FL/Fl.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Scroll.H>
#include <FL/fl_draw.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Check_Button.H>

using Control = Fl_Widget;

#include "LBFL/BorderLayout.hpp"
#include "LBFL/GridLayout.hpp"
#include "LBFL/FileChooser.hpp"
#include "LBFL/Arg.hpp"
#include "LBFL/TypedArg.hpp"
#include "LBFL/TextDisplay.hpp"

using namespace QutBio;
using namespace std;
using namespace LBFL;

#include "AlphabetChooser.hpp"
#include "Page.hpp"
