#pragma once

#include "Alphabet.hpp"
#include "Args.hpp"
#include "DistanceType.hpp"
#include "SimilarityMatrix.hpp"

namespace QutBio {
	class AlphabetHelper {
	public:

		/** Deciphers the arguments used to define an alphabet and distance measure.
		 * Possibilities are:
		 *
		 *	alphabet in { DNA, RNA, AA, DEFAULT }:
		 *	-	the corresponding standard alphabet will be returned.
		 *	-	For AA, you need to supply a matrixId or a matrixFile.
		 *	-	distance returned will be UngappedEdit regardless of other values.
		 *			If you want to do something more elaborate, set up a Similarity Matrix in a
		 *			file and use a custom alphabet instead.
		 *
		 *	alphabet not in the list above:
		 *	-	Then the string following the --alphabet argument will become the symbols
		 *		of a new custom alphabet.
		 *
		 *	matrixId:
		 *		if this is present, we know we're using one of the BLOSUM matrices, so
		 *			the alphabet returned will be AA,
		 *			the appropriate BLOSUM matrix will be returned.
		 *
		 *			dist == "HalperinEtAl" or dist == "BlosumDistance" or dist == "UngappedEdit":
		 *				the value is recorded without any special further action.
		 *			dist not present: 
		 *				BlosumDistance will be returned.
		 *			
		 *		Otherwise:
		 *			the argument matrixFile must be present.
		 *			The custom similarity matrix is loaded from the file and will be returned.
		 *			The value of distance returned will be BlosumDistance. If you have actual
		 *				distances, negate them to produce similarities, so that BlosumDistance will
		 *				reverse the negation and produce positive distances. It is what it is.
		 */

		static bool GetAlphabetAndMatrix(
			Args & arguments, 
			pAlphabet & alphabet, 
			pSimilarityMatrix & matrix,
			pDistanceType &distance,
			ostream & err
		) {
			bool OK = true;

			if (arguments.IsDefined("alphabet")) {
				string alphaName;

				if (!arguments.Get("alphabet", alphaName)) {
					err << "Unable to parse argument 'alphabet'." << endl;
					OK = false;
				}

				alphabet = Alphabets::ByName( alphaName );

				if (alphabet != Alphabets::AA()) {
					matrix = 0;
					distance = DistanceType::UngappedEdit();
					return OK;
				}
			}
			else {
				throw Exception("Argument --alphabet required.", FileAndLine);
			}
			
			if (arguments.IsDefined("matrixId")) {
				int matrixId = 0;

				if (!(arguments.Get("matrixId", matrixId))) {
					err << "Argument 'matrixId' not valid." << endl;
					OK = false;
				}

				vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for (auto x : matrices) {
					if (x == matrixId) { found = true; }
				}

				if (!found) {
					matrix = 0;
					err << "Matrix id not recognised." << endl;
					OK = false;
				}
				else {
					matrix = SimilarityMatrix::GetBlosum( matrixId );

					if ( arguments.IsDefined( "distance" ) ) {
						if (!arguments.Get("distance", DistanceType::Values(), distance)) {
							err << "Argument 'distance' is not valid." << endl;
							OK = false;
						}
					}
				}
			}

			else if (arguments.IsDefined("matrixFile")) {
				string matrixFile;
				bool isCaseSensitive = true;
				arguments.Get("matrixFile", matrixFile);
				
				if (!arguments.GetOptionalArgument("isCaseSensitive", isCaseSensitive)) {
					err << "Argument 'isCaseSensitive', if supplied, must be true or false." << endl;
					OK = false;
				}
				
				matrix = SimilarityMatrix::GetMatrix( alphabet, DistanceType::Custom(), -1, matrixFile);
				distance = DistanceType::BlosumDistance();
			}

			else {
				err << "Must have either 'alphabet', 'matrixId', or 'matrixFile' defined in arguments." << endl;
				OK = false;
			}

			return OK;
		}
	};
}
