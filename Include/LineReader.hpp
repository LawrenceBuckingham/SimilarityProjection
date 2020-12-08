#pragma once

#include <istream>
#include <string>

using namespace std;

namespace QutBio {
	class LineReader {
	private:
		string currentLine;
		istream& reader;

	public:
		LineReader( istream& reader ) : reader(reader) {
			Advance();
		}

		virtual ~LineReader() {}

		const string& CurrentLine() const { return currentLine; }

		bool Advance() {
			if ( ! Ok() ) return false;
			std::getline( reader, currentLine );
			return true;
		}

		bool Ok () {
			return !(reader.fail() || reader.bad() || reader.eof());
		}
	};
}