#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <ios>
#include <string>
#include <vector>
#include <limits>
#include <climits>
#include <istream>
#include <cstdint>

#include "Exception.hpp"
#include "Delegates.hpp"
#include "Util.hpp"

namespace QutBio {
	using namespace std;

	class CsvWriter;
	class CsvReader;

	struct ICsvWriter
	{
		virtual ~ICsvWriter() {}
		virtual void Write(CsvWriter &w) const = 0;
	};

	struct ICsvReader
	{
		virtual ~ICsvReader() {}
		virtual void Read(CsvReader &r) = 0;
	};

	/// <summary> CsvIO provides low-level I/O 
	///		functions to parse and generate Comma Separated Value (.CSV) text 
	///		streams as implicitly defined by the behaviour of Microsoft Excel.
	///	<para>
	///		The definition can be found in RFC 4180 - Common Format and MIME 
	///		Type for CSV Files - October 2005:<br />
	///		<a href="http://tools.ietf.org/html/rfc4180">http://tools.ietf.org/html/rfc4180</a>  
	///	</para>
	///	<para>
	///		This is a JavaScript port of source code originally implemented as part of the 
	///		CIEAM HAMISH Reliability Basics product.
	///	</para>
	/// </summary>

	class CsvReader
	{

#pragma region
		/// Internal representation of the separator character.
		int separator;

		/// Internal representation of the escape (quote) character.
		int quoteChar;

		static const int NUL = 0;
		static const int EOF_ = -1;
		static const int CR = 0x0D;
		static const int LF = 0x0A;
		static const int CRLF = -2;
		static const int EOS = -3;

		int currentChar = 0;
		int nextChar = 0;
		istream & reader;
		bool parsingString = false;
		vector<string> * currentRow;
		string * currentField = 0;
		int position = -1;
		int lineNumber = -1;
#pragma endregion Variables and constants used by Parse function.

	public:

		/// <summary> The default field separator: a comma.
		/// </summary>
		static const char DefaultSeparator = ',';

		/// <summary> The default escape character: a double quote.
		/// </summary>
		static const char DefaultQuoteChar = '"';


		/// <summary> Initialises a CsvIo object with default separator and quote character.
		/// </summary>

		CsvReader(istream & reader) : reader(reader)
		{
			ConstructorImpl(DefaultSeparator, DefaultQuoteChar);
		}

		/// <summary> Initialises a CsvIo object with explicit separator and default quote character.
		/// </summary>

		CsvReader(istream & reader, char separator) : reader(reader)
		{
			ConstructorImpl(separator, DefaultQuoteChar);
		}

		/// <summary> Initialises a CsvIO object, setting the field separator and quote character.
		///		<para>Records are always separated by either carriage return (CR) or line feed (LF) characters.</para>
		///		<para><strong>NB:</strong> This is available only if the assembly is referenced directly.</para>
		/// </summary>
		/// <param name="separator">
		///		The character that will be used to separate fields within a record.
		///		<para>
		///			This may not be carriage-return, line-feed or the same as the "quote" character.
		///		</para>
		/// </param>
		/// <param name="quoteChar">
		///		The escape (quote) character that will be used to mark up fields which contain
		///		embedded carriage-return, line-feed, separator or quote characters.
		/// </param>
		/// <exception cref="Exception">
		///		Thrown if invalid separator or quote char are supplied.
		/// </exception>

		CsvReader(istream & reader, char separator, char quoteChar) : reader(reader)
		{
			ConstructorImpl(separator, quoteChar);
		}

		/* Excerpt from RFC 4180       Common Format and MIME Type for CSV Files    October 2005
		 *		http://tools.ietf.org/html/rfc4180
		 *
		 *		The grammar, modified from RFC by removing the optional header, is:
		 *
		 *		file = record *(CRLF record) [CRLF]
		 *		record = field *(COMMA field)
		 *		field = (escaped / non-escaped)
		 *		escaped = DQUOTE *(TEXTDATA / COMMA / CR / LF / 2DQUOTE) DQUOTE
		 *		non-escaped = *TEXTDATA
		 *		COMMA = %x2C
		 *		CR = %x0D ;as per section 6.1 of RFC 2234 [2]
		 *		DQUOTE =  %x22 ;as per section 6.1 of RFC 2234 [2]
		 *		LF = %x0A ;as per section 6.1 of RFC 2234 [2]
		 *		CRLF = CR LF ;as per section 6.1 of RFC 2234 [2]
		 *		TEXTDATA =  %x20-21 / %x23-2B / %x2D-7E
		 */

		 /// <summary> Returns true iff the current character is an end-of-line marker (CR, LF or CRLF).
		 /// </summary>
		bool IsEOL() { return currentChar == CR || currentChar == LF || currentChar == CRLF; }

		/// <summary> Returns true iff the end-of-file has been encountered.
		/// </summary>

		bool IsEOF() { return currentChar == EOF_; }

		/// <summary> Returns true iff the parser is at the beginning of the file.
		/// </summary>

		bool IsBOF() { return currentChar == NUL; }

		/// <summary> Returns true iff the parser is at a separator (by default, COMMA).
		/// </summary>

		bool IsCOMMA() { return currentChar == separator; }

		/// <summary> Returns true iff the parse is at the closing quote that marks the end of a string.
		/// </summary>

		bool IsEOS() { return currentChar == EOS; }

		/// <summary> Parses a CSV-formatted dataset from the supplied text reader.
		///		<para><strong>NB:</strong> This is available only if the assembly is referenced directly.</para>
		/// </summary>
		/// <param name="reader">
		///		A TextReader from which input is obtained.
		/// </param>
		/// <returns>
		///		Returns a list of records, one for each row in the dataset. Each of record 
		///		is an array of string containing one element for each field in the row. 
		/// </returns>

		void Read(vector<vector<string>> & rows, int observations = INT_MAX)
		{
			ReadFile(rows, observations);
		}

		/// <summary> Parses a CSV-formatted dataset from the supplied text reader and invokes the .
		///		<para><strong>NB:</strong> This is available only if the assembly is referenced directly.</para>
		/// </summary>
		/// <param name="reader">
		///		A TextReader from which input is obtained.
		/// </param>
		/// <returns>
		///		Returns a list of records, one for each row in the dataset. Each of record 
		///		is an array of string containing one element for each field in the row. 
		/// </returns>

		template<typename T>
		void Parse(vector<T> & records, int observations = INT_MAX)
		{
			ParseFile(records, observations);
		}

		/// <summary> Parses a CSV-formatted dataset from the supplied text reader and invokes the .
		///		<para><strong>NB:</strong> This is available only if the assembly is referenced directly.</para>
		/// </summary>
		/// <param name="reader">
		///		A TextReader from which input is obtained.
		/// </param>
		/// <returns>
		///		Returns a list of records, one for each row in the dataset. Each of record 
		///		is an array of string containing one element for each field in the row. 
		/// </returns>

		template<typename T> void Stream(Func1<T &, bool> process, Action loadComplete, int observations = INT_MAX)
		{
			StreamFile(process, loadComplete, observations);
		}

		/// <summary> Reads an entire file, which is defined by a production of the form
		///		File = record *(CRLF record) [CRLF].
		///	<para> The function supplied in argument process is called once for each non-empty
		///		record in the stream. processing continues until all records have been read, 
		///		or until process returns false.
		///	</para>
		///	<para>
		///		After record processing is finished, the loadComplete callback is invoked.
		///	</para>
		///	<para>
		///		On entry, one character from the text reader has been consumed and stored in
		///		nextChar while currentChar is uninitialised.
		///	</para>
		///	<para>
		///		On exit, all available characters in the text reader have been consumed and 
		///		currentChar is EOF_.
		///	</para>
		/// </summary>

		void StreamRecords(
			Func1<vector<string> &, bool> process,
			Action loadComplete,
			const unsigned int observations = numeric_limits<unsigned int>::max()
		)
		{
			size_t count = 0;
			bool keep_going = true;

			do
			{
				vector<string> currentRow;
				ReadRecord(currentRow);

				if (currentRow.size() > 1 || currentRow[0].length() > 0)
				{
					keep_going = process(currentRow);
					count++;
				}
			}
			while (keep_going && IsEOL() && (count < observations));

			loadComplete();
		}

		/// <summary> Reads a record, of the form 
		///		Record = Field *(COMMA field)
		///	<para>
		///		On entry, it is assumed that currentChar is either uninitialised or retains
		///		the end-of-line marker which terminated the previous record. Either way, this 
		///		value is discarded. 
		///		As an attempt to continue the parse in the event of malformed input, we skip 
		///		characters until this assumption becomes true.
		///	</para>
		///	<para>
		///		On exit, currentChar is either CRLF or EOF_, currentRow has been set back to null,
		///		and a new row has been added to the rows collection with the fields garnered from 
		///		this row.
		///	</para>
		/// </summary>

		void ReadRecord(vector<string> & currentRow)
		{
			while (!(IsBOF() || IsEOF() || IsEOL())) ReadChar();

			if (IsEOF()) return;

			this->currentRow = &currentRow;

			do ReadField(); while (currentChar == separator);
		}

		/// <summary> Reads a single field, implementing the three productions:
		///		field = (escaped / non-escaped)
		/// 	escaped = DQUOTE *(TEXTDATA / COMMA / CR / LF / 2DQUOTE) DQUOTE
		/// 	non-escaped = *TEXTDATA
		/// <para>
		///		On entry, currentChar should be either NUL, EOLN or COMMA. That is,
		///		we are about to process a newly encountered field. If this is not 
		///		the case, we skip forward until EOF_ or one of these is found to e
		/// </para>
		/// </summary>

		void ReadField(string & field)
		{
			while (!(IsBOF() || IsEOF() || IsEOL() || IsCOMMA())) ReadChar();

			if (IsEOF()) return;

			currentField = &field;
			currentField->clear();

			if (nextChar == quoteChar) ReadEscaped(); else ReadNonEscaped();

			currentField = 0;
		}

	private:

		void ConstructorImpl(char separator, char quoteChar)
		{
			if (separator == CR || separator == LF || separator == quoteChar)
			{
				throw Exception("separator may not be CR, LF or the escape (quote) character.", __FILE__, __LINE__);
			}

			if (quoteChar == CR || quoteChar == LF)
			{
				throw Exception("quote character may not be CR or LF.", __FILE__, __LINE__);
			}

			this->currentChar = NUL;
			this->nextChar = EOF_;

			this->separator = (int)separator;
			this->quoteChar = (int)quoteChar;

			LowLevelRead();
		}

		/// <summary> Consumes a character from the input stream, stores it in nextChar and updates the position.
		/// </summary>

		void LowLevelRead()
		{
			nextChar = reader.get();

			if (nextChar != EOF_) position++;
		}

		/// <summary> Copies the buffered nextChar value to currentChar,
		///		with appropriate processing for digraphs such as CRLF and 2DQUOTE.
		///		<para>
		///			When this returns, currentChar contains the "most recently read symbol",
		///			and nextChar contains the first character of the next symbol to be read.
		///			This will _usually_ be the next character, but it may be the first character
		///			of a digraph.
		///		</para>
		/// </summary>

		void ReadChar()
		{
			currentChar = nextChar;

			if (nextChar == CR)
			{
#pragma region
				LowLevelRead();

				if (nextChar == LF)
				{
					LowLevelRead();
					currentChar = CRLF;
				}

				lineNumber++;
				position = 0;
#pragma endregion convert CR LF to CRLF
			}

			else if (nextChar == LF)
			{
#pragma region Process EOL
				LowLevelRead();

				lineNumber++;
				position = 0;
#pragma endregion Process EOL
			}
			else if (parsingString && nextChar == quoteChar)
			{
#pragma region
				nextChar = reader.get();

				if (nextChar == quoteChar)
				{
					nextChar = reader.get();
				}
				else
				{
					currentChar = EOS;
				}
#pragma endregion Convert DQUOTE DQUOTE to either single quote or end of string.
			}
			else
			{
#pragma region
				nextChar = reader.get();
#pragma endregion No other conversion necessary... get next char.
			}
		}


		/// <summary> Reads an entire file, which is defined by a production of the form
		///		File = record *(CRLF record) [CRLF]
		/// <para>
		///		On entry, one character from the text reader has been consumed and stored in
		///		nextChar while currentChar is uninitialised.
		/// </para>
		/// <para>
		///		On exit, all available characters in the text reader have been consumed and 
		///		currentChar is EOF_.
		/// </para>
		/// </summary>

		void ReadFile(vector<vector<string>> & rows, const unsigned int observations)
		{
			do
			{
				vector<string> currentRow;
				ReadRecord(currentRow);

				if (currentRow.size() > 1 || (currentRow.size() == 1 && currentRow[0].length() > 0))
				{
					rows.push_back(currentRow);
				}
			}
			while (IsEOL() && rows.size() < observations);
		}

		/// <summary> Reads an entire file, which is defined by a production of the form
		///		File = record *(CRLF record) [CRLF]
		/// <para>
		///		On entry, one character from the text reader has been consumed and stored in
		///		nextChar while currentChar is uninitialised.
		/// </para>
		/// <para>
		///		On exit, all available characters in the text reader have been consumed and 
		///		currentChar is EOF_.
		/// </para>
		/// </summary>

		template<typename T> void ParseFile(vector<T> & records, const unsigned int observations)
		{
			do
			{
				vector<string> currentRow;
				ReadRecord(currentRow);

				if (currentRow.size() > 1 || currentRow[0].length() > 0)
				{
					T record;
					record.Parse(currentRow);
					records.push_back(record);
				}
			}
			while (IsEOL() && records.size() < observations);
		}

		/// <summary> Reads an entire file, which is defined by a production of the form
		///		File = record *(CRLF record) [CRLF]
		/// <para>
		///		On entry, one character from the text reader has been consumed and stored in
		///		nextChar while currentChar is uninitialised.
		/// </para>
		/// <para>
		///		On exit, all available characters in the text reader have been consumed and 
		///		currentChar is EOF_.
		/// </para>
		/// </summary>

		template<typename T> void StreamFile(
			Func1<T &, bool> process,
			Action loadComplete,
			const unsigned int observations
		)
		{
			size_t count = 0;
			bool keep_going = true;

			do
			{
				vector<string> currentRow;
				ReadRecord(currentRow);

				if (currentRow.size() > 1 || currentRow[0].length() > 0)
				{
					T record;
					record.Parse(currentRow);
					keep_going = process(record);
					count++;
				}
			}
			while (keep_going && IsEOL() && (count < observations));

			loadComplete();
		}

		/// <summary> Reads a single field, implementing the three productions:
		///		field = (escaped / non-escaped)
		/// 	escaped = DQUOTE *(TEXTDATA / COMMA / CR / LF / 2DQUOTE) DQUOTE
		/// 	non-escaped = *TEXTDATA
		/// <para>
		///		On entry, currentChar should be either NUL, EOLN or COMMA. That is,
		///		we are about to process a newly encountered field. If this is not 
		///		the case, we skip forward until EOF_ or one of these is found to e
		/// </para>
		/// </summary>

		void ReadField()
		{
			while (!(IsBOF() || IsEOF() || IsEOL() || IsCOMMA())) ReadChar();

			if (IsEOF()) return;

			string f;
			currentField = &f;
			currentField->clear();

			if (nextChar == quoteChar) ReadEscaped(); else ReadNonEscaped();

			currentRow->push_back(f);
			currentField = 0;
		}

		/// <summary> Reads an unquoted field, implementing the production:
		/// 	non-escaped = *TEXTDATA
		/// 	<para>
		/// 		On entry, current char is: beginning of file, separator or end of line 
		/// 		(the previous line). 
		/// 	</para>
		/// 	<para>
		/// 		On exit, parse has advanced to the end of the current record or the next 
		/// 		separator or the end of file.
		/// 	</para>
		/// </summary>
		void ReadNonEscaped()
		{
			for (;; )
			{
				ReadChar();

				if (IsEOF() || IsEOL() || IsCOMMA()) break;

				currentField->push_back((char)(currentChar));
			}
		}

		/// <summary> Reads an escaped (quoted) field, implementing the production:
		///		escaped = DQUOTE *(TEXTDATA / COMMA / CR / LF / 2DQUOTE) DQUOTE
		///		<para>
		///			On entry, it is assumed that the parser is looking at the opening quote
		///			of an escaped field. This value is discarded and we set parsingString to 
		///			indicate that special treatment is required for quote characters. The field 
		///			is terminated by a quote char followed by anything other than a quote char, 
		///			including EOF_.
		///		</para>
		///		<para>
		///			According to the RFC, there can be no other content in the field 
		///			outside (particularly, _after_) the quotes.
		///			Therefore, once the string terminator is found, we advance until EOF_, EOL or COMMA
		///			is seen.
		///		</para>
		/// </summary>

		void ReadEscaped()
		{
			try
			{
				ReadChar();

				if (!(currentChar == quoteChar))
				{
					throw Exception("Expecting quote character at beginning of quoted field.", __FILE__, __LINE__);
				}

				parsingString = true;

				for (;; )
				{
					ReadChar();

					if (IsEOF() || IsEOS()) break;

					currentField->push_back(IsEOL() ? '\n' : (char)(currentChar));
				}

				parsingString = false;

				while (!(IsEOF() || IsCOMMA() || IsEOL())) ReadChar();
			}
			catch (Exception & e)
			{
				parsingString = false;
				throw e;
			}
		}

	private:
		string scanBuff;

	public:
		CsvReader &operator>>(string & s)
		{
			ReadField(s);
			return *this;
		}

		CsvReader & operator>>(int &i)
		{
			ReadField(scanBuff);
			i = Int::Parse(scanBuff);
			return *this;
		}

		CsvReader & operator>>(uint &i)
		{
			ReadField(scanBuff);
			i = Uint::Parse(scanBuff);
			return *this;
		}

		CsvReader & operator>>( uint64_t &i ) {
			ReadField( scanBuff );
			i = Uint64::Parse( scanBuff );
			return *this;
		}


		CsvReader & operator>>( int64_t &i ) {
			ReadField( scanBuff );
			i = Int64::Parse( scanBuff );
			return *this;
		}

		CsvReader & operator>>(double &i)
		{
			ReadField(scanBuff);
			i = Double::Parse(scanBuff);
			return *this;
		}

		CsvReader & operator>>(bool &i)
		{
			ReadField(scanBuff);
			i = Bool::Parse(scanBuff);
			return *this;
		}

		CsvReader & operator>>(ICsvReader &obj)
		{
			obj.Read(*this);
			return *this;
		}

		CsvReader & operator>>(vector<string> &record)
		{
			ReadRecord(record);
			return *this;
		}

		CsvReader & operator>>(vector<vector<string>> &table)
		{
			Read(table);
			return *this;
		}
	};

	/// <summary> CsvIO provides low-level I/O 
	///		functions to parse and generate Comma Separated Value (.CSV) text 
	///		streams as implicitly defined by the behaviour of Microsoft Excel.
	///	<para>
	///		The definition can be found in RFC 4180 - Common Format and MIME 
	///		Type for CSV Files - October 2005:<br />
	///		<a href="http://tools.ietf.org/html/rfc4180">http://tools.ietf.org/html/rfc4180</a>  
	///	</para>
	///	<para>
	///		This is a JavaScript port of source code originally implemented as part of the 
	///		CIEAM HAMISH Reliability Basics product.
	///	</para>
	/// </summary>

	class CsvWriter
	{
	private:
		bool isBol = true;

		/// <summary> Internal representation of the separator character.
		/// </summary>
		int separator;

		/// <summary> Internal representation of the escape (quote) character.
		/// </summary>
		int quoteChar;
		void ConstructorImpl(char separator, char quoteChar)
		{
			if (separator == CR || separator == LF || separator == quoteChar)
			{
				throw new Exception("separator may not be CR, LF or the escape (quote) character.", __FILE__, __LINE__);
			}

			if (quoteChar == CR || quoteChar == LF)
			{
				throw new Exception("quote character may not be CR or LF.", __FILE__, __LINE__);
			}

			this->separator = (int)separator;
			this->quoteChar = (int)quoteChar;
		}

		/* Excerpt from RFC 4180       Common Format and MIME Type for CSV Files    October 2005
		 *		http://tools.ietf.org/html/rfc4180
		 *
		 *		The grammar, modified from RFC by removing the optional header, is:
		 *
		 *		file = record *(CRLF record) [CRLF]
		 *		record = field *(COMMA field)
		 *		field = (escaped / non-escaped)
		 *		escaped = DQUOTE *(TEXTDATA / COMMA / CR / LF / 2DQUOTE) DQUOTE
		 *		non-escaped = *TEXTDATA
		 *		COMMA = %x2C
		 *		CR = %x0D ;as per section 6.1 of RFC 2234 [2]
		 *		DQUOTE =  %x22 ;as per section 6.1 of RFC 2234 [2]
		 *		LF = %x0A ;as per section 6.1 of RFC 2234 [2]
		 *		CRLF = CR LF ;as per section 6.1 of RFC 2234 [2]
		 *		TEXTDATA =  %x20-21 / %x23-2B / %x2D-7E
		 */

#pragma region
		static const int NUL = 0;
		static const int EOF_ = -1;
		static const int CR = 0x0D;
		static const int LF = 0x0A;
		static const int CRLF = -2;
		static const int EOS = -3;

		ostream & writer;

#pragma endregion Variables and constants used by Parse function.

	public:
		/// <summary> Initialises a CsvIo object with default separator and quote character.
		/// </summary>
		/// <summary> The default field separator: a comma.
		/// </summary>
		static const char DefaultSeparator = ',';

		/// <summary> The default escape character: a double quote.
		/// </summary>
		static const char DefaultQuoteChar = '"';

		CsvWriter(ostream & writer) : writer(writer)
		{
			ConstructorImpl(DefaultSeparator, DefaultQuoteChar);
		}

		/// <summary> Initialises a CsvIo object with explicit separator and default quote character.
		/// </summary>

		CsvWriter(ostream & writer, char separator) : writer(writer)
		{
			ConstructorImpl(separator, DefaultQuoteChar);
		}

		/// <summary> Initialises a CsvIO object, setting the field separator and quote character.
		///		<para>Records are always separated by either carriage return (CR) or line feed (LF) characters.</para>
		///		<para><strong>NB:</strong> This is available only if the assembly is referenced directly.</para>
		/// </summary>
		/// <param name="separator">
		///		The character that will be used to separate fields within a record.
		///		<para>
		///			This may not be carriage-return, line-feed or the same as the "quote" character.
		///		</para>
		/// </param>
		/// <param name="quoteChar">
		///		The escape (quote) character that will be used to mark up fields which contain
		///		embedded carriage-return, line-feed, separator or quote characters.
		/// </param>
		/// <exception cref="Exception">
		///		Thrown if invalid separator or quote char are supplied.
		/// </exception>

		CsvWriter(ostream & writer, char separator, char quoteChar) : writer(writer)
		{
			ConstructorImpl(separator, quoteChar);
		}

		/// <summary> Writes a collection of records to a stream in CSV format, encoding 
		/// fields as required.
		///		<para><strong>NB:</strong> This is available only if the assembly is referenced directly.</para>
		/// </summary>
		/// <param name="records">
		///		A list of records, each of which is a list of strings which are to be serialised in 
		///		CSV format.
		/// </param>
		/// <param name="writer">
		///		A TextWriter to which the resulting character stream will be written.
		/// </param>
		void Write(const vector<vector<string>> & records)
		{
			for (auto record = records.begin(); record != records.end(); record++)
			{
				Write(*record);
			}
		}

		void Write(const vector<string> & record)
		{
			WriteFields(record);
			Ln();
		}

		CsvWriter & Ln()
		{
			writer.put('\n');
			isBol = true;
			return *this;
		}

		CsvWriter & Sep()
		{
			writer.put((char)separator);
			isBol = false;
			return *this;
		}

		template< typename T>
		CsvWriter & WriteFields(const T & fields)
		{
			bool deja = false;

			for (auto &field : fields)
			{
				if (deja) Sep(); else deja = true;

				Write(field);
			}

			return *this;
		}

		CsvWriter & Write(const ICsvWriter & csvWriter)
		{
			// TODO: figure out the proper way to overload the Write template function.
			csvWriter.Write(*this);
			return *this;
		}

		template <typename T>
		CsvWriter & Write(const T & t)
		{
			string s = Util::ToString(t);
			Write(s);
			return *this;
		}

		template <typename T>
		CsvWriter & WriteSep(const T & t)
		{
			return Write(t).Sep();
		}

		template <typename T>
		CsvWriter & WriteLn(const T & t)
		{
			return Write(t).Ln();
		}

		CsvWriter & Write(const string & field)
		{
			char sepC = (char)separator;
			char quoteC = (char)quoteChar;
			char crC = (char)CR;
			char lfC = (char)LF;

			char escapeWorthyCharacters[] = { sepC, quoteC, crC, lfC };

			bool escaped = field.find_first_of(escapeWorthyCharacters, 0, 4) != string::npos;

			if ( !isBol ) {
				writer.put(sepC);
			}

			if (escaped) writer.put(quoteC);

			for (size_t j = 0; j < field.length(); j++)
			{
				char c = field[j];

				if (c == quoteC) writer.put(c);

				writer.put(c);
			}

			if (escaped) writer.put(quoteC);

			isBol = false;
			return *this;
		}

	public:
		/// Serialise a single char.
		///	@param ch A character. Two characters are treated specially: 
		///		',' causes the field separator to be written;
		///		'\n' causes the record separator to be written;
		///		Anything else is treated as a data field of length 1.
		///		You can serialise a filed consisting of a single occurrence of 
		///		these special chars by putting it in a string.
		///	@returns *this
		CsvWriter &operator<<(char ch)
		{
			if (ch == ',') Sep();
			else if (ch == '\n') Ln();
			else
			{
				char s[] = { ch, 0 };
				Write(s);
			}
			return *this;
		}

		CsvWriter &operator<<(const string & s)
		{
			Write(s);
			return *this;
		}

		CsvWriter &operator<<(const char *s)
		{
			Write(s);
			return *this;
		}

		CsvWriter & operator<<(int &i)
		{
			Write(i);
			return *this;
		}

		CsvWriter & operator<<(uint i)
		{
			Write(i);
			return *this;
		}

		CsvWriter & operator<<(uint64_t i)
		{
			Write(i);
			return *this;
		}

		CsvWriter & operator<<(double i)
		{
			Write(i);
			return *this;
		}

		CsvWriter & operator<<(bool i)
		{
			Write<bool>(i);
			return *this;
		}

		CsvWriter & operator<<(const ICsvWriter &obj)
		{
			obj.Write(*this);
			return *this;
		}

		CsvWriter & operator<<(const vector<string> &record)
		{
			Write(record);
			return *this;
		}

		CsvWriter & operator<<(const vector<vector<string>> &table)
		{
			this->Write(table);
			return *this;
		}
	};
}

