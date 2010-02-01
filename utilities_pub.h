/*    Utility functions for testing file existence, reading/parsing command-line 
 * options, etc.
 *
 */

#ifndef _UTILITIES_PUB_H_
#define _UTILITIES_PUB_H_

#include <string>
#include <vector>
using namespace std;

/* constants for use parameter "restriction" when calling NotANumber(): */
#define kAnyInt          0
#define kPosInt          1
#define kAnyReal         2
#define kPosReal         3


// String-processing functions

// Splits a string and returns substrings as elements of tokens.
void SplitString( const string& str, vector<string>& tokens, const string& delimiters = "\t " );

// Removes remainder of string after first occurance of delimiter.
void ChopComment( string& inputString, char delimiter = '#' );

// Removes leading and trailing whitespace ("\t ")from a string
void TrimWhitespace( string& stringToModify );


bool FileExists(const char * filename);


void CommandLineError( char errorString[] );


bool NotANumber( char theString[], int index, int restriction );


#endif /* _UTILITIES_PUB_H_ */