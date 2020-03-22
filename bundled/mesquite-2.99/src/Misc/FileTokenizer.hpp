/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#ifndef FILE_TOKENIZER_HPP
#define FILE_TOKENIZER_HPP

#include "Mesquite.hpp"
#include <cstdio>
#include <sys/types.h>

namespace MESQUITE_NS
{

  class MsqError;

/** 
 * \class  FileTokenizer
 * \brief  Parse a file as space-separated tokens
 * \author Jason Kraftcheck
 * \date   30 Sept 2004
 *
 * Read a file, separating it into space-separated tokens.
 * This is provided in place of using the standard C or C++
 * file parsing routines because it counts lines, which is
 * useful for error reporting.  Also provides some useful
 * utility methods for parsing VTK files (which is the 
 * intended use of this implementation.)
 *
 * Uses raw reads/writes, implementing internal buffering.
 * Token size may not exceed buffer size.
 */

class FileTokenizer 
{
  public:
  
      /** \brief constructor 
       * 
       * \param file_ptr The file to read from.
       */
    FileTokenizer( std::FILE* file_ptr );
    
      /** \brief destructor : closes file.
       *
       * The destructor closes the passed file handle.   This
       * is done as a convenience feature.  If the caller
       * creates an instance of this object on the stack, the
       * file will automatically be closed when the caller
       * returns.
       */
    ~FileTokenizer();
    
      /** \brief get next token
       *
       * Get the next whitesapce-deliminated token from the file.
       * NOTE: The returned string is only valid until the next
       *       call to any of the functions in this class that
       *       read from the file.
       *
       * \return A pointer to the buffer space containing the string,
       *         or NULL if an error occured.
       */
    const char* get_string( MsqError& err );
    
      /** \brief check for newline
       *
       * Consume whitespace upto and including the next newline.
       * If a non-space character is found before a newline, 
       * the function will stop, set the error message, and
       * return false.
       * 
       * \return True if a newline was found before any non-space
       *         character.  False otherwise.
       */
    bool get_newline( MsqError& err );
    
    
      /** \brief Parse a sequence of double values.
       *
       * Read the specified number of space-deliminated doubles.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_doubles( size_t count, double* array, MsqError& err );
     
    
      /** \brief Parse a sequence of float values.
       *
       * Read the specified number of space-deliminated doubles.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_floats( size_t count, float* array, MsqError& err );
   
      /** \brief Parse a sequence of integer values.
       *
       * Read the specified number of space-deliminated ints.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_integers( size_t count, int* array, MsqError& err );
   
      /** \brief Parse a sequence of integer values.
       *
       * Read the specified number of space-deliminated ints.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_long_ints( size_t count, long* array, MsqError& err );
   
      /** \brief Parse a sequence of integer values.
       *
       * Read the specified number of space-deliminated ints.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_short_ints( size_t count, short* array, MsqError& err );
   
      /** \brief Parse a sequence of integer values.
       *
       * Read the specified number of space-deliminated ints.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_bytes( size_t count, unsigned char* array, MsqError& err );
    
      /** \brief Parse a sequence of bit or boolean values.
       *
       * Read the specified number of space-deliminated values.
       *
       * \param count   The number of values to read.
       * \param array   The memory at which to store the values.
       * \return true if successful, false otherwise.
       */
    bool get_booleans( size_t count, bool* array, MsqError& err );
  
      /** 
       * Check for end-of-file condition.
       */
    bool eof() const;
    
      /** 
       * Get the line number the last token was read from.
       */
    int line_number() const { return lineNumber; }
    
      /** 
       * Put current token back in buffer.  Can only unget one token.
       */
    void unget_token();
    
      /**
       * Match current token to passed string.  If token
       * doesn't match, set error message.
       */
    bool match_token( const char* string, MsqError& err );
    
      /**
       * Match the current token to one of an array of strings.  
       * Sets the error message if the current token doesn't
       * match any of the input strings.
       *
       * \param  string_list A NULL-terminated array of strings.
       * \return One greater than the index of the matched
       *         string, or zero if no match.
       */
    int match_token( const char* const* string_list, MsqError& err );
  
  private:
  
      /** Internal implementation of get_doubles() */
    bool get_double_internal( double& result, MsqError& err );
      /** Internal implementation of get_long_ints() */
    bool get_long_int_internal( long& result, MsqError& err );
      /** Internal implementation of get_Booleans() */
    bool get_boolean_internal( bool& result, MsqError& err );
  
      /** Internal implementation of get_floats() */
    bool get_float_internal( float& result, MsqError& err );
      /** Internal implementation of get_integers() */
    bool get_integer_internal( int& result, MsqError& err );
      /** Internal implementation of get_short_ints() */
    bool get_short_int_internal( short& result, MsqError& err );
      /** Internal implementation of get_bytes() */
    bool get_byte_internal( unsigned char& result, MsqError& err );
  
      /** Pointer to standard C FILE struct */
    std::FILE* filePtr;
    
      /** Input buffer */
    char buffer[512];
    
      /** One past the end of the last token returned */
    char *nextToken;
      /** One past the last used byte of the buffer */
    char *bufferEnd;
    
      /** Line number of last returned token */
    int lineNumber;
    
      /** The whitespace character marking the end of the 
       *  last returned token.  Saved here because if it
       *  is a newline, the line count will need to be
       *  incremented when the next token is returned.
       */
    char lastChar;
};

} // namespace Mesquite

#endif
