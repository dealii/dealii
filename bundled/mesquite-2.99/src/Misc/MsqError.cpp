/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
   
  ***************************************************************** */
/*!
  \file   MsqError.cpp
  \brief  Used to hold the error state and return it to the application.
  \author Jason Kraftcheck
  \date   2004-09-17
*/

#include "MsqError.hpp"
#include "Mesquite.hpp"

#include <ostream>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include <cstring>
using std::strncpy;

namespace MESQUITE_NS {

const char* MsqError::error_message() const
{
  static const char* const error_messages[] = {
   "No Error",
   "<unknown>",
   "Out of memory",
   "Invalid argument",
   "Data not initialized",
   "Invalid state",
   "File access error",
   "File format error",
   "Syntax error",
   "I/O error",
   "Invalid mesh",
   "No storage mode for PatchData",
   "Not implemented",
   "Internal error",
   "Interrupted",
   "Duplicate tag name",
   "Tag not found",
   "Unsupported element type",
   "Parallel Error - error occurred on at least one processor",
   "barruer violated when processing barrier Target Metric",
   "Invalid Error Code"
  };
  
    /* If this is ever false, it should be caught by a unit test. 
       Do an assert here so the unit test fails.
       This asserts that all error codes have a string in the above list. */
  assert( sizeof(error_messages) == sizeof(char*) * (LAST_ERROR_CODE+1) );
  
  if (!errorMessage.empty())
    return errorMessage.c_str();
  
  if (errorCode >= 0 && errorCode < LAST_ERROR_CODE)
    return error_messages[errorCode];
  
  return error_messages[LAST_ERROR_CODE];
}

MsqError::~MsqError() {}

bool MsqError::Setter::set( const std::string& msg, ErrorCode num )
{
  return mErr.set_error( num, msg.c_str() ) 
      && mErr.push( functionName, fileName, lineNumber );
}

bool MsqError::Setter::set( const char* msg, ErrorCode num )
{
  return mErr.set_error( num, msg )
      && mErr.push( functionName, fileName, lineNumber );
}

bool MsqError::Setter::set( ErrorCode num )
{
  return mErr.set_error( num )
      && mErr.push( functionName, fileName, lineNumber );
}

bool MsqError::Setter::set( ErrorCode num, const char* format, ... )
{
  char buffer[1024];
  
#if defined(HAVE_VSNPRINTF)
  va_list args;
  va_start( args, format );
  vsnprintf( buffer, sizeof(buffer), format, args );
  va_end( args );
#elif defined(HAVE__VSNPRINTF)
  va_list args;
  va_start( args, format );
  _vsnprintf( buffer, sizeof(buffer), format, args );
  va_end( args );
#elif defined(HAVE_VSPRINTF)
  va_list args;
  va_start( args, format );
  vsprintf( buffer, format, args );
  va_end( args );
#else
  strncpy( buffer, format, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';
#endif

  return mErr.set_error( num, buffer )
      && mErr.push( functionName, fileName, lineNumber );
}

bool MsqError::push( const char* function, const char* file, int line )
{
  stackTrace.push_back( Trace(function, file, line) );
  return true;
}

bool MsqError::set_error( ErrorCode num, const char* msg )
{
  errorCode = num;
  stackTrace.clear();
  
  if (msg)
    errorMessage = msg;
  else
	  // MS VC6 doesn't have string::clear()!
	errorMessage.resize(0);
    
  return num != NO_ERROR;
}

void MsqError::clear()
{
  errorCode = NO_ERROR;
	// MS VC6 doesn't have string::clear()!
  errorMessage.resize(0);
  stackTrace.clear();
}

std::ostream& operator<<( std::ostream& str, const MsqError::Trace& tr ) 
{
  return (str << tr.function << " at " << tr.file << ":" << tr.line);
}

std::ostream& operator<<( std::ostream& str, const MsqError& err ) 
{
  str << "MESQUITE ERROR " << (int)err.error_code() << " : " 
      << err.error_message() << std::endl;

  MsqError::StackTrace::const_iterator iter = err.stack().begin();
  const MsqError::StackTrace::const_iterator end = err.stack().end();
  if (iter != end)
  {
    str << "  at " << *iter << std::endl;
    ++iter;
  }
  for ( ; iter != end; ++iter)
    str << "  in " << *iter << std::endl;
  
  return str;
}

MsqPrintError::~MsqPrintError()
  { if (error()) outputStream << *this << std::endl; }

}
