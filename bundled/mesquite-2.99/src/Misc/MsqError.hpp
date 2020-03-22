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

#ifndef MSQ_ERROR_HPP
#define MSQ_ERROR_HPP

#include "Mesquite.hpp"

#include <iosfwd>
#include <list>
#include <string>

namespace MESQUITE_NS
{

/**\defgroup error Mesquite internal error handling
 */
/*@{*/

#ifndef MSQ_FUNCTION
#  define MSQ_FUNCTION ""
#endif


/**\brief Mesquite's Error Checking macro.
 *
 * Check the status of an MsqError.  Returns true as pushes the current
 * file/line onto the stack trace if the error flag is true.  Returns
 * false otherwise.  
 */
#define MSQ_CHKERR(err)  (err.error() && err.push(MSQ_FUNCTION,__FILE__,__LINE__))

/**\brief If passed error is true, return from a void function.
 *
 * Shorthand for if (MSQ_CHKERR(err)) return
 */
#define MSQ_ERRRTN(err)  if (MSQ_CHKERR(err)) return

/**\brief Return zero/NULL on error.
 */
#define MSQ_ERRZERO(err)  if (MSQ_CHKERR(err)) return 0

/**\brief Return false on error.
 */
#define MSQ_ERRFALSE(err)  if (MSQ_CHKERR(err)) return false

/**\brief Macro to set error - use err.clear() to clear.
 *
 * Set the error object to the specified error code and optional error message,
 * and push the current file/line onto the stack trace.  
 * Examples:
 *  - MSQ_SETERR( err )( MsqError::INVALID_ARG );
 *  - MSQ_SETERR( err )( "foo cannot be zero", MsqError::INVALID_ARG );
 *  - MSQ_SETERR( err )( MsqError::INVALID_ARG, "foo = %d", foo );
*/
#define MSQ_SETERR(err) \
Mesquite::MsqError::setter(err, MSQ_FUNCTION,__FILE__,__LINE__).set

/**
 *\class MsqError
 *\brief  Used to hold the error state and return it to the application.
 *\author Jason Kraftcheck
 *\date   2004-09-17
 *
 * Used to hold error state and related information.
 * Internal Mesquite code should access this object via
 * the MSQ_SETERR() and MSQ_CHKERR() macros. 
 *
 * For applications, the cast-to-bool operator and << operator
 * are provided for convenient, if simple access to this data.
 * E.g.:  if (err) cout << err << endl;
 * 
 * There are two options for an application to gain more detailed
 * access to the error data.  The application may either access
 * the data stored in this class via the provided methods or
 * subclass MsqError, overriding set_error() and push()
 * to handle the error data as it is generated.
*/
class MsqError {
public:

    /* NOTE: If you add an error to this list, make sure
       a) you add it *before* LAST_ERROR_CODE
       b) you add the corresponding string in MsqError.cpp
     */
    /** \brief Error codes
     */
  enum ErrorCode {
    NO_ERROR          =  0,/**< no error */
    UNKNOWN_ERROR,         /**< unknown error occured */
    OUT_OF_MEMORY,         /**< unable to allocate the necessary memory */
    INVALID_ARG ,          /**< invalid function argument passed */
    NOT_INITIALIZED,       /**< object not initialized */
    INVALID_STATE,         /**< object is in an invalid state */
    FILE_ACCESS,           /**< File cannot be opened/created. */
    FILE_FORMAT,           /**< Wrong file type */
    PARSE_ERROR,           /**< Error parsing input (or input file) */
    IO_ERROR,              /**< An I/O error occured (e.g. read from file failed.) */
    INVALID_MESH,          /**< The mesh is invalid */
    NO_PD_STORAGE_MODE,    /**< no storage mode chosen within PatchData */
    NOT_IMPLEMENTED,       /**< requested functionality is not (yet) implemented */
    INTERNAL_ERROR,        /**< A bug in Mesquite */
    INTERRUPTED,           /**< Application or user interrupted operation */
    TAG_ALREADY_EXISTS,    /**< Attempt to create tag that already exists */
    TAG_NOT_FOUND,         /**< Specified tag does not exist */
    UNSUPPORTED_ELEMENT,   /**< the element type is not supported. */
    PARALLEL_ERROR,        /**< an error occurred in parallel > */
    BARRIER_VIOLATED,      /**< barruer violated when processing barrier Target Metric */
    LAST_ERROR_CODE
  };
  
  //! \brief resets error object to non-active state (no error).
  MESQUITE_EXPORT void clear(); 
  
  //! \brief Check if an error has occured
  inline bool error() const       { return NO_ERROR != errorCode; }
  //! \brief Check if an error has occured
  inline operator bool() const    { return NO_ERROR != errorCode; }

  //! \brief Initialize to cleared state.
  MESQUITE_EXPORT MsqError() : errorCode(NO_ERROR) { } 
  
  //! Destructor - empty but must declar virtual destrucor if virtual functions.
  MESQUITE_EXPORT virtual ~MsqError();

 /* ************************************************************ *
  *            Low-level access to error data
  * ************************************************************ */  

  //! Get error code
  inline ErrorCode error_code() const    { return errorCode; }

  //!\class Trace
  //!\brief One line of stack trace data
  struct MESQUITE_EXPORT Trace {
    std::string function;
    std::string file;
    int line;
    
    Trace( const char* fun, const char* fil, int lin )
     : function(fun), file(fil), line(lin) {}
  };
  
  //! Get error message
  MESQUITE_EXPORT const char* error_message() const;
  
  //! Container type used to store stack trace.
  //! Return type for stack()
  typedef std::list<Trace> StackTrace;

  //! Get stack trace
  inline const StackTrace& stack() const { return stackTrace; }


 /* ************************************************************ *
  *                        Set error data
  * ************************************************************ */  
  
  //! Add to back-trace of call stack.  Called by MSQ_CHKERR.
  //! Must always return true.
  MESQUITE_EXPORT virtual bool push( const char* function, const char* file, int line );
  
  //! Initialize the error object with the passed data.
  MESQUITE_EXPORT virtual bool set_error( ErrorCode num, const char* msg = 0 );
  
  //!\class setter
  //! Used for implementing pre-processor macros for internal use
  //! in Mesquite.
  class MESQUITE_EXPORT Setter {
    public:
      Setter(MsqError& err, const char* function, const char* file, int line) 
        : mErr(err), functionName(function), fileName(file), lineNumber(line) {}
      
      bool set( ErrorCode num );
      bool set( const char* message, ErrorCode num );
      bool set( const std::string& message, ErrorCode num );
      bool set( ErrorCode num, const char* format, ... )
      #ifdef __GNUC__
        __attribute__ ((format (printf, 3, 4)))
      #endif
        ; // ending semicolon for set( ErrorCode num, const char* format, ... )
    private:
      MsqError& mErr;
      const char* functionName;
      const char* fileName;
      int lineNumber;
  };

  static inline Setter setter( MsqError& err, const char* function, const char* file, int line )
  { return Setter( err, function, file, line ); }

private:

  ErrorCode errorCode;
  std::string errorMessage;
  StackTrace stackTrace;
};


  //! Print message and stack trace
MESQUITE_EXPORT std::ostream& operator<<( std::ostream&, const MsqError& );
  //! Print MsqError::Trace
MESQUITE_EXPORT std::ostream& operator<<( std::ostream&, const MsqError::Trace& );

/**
 *\class MsqPrintError
 *\brief  Utility class for printing error data - used in Mesquite tests.
 *\author Jason Kraftcheck
 *\date   2004-10-11
 *
 * A subclass of MsqError.  Behaves the same as MsqError, except that
 * it will automatically print itself to the specified ostream upon
 * destruction if an error occured.  For objections of this type 
 * declared on the stack (not new'd), this means that the error will
 * be printed when the function returns (if an error occured.)
*/
class MsqPrintError : public MsqError
{
  public:
      //!\brief Initialize with ostream to print error data to.
    MESQUITE_EXPORT MsqPrintError( std::ostream& stream )
      : outputStream(stream) {}
    
      //!\brief On destruction, conditionally prints error data.
    MESQUITE_EXPORT virtual ~MsqPrintError( );
    
  private:
      
    std::ostream& outputStream;
};

/*@}*/

} // namespace


#endif


