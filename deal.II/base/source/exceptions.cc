//----------------------------  exceptions.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  exceptions.cc  ---------------------------


//TODO:[WB] (compiler) replace s.c_str() by s when that is possible


#include <base/exceptions.h>
#include <string>
#include <cstdlib>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif



ExceptionBase::ExceptionBase () :
		file(""), line(0), function(""), cond(""), exc("")
{};



ExceptionBase::ExceptionBase (const char* f, const int l, const char *func,
			      const char* c, const char *e) :
		file(f), line(l), function(func), cond(c), exc(e)
{};



ExceptionBase::~ExceptionBase () throw ()
{};



void ExceptionBase::SetFields (const char* f,
			       const int l,
			       const char *func,
			       const char *c,
			       const char *e)
{
  file = f;
  line = l;
  function = func;
  cond = c;
  exc  = e;
};



void ExceptionBase::PrintExcData (std::ostream &out) const
{
  out << "An error occurred in line <" << line
      << "> of file <" << file
      << "> in function" << std::endl
      << "    " << function << std::endl
      << "The violated condition was: "<< std::endl
      << "    " << cond << std::endl
      << "The name and call sequence of the exception was:" << std::endl
      << "    " << exc  << std::endl
      << "Additional Information: " << std::endl;
};



void ExceptionBase::PrintInfo (std::ostream &out) const
{
  out << "(none)" << std::endl;
};



const char * ExceptionBase::what () const throw ()
{
				   // have a place where to store the
				   // description of the exception as a char *
				   //
				   // this thing obviously is not multi-threading
				   // safe, but we don't care about that for now
				   //
				   // we need to make this object static, since
				   // we want to return the data stored in it
				   // and therefore need a liftime which is
				   // longer than the execution time of this
				   // function
  static std::string description;
				   // convert the messages printed by the
				   // exceptions into a std::string
#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream converter;
#else
  std::ostrstream converter;
#endif

  converter << "--------------------------------------------------------"
	    << std::endl;
				   // put general info into the std::string
  PrintExcData (converter);
				   // put in exception specific data
  PrintInfo (converter);
  
  converter << "--------------------------------------------------------"
	    << std::endl
	    << std::ends;

  description = converter.str();

  return description.c_str();
};



namespace deal_II_exceptions
{

  std::string additional_assert_output;

  void set_additional_assert_output (const char * const p)
  {
    additional_assert_output = p;
  };

  
  
  namespace internals 
  {
    
				     /**
				      * Number of exceptions dealt
				      * with so far. Zero at program
				      * start. Messages are only
				      * displayed if the value is
				      * zero.
				      */
    static unsigned int n_treated_exceptions;
  

    void issue_error_assert (const char *file,
			     int         line,
			     const char *function,
			     const char *cond,
			     const char *exc_name,
			     ExceptionBase &         e)
    {
				       // fill the fields of the
				       // exception object
      e.SetFields (file, line, function, cond, exc_name);
      
				       // if no other exception has
				       // been displayed before, show
				       // this one
      if (n_treated_exceptions == 0)
	{
	  std::cerr << "--------------------------------------------------------"
		    << std::endl;
					   // print out general data
	  e.PrintExcData (std::cerr);
					   // print out exception
					   // specific data
	  e.PrintInfo (std::cerr);
	  std::cerr << "--------------------------------------------------------"
		    << std::endl;

					   // if there is more to say,
					   // do so
	  if (!additional_assert_output.empty())
	    std::cerr << additional_assert_output << std::endl;
	}
      else
	{
					   // if this is the first
					   // follow-up message,
					   // display a message that
					   // further exceptions are
					   // suppressed
	  if (n_treated_exceptions == 1)
	    std::cerr << "******** More assertions fail but messages are suppressed! ********"
		      << std::endl;
	};
      
				       // increase number of treated
				       // exceptions by one
      n_treated_exceptions++;
    
    
				       // abort the program now since
				       // something has gone horribly
				       // wrong. however, there is one
				       // case where we do not want to
				       // do that, namely when another
				       // exception, possibly thrown
				       // by AssertThrow is active,
				       // since in that case we will
				       // not come to see the original
				       // exception. in that case
				       // indicate that the program is
				       // not aborted due to this
				       // reason.
      if (std::uncaught_exception() == true)
	{
					   // only display message once
	  if (n_treated_exceptions == 0)
	    std::cerr << "******** Program is not aborted since another exception is active! ********"
		      << std::endl;
	}
      else
	std::abort ();
    };



    void abort ()
    {
      std::abort ();
    };
    
  };
  
};
