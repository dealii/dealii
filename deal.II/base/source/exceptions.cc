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

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


unsigned int ExceptionBase::n_treated_exceptions;


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
