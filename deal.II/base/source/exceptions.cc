//----------------------------  exceptions.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  exceptions.cc  ---------------------------


#include <strstream>


ExceptionBase::ExceptionBase () :
		file(""), line(0), function(""), cond(""), exc("")
{};


ExceptionBase::ExceptionBase (const char* f, const int l, const char *func,
			      const char* c, const char *e) :
		file(f), line(l), function(func), cond(c), exc(e)
{};


void ExceptionBase::SetFields (const char* f,
			       const int l,
			       const char *func,
			       const char *c,
			       const char *e) {
  file = f;
  line = l;
  function = func;
  cond = c;
  exc  = e;
};


void ExceptionBase::PrintExcData (ostream &out) const {
  out << "An error occurred in line <" << line
      << "> of file <" << file
      << "> in function" << endl
      << "    " << function << endl
      << "The violated condition was: "<< endl
      << "    " << cond << endl
      << "The name and call sequence of the exception was:" << endl
      << "    " << exc  << endl
      << "Additional Information: " << endl;
};


void ExceptionBase::PrintInfo (ostream &out) const {
  out << "(none)" << endl;
};


const char * ExceptionBase::what () const {
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
  static string description;
				   // convert the messages printed by the
				   // exceptions into a string
  ostrstream converter;

  converter << "--------------------------------------------------------"
	    << endl;
				   // put general info into the string
  PrintExcData (converter);
				   // put in exception specific data
  PrintInfo (converter);
  
  converter << "--------------------------------------------------------"
	    << endl
	    << ends;

  description = converter.str();

  return description.c_str();
};

