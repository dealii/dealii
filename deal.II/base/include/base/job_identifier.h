//----------------------------  job_identifier.h  ---------------------------
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
//----------------------------  job_identifier.h  ---------------------------
#ifndef __deal2__job_identifier_h
#define __deal2__job_identifier_h


#include <string>

/**
 * Identification of a program run. #JobIdentifier# determines the
 * start time of a program run and stores it as a program
 * identifier. There exists a library object #dealjobid# of this
 * class. This object can be accessed by all output functions to
 * provide an id for the current job.
 */
class JobIdentifier
{
public:
  JobIdentifier();
  ~JobIdentifier();
  
  static const char* program_id();
  const string operator () ()const;
private:
  string id;
};

/*------------------------------ Inline functions ------------------------------*/

extern JobIdentifier dealjobid;

#endif
