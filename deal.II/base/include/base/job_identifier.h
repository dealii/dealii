//----------------------------  job_identifier.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
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
 * Identification of a program run. @p{JobIdentifier} determines the
 * start time of a program run and stores it as a program
 * identifier. There exists a library object @p{dealjobid} of this
 * class. This object can be accessed by all output functions to
 * provide an id for the current job.
 */
class JobIdentifier
{
  public:
				     /**
				      * Constructor. Set program
				      * identifier to value of
				      * @p{program_id} concatenated
				      * with the present time.
				      */
    JobIdentifier();

				     /**
				      * ???
				      */
    static const char* program_id();

				     /**
				      * Return the value of @p{id}.
				      */
    const std::string operator () ()const;
    
  private:
				     /**
				      * String holding the identifier
				      * of the presently running
				      * program.
				      */
    std::string id;
};


/*------------------------------ Inline functions ------------------------------*/


/**
 * Global object to identify the presently running program.
 */
extern JobIdentifier dealjobid;

#endif
