// $Id$

#ifndef __jobidentifier_H
#define __jobidentifier_H

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
