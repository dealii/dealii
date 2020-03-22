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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-

#include "Mesquite.hpp"
#include "mesquite_version.h"

const char* Mesquite::version_string(bool)
{
  return MSQ_VERSION_STRING;
}

unsigned int Mesquite::major_version_number()
{
  return MSQ_VERSION_MAJOR;
}

unsigned int Mesquite::minor_version_number()
{
  return MSQ_VERSION_MINOR;
}

unsigned int Mesquite::patch_version_number()
{
#ifdef MSQ_VERSION_PATCH
  return MSQ_VERSION_PATCH;
#else
  return 0;
#endif
}

Mesquite::ReleaseType Mesquite::release_type()
{  
#if MSQ_VERSION_MINOR == 99
  return Mesquite::ALPHA;
#elif defined(MSQ_VERSION_PATCH)
  return Mesquite::RELEASE;
#else
  return Mesquite::BETA;
#endif
}

