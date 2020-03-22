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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE:  7-Nov-02 at 16:22:26
//  LAST-MOD:  8-Nov-02 at 10:27:00 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file QualityImprover.cpp

Implements a couple of default virtual functions of the virtual class

 \author Thomas Leurent
 */
// DESCRIP-END.
//

#include "QualityImprover.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

    /*! The default constructor initialises a few member variables
        to default values.
        This can be reused by concrete class constructor. */    
QualityImprover::QualityImprover()
{
  defaultOuterCriterion = new TerminationCriterion;
  defaultInnerCriterion = new TerminationCriterion;
  defaultOuterCriterion->add_iteration_limit( 1 );
  outerTerminationCriterion = defaultOuterCriterion;
  innerTerminationCriterion = defaultInnerCriterion;
}
 
QualityImprover::~QualityImprover()
{
	delete defaultOuterCriterion;
	delete defaultInnerCriterion;
}
    
void QualityImprover::initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                        const Settings* settings,
                                        MsqError& err )
{
  innerTerminationCriterion->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
  outerTerminationCriterion->initialize_queue( mesh_and_domain, settings, err ); MSQ_ERRRTN(err);
}

} // namespace MESQUITE_NS

