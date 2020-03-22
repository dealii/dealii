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
/*!
  \file   QualityImprover.hpp
  \brief  

  The Quality Improver Class is the base class for all the algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_QualityImprover_hpp 
#define Mesquite_QualityImprover_hpp

#include <string>

#include "Mesquite.hpp"
#include "Instruction.hpp"
#include "TerminationCriterion.hpp"

namespace MESQUITE_NS
{
  class PatchSet;
  class MsqError;
  class Mesh;
  class MeshDomain;
  class Settings;

  
  /*! \class QualityImprover
    \brief Base class for all quality improvers.
    Mote that the PatchData settings are inherited from the PathDataUser class. 

  */ 
  class MESQUITE_EXPORT QualityImprover : public Instruction
  {
  public:

    // Constructor is protected ... see below.
    
     // virtual destructor ensures use of polymorphism during destruction
    virtual ~QualityImprover();

      //!Sets in the termination criterion for the concrete solver's
      //! optimization.
    void set_inner_termination_criterion(TerminationCriterion* crit)
      {
        innerTerminationCriterion=crit;
      }
      //!Sets in the termination criterion for the outer loop over 
      //! patches.
    void set_outer_termination_criterion(TerminationCriterion* crit)
      {
        outerTerminationCriterion=crit;
      }
      
    virtual PatchSet* get_patch_set() = 0;
    
    virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );

  protected:

    /*! The default constructor initialises a few member variables
        to default values.
        This can be reused by concrete class constructor. */    
    QualityImprover();
    
      //!return the outer termination criterion pointer 
    TerminationCriterion* get_outer_termination_criterion()
      { return outerTerminationCriterion; }
      //!return the inner termination criterion pointer       
    TerminationCriterion* get_inner_termination_criterion()
      { return innerTerminationCriterion; } 
    
  private:
  
    TerminationCriterion* innerTerminationCriterion;
    TerminationCriterion* outerTerminationCriterion;
      //default TerminationCriterion for outer loop will be set in constructor
    TerminationCriterion* defaultOuterCriterion;
      //default TerminationCriterion for inner loop set by concrete improver
    TerminationCriterion* defaultInnerCriterion;
  };

}

#endif
