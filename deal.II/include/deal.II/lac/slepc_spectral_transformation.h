//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//    Author: Toby D. Young, Polish Academy of Sciences, 2009
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__slepc_spectral_transformation_h
#define __deal2__slepc_spectral_transformation_h


#include <base/config.h>

#ifdef DEAL_II_USE_SLEPC

#  include <base/std_cxx1x/shared_ptr.h>
#  include <lac/exceptions.h>

#  include <petscksp.h>
#  include <slepceps.h>

DEAL_II_NAMESPACE_OPEN


namespace SLEPcWrappers
{

/**
 * Base class for spectral transformation classes using the SLEPc
 * solvers which are selected based on flags passed to the spectral
 * transformation.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2009
 **/
  class TransformationBase
    {
    public:

                                   /**
                                    * Constructor. 
                                    */
      TransformationBase ();

                                   /**
                                    * Destructor.
                                    */
      virtual ~TransformationBase ();

                                   /**
				    * Record the EPS object that is associated
				    * to the spectral transformation
                                    */
      void set_context (EPS &eps);

    protected:

      virtual void set_transformation_type (ST &st) const = 0;

    private:

                                   /**
                                    * Objects of this type are
                                    * explicitly created, but are
                                    * destroyed when the surrounding
                                    * solver object goes out of scope,
                                    * or when we assign a new value to
                                    * the pointer to this object. The
                                    * respective Destroy functions are
                                    * therefore written into the
                                    * destructor of this object, even
                                    * though the object does not have
                                    * a constructor.
                                    */
      struct TransformationData
      {

                                   /**
                           	    * Destructor.
                              	    */
	~TransformationData ();

                                   /**
                                    * Objects for Eigenvalue Problem
                                    * Solver.
                                    */
	ST st;
      };

      std_cxx1x::shared_ptr<TransformationData> transformation_data;
    };

/**
 * An implementation of the transformation interface using the SLEPc
 * Shift.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2009
 */
  class TransformationShift : public TransformationBase
    {
    public:

                                   /**
                                    * Standardized data struct to
                                    * pipe additional data to the
                                    * solver.
                                    */
      struct AdditionalData
      {

                                   /**
                                    * Constructor. By default, set the
                                    * shift parameter to zero.
                                    */
	AdditionalData (const double shift_parameter = 0);

                                   /**
                                    * Shift parameter.
                                    */
	const double shift_parameter;
      };


                                   /**
                                    * Constructor.
                                    */
      TransformationShift (const AdditionalData &data = AdditionalData());


    protected:

                                   /**
                                    * Store a copy of the flags for this
                                    * particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Spectral
                                    * Transformation context object,
                                    * and sets the type of spectral
                                    * transformationthat is
                                    * appropriate for this class.
                                    */
      virtual void set_transformation_type (ST &st) const;
    };

/**
 * An implementation of the transformation interface using the SLEPc
 * Shift and Invert.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2009
 */
  class TransformationShiftInvert : public TransformationBase
    {
    public:

                                   /**
                                    * Standardized data struct to
                                    * pipe additional data to the
                                    * solver.
                                    */
      struct AdditionalData
      {
                                   /**
                                    * Constructor. By default, set the
                                    * shift parameter to zero.
                                    */
	AdditionalData (const double shift_parameter = 0);

                                   /**
                                    * Shift parameter.
                                    */
	const double shift_parameter;
      };


                                   /**
                                    * Constructor.
                                    */
      TransformationShiftInvert (const AdditionalData &data = AdditionalData());

    protected:

                                   /**
                                    * Store a copy of the flags for this
                                    * particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Spectral
                                    * Transformation context object,
                                    * and sets the type of spectral
                                    * transformationthat is
                                    * appropriate for this class.
                                    */
      virtual void set_transformation_type (ST &st) const;
    };

/**
 * An implementation of the transformation interface using the SLEPc
 * Spectrum Folding.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2009
 */
  class TransformationSpectrumFolding : public TransformationBase
    {
    public:

                                   /**
                                    * Standardized data struct to
                                    * pipe additional data to the
                                    * solver.
                                    */
      struct AdditionalData
      {
                                   /**
                                    * Constructor. By default, set the
                                    * shift parameter to zero.
                                    */
	AdditionalData (const double shift_parameter = 0);

                                   /**
                                    * Shift parameter.
                                    */
	const double shift_parameter;
      };


                                   /**
                                    * Constructor.
                                    */
      TransformationSpectrumFolding (const AdditionalData &data = AdditionalData());

    protected:

                                   /**
                                    * Store a copy of the flags for this
                                    * particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Spectral
                                    * Transformation context object,
                                    * and sets the type of spectral
                                    * transformationthat is
                                    * appropriate for this class.
                                    */
      virtual void set_transformation_type (ST &st) const;
    };

/**
 * An implementation of the transformation interface using the SLEPc
 * Cayley.
 *
 * @ingroup SLEPcWrappers
 * @author Toby D. Young 2009
 */
  class TransformationCayley : public TransformationBase
    {
    public:

                                   /**
                                    * Standardized data struct to pipe
                                    * additional data to the solver.
                                    */
      struct AdditionalData
      {
                                   /**
                                    * Constructor. Requires two shift
                                    * parameters
                                    */
	AdditionalData (const double shift_parameter     = 0,
			const double antishift_parameter = 0);

                                   /**
                                    * Shift and antishift parameter.
                                    */
	const double shift_parameter;
	const double antishift_parameter;
      };


                                   /**
                                    * Constructor.
                                    */
      TransformationCayley (const double shift,
			    const double antishift);

    protected:

                                   /**
                                    * Store a copy of the flags for this
                                    * particular solver.
                                    */
      const AdditionalData additional_data;

                                   /**
                                    * Function that takes a Spectral
                                    * Transformation context object,
                                    * and sets the type of spectral
                                    * transformationthat is
                                    * appropriate for this class.
                                    */
      virtual void set_transformation_type (ST &st) const;
    };

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_SLEPC

/*--------------------   slepc_spectral_transformation.h   ------------------*/

#endif

/*--------------------   slepc_spectral_transformation.h   ------------------*/
