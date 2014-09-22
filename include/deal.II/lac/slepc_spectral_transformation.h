// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef __deal2__slepc_spectral_transformation_h
#define __deal2__slepc_spectral_transformation_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SLEPC

#  include <deal.II/base/std_cxx11/shared_ptr.h>
#  include <deal.II/lac/exceptions.h>

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
   * <code>SLEPcWrappers::TransformationXXX</code>, where
   * <code>XXX</code> is your favourite transformation type, can then
   * be implemented in application codes in the following way for
   * <code>XXX=INVERT</code> with the solver object
   * <code>eigensolver</code>:
   * @code
   *  // Set a transformation, this one shifts the eigenspectrum by 3.142..
   *  SLEPcWrappers::TransformationShift::AdditionalData additional_data (3.142);
   *  SLEPcWrappers::TransformationShift shift (additional_data);
   *  eigensolver.set_transformation (shift);
   * @endcode
   * and later calling the <code>solve()</code> function as usual:
   * @code
   *  SolverControl solver_control (1000, 1e-9);
   *  SolverArnoldi system (solver_control, mpi_communicator);
   *  eigensolver.solve (A, B, lambda, x, size_of_spectrum);
   * @endcode
   *
   * @note These options can also be set at the commandline.
   *
   * @ingroup SLEPcWrappers
   * @author Toby D. Young 2009, 2013
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

    std_cxx11::shared_ptr<TransformationData> transformation_data;
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
   * Spectrum Folding. This transformation type has been removed in
   * SLEPc 3.5.0 and thus cannot be used in the newer versions.
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

#endif // DEAL_II_WITH_SLEPC

/*--------------------   slepc_spectral_transformation.h   ------------------*/

#endif

/*--------------------   slepc_spectral_transformation.h   ------------------*/
