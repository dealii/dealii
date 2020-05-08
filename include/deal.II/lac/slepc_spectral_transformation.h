// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_slepc_spectral_transformation_h
#  define dealii_slepc_spectral_transformation_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_SLEPC

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_solver.h>

#    include <petscksp.h>

#    include <slepceps.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#    ifndef DOXYGEN
namespace PETScWrappers
{
  // forward declarations
  class SolverBase;
} // namespace PETScWrappers
#    endif

namespace SLEPcWrappers
{
  /**
   * Base class for spectral transformation classes using the SLEPc solvers
   * which are selected based on flags passed to the spectral transformation.
   *
   * <code>SLEPcWrappers::TransformationXXX</code>, where <code>XXX</code> is
   * your favourite transformation type, can then be implemented in
   * application codes in the following way for <code>XXX=INVERT</code> with
   * the solver object <code>eigensolver</code>:
   * @code
   * // Set a transformation, this one shifts the eigenspectrum by 3.142..
   * SLEPcWrappers::TransformationShift::AdditionalData
   *   additional_data(3.142);
   * SLEPcWrappers::TransformationShift shift(mpi_communicator,additional_data);
   * eigensolver.set_transformation(shift);
   * @endcode
   * and later calling the <code>solve()</code> function as usual:
   * @code
   * SolverControl solver_control (1000, 1e-9);
   * SolverArnoldi system (solver_control, mpi_communicator);
   * eigensolver.solve (A, B, lambda, x, size_of_spectrum);
   * @endcode
   *
   * @note These options can also be set at the command line.
   *
   * @ingroup SLEPcWrappers
   * @author Toby D. Young 2009, 2013; and Denis Davydov 2015.
   */
  class TransformationBase
  {
  protected:
    /**
     * Constructor.
     */
    TransformationBase(const MPI_Comm &mpi_communicator);

  public:
    /**
     * Destructor.
     */
    virtual ~TransformationBase();

    /**
     * Set a flag to indicate how the transformed matrices are being stored in
     * the spectral transformations.
     *
     * The possible values are given by the enumerator STMatMode in the SLEPc
     * library
     * http://www.grycap.upv.es/slepc/documentation/current/docs/manualpages/ST/STMatMode.html
     */
    void
    set_matrix_mode(const STMatMode mode);

    /**
     * Set solver to be used when solving a system of linear algebraic
     * equations inside the eigensolver.
     */
    void
    set_solver(const PETScWrappers::SolverBase &solver);

  protected:
    /**
     * SLEPc spectral transformation object.
     */
    ST st;

    // Make the solver class a friend, since it needs to set spectral
    // transformation object.
    friend class SolverBase;
  };

  /**
   * An implementation of the transformation interface using the SLEPc Shift.
   *
   * @ingroup SLEPcWrappers
   * @author Toby D. Young 2009
   */
  class TransformationShift : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the shift parameter to zero.
       */
      AdditionalData(const double shift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;
    };


    /**
     * Constructor.
     */
    TransformationShift(const MPI_Comm &      mpi_communicator,
                        const AdditionalData &data = AdditionalData());


  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the transformation interface using the SLEPc Shift
   * and Invert.
   *
   * @ingroup SLEPcWrappers
   * @author Toby D. Young 2009
   */
  class TransformationShiftInvert : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the shift parameter to zero.
       */
      AdditionalData(const double shift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;
    };


    /**
     * Constructor.
     */
    TransformationShiftInvert(const MPI_Comm &      mpi_communicator,
                              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    // Make the solver class a friend, since it may need to set target
    // equal the provided shift value.
    friend class SolverBase;
  };

  /**
   * An implementation of the transformation interface using the SLEPc
   * Spectrum Folding. This transformation type has been removed in SLEPc
   * 3.5.0 and thus cannot be used in the newer versions.
   *
   * @ingroup SLEPcWrappers
   * @author Toby D. Young 2009
   */
  class TransformationSpectrumFolding : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the shift parameter to zero.
       */
      AdditionalData(const double shift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;
    };


    /**
     * Constructor.
     */
    TransformationSpectrumFolding(
      const MPI_Comm &      mpi_communicator,
      const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the transformation interface using the SLEPc Cayley.
   *
   * @ingroup SLEPcWrappers
   * @author Toby D. Young 2009
   */
  class TransformationCayley : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. Requires two shift parameters
       */
      AdditionalData(const double shift_parameter     = 0,
                     const double antishift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;

      /**
       * Antishift parameter.
       */
      const double antishift_parameter;
    };


    /**
     * Constructor.
     */
    TransformationCayley(const MPI_Comm &      mpi_communicator,
                         const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };

} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_SLEPC

/*--------------------   slepc_spectral_transformation.h   ------------------*/

#endif

/*--------------------   slepc_spectral_transformation.h   ------------------*/
