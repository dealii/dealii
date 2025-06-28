// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <cassert>
#include <iostream>

#include "../tests.h"

// These tests were used while building extensions for
// PetscScalar=complex. They check that std::complex<double> and
// PetscScalar complex are compatible without explicitly casting.
//
// These tests look archaic - they are the original tests used to
// figure out how PetscScalar complex works...

// Multiply a PETSc complex number by a complex number.
void
multiply_petsc_complex_by_a_number()
{
  deallog << "Check PetscScalar multiply" << std::endl;

  // make alpha
  const PetscScalar alpha = 1.0 + 2.0 * PETSC_i;

  // make beta
  const PetscScalar beta = 2.0 + 1.0 * PETSC_i;

  // multiply them together
  const PetscScalar gamma = alpha * beta;

  // Check real access
  AssertThrow(PetscRealPart(gamma) == 0., ExcInternalError());
  AssertThrow(PetscImaginaryPart(gamma) == 5., ExcInternalError());

  // This is legal too!
  AssertThrow(gamma == std::complex<double>(0., 5), ExcInternalError());

  deallog << "OK" << std::endl;
}

// Divide a PETSc complex number by a real number.
void
divide_petsc_complex_by_a_number()
{
  deallog << "Check PetscScalar division" << std::endl;

  // make alpha
  PetscScalar alpha = 1.0 + 2.0 * PETSC_i;

  // make beta <- alpha/2
  const PetscScalar beta = alpha / 2.;

  // Check real access
  AssertThrow(PetscRealPart(alpha) == 1., ExcInternalError());
  AssertThrow(PetscRealPart(beta) == 0.5, ExcInternalError());

  // operate on alpha (now alpha==beta)
  alpha /= 2.;

  // Check complex access
  AssertThrow(PetscImaginaryPart(alpha) == 1., ExcInternalError());
  AssertThrow(PetscImaginaryPart(beta) == 1., ExcInternalError());

  // This is legal too!
  AssertThrow(alpha == beta, ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a std::complex number from an PETSc complex number
void
make_std_complex_from_petsc_complex()
{
  deallog << "Check std initialised from PetscScalar" << std::endl;
  const PetscScalar          alpha = 1.0 + 2.0 * PETSC_i;
  const std::complex<double> beta  = alpha;

  // These should be the same of course
  AssertThrow(alpha == beta, ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a PETSc complex number from an std::complex number
void
make_petsc_complex_from_std_complex()
{
  deallog << "Check PetscScalar initialised from std" << std::endl;
  const std::complex<double> beta(1, 2);
  const PetscScalar          alpha = beta;

  // These should be the same of course
  AssertThrow(alpha == beta, ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a PETSc complex number directly.
void
make_petsc_complex()
{
  deallog << "Check a nonzero PetscScalar" << std::endl;

  // init
  const PetscScalar alpha = 1.0 + 2.0 * PETSC_i;

  // Test if PetscScalar is an std::complex.
  AssertThrow(alpha == std::complex<double>(1., 2.), ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a PETSc complex number directly only, check that it is
// initialised to 0+0i.
void
init_petsc_complex()
{
  deallog << "Check PetscScalar initialization" << std::endl;

  // Initialise (no argument) to zero.
  const PetscScalar alpha;
  AssertThrow(alpha == std::complex<double>(0., 0.), ExcInternalError());

  // Initialise (real argument) to zero.
  const PetscScalar beta = 0.;
  AssertThrow(beta == alpha, ExcInternalError());

  // Initialise (real+complex argument) to zero.
  const PetscScalar gamma = 0. + 0. * PETSC_i;
  AssertThrow(gamma == beta, ExcInternalError());

  // If alpha==beta==gamma, then:
  AssertThrow(alpha == gamma, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        // initialization (zero and nonzero)
        init_petsc_complex();
        make_petsc_complex();

        // initialization from std::complex (and vice versa)
        make_petsc_complex_from_std_complex();
        make_std_complex_from_petsc_complex();

        // some operators
        divide_petsc_complex_by_a_number();
        multiply_petsc_complex_by_a_number();
      }
    }

  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
