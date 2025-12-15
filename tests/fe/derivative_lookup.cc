// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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

#include <deal.II/base/config.h>



/**
 * Test for the FE_Hermite method function 
 * get_dofs_corresponding_to_outward_normal_derivatives, which produces
 * a 2D table listing all the local DoF indices of DoFs with a non-zero 
 * value for a specified normal derivative on each element face. For
 * instance, if a @p derivative_order of 0 is passed to the function, it
 * will list all local DoFs that are non-zero on the boundaries. 
 * 
 * This function is necessary for Hermite elements since all DoFs are
 * allocated at nodes, so existing functions for applying boundary
 * conditions apply too many constraints if used with Hermite elements.
 */

  Table<2, unsigned int>
  get_dofs_corresponding_to_outward_normal_derivatives(
    const unsigned int derivative_order) const;