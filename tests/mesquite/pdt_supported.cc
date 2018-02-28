/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

/**
 * Check if Mesquite has the capabilities to support distributed triangulations.
 *
 * This should already be picked up correctly if we use an external version of
 * Mesquite, so this test primarily targets the bundled version in which we've
 * patched this bug.
 */


#include <Mesquite_all_headers.hpp>

#include "../tests.h"

int
main()
{
  initlog();

#if defined(DEAL_II_WITH_MESQUITE) ||         \
  (defined(DEAL_II_TRILINOS_WITH_MESQUITE) && \
   defined(DEAL_II_TRILINOS_MESQUITE_SUPPORTS_PDT))
  Mesquite2::Mesh *           mesh = nullptr;
  Mesquite2::ParallelMeshImpl parallel_mesh(mesh);
#endif

  deallog << "OK" << std::endl;
}
