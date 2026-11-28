// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check polynomial L2 projection with FE_Q_Hierarchical<3> of degrees one and
// two on seven locally refined hyper-ball meshes with nonstandard face
// orientations. Require exact projection for polynomials in the FE space and
// finite errors and coefficient norms for higher-degree polynomials.
#include "project_common.h"

template <int dim>
void
test()
{
  for (unsigned int p = 1; p < 6 - dim; ++p)
    test_with_wrong_face_orientation(FE_Q_Hierarchical<dim>(p), p, 0, true);
}
