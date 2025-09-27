// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests FE_Nedelec finite elements of orders 0...4 on a two-cells 3D mesh
// with a shared face orientation 3. See the notes in the header file for
// more detail.

#include "fe_nedelec_face_edge_orientation.h"

int
main()
{
  unsigned char combined_face_orientation = 3;

  initlog();
  run<3, 0>(combined_face_orientation);
  run<3, 1>(combined_face_orientation);
  run<3, 2>(combined_face_orientation);
  run<3, 3>(combined_face_orientation);
  run<3, 4>(combined_face_orientation);
}
