// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

// check VectorView::checkReadOnlyConstructor

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>
#include <cmath>
#include <fstream>
#include <iomanip>

template<typename number>
  void
  checkReadOnlyConstructor(const Vector<number> &V)
  {
    deallog << "Read-only constructor" << std::endl;
    VectorView<number> VV(V.size(), V.begin());

    deallog << "Printing Vector<number>" << std::endl;
    for (unsigned int i = 0; i < V.size(); ++i)
      deallog << V(i) << '\t';
    deallog << std::endl;

    deallog << "Printing VectorView<number> pointing to Vector<number>"
        << std::endl;
    for (unsigned int i = 0; i < VV.size(); ++i)
      deallog << VV(i) << '\t';
    deallog << std::endl;

    /* deallog << "Incrementing Vector<number> elements using Read-only handle of VectorView<number>" << std::endl;
     deallog << "Function fails beyond this point" << std::endl;
     for (unsigned int i=0; i<VV.size(); ++i)
     VV(i)=VV(i)+1; */
  }

int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Vector<double> V1(5);
  V1(0) = 1;
  V1(1) = 2;
  V1(2) = 3;
  V1(3) = 4;
  V1(4) = 5;

  const Vector<double> V2(V1);

  checkReadOnlyConstructor<double>(V2);
}

