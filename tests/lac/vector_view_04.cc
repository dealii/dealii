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

// check VectorView::checkReinit2

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>
#include <cmath>
#include <fstream>
#include <iomanip>

template<typename number>
  void
  checkReinit2(Vector<number> &V)
  {
    deallog
        << "Reinit a ReadWrite VectorView<number> with Vector<number> and const size"
        << std::endl;

    deallog
        << "Creating dummy Vector<number> of size V.size() and filling with zeros"
        << std::endl;

    Vector<number> _V(V.size());
    for (unsigned int i = 0; i < _V.size(); i++)
      _V(i) = 0;

    deallog << "Creating VectorView<number> pointing to dummy Vector<number>"
        << std::endl;
    VectorView<number> VV(_V.size(), _V.begin());

    deallog << "Printing dummy Vector<number>" << std::endl;
    for (unsigned int i = 0; i < _V.size(); ++i)
      deallog << _V(i) << '\t';
    deallog << std::endl;

    deallog << "Printing VectorView<number> pointing to dummy Vector<number>"
        << std::endl;
    for (unsigned int i = 0; i < VV.size(); ++i)
      deallog << VV(i) << '\t';
    deallog << std::endl;

    deallog << "Reinit VectorView<number> to half of Vector<number>"
        << std::endl;
    VV.reinit(V.size() / 2, V.begin());

    deallog << "Printing Vector<number>" << std::endl;
    for (unsigned int i = 0; i < V.size(); ++i)
      deallog << V(i) << '\t';
    deallog << std::endl;

    deallog << "Printing VectorView<number> pointing to half of Vector<number>"
        << std::endl;
    for (unsigned int i = 0; i < VV.size(); ++i)
      deallog << VV(i) << '\t';
    deallog << std::endl;

    deallog
        << "Incrementing Vector<number> elements using Read-write handle of VectorView<number>"
        << std::endl;
    for (unsigned int i = 0; i < VV.size(); ++i)
      VV(i) = VV(i) + 1;

    deallog << "Printing modified Vector<number>" << std::endl;
    for (unsigned int i = 0; i < V.size(); ++i)
      deallog << V(i) << '\t';
    deallog << std::endl;
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

  checkReinit2<double>(V1);
}

