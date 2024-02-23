// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

#include "../tests.h"

class Test
{
public:
  void
  generate_all_dim(const char *myname)
  {
    generate<1, 1>(myname);
    generate<2, 2>(myname);
    generate<3, 3>(myname);
  }


  template <int dim, int spacedim>
  void
  generate(const char *myname)
  {
    std::string name = myname;
    deallog << "Read " << name << std::endl;

    std::unique_ptr<FiniteElement<dim, spacedim>> fe =
      FETools::get_fe_by_name<dim, spacedim>(name);

    deallog << "Generated :" << std::endl;
    deallog << fe->get_name() << std::endl;
  }


  void
  generate_all_codim(const char *myname)
  {
    generate<1, 1>(myname);
    generate<1, 2>(myname);
    generate<2, 2>(myname);
    generate<1, 3>(myname);
    generate<2, 3>(myname);
    generate<3, 3>(myname);
  }
};

int
main()
{
  initlog();


  Test gen;
  // For some of the finite element types, their
  // implementation or instantiation for codimension != 0 is still
  // missing. In these case we just test for spacedim==dim

  // generate_all_dim
  gen.generate_all_dim("FE_Q_Hierarchical(1)");
  gen.generate<2, 2>("FE_ABF(0)");
  //  gen.generate<3,3>("FE_ABF(1)");
  gen.generate<2, 2>("FE_BDM(1)");
  gen.generate<3, 3>("FE_BDM(1)");
  gen.generate<2, 2>("FE_DGBDM(1)");
  gen.generate<3, 3>("FE_DGBDM(1)");
  gen.generate<2, 2>("FE_DGNedelec(1)");
  gen.generate<3, 3>("FE_DGNedelec(1)");
  gen.generate<2, 2>("FE_DGRaviartThomas(1)");
  gen.generate<3, 3>("FE_DGRaviartThomas(1)");
  gen.generate<2, 2>("FE_RaviartThomas(1)");
  gen.generate<3, 3>("FE_RaviartThomas(1)");
  gen.generate<2, 2>("FE_RaviartThomasNodal(1)");
  gen.generate<3, 3>("FE_RaviartThomasNodal(1)");
  gen.generate<2, 2>("FE_Nedelec(1)");
  gen.generate<3, 3>("FE_Nedelec(1)");
  gen.generate<2, 2>("FE_NedelecSZ(1)");
  gen.generate<3, 3>("FE_NedelecSZ(1)");
  gen.generate_all_dim("FE_DGPNonparametric(1)");
  gen.generate_all_dim("FE_DGPMonomial(1)");
  gen.generate_all_dim("FE_FaceQ(1)");
  gen.generate_all_dim("FE_FaceP(1)");
  gen.generate<2, 2>("FE_RannacherTurek(0)");
  gen.generate_all_dim("FE_Q_Hierarchical(2)");
  gen.generate<2, 2>("FE_ABF(2)");
  //  gen.generate<3,3>("FE_ABF(2)");
  gen.generate<2, 2>("FE_BDM(2)");
  gen.generate<3, 3>("FE_BDM(2)");
  gen.generate<2, 2>("FE_DGBDM(2)");
  gen.generate<3, 3>("FE_DGBDM(2)");
  gen.generate<2, 2>("FE_DGNedelec(2)");
  gen.generate<3, 3>("FE_DGNedelec(2)");
  gen.generate<2, 2>("FE_DGRaviartThomas(2)");
  gen.generate<3, 3>("FE_DGRaviartThomas(2)");
  gen.generate<2, 2>("FE_RaviartThomas(2)");
  gen.generate<3, 3>("FE_RaviartThomas(2)");
  gen.generate<2, 2>("FE_RaviartThomasNodal(2)");
  gen.generate<3, 3>("FE_RaviartThomasNodal(2)");
  gen.generate<2, 2>("FE_Nedelec(2)");
  gen.generate<3, 3>("FE_Nedelec(2)");
  gen.generate<2, 2>("FE_NedelecSZ(2)");
  gen.generate<3, 3>("FE_NedelecSZ(2)");
  gen.generate_all_dim("FE_DGPNonparametric(2)");
  gen.generate_all_dim("FE_DGPMonomial(2)");
  gen.generate_all_dim("FE_FaceQ(2)");
  gen.generate_all_dim("FE_FaceP(2)");

  // generate_all_codim
  gen.generate_all_codim("FE_Bernstein(1)");
  gen.generate_all_codim("FE_DGP(1)");
  gen.generate_all_codim("FE_DGQ(1)");
  gen.generate_all_codim("FE_DGQLegendre(1)");
  gen.generate_all_codim("FE_Nothing()");
  gen.generate_all_codim("FE_DGQArbitraryNodes(QGauss(2))");
  gen.generate_all_codim("FE_Q_Bubbles(1)");
  gen.generate_all_codim("FE_Q_DG0(1)");
  gen.generate_all_codim("FE_Q_iso_Q1(1)");
  gen.generate_all_codim("FE_Q(1)");
  gen.generate_all_codim("FE_Bernstein(2)");
  gen.generate_all_codim("FE_DGP(2)");
  gen.generate_all_codim("FE_DGQ(2)");
  gen.generate_all_codim("FE_DGQArbitraryNodes(QGauss(3))");
  gen.generate_all_codim("FE_Q_Bubbles(2)");
  gen.generate_all_codim("FE_Q_DG0(2)");
  gen.generate_all_codim("FE_Q_iso_Q1(2)");
  gen.generate_all_codim("FE_Q(2)");
  gen.generate_all_codim("FE_Bernstein(2)");
  gen.generate_all_codim("FE_DGQHermite(3)");

  // Simplex finite elements
  gen.generate_all_codim("FE_SimplexP(1)");
  gen.generate_all_codim("FE_SimplexP(2)");
  gen.generate_all_codim("FE_SimplexDGP(1)");
  gen.generate_all_codim("FE_SimplexDGP(2)");
  gen.generate_all_codim("FE_SimplexP_Bubbles(2)");

  // Mixed grids finite elements
  gen.generate<3, 3>("FE_PyramidP(1)");
  gen.generate<3, 3>("FE_PyramidDGP(1)");
  gen.generate<3, 3>("FE_WedgeP(1)");
  gen.generate<3, 3>("FE_WedgeDGP(1)");

  // systems
  gen.generate_all_dim("FESystem[FE_Q_Hierarchical(1)^2-FE_Q_Hierarchical(1)]");
  gen.generate_all_dim(
    "FESystem[FE_DGPNonparametric(1)^2-FE_Q_Hierarchical(1)]");
  gen.generate_all_dim("FESystem[FE_DGQ(1)^2-FE_Q_Hierarchical(1)]");
  gen.generate_all_dim("FESystem[FE_Q(1)^2-FE_Q_Hierarchical(1)]");

  gen.generate_all_dim(
    "FESystem[FE_Q_Hierarchical(1)^2-FE_DGPNonparametric(1)]");
  gen.generate_all_dim(
    "FESystem[FE_DGPNonparametric(1)^2-FE_DGPNonparametric(1)]");
  gen.generate_all_dim("FESystem[FE_DGQ(1)^2-FE_DGPNonparametric(1)]");
  gen.generate_all_dim("FESystem[FE_Q(1)^2-FE_DGPNonparametric(1)]");

  gen.generate_all_dim("FESystem[FE_Q_Hierarchical(1)^2-FE_DGQ(1)]");
  gen.generate_all_dim("FESystem[FE_DGPNonparametric(1)^2-FE_DGQ(1)]");
  gen.generate_all_codim("      FESystem[FE_DGQ(1)^2     -FE_DGQ(1)]");
  gen.generate_all_codim("FESystem[FE_Q(1)^2-FE_DGQ(1)]");
  gen.generate_all_codim("FESystem[FE_Q(1)^2-FE_DGQArbitraryNodes(QGauss(2))]");
  gen.generate_all_codim("FESystem[FE_DGQArbitraryNodes(QGauss(2))^2-FE_Q(1)]");

  gen.generate_all_dim("FESystem[FE_Q_Hierarchical(1)^2-FE_Q(1)]");
  gen.generate_all_dim("FESystem[FE_DGPNonparametric(1)^2-FE_Q(1)]");
  gen.generate_all_codim("FESystem[FE_DGQ(1)^2-FE_Q(1)]");
  gen.generate_all_codim("FESystem[FE_Q(1)^2-FE_Q(1)]");

  return 0;
}
