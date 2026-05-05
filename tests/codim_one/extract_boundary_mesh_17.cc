#include "deal.II/base/exception_macros.h"
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <iostream>

#include "../tests.h"

// Check the orientation is correct by using the divergence theorem

// int_{d Omega}  x \cdot n dS = dim * |Omega| = dim * 1 = dim

using namespace dealii;

template <int dim, bool simplex>
void
run(const std::string &label)
{
  Triangulation<dim> hex_tria, vol_tria;
  GridGenerator::hyper_cube(hex_tria);
  hex_tria.refine_global(2);
  if constexpr (simplex)
    GridGenerator::convert_hypercube_to_simplex_mesh(hex_tria, vol_tria);
  else
    vol_tria.copy_triangulation(hex_tria);

  Triangulation<dim - 1, dim> surf;
  GridGenerator::extract_boundary_mesh(vol_tria, surf);


  std::unique_ptr<FiniteElement<dim - 1, dim>> fe;
  std::unique_ptr<Quadrature<dim - 1>>         quad;
  if (simplex)
    {
      if constexpr (dim == 2)
        fe = std::make_unique<FE_Q<dim - 1, dim>>(1);
      else
        fe = std::make_unique<FE_SimplexP<dim - 1, dim>>(1);
      if constexpr (dim == 2)
        quad = std::make_unique<QGauss<dim - 1>>(2);
      else
        quad = std::make_unique<QGaussSimplex<dim - 1>>(2);
    }
  else
    {
      fe   = std::make_unique<FE_Q<dim - 1, dim>>(1);
      quad = std::make_unique<QGauss<dim - 1>>(2);
    }
  MappingFE<dim - 1, dim> mapping(*fe);
  FEValues<dim - 1, dim>  fev(mapping,
                             *fe,
                             *quad,
                             update_quadrature_points | update_normal_vectors |
                               update_JxW_values);

  double integral = 0.0;
  for (const auto &cell : surf.active_cell_iterators())
    {
      fev.reinit(cell);
      for (unsigned int q = 0; q < fev.n_quadrature_points; ++q)
        integral +=
          (fev.quadrature_point(q) * fev.normal_vector(q)) * fev.JxW(q);
    }

  AssertThrow(std::fabs(integral - dim) < 1e-14, ExcInternalError());
  deallog << "Divergence theorem for " << label << ": ok " << std::endl;
}

int
main()
{
  initlog();
  deallog.depth_console(0);

  run<2, false>("2D hex   ");
  run<2, true>("2D simplex");
  run<3, false>("3D hex   ");
  run<3, true>("3D simplex");
}
