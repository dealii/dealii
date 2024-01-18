#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

#include <vector>

#include "../tests.h"

// Test that Jacobians are correctly scaled on faces. In particular - this
// didn't used to work with MappingFEField since it did not apply the sqrt(2)
// factor in 2D correctly.


template <int dim, int spacedim = dim>
void
test(const Triangulation<dim, spacedim> &tria,
     const Mapping<dim, spacedim>       &mapping,
     const Quadrature<dim - 1>          &face_quadrature)
{
  Assert(tria.get_reference_cells().size() == 1, ExcNotImplemented());
  const ReferenceCell reference_cell = tria.get_reference_cells().front();
  Assert(reference_cell == ReferenceCells::get_simplex<dim>() ||
           reference_cell == ReferenceCells::get_hypercube<dim>(),
         ExcNotImplemented());

  FE_Nothing<dim, spacedim>   fe_nothing(reference_cell);
  FEFaceValues<dim, spacedim> face_values(mapping,
                                          fe_nothing,
                                          face_quadrature,
                                          update_JxW_values);
  for (const auto &cell : tria.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      {
        face_values.reinit(cell, face);
        double measure = 0.0;
        for (unsigned int q = 0; q < face_quadrature.size(); ++q)
          measure += face_values.JxW(q);
        Assert(std::abs(face->measure() - measure) < 1e-10 * measure,
               ExcMessage("face measures should be equal"));
        deallog << "face measure = " << face->measure() << ", " << measure
                << std::endl;
      }
}

int
main(int argc, char **argv)
{
  initlog();

  {
    Triangulation<2> tria_hex;
    GridGenerator::hyper_cube(tria_hex);
    Triangulation<2> tria_simplex;
    GridGenerator::convert_hypercube_to_simplex_mesh(tria_hex, tria_simplex);

    // easy test - should work with MappingFE
    deallog << "Test with MappingFE" << std::endl;
    test(
      tria_simplex,
      ReferenceCells::get_simplex<2>().template get_default_linear_mapping<2>(),
      QGaussSimplex<1>(2));

    // harder test - at the time this test was initially written this did not
    // work
    deallog << "Test with MappingFEField" << std::endl;
    FESystem<2>   fe(FE_SimplexP<2>(2), 2);
    DoFHandler<2> dof_handler(tria_simplex);
    dof_handler.distribute_dofs(fe);

    Vector<double> position(dof_handler.n_dofs());
    VectorTools::interpolate(dof_handler,
                             Functions::IdentityFunction<2>(),
                             position);
    MappingFEField<2, 2, Vector<double>> mapping(dof_handler, position);
    test(tria_simplex, mapping, QGaussSimplex<1>(2));
  }

  {
    Triangulation<3> tria_hex;
    GridGenerator::hyper_cube(tria_hex);
    Triangulation<3> tria_simplex;
    GridGenerator::convert_hypercube_to_simplex_mesh(tria_hex, tria_simplex);

    // easy test - should work with MappingFE
    deallog << "Test with MappingFE" << std::endl;
    test(
      tria_simplex,
      ReferenceCells::get_simplex<3>().template get_default_linear_mapping<3>(),
      QGaussSimplex<2>(2));

    // harder test - at the time this test was initially written this did not
    // work
    deallog << "Test with MappingFEField" << std::endl;
    FESystem<3>   fe(FE_SimplexP<3>(2), 3);
    DoFHandler<3> dof_handler(tria_simplex);
    dof_handler.distribute_dofs(fe);

    Vector<double> position(dof_handler.n_dofs());
    VectorTools::interpolate(dof_handler,
                             Functions::IdentityFunction<3>(),
                             position);
    MappingFEField<3, 3, Vector<double>> mapping(dof_handler, position);
    test(tria_simplex, mapping, QGaussSimplex<2>(2));
  }
}
