#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_nothing.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Verify that we can set up a mesh consisting entirely of FE_Nothing FEs. While
// this by itself would be useless, its possible in other contexts to get
// DoFHandlers with zero DoFs on them with uncommon parallel data distributions.
// This previously crashed with a segmentation fault inside
// DoFHandlerImplementation::Policy::Sequential::renumber_dofs().

int
main()
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);
  FE_Nothing<2> fe;
  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  DoFRenumbering::random(dof_handler);
  deallog << "Number of DoFs = " << dof_handler.n_dofs() << std::endl;
}
