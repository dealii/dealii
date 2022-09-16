
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



int
main()
{
  initlog();

  const unsigned int dim       = 2;
  const std::string  inputFile = SOURCE_DIR "/merge_01.inp";

  Triangulation<dim> triangulation_(Triangulation<dim>::maximum_smoothing);

  GridIn<dim> gridIn;
  gridIn.attach_triangulation(triangulation_);
  std::ifstream gridInputStream(inputFile);
  gridIn.read(gridInputStream, GridIn<dim>::abaqus);

  FE_Q<dim>       ppFE_(1);
  DoFHandler<dim> ppDoFHandler_(triangulation_);
  ppDoFHandler_.distribute_dofs(ppFE_);

  Vector<double> materialIDs_(
    ppDoFHandler_.get_triangulation().n_active_cells());


  const std::string               fileNameBase     = "mwe-stresses";
  const DataOutBase::OutputFormat dataOutputFormat = DataOutBase::vtk;

  MappingQ1<dim> qMapping;

  DataOut<dim> dataOutPrimary, dataOutSecondary;
  dataOutPrimary.attach_dof_handler(ppDoFHandler_);
  dataOutSecondary.attach_dof_handler(ppDoFHandler_);

  dataOutPrimary.set_cell_selection(
    [](const typename Triangulation<dim>::cell_iterator &cell) {
      return (cell->is_active() && cell->material_id() == 1);
    });
  dataOutSecondary.set_cell_selection(
    [](const typename Triangulation<dim>::cell_iterator &cell) {
      return (cell->is_active() && cell->material_id() == 2);
    });

  {
    // first add the material-ID data
    dataOutPrimary.add_data_vector(materialIDs_, "abc");
    dataOutPrimary.build_patches(qMapping,
                                 2); //, DataOut<dim>::curved_inner_cells);
  }
  {
    // first add the material-ID data
    dataOutSecondary.add_data_vector(materialIDs_, "abc");
    dataOutSecondary.build_patches(qMapping, 2);
  }

  // error here:
  dataOutPrimary.merge_patches(dataOutSecondary);

  // now write to file
  dataOutPrimary.write(deallog.get_file_stream(), dataOutputFormat);
}
