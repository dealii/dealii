// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat

// Test the different FEValuesBase::get_function_values

#include <fstream>

#include <base/quadrature_lib.h>
#include <base/logstream.h>

#include <lac/vector.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>

#include <grid/grid_generator.h>


// Call this function with a system consisting of several copies of
// the SAME element
template<int dim>
void vector_values(const FiniteElement<dim>& fe)
{
  Assert(fe.n_base_elements() == 1, ExcNotImplemented());
  deallog.push(fe.get_name());
  
  QTrapez<dim> quadrature;
  std::vector<unsigned int> renumbering(fe.dofs_per_cell);
  std::vector<std::vector<unsigned int> > component_start;
  FETools::compute_component_wise(fe, renumbering, component_start);

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof);
  
  Vector<float> v(dof.n_dofs());
  for (unsigned int i=0;i<v.size();++i)
    v(i) = i;
  
  FEValues<dim> feval(fe, quadrature, update_values);
  std::vector<Vector<double> > local(quadrature.n_quadrature_points,
				     Vector<double>(fe.n_components()));
  
  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  const typename DoFHandler<dim>::active_cell_iterator end = dof.end();

  unsigned int cell_no = 0;
  while (cell != end)
    {
      deallog << "Cell " << cell_no++ << std::endl;
      feval.reinit(cell);
      feval.get_function_values(v, local);
      for (unsigned int c=0;c<fe.n_components();++c)
	{
	  deallog << "Component " << c;
	  for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
	    deallog << '\t' << (int) local[k](c);
	  deallog << std::endl;
	}
      ++cell;
    }
  deallog.pop();
}

template<int dim>
void test_vectors()
{
  FE_Q<dim> q1(1);

  FESystem<dim> q1_3(q1,3);
  vector_values(q1_3);
}


int main()
{
  std::ofstream logfile("function.output");
  deallog.attach(logfile);
//  deallog.depth_console(0);

  test_vectors<2>();
  
  return 0;
}
