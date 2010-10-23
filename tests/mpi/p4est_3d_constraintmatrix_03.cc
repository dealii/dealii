//---------------------------------------------------------------------------
//    $Id: 2d_coarse_grid_01.cc 17444 2008-10-31 19:35:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check ConstraintMatrix.distribute() for a distributed mesh
// with Trilinos; manual check of the graphical output...
// Mesh: shell with random refinement

#include "../tests.h"
#include "coarse_grid_common.h"
#include <base/logstream.h>
#include <base/tensor.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <grid/filtered_iterator.h>

#include <grid/tria.h>
#include <distributed/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <distributed/solution_transfer.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_boundary_lib.h>
#include <numerics/data_out.h>
#include <lac/trilinos_block_vector.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>

#include <fstream>
#include <sstream>

const double R0      = 0.5;//6371000.-2890000.;
const double R1      = 1.0;//6371000.-  35000.;
const double T0      = 1.0;
const double T1 = 2.0;



template<int dim>
class FilteredDataOut : public DataOut<dim>
{
  public:
    FilteredDataOut (const unsigned int subdomain_id)
		    :
		    subdomain_id (subdomain_id)
      {}

    virtual typename DoFHandler<dim>::cell_iterator
    first_cell ()
      {
	typename DoFHandler<dim>::active_cell_iterator
	  cell = this->dofs->begin_active();
	while ((cell != this->dofs->end()) &&
	       (cell->subdomain_id() != subdomain_id))
	  ++cell;

	return cell;
      }

    virtual typename DoFHandler<dim>::cell_iterator
    next_cell (const typename DoFHandler<dim>::cell_iterator &old_cell)
      {
	if (old_cell != this->dofs->end())
	  {
	    const IteratorFilters::SubdomainEqualTo
	      predicate(subdomain_id);

	    return
	      ++(FilteredIterator
		 <typename DoFHandler<dim>::active_cell_iterator>
		 (predicate,old_cell));
	  }
	else
	  return old_cell;
      }

  private:
    const unsigned int subdomain_id;
};

template <int dim>
class TemperatureInitialValues : public Function<dim>
{
  public:
    TemperatureInitialValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &value) const;
};



template <int dim>
double
TemperatureInitialValues<dim>::value (const Point<dim>  &p,
				      const unsigned int) const
{
  const double r = p.norm();
  const double h = R1-R0;

  const double s = (r-R0)/h;
// see http://www.wolframalpha.com/input/?i=plot+(sqrt(x^2%2By^2)*0.95%2B0.05*sin(6*atan2(x,y))),+x%3D-1+to+1,+y%3D-1+to+1

  double s_mod = s*0.95 + 0.05*sin(6.0*atan2(p(0),p(1)));
//alternative:    http://www.wolframalpha.com/input/?i=plot+atan((sqrt(x^2%2By^2)*0.95%2B0.05*sin(6*atan2(x,y))-0.5)*10)/pi%2B0.5,+x%3D-1+to+1,+y%3D-1+to+1
//    s_mod = atan((s_mod-0.5)*10.0)/dealii::numbers::PI+0.5;

//    return T1+(T0-T1)*(1.0-s_mod);

				   //old:
//    return T1+(T0-T1)*((1-s)*(1-s)); //old
//    return T1+(T0-T1)*((1-s)); //linear
  return p(0)*T1+p(1)*(T0-T1); //simple
}


template <int dim>
void
TemperatureInitialValues<dim>::vector_value (const Point<dim> &p,
					     Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = TemperatureInitialValues<dim>::value (p, c);
}


template<int dim>
void test()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);


//  GridGenerator::hyper_cube(tr);
  
  
  GridGenerator::hyper_shell (tr,
			      Point<dim>(),
			      R0,
			      R1,
			      12,
			      true);
  static HyperShellBoundary<dim> boundary;
//  tr.set_boundary (0, boundary);
//tr.set_boundary (1, boundary);
  
  tr.refine_global (1);
  if (1)
    for (unsigned int step=0; step<5;++step)
      {
	typename Triangulation<dim>::active_cell_iterator
	  cell = tr.begin_active(),
	  endc = tr.end();
      
	for (; cell!=endc; ++cell)
	  if (std::rand()%42==1)
	    cell->set_refine_flag ();
      
	tr.execute_coarsening_and_refinement ();
      }
 
  DoFHandler<dim> dofh(tr);

  static FE_Q<dim> fe(3);
  
  dofh.distribute_dofs (fe);

  IndexSet owned_set = dofh.locally_owned_dofs();
  
  IndexSet dof_set;
  DoFTools::extract_locally_active_dofs (dofh, dof_set);

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);

  VectorTools::interpolate(* static_cast<DoFHandler<dim>* >(&dofh),
			   TemperatureInitialValues<dim>(),
			   x);
  x.compress();
  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);
  x_rel = x;
  
  for (unsigned int steps=0;steps<7;++steps)
    { 
      {
	typename Triangulation<dim>::active_cell_iterator
	  cell = tr.begin_active(),
	  endc = tr.end();
	
	for (; cell!=endc; ++cell)
	  if (!cell->is_artificial() && !cell->is_ghost())
	    {
	      if (std::rand()%12==1)
		cell->set_refine_flag ();
	      else if (std::rand()%7==1)
		cell->set_coarsen_flag ();
	    }	
      }
      for (typename Triangulation<dim>::cell_iterator
	     cell = tr.begin();
	   cell != tr.end(); ++cell)
	{
	  if (!cell->has_children())
	    continue;
	  
	  bool coarsen_me = false;
	  for (unsigned int i=0;i<cell->n_children();++i)
	    if (cell->child(i)->coarsen_flag_set())
	      {
		coarsen_me = true;
		break;
	      }
	  if (coarsen_me)
	    for (unsigned int i=0;i<cell->n_children();++i)
	      {
		if (cell->child(i)->is_artificial())
		  ;//		  std::cout << "art" << std::endl;
		else if (cell->child(i)->is_ghost())
		  ;//		  std::cout << "ghost" << std::endl;
		else
		  {

		    cell->child(i)->clear_refine_flag();
		    cell->child(i)->set_coarsen_flag();
		  }


	      }

	}
      
      parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::Vector>
	trans(dofh);
      tr.prepare_coarsening_and_refinement();
	

      trans.prepare_for_coarsening_and_refinement(x_rel);
	
      tr.execute_coarsening_and_refinement ();
	

      static FE_Q<dim> fe(1);
  
      dofh.distribute_dofs (fe);

      owned_set = dofh.locally_owned_dofs();
  
      DoFTools::extract_locally_active_dofs (dofh, dof_set);

      DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

      x.reinit(owned_set, MPI_COMM_WORLD);

      trans.interpolate(x);
      
      x_rel.reinit(relevant_set, MPI_COMM_WORLD);
      x_rel = 0;
      x_rel.compress();
  
      ConstraintMatrix cm(relevant_set);
      DoFTools::make_hanging_node_constraints (* static_cast<DoFHandler<dim>* >(&dofh), cm);
/*  std::vector<bool> velocity_mask (dim+1, true);
  
    velocity_mask[dim] = false;
				    
    VectorTools::interpolate_boundary_values (static_cast<const DoFHandler<dim>&>(dofh),
    0,
    ZeroFunction<dim>(1),
    cm,
    velocity_mask);				   */

      cm.close ();

      cm.distribute(x);
      x_rel = x;
    }


  TrilinosWrappers::MPI::Vector x_ref;
  x_ref.reinit(owned_set, MPI_COMM_WORLD);
  
  VectorTools::interpolate(* static_cast<DoFHandler<dim>* >(&dofh),
			   TemperatureInitialValues<dim>(),
			   x_ref);
  x_ref.compress();
  
  x_ref -= x;
  double err = x_ref.linfty_norm();
  if (err>1.0e-12)
    if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
      deallog << "err:" << err << std::endl;
  
//	x_rel=x_ref; //uncomment to output error
  
  std::vector<std::string> solution_names(1,"T");
  
  FilteredDataOut<dim> data_out (tr.locally_owned_subdomain());
  data_out.attach_dof_handler (dofh);
  
  data_out.add_data_vector(x_rel, solution_names);
  
  data_out.build_patches (1);
  const std::string filename = ("p4est_3d_constraintmatrix_03/solution." +
				Utilities::int_to_string
				(tr.locally_owned_subdomain(), 4) +
				".vtu");
  std::ofstream output (filename.c_str());
  data_out.write_vtu (output);
  
  tr.set_boundary (0);
  tr.set_boundary (1);

  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Init (&argc,&argv);
#else
  (void)argc;
  (void)argv;
#endif

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("p4est_3d_constraintmatrix_03").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    test<3>();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Finalize();
#endif
}
