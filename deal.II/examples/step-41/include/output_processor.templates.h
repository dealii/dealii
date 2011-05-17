#include "../include/output_processor.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe.h>

template <int dim, typename VECTOR>
OutputProcessor<dim,VECTOR>::OutputProcessor (const unsigned int n_mpi, 
                                              const unsigned int this_mpi) :
  n_mpi_processes(n_mpi),
  this_mpi_process(this_mpi),
  data_out(this_mpi)
{
  initialized = false;
}

template <int dim, typename VECTOR>
void OutputProcessor<dim,VECTOR>::declare_parameters (ParameterHandler &prm)
{
  prm.declare_entry ("Solution names", "u", Patterns::Anything(), 
		     "This is used to output the solution to a file.");
  
  prm.declare_entry ("Output partitioning", "false", Patterns::Bool());
  
  prm.enter_subsection("Solution Out Format");
  DataOut<dim>::declare_parameters(prm);
  prm.set("Output format", "vtk");
  prm.leave_subsection();
}
template <int dim, typename VECTOR>
void OutputProcessor<dim,VECTOR>::parse_parameters (ParameterHandler &prm)
{
  output_partitioning = prm.get_bool ("Output partitioning");

  headers = Utilities::split_string_list(prm.get ("Solution names"));

  prm.enter_subsection("Solution Out Format");
  data_out.parse_parameters(prm);
  prm.leave_subsection();

  initialized = true;
}


template <int dim, typename VECTOR>
void OutputProcessor<dim,VECTOR>::prepare_data_output(const DoFHandler<dim> &dh,
						    const std::string &filename) {
  AssertThrow(initialized, ExcNotInitialized());
  if(data_out.default_suffix() != "") {
    
      // If the output is needed and we have many processes, just output
      // the one we need *in intermediate format*.
      std::string fname = filename;
      if(n_mpi_processes > 1) {
	  fname += ("_" + Utilities::int_to_string(this_mpi_process, 2) + 
		    data_out.default_suffix(DataOutBase::deal_II_intermediate)) ;
      } else {
	  fname += data_out.default_suffix();
      }
      
      deallog << "Will write on file: " << fname.c_str() << std::endl;
      output_file.open(fname.c_str());
      AssertThrow(output_file, ExcIO());
      data_out.attach_dof_handler (dh);
      
      if(n_mpi_processes > 1) {
	  // Output the partitioning
	  if(output_partitioning) {
	      std::vector<unsigned int> partition_int (dh.get_tria().n_active_cells());
	      GridTools::get_subdomain_association (dh.get_tria(), partition_int);
	      Vector<double> partitioning(partition_int.begin(),
					  partition_int.end());
	      static Vector<double> static_partitioning;
	      static_partitioning.swap(partitioning);
	      data_out.add_data_vector (static_partitioning, "partitioning");
	  }
      }
  }
}


template <int dim, typename VECTOR>
void OutputProcessor<dim,VECTOR>::add_data_vector(const VECTOR &data_vector, 
						const std::string &desc)
{
  AssertThrow(initialized, ExcNotInitialized());
  deallog.push("AddingData");
  std::vector<std::string> dd = Utilities::split_string_list(desc);
  if(data_out.default_suffix() != "") {
    if (dd.size() ==1 ) 
      data_out.add_data_vector (data_vector, desc);
    else 
      data_out.add_data_vector (data_vector, dd);
    deallog << "Added data: " << desc << std::endl;
  }
  deallog.pop();
}


template <int dim, typename VECTOR>
void OutputProcessor<dim,VECTOR>::write_data_and_clear() {
  AssertThrow(initialized, ExcNotInitialized());
  AssertThrow(output_file, ExcIO());
  deallog.push("WritingData");
  if(data_out.default_suffix() != "") {
      data_out.build_patches();
    
      if(n_mpi_processes > 1) {
          data_out.write_deal_II_intermediate(output_file);
      } else {
          data_out.write(output_file);
      }
      deallog << "Wrote output file." << std::endl;
      data_out.clear();
      output_file.close();
      deallog << "Reset output." << std::endl;
  }
  deallog.pop();
}

template <int dim, typename VECTOR>
void OutputProcessor<dim,VECTOR>::dump_vector (const VECTOR &, 
					     const std::string & )
{
  Assert(false, ExcNotImplemented());
  // Only specializations exist.
}
