#include "../include/vector_space.h"
#include "../include/utilities.h"
#include <base/logstream.h>
#include <grid/grid_generator.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <fe/mapping_cartesian.h>
#include <fe/mapping_q.h>
#include <fe/fe.h>
#include <base/utilities.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>
#include <numerics/solution_transfer.h>
#include <numerics/vectors.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

using namespace std;

template <int dim> 
VectorSpace<dim>::VectorSpace() :
    fe(0, "Vector Space FE"),
    dh(0, "Vector Space This DH"),
    last_dh(0, "Vector Space Last DH"),
    coarse(0, "Vector Space Coarse Grid"),
    tria(0, "Vector Space This Grid"),
    last_tria(0, "Vector Space Last Grid"),
    mapping(0, "Vetor Space Mapping")
{}

template <int dim>
VectorSpace<dim>::VectorSpace(ParameterHandler &prm, 
			      Triangulation<dim> &dd,
			      const unsigned int n_mpi,
			      const unsigned int this_mpi) :
    fe(0, "Vector Space FE"),
    dh(0, "Vector Space This DH"),
    last_dh(0, "Vector Space Last DH"),
    coarse(0, "Vector Space Coarse Grid"),
    tria(0, "Vector Space This Grid"),
    last_tria(0, "Vector Space Last Grid"),
    mapping(0, "Vetor Space Mapping"),
    n_mpi_processes(n_mpi),
    this_mpi_process(this_mpi)
{
  reinit(prm, dd, n_mpi, this_mpi);
}

template <int dim> 
VectorSpace<dim>::~VectorSpace()
{
  // Get rid of the dof handler object.
  smart_delete(dh);
  smart_delete(last_dh);
  smart_delete(fe);
  coarse = 0;
  smart_delete(tria);
  smart_delete(last_tria);
  smart_delete(mapping);
}

template <int dim>
void VectorSpace<dim>::save_step() {
  deallog.push("Save");
  // Get rid of old dh and triangulation
  last_dh->clear();
  last_tria->clear();
  
  deallog << "Saving triangulation and dofs..." << std::endl;
  // Copy old triangulation (which is the one just used)
  last_tria->copy_triangulation(*tria);
  
  // swap back and forth, so that we have the same dofs...
  dh.swap(last_dh);
  redistribute_dofs();
  dh.swap(last_dh);
  deallog.pop();
}


template <int dim>
void VectorSpace<dim>::save_step_and_reinit() {
  deallog.push("Save");
  deallog << "Swapping tria and dh with old_tria and old_dh." << std::endl;
  // Swap dofhandler.
  last_dh.swap(dh);
  // Swap triangulation.
  last_tria.swap(tria);
  // Now call reinit
  reinit();
  deallog.pop();
}


template <int dim>
void VectorSpace<dim>::refine_globally() 
{
    deallog.push("GlobalRefinement");
    tria->refine_global (1);
    redistribute_dofs();
    deallog.pop();
}

template <int dim>
void VectorSpace<dim>::flag_for_refinement(Vector<float> &error) 
{
    deallog.push("Refining");
    deallog << "Strategy: " << refinement_strategy << std::endl;
    if(refinement_strategy == "global") {
      tria->set_all_refine_flags();
    } else if(refinement_strategy == "fixed_number") {
      if(max_cells > 0)
	GridRefinement::refine_and_coarsen_fixed_number(*tria, error, 
							top_fraction, bottom_fraction, max_cells);
      else
	GridRefinement::refine_and_coarsen_fixed_number(*tria, error, 
							top_fraction, bottom_fraction);
    } else if(refinement_strategy == "fixed_fraction") {
      if(max_cells > 0)
	GridRefinement::refine_and_coarsen_fixed_fraction(*tria, error, 
							  top_fraction, bottom_fraction, max_cells);
      else
	GridRefinement::refine_and_coarsen_fixed_fraction(*tria, error, 
							  top_fraction, bottom_fraction);
    } else if(refinement_strategy == "optimize") {
      GridRefinement::refine_and_coarsen_optimize(*tria, error);
    } else {
      Assert(false, ExcInternalError());
    }
    deallog.pop();
}



template <int dim>
void VectorSpace<dim>::reorder(int level,
			       std::vector<unsigned int> target_comps) {

  if(level<0) level = dh->get_tria().n_levels()-1;

  // Renumber by subdomain association only if needed
  if(n_mpi_processes > 1) {
    deallog << "Renumbering subdomain wise." << std::endl;
    DoFRenumbering::subdomain_wise(static_cast<DoFHandler<dim>&>(*dh));
  }
  
  for(unsigned int i=0; i<ordering.size(); ++i) {
    if(ordering[i] == "cuth") {
      deallog << "Renumbering with Cuthill McKee algorithm." << std::endl;
      DoFRenumbering::Cuthill_McKee(static_cast<DoFHandler<dim>&>(*dh));
    } else if(ordering[i] == "comp") {
      deallog << "Renumbering component wise." << std::endl;
      DoFRenumbering::component_wise(static_cast<DoFHandler<dim>&>(*dh), target_comps);
    } else if(ordering[i] == "upwind") {
      DoFRenumbering::downstream(static_cast<DoFHandler<dim>&>(*dh), wind);
    } else if(ordering[i] == "none") {}
    else {
      AssertThrow(false, ExcMessage("This Reordering not implemented"));
    }
  }
    
}

template <int dim>
void VectorSpace<dim>::measure_mesh() {
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  endc = tria->end();
  h_max = 0;
  h_min = 1000;
  for(cell = tria->begin_active(); cell != endc; ++cell) {
    h_max = std::max(h_max, cell->diameter());
    h_min = std::min(h_min, cell->diameter());
  }
  deallog << "Max diameter of a cell: " << h_max << std::endl
	  << "Min diameter of a cell: " << h_min << std::endl;
}

template <int dim>
void VectorSpace<dim>::reinit(ParameterHandler &prm, 
			      Triangulation<dim> &dd,
			      const unsigned int n_mpi,
			      const unsigned int this_mpi,
			      const std::string space_name) 
{
  deallog.push("VECTORSPACE");
  n_mpi_processes = n_mpi;
  this_mpi_process = this_mpi;
  
  smart_delete(dh);
  smart_delete(last_dh);
  smart_delete(fe);
  smart_delete(mapping);
  smart_delete(coarse);
  smart_delete(last_tria);

  // Parse Parameters
  parse_parameters(prm, space_name);

  // generate Finite element
  FiniteElement<dim> * fe_p=0;
    
  try {
      fe_p = FETools::get_fe_from_name<dim>(fe_name);
  } catch (...) {

    // Did not recognize the finite element to use. Try your own
    // one...
    // fe_p = new FESystem<dim>(FE_BGH<dim>(2),1, FE_DGP<dim>(1), 1);
    throw;
  }
  fe = fe_p;
  deallog << "Finite Element Space: " << fe->get_name() << endl;

  if(map_type == 0) {
    mapping= new MappingCartesian<dim>();
  } else {
    mapping = new MappingQ<dim>(map_type);
  }
    
  // Store a pointer to the triangulation
  coarse = &dd;

  // Now generate an empty triangulation.
  last_tria = new Triangulation<dim>();
  
  // And a copy of the pristine version of the triangulation.
  tria = new Triangulation<dim>();
  tria->copy_triangulation(*coarse);

  // And we initialize the current triangulation.
  initialize_mesh();

  // Generate new DoFHandlers
  dh = new MGDoFHandler<dim>(*tria);
  last_dh = new MGDoFHandler<dim>(*last_tria);
  
  // Now generates the boundary maps
  std::vector<std::vector<std::string> > bcs;
  bcs.push_back(d_bc);
  bcs.push_back(n_bc);
  bcs.push_back(o_bc);

  std::vector<std::map<char, std::vector<bool> > > bc_maps(3);
    
  for(unsigned int bcnumber=0; bcnumber<3; ++bcnumber) {
    // Get some aliases for the bcs
    std::vector<std::string> & bc = bcs[bcnumber];
    std::map<char, std::vector<bool> > & bc_map = bc_maps[bcnumber];
      
    if(bc.size()) 
      for(unsigned int i=0; i < bc.size(); ++i) {
	std::vector<std::string> id_and_comps = 
	  Utilities::split_string_list(bc[i], ':');
	AssertThrow(id_and_comps.size() == 2,
		    ExcMessage("Wrong Format for boundary indicator map."));
	std::vector<int> ids = 
	  Utilities::string_to_int(Utilities::split_string_list(id_and_comps[0]));
	std::vector<int> comps = 
	  Utilities::string_to_int(Utilities::split_string_list(id_and_comps[1]));
	  
	unsigned int n_c = fe->n_components();
	std::vector<bool> filter(n_c, false);
	// Now check that the components make sense
	for(unsigned int i=0; i<comps.size(); ++i) {
	  AssertThrow((unsigned int) comps[i] < n_c,
		      ExcIndexRange(comps[i], 0, n_c));
	  filter[ comps[i] ] = true;
	}
	// And now save these components. Merge them if the map is
	// already stored.
	for(unsigned int i = 0; i<ids.size(); ++i) {
	  // For each of the id, save the map just generated
	  if(bc_map.find(ids[i]) != bc_map.end()) {
	    // For each of the id, add the components that appear in
	    // the map just generated.
	    Assert(bc_map[ ids[i] ].size() == (unsigned int) n_c,
		   ExcDimensionMismatch(bc_map[ ids[i] ].size(), n_c));
	    for(unsigned int j=0; j<n_c; ++j) {
	      bc_map[ ids[i] ][j] = filter[j];
	    }
	  } else {
	    bc_map[ids[i]] = filter;
	  }
	}
      }
  }
  dirichlet_bc = bc_maps[0];
  neumann_bc = bc_maps[1];
  other_bc = bc_maps[2];

  deallog.pop();
}



template <int dim>
void VectorSpace<dim>::reinit() {
  // In this case, we only get rid of the current dofhandler
  // and triangulation, assuming they have already been saved.
  dh->clear();
  tria->clear();
  // And we restore the initial copy of the triangulation.
  tria->copy_triangulation(*coarse);
}

template <int dim>
void VectorSpace<dim>::initialize_mesh() 
{
  deallog.push("InitializeMesh");
  tria->refine_global (global_refinement);
  deallog << "Active Cells: "
	  << tria->n_active_cells()
	  << endl;
  deallog.pop();
}

template <int dim>
void VectorSpace<dim>::redistribute_dofs(std::vector<unsigned int> target_components ) {
    // Repartition the triangulation
    if(n_mpi_processes > 1) 
      GridTools::partition_triangulation (n_mpi_processes, *tria);
    
    // Distort the mesh if necessary
    if(distortion > 0.) tria->distort_random(distortion);
    
    // Measure it
    measure_mesh();
    
    // And Distribute degrees of freedom
    dh->distribute_dofs(*fe);
    
    // Now count what you just distributed
    count_dofs(target_components);
    
    reorder(-1, target_components);
    hang.clear();
    DoFTools::make_hanging_node_constraints (*dh, hang);
    hang.close();
    deallog << "Number of constrained DOFS: " << hang.n_constraints() << std::endl;
}

template <int dim>
void VectorSpace<dim>::declare_parameters(ParameterHandler &prm, const std::string space_name) 
{

   prm.enter_subsection(space_name);

    prm.declare_entry ("Finite element space", "FE_Q(1)", 
		       Patterns::Anything(),
		       "The finite element space to use. For vector "
		       "finite elements use the notation "
		       "FESystem[FE_Q(2)^2-FE_DGP(1)] (e.g. Navier-Stokes). ");

    prm.declare_entry ("Mapping degree", "1", Patterns::Integer(),
		       "Degree of the mapping. If 0 is used, then a Cartesian mapping is assumed.");
    prm.declare_entry ("Dof ordering", "cuth, comp", Patterns::Anything(),
		       "Ordering of the degrees of freedom: none, comp, cuth, upwind.");
    prm.declare_entry ("Wind direction", ".01, .01, 1", Patterns::Anything(),
		       "Direction of the wind for upwind ordering of the mesh. ");

    prm.declare_entry ("Dirichlet boundary map", "1:0", Patterns::Anything(),
		       "Boundary indicator, followed by semicolomn and a list"
		       " of components to which this boundary conditions apply. "
		       "More boundary indicators can be separated by semicolumn. "
		       "1:0,1,4 ; 2,4:0,2");
    prm.declare_entry ("Neumann boundary map", "2:0", Patterns::Anything(),
		       "Boundary indicators, followed by semicolomn and a list of "
		       "components to which this boundary conditions apply. More "
		       "boundary indicators can be separated by semicolumn. "
		       "1:0,1,4 ; 2,4:0,2");
    prm.declare_entry ("Other boundary map", "3:0", Patterns::Anything(),
		       "Boundary indicator, followed by semicolomn and a list of "
		       "components to which this boundary conditions apply. More "
		       "boundary indicators can be separated by semicolumn. "
		       "1:0,1,4 ; 2,4:0,2");
    
    prm.enter_subsection("Grid Parameters");
    
    prm.declare_entry ("Global refinement", "4", Patterns::Integer());
    prm.declare_entry ("Distortion coefficient", "0", Patterns::Double(),
		       "If this number is greater than zero, the mesh is distorted"
		       " upon refinement in order to disrupt its structureness.");
    
    prm.declare_entry("Refinement strategy", 
		      "fixed_number", Patterns::Selection("fixed_number|fixed_fraction|optimize|global"),
		      "fixed_number: the Top/Bottom threshold fraction of cells are flagged for "
		      "refinement/coarsening. "
		      "fixed_fraction: the cells whose error is Top/Bottom fraction of the total "
		      "are refined/coarsened. optmized: try to reach optimal error distribution, "
		      "assuming error is divided by 4 upon refining. global: refine all cells.");
    prm.declare_entry("Bottom fraction", ".3", Patterns::Double());
    prm.declare_entry("Top fraction", ".3", Patterns::Double());
    prm.declare_entry("Max number of cells", "0", Patterns::Integer(),
		      "A number of zero means no limit. ");
    prm.leave_subsection();
  
    
    prm.leave_subsection();
}

    template <int dim>
void VectorSpace<dim>::parse_parameters(ParameterHandler &prm, const std::string space_name) 
{
    prm.enter_subsection(space_name);

    fe_name = prm.get ("Finite element space"); 
    map_type = prm.get_integer("Mapping degree");
    std::string all_ordering = prm.get ("Dof ordering");

    d_bc =  Utilities::split_string_list(prm.get ("Dirichlet boundary map"), ';');
    n_bc =  Utilities::split_string_list(prm.get ("Neumann boundary map"), ';');
    o_bc =  Utilities::split_string_list(prm.get ("Other boundary map"), ';');
    
    std::vector<std::string> wind_str = 
	Utilities::split_string_list(prm.get ("Wind direction") );
    for(unsigned int i=0; (i<wind_str.size()) && (i<dim); ++i) 
	    sscanf(wind_str[i].c_str(), "%lf", &wind[i]);
				     
    prm.enter_subsection("Grid Parameters");
    
    global_refinement		= prm.get_integer ("Global refinement");
    distortion			= prm.get_double("Distortion coefficient");
    
    refinement_strategy		= prm.get("Refinement strategy");
    enable_local_refinement	= !(refinement_strategy == "global");
    
    bottom_fraction		= prm.get_double("Bottom fraction");
    top_fraction		= prm.get_double("Top fraction");
    max_cells			= prm.get_integer("Max number of cells");
    prm.leave_subsection();

    prm.leave_subsection();

    ordering = Utilities::split_string_list(all_ordering);
}


template <int dim>
void VectorSpace<dim>::interpolate_dirichlet_bc(const Function<dim> & f,  
						std::map<unsigned int, double> & bvalues) 
{
  deallog.push("DBC");
  std::map<char, std::vector<bool> >::iterator 
    dmap = dirichlet_bc.begin(),
    dmapend = dirichlet_bc.end();
  
  unsigned int last_counted = 0;
  for(; dmap!=dmapend; ++dmap) {
    char id = dmap->first;
    std::vector<bool> &filter = dmap->second;
    
        deallog << (int) id << " :";
    const unsigned int n_components = fe->n_components();
    for(unsigned int i=0; i < n_components; ++i) {
      if(filter[i]) ;//deallog << i << ", ";
    }
    VectorTools::interpolate_boundary_values(*mapping, (const DoFHandler<dim>&) *dh, 
					     // Dirichlet boundary only...
					     id, f, bvalues, filter);
    deallog << " #: " << bvalues.size() - last_counted << std::endl;
    last_counted = bvalues.size();
  }
  deallog << "Total DBC dofs: " << bvalues.size() << std::endl;
  deallog.pop();
}



template <int dim>
void VectorSpace<dim>::count_dofs(std::vector<unsigned int> target_components) {
  // Count dofs per processor
  local_dofs_per_process.resize(n_mpi_processes);
  
  // Output dofs per processor
  deallog << "Number of DOFS: " << dh->n_dofs() 
	  << ", (by partition:";
  for (unsigned int p=0; p<n_mpi_processes; ++p) {
    unsigned int nl = DoFTools::
      count_dofs_with_subdomain_association (*dh, p);
    local_dofs_per_process[p] = nl;
    deallog << (p==0 ? ' ' : '+') << nl;
  }    
  deallog << ")" << std::endl;
  
  // Count number of blocks
  number_of_blocks = ( target_components.size() ? 
		       target_components[target_components.size()-1]+1 : 1 );
  
  deallog << "Number of blocks: " << number_of_blocks << std::endl;
  
  if( number_of_blocks > 1 ) {
    // Count dofs per blocks.
    dofs_per_block.resize( number_of_blocks );

    DoFTools::count_dofs_per_component (*dh, dofs_per_block, false, target_components);
    
    while(dofs_per_block.size() > number_of_blocks) {
      deallog << "Old Deal.II!!!!"<< std::endl;
      dofs_per_block.pop_back();
    }
      
    // Output the dofs per block
    deallog << "Dofs by Block :";
    for (unsigned int p=0; p<number_of_blocks; ++p) {
      deallog << (p==0 ? ' ' : '+') << dofs_per_block[p];
    }    
    deallog << std::endl;
    
    // Now count dofs per block per processor
    dofs_per_block_per_process.resize(number_of_blocks,
				      std::vector<unsigned int>(n_mpi_processes, 0));
    dofs_per_process_per_block.resize(n_mpi_processes,
				      std::vector<unsigned int>(number_of_blocks, 0));

    // We now have a vector that contains all subdomain ids for each dof
    std::vector<unsigned int> subdomain_association(dh->n_dofs());
    DoFTools::get_subdomain_association (*dh, subdomain_association);
    
    // Now cycle on each block to find out who belongs here
    for(unsigned int block = 0; block < number_of_blocks; ++block) {
      
      // Select dofs belonging to this block
      deallog << "Block " << block << ":";
      std::vector<bool> mask(fe->n_components(), false);
      for(unsigned int i=0; i<mask.size(); ++i) {
	if(target_components[i] == block) mask[i] = true;
	deallog << ( i==0 ? ' ' : ',' ) << mask[i];
      }
      deallog << ": ";
      
      std::vector<bool> this_block_dofs(dh->n_dofs());
      DoFTools::extract_dofs(*dh, mask, this_block_dofs);
      
      // Reset counter
      for(unsigned int j=0; j<n_mpi_processes; ++j) {
	dofs_per_block_per_process[block][j] = 0;
	dofs_per_process_per_block[j][block] = 0;
      }
      
      // Count dofs on this block
      for(unsigned int i=0; i< dh->n_dofs(); ++i)
	if(this_block_dofs[i]) {
	  ++dofs_per_block_per_process[block][subdomain_association[i]];
	  ++dofs_per_process_per_block[subdomain_association[i]][block];
	}
      unsigned int check_pb=0;
      unsigned int check_pp=0;
      
      // Output the dofs per block per processor
      for(unsigned int j=0; j<n_mpi_processes; ++j) {
	deallog << (j==0 ? ' ' : '+') << dofs_per_block_per_process[block][j];
	check_pb += dofs_per_block_per_process[block][j];
	check_pp += dofs_per_process_per_block[j][block];
      }    
      deallog << std::endl;
      Assert(check_pb == dofs_per_block[block],
	     ExcDimensionMismatch(check_pb, dofs_per_block[block]));
      Assert(check_pp == dofs_per_block[block],
	     ExcDimensionMismatch(check_pp, dofs_per_block[block]));
    }
  } else {
    //
    dofs_per_block.resize(1);
    dofs_per_block[0] = dh->n_dofs();
    // 
    dofs_per_block_per_process.resize(1);
    dofs_per_block_per_process[0] = local_dofs_per_process;
    dofs_per_process_per_block.resize(n_mpi_processes);
    for(unsigned int j=0; j<n_mpi_processes; ++j)
      dofs_per_process_per_block[j] = 
	std::vector<unsigned int>(1, local_dofs_per_process[j]);
      
  }
  
  // Save the number of cells.
  number_of_cells = tria->n_active_cells();
  deallog << "Number of CELLS: " << number_of_cells << endl;
  
  // Now save the number of local dofs...
  n_local_dofs = local_dofs_per_process[this_mpi_process];
}
