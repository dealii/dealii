/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/mg_dof.h>
#include <grid/dof_levels.h>
#include <grid/mg_dof_accessor.h>
#include <grid/dof_constraints.h>
#include <grid/tria_levels.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>
#include <lac/dsmatrix.h>

#include <algorithm>



/* ------------------------ MGVertexDoFs ----------------------------------- */

template <int dim>
MGDoFHandler<dim>::MGVertexDoFs::MGVertexDoFs () :
		coarsest_level (1<<30),
		finest_level (0),
		indices (0)
{};



template <int dim>
void MGDoFHandler<dim>::MGVertexDoFs::init (const unsigned int cl,
					    const unsigned int fl,
					    const unsigned int dofs_per_vertex) {
  Assert (indices == 0, ExcInternalError());

  coarsest_level = cl;
  finest_level   = fl;

  const unsigned int n_levels = finest_level-coarsest_level + 1;
  
  indices = new int[n_levels * dofs_per_vertex];
  Assert (indices != 0, ExcNoMemory ());

  for (unsigned int i=0; i<n_levels*dofs_per_vertex; ++i)
    indices[i] = -1;
};



template <int dim>
MGDoFHandler<dim>::MGVertexDoFs::~MGVertexDoFs () {
  Assert (indices != 0, ExcInternalError ());
  
  delete[] indices;
};



template <int dim>
inline
unsigned int MGDoFHandler<dim>::MGVertexDoFs::get_coarsest_level () const {
  return coarsest_level;
};



template <int dim>
inline
unsigned int MGDoFHandler<dim>::MGVertexDoFs::get_finest_level () const {
  return finest_level;
};






/* ------------------------ MGDoFHandler ------------------------------------- */

template <int dim>
MGDoFHandler<dim>::MGDoFHandler (Triangulation<dim> *tria) :
		DoFHandler<dim> (tria)
{};



template <int dim>
MGDoFHandler<dim>::~MGDoFHandler () {};






#if deal_II_dimension == 1

template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::begin_raw (const unsigned int level) const {
  return begin_raw_line (level);
};



template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::begin (const unsigned int level) const {
  return begin_line (level);
};



template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::begin_active (const unsigned int level) const {
  return begin_active_line (level);
};



template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::end () const {
  return end_line ();
};




template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::last_raw () const {
  return last_raw_line ();
};



template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::last_raw (const unsigned int level) const {
  return last_raw_line (level);
};



template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::last () const {
  return last_line ();
};



template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::last (const unsigned int level) const {
  return last_line (level);
};



template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::last_active () const {
  return last_active_line ();
};



template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::last_active (const unsigned int level) const {
  return last_active_line (level);
};



template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::begin_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::face_iterator
MGDoFHandler<1>::begin_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::active_face_iterator
MGDoFHandler<1>::begin_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::end_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::last_raw_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::last_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::face_iterator
MGDoFHandler<1>::last_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::face_iterator
MGDoFHandler<1>::last_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::active_face_iterator
MGDoFHandler<1>::last_active_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
MGDoFDimensionInfo<1>::active_face_iterator
MGDoFHandler<1>::last_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};




template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::begin_raw_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::begin_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::begin_active_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::end_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::last_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::last_active_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::last_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::last_active_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};

#endif



#if deal_II_dimension == 2

template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::begin_raw (const unsigned int level) const {
  return begin_raw_quad (level);
};



template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::begin (const unsigned int level) const {
  return begin_quad (level);
};



template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::begin_active (const unsigned int level) const {
  return begin_active_quad (level);
};



template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::end () const {
  return end_quad ();
};



template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::last_raw () const {
  return last_raw_quad ();
};



template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::last_raw (const unsigned int level) const {
  return last_raw_quad (level);
};



template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::last () const {
  return last_quad ();
};



template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::last (const unsigned int level) const {
  return last_quad (level);
};



template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::last_active () const {
  return last_active_quad ();
};



template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::last_active (const unsigned int level) const {
  return last_active_quad (level);
};


template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::begin_raw_face (const unsigned int level) const {
  return begin_raw_line (level);
};



template <>
MGDoFDimensionInfo<2>::face_iterator
MGDoFHandler<2>::begin_face (const unsigned int level) const {
  return begin_line (level);
};



template <>
MGDoFDimensionInfo<2>::active_face_iterator
MGDoFHandler<2>::begin_active_face (const unsigned int level) const {
  return begin_active_line (level);
};



template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::end_face () const {
  return end_line ();
};



template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::last_raw_face () const {
  return last_raw_line ();
};



template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::last_raw_face (const unsigned int level) const {
  return last_raw_line (level);
};



template <>
MGDoFDimensionInfo<2>::face_iterator
MGDoFHandler<2>::last_face () const {
  return last_line ();
};



template <>
MGDoFDimensionInfo<2>::face_iterator
MGDoFHandler<2>::last_face (const unsigned int level) const {
  return last_line (level);
};



template <>
MGDoFDimensionInfo<2>::active_face_iterator
MGDoFHandler<2>::last_active_face () const {
  return last_active_line ();
};



template <>
MGDoFDimensionInfo<2>::active_face_iterator
MGDoFHandler<2>::last_active_face (const unsigned int level) const {
  return last_active_line (level);
};

#endif





template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::begin_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->begin_raw_line(level)->level(),
			    tria->begin_raw_line(level)->index(),
			    this);
};



template <int dim>
typename MGDoFHandler<dim>::line_iterator
MGDoFHandler<dim>::begin_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->begin_line(level)->level(),
			tria->begin_line(level)->index(),
			this);
};



template <int dim>
typename MGDoFHandler<dim>::active_line_iterator
MGDoFHandler<dim>::begin_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->begin_active_line(level)->level(),
			       tria->begin_active_line(level)->index(),
			       this);
};



template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::begin_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->begin_raw_quad(level)->level(),
			    tria->begin_raw_quad(level)->index(),
			    this);
};



template <int dim>
typename MGDoFHandler<dim>::quad_iterator
MGDoFHandler<dim>::begin_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->begin_quad(level)->level(),
			tria->begin_quad(level)->index(),
			this);
};



template <int dim>
typename MGDoFHandler<dim>::active_quad_iterator
MGDoFHandler<dim>::begin_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->begin_active_quad(level)->level(),
			       tria->begin_active_quad(level)->index(),
			       this);
};



template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::end_line () const {
  return raw_line_iterator (tria, -1, -1, this);
};



template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::end_quad () const {
  return raw_quad_iterator (tria, -1, -1, this);
};



template <int dim>
typename MGDoFDimensionInfo<dim>::raw_cell_iterator
MGDoFHandler<dim>::end_raw (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  end() :
	  begin_raw (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::cell_iterator
MGDoFHandler<dim>::end (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::active_cell_iterator
MGDoFHandler<dim>::end_active (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::raw_face_iterator
MGDoFHandler<dim>::end_raw_face (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  end_face() :
	  begin_raw_face (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::face_iterator
MGDoFHandler<dim>::end_face (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  face_iterator(end_face()) :
	  begin_face (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::active_face_iterator
MGDoFHandler<dim>::end_active_face (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  active_face_iterator(end_face()) :
	  begin_active_face (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::raw_line_iterator
MGDoFHandler<dim>::end_raw_line (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  end_line() :
	  begin_raw_line (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::line_iterator
MGDoFHandler<dim>::end_line (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  line_iterator(end_line()) :
	  begin_line (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::active_line_iterator
MGDoFHandler<dim>::end_active_line (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  active_line_iterator(end_line()) :
	  begin_active_line (level+1));
};




template <int dim>
typename MGDoFDimensionInfo<dim>::raw_quad_iterator
MGDoFHandler<dim>::end_raw_quad (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  end_quad() :
	  begin_raw_quad (level+1));
};



template <int dim>
typename MGDoFDimensionInfo<dim>::quad_iterator
MGDoFHandler<dim>::end_quad (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  quad_iterator(end_quad()) :
	  begin_quad (level+1));
};




template <int dim>
typename MGDoFDimensionInfo<dim>::active_quad_iterator
MGDoFHandler<dim>::end_active_quad (const unsigned int level) const {
  return (level == mg_levels.size()-1 ?
	  active_quad_iterator(end_quad()) :
	  begin_active_quad (level+1));
};



template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::last_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->last_raw_line(level)->level(),
			    tria->last_raw_line(level)->index(),
			    this);
};



template <int dim>
typename MGDoFHandler<dim>::line_iterator
MGDoFHandler<dim>::last_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->last_line(level)->level(),
			tria->last_line(level)->index(),
			this);
};



template <int dim>
typename MGDoFHandler<dim>::active_line_iterator
MGDoFHandler<dim>::last_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->last_active_line(level)->level(),
			       tria->last_active_line(level)->index(),
			       this);
};



template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::last_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->last_raw_quad(level)->level(),
			    tria->last_raw_quad(level)->index(),
			    this);
};




template <int dim>
typename MGDoFHandler<dim>::quad_iterator
MGDoFHandler<dim>::last_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->last_quad(level)->level(),
			tria->last_quad(level)->index(),
			this);
};




template <int dim>
typename MGDoFHandler<dim>::active_quad_iterator
MGDoFHandler<dim>::last_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->last_active_quad(level)->level(),
			       tria->last_active_quad(level)->index(),
			       this);
};




template <int dim>
typename MGDoFHandler<dim>::raw_line_iterator
MGDoFHandler<dim>::last_raw_line () const {
  return last_raw_line (mg_levels.size()-1);
};



template <int dim>
typename MGDoFHandler<dim>::raw_quad_iterator
MGDoFHandler<dim>::last_raw_quad () const {
  return last_raw_quad (mg_levels.size()-1);
};



template <int dim>
typename MGDoFHandler<dim>::line_iterator
MGDoFHandler<dim>::last_line () const {
  return last_line (mg_levels.size()-1);
};



template <int dim>
typename MGDoFHandler<dim>::quad_iterator
MGDoFHandler<dim>::last_quad () const {
  return last_quad (mg_levels.size()-1);
};



template <int dim>
typename MGDoFHandler<dim>::active_line_iterator
MGDoFHandler<dim>::last_active_line () const {
  return last_active_line (mg_levels.size()-1);
};




template <int dim>
typename MGDoFHandler<dim>::active_quad_iterator
MGDoFHandler<dim>::last_active_quad () const {
  return last_active_quad (mg_levels.size()-1);
};






//------------------------------------------------------------------







template <int dim>
void MGDoFHandler<dim>::distribute_dofs (const FiniteElement<dim> &fe) {
				   // first distribute global dofs
  DoFHandler<dim>::distribute_dofs (fe);


				   // reserve space for the MG dof numbers
  reserve_space ();

				   // clear user flags because we will
				   // need them
  tria->clear_user_flags ();

				   // now distribute indices on each level
				   // separately
  for (unsigned int level=0; level<tria->n_levels(); ++level)
    {
      unsigned int next_free_dof = 0;   
      active_cell_iterator cell = begin(level),
			   endc = end(level);

      for (; cell != endc; ++cell) 
	next_free_dof = distribute_dofs_on_cell (cell, next_free_dof);
  
      mg_used_dofs[level] = next_free_dof;
    };
};



#if deal_II_dimension == 1

template <>
unsigned int
MGDoFHandler<1>::distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int          next_free_dof) {

				   // distribute dofs of vertices
  if (selected_fe->dofs_per_vertex > 0)
    for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
      {
	cell_iterator neighbor = cell->neighbor(v);
	
	if (neighbor.state() == valid)
	  {
					     // has neighbor already been processed?
	    if (neighbor->user_flag_set() &&
		(neighbor->level() == cell->level()))
					       // copy dofs if the neighbor is on
					       // the same level (only then are
					       // mg dofs the same)
	      {
		if (v==0) 
		  for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
		    cell->set_mg_vertex_dof_index (0, d,
						   neighbor->mg_vertex_dof_index (1, d));
		else
		  for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
		    cell->set_mg_vertex_dof_index (1, d,
						   neighbor->mg_vertex_dof_index (0, d));
		
						 // next neighbor
		continue;
	      };
	  };
	
					 // otherwise: create dofs newly
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  cell->set_mg_vertex_dof_index (v, d, next_free_dof++);
      };
  
				   // dofs of line
  if (selected_fe->dofs_per_line > 0)
    for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
      cell->set_mg_dof_index (d, next_free_dof++);

				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif


#if deal_II_dimension == 2

template <>
unsigned int
MGDoFHandler<2>::distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int          next_free_dof) {
  if (selected_fe->dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->mg_vertex_dof_index(vertex, 0) == -1)
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  cell->set_mg_vertex_dof_index (vertex, d, next_free_dof++);
    
  				   // for the four sides
  if (selected_fe->dofs_per_line > 0)
    for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
      {
	line_iterator line = cell->line(side);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (line->mg_dof_index(0) == -1)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
	    line->set_mg_dof_index (d, next_free_dof++);	    
      };
  

      				       // dofs of quad
  if (selected_fe->dofs_per_quad > 0)
    for (unsigned int d=0; d<selected_fe->dofs_per_quad; ++d)
      cell->set_mg_dof_index (d, next_free_dof++);

  
				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif




template <int dim>
unsigned int MGDoFHandler<dim>::n_dofs (const unsigned int level) const {
  Assert (level < mg_used_dofs.size(), ExcInvalidLevel(level));
  
  return mg_used_dofs[level];
};



template <int dim>
void MGDoFHandler<dim>::make_sparsity_pattern (const unsigned int  level,
					       dSMatrixStruct     &sparsity) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (sparsity.n_rows() == n_dofs(level),
	  ExcDifferentDimensions (sparsity.n_rows(), n_dofs(level)));
  Assert (sparsity.n_cols() == n_dofs(level),
	  ExcDifferentDimensions (sparsity.n_cols(), n_dofs(level)));

  const unsigned int n_dofs = selected_fe->total_dofs;
  vector<int> dofs_on_this_cell(n_dofs);
  cell_iterator cell = begin_active(level),
		endc = end(level);
  for (; cell!=endc; ++cell) 
    {
      cell->get_mg_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<n_dofs; ++i)
	for (unsigned int j=0; j<n_dofs; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    };
};




template <int dim>
void MGDoFHandler<dim>::renumber_dofs (const unsigned int      level,
				       const RenumberingMethod method,
				       const vector<int>      &starting_points) {
				   // make the connection graph
  dSMatrixStruct sparsity (n_dofs(level), max_couplings_between_dofs());
  make_sparsity_pattern (level, sparsity);
    
  int n_dofs = sparsity.n_rows();
				   // store the new dof numbers; -1 means
				   // that no new number was chosen yet
				   //
				   // the commented line is what would be the
				   // correct way to do, but gcc2.8 chokes
				   // over that. The other lines are a
				   // workaround
//  vector<int> new_number(n_dofs, -1);
  vector<int> new_number;
  new_number.resize (n_dofs, -1);
  
				   // store the indices of the dofs renumbered
				   // in the last round. Default to starting
				   // points
  vector<int> last_round_dofs (starting_points);
  
				   // delete disallowed elements
  for (unsigned int i=0; i<last_round_dofs.size(); ++i)
    if ((last_round_dofs[i]<0) || (last_round_dofs[i]>=n_dofs))
      last_round_dofs[i] = -1;
  
  remove_if (last_round_dofs.begin(), last_round_dofs.end(),
	     bind2nd(equal_to<int>(), -1));
  
				   // now if no valid points remain:
				   // find dof with lowest coordination
				   // number
  
  if (last_round_dofs.size() == 0)
    {
      int          starting_point   = -1;
      unsigned int min_coordination = n_dofs;
      for (int row=0; row<n_dofs; ++row) 
	{
	  unsigned int j;
	  for (j=sparsity.get_rowstart_indices()[row];
	       j<sparsity.get_rowstart_indices()[row+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == -1)
	      break;
					   // post-condition after loop:
					   // coordination is now
					   // j-rowstart[row]
	  if (j-sparsity.get_rowstart_indices()[row] <  min_coordination)
	    {
	      min_coordination = j-sparsity.get_rowstart_indices()[row];
	      starting_point   = row;
	    };
	};
				       // initialize the first dof
      last_round_dofs.push_back (starting_point);
    };
  

				   // store next free dof index
  int         next_free_number = 0;

				   // enumerate the first round dofs
  for (unsigned int i=0; i!=last_round_dofs.size(); ++i)
    new_number[last_round_dofs[i]] = next_free_number++;

  bool all_dofs_renumbered = false;

				   // now do as many steps as needed to
				   // renumber all dofs
  while (!all_dofs_renumbered) 
    {
				       // store the indices of the dofs to be
				       // renumbered in the next round
      vector<int> next_round_dofs;

				       // find all neighbors of the
				       // dofs numbered in the last
				       // round
      for (unsigned int i=0; i<last_round_dofs.size(); ++i)
	for (unsigned int j=sparsity.get_rowstart_indices()[last_round_dofs[i]];
	     j<sparsity.get_rowstart_indices()[last_round_dofs[i]+1]; ++j)
	  if (sparsity.get_column_numbers()[j] == -1)
	    break;
	  else
	    next_round_dofs.push_back (sparsity.get_column_numbers()[j]);
      
				       // sort dof numbers
      sort (next_round_dofs.begin(), next_round_dofs.end());

				       // delete multiple entries
      vector<int>::iterator end_sorted;
      end_sorted = unique (next_round_dofs.begin(), next_round_dofs.end());
      next_round_dofs.erase (end_sorted, next_round_dofs.end());

				       // eliminate dofs which are
				       // already numbered
      for (int s=next_round_dofs.size()-1; s>=0; --s)
	if (new_number[next_round_dofs[s]] != -1)
	  next_round_dofs.erase (&next_round_dofs[s]);

				       // check whether there are any new
				       // dofs
      all_dofs_renumbered = (next_round_dofs.size() == 0);
      if (all_dofs_renumbered)
					 // end loop if possible
	continue;
      

				       // store for each coordination
				       // number the dofs with these
				       // coordination number
      multimap<unsigned int, int> dofs_by_coordination;
      
				       // find coordination number for
				       // each of these dofs
      for (vector<int>::iterator s=next_round_dofs.begin();
	   s!=next_round_dofs.end(); ++s) 
	{
	  unsigned int coordination = 0;
	  for (unsigned int j=sparsity.get_rowstart_indices()[*s];
	       j<sparsity.get_rowstart_indices()[*s+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == -1)
	      break;
	    else
	      ++coordination;

					   // insert this dof at its
					   // coordination number
	  const pair<const unsigned int, int> new_entry (coordination, *s);
	  dofs_by_coordination.insert (new_entry);
	};
      
				       ////
      multimap<unsigned int, int>::iterator i;
      for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i) 
	new_number[i->second] = next_free_number++;

				       // after that: copy this round's
				       // dofs for the next round
      last_round_dofs = next_round_dofs;
    };

#ifdef DEBUG
				   //  test for all indices numbered
  if (find (new_number.begin(), new_number.end(), -1) != new_number.end())
    Assert (false, ExcRenumberingIncomplete());
  Assert (next_free_number == n_dofs,
	  ExcRenumberingIncomplete());
#endif

  switch (method) 
    {
      case Cuthill_McKee:
	    break;
      case reverse_Cuthill_McKee:
      {
	for (vector<int>::iterator i=new_number.begin(); i!=new_number.end(); ++i)
	  *i = n_dofs-*i;
	break;
      };
      default:
	    Assert (false, ExcNotImplemented());
    };

				   // actually perform renumbering;
				   // this is dimension specific and
				   // thus needs an own function
  do_renumbering (level, new_number);
};




#if deal_II_dimension == 1

template <>
void MGDoFHandler<1>::do_renumbering (const unsigned int level,
				      const vector<int> &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level), ExcRenumberingIncomplete());
  
				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // ways it's not only more correct but also
				   // faster
  for (vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	i->set_index (level, d, selected_fe->dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						selected_fe->dofs_per_vertex)]);

  for (vector<int>::iterator i=mg_levels[level]->line_dofs.begin();
       i!=mg_levels[level]->line_dofs.end(); ++i) 
    {
      Assert (*i != -1, ExcInternalError());
      *i = new_numbers[*i];
    };
};

#endif




#if deal_II_dimension == 2

template <>
void MGDoFHandler<2>::do_renumbering (const unsigned int  level,
				      const vector<int>  &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level), ExcRenumberingIncomplete());
  
  for (vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	i->set_index (level, d, selected_fe->dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						selected_fe->dofs_per_vertex)]);
  
  for (vector<int>::iterator i=mg_levels[level]->line_dofs.begin();
       i!=mg_levels[level]->line_dofs.end(); ++i)
    {
      Assert (*i != -1, ExcInternalError());
      *i = new_numbers[*i];
    };
  for (vector<int>::iterator i=mg_levels[level]->quad_dofs.begin();
       i!=mg_levels[level]->quad_dofs.end(); ++i)
    {
      Assert (*i != -1, ExcInternalError());
      *i = new_numbers[*i];
    };
};

#endif




#if deal_II_dimension == 1

template <>
void MGDoFHandler<1>::reserve_space () {
  const unsigned int dim = 1;
  
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());

				   //////////////////////////
				   // DESTRUCTION
  
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<mg_levels.size(); ++i)
    delete mg_levels[i];
  mg_levels.resize (0);

				   // also delete vector of vertex indices
				   // this calls the destructor which
				   // must free the space
  mg_vertex_dofs.resize (0);

				   ////////////////////////////
				   // CONSTRUCTION
  
				   // first allocate space for the
				   // lines on each level
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      mg_levels.push_back (new DoFLevel<1>);

      mg_levels.back()->line_dofs = vector<int>(tria->levels[i]->lines.lines.size() *
						selected_fe->dofs_per_line,
						-1);
    };

				   // now allocate space for the
				   // vertices. To this end, we need
				   // to construct as many objects as
				   // there are vertices and let them
				   // allocate enough space for their
				   // vertex indices on the levels they
				   // live on. We need therefore to
				   // count to how many levels a cell
				   // belongs to, which we do by looping
				   // over all cells and storing the
				   // maximum and minimum level each
				   // vertex we pass by  belongs to
  mg_vertex_dofs.resize (tria->vertices.size());

  vector<unsigned int> min_level (tria->vertices.size(), tria->n_levels());
  vector<unsigned int> max_level (tria->vertices.size(), 0);

  Triangulation<dim>::cell_iterator cell = tria->begin(),
				    endc = tria->end();
  for (; cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	 ++vertex)
      {
	const unsigned int vertex_index = cell->vertex_index(vertex);
	if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
	  min_level[vertex_index] = cell->level();
	if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
	  max_level[vertex_index] = cell->level();
      };

				   // now allocate the needed space
  for (unsigned int vertex=0; vertex<tria->vertices.size(); ++vertex)
    {
      Assert (min_level[vertex] < tria->n_levels(),   ExcInternalError());
      Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());

      mg_vertex_dofs[vertex].init (min_level[vertex],
				   max_level[vertex],
				   selected_fe->dofs_per_vertex);
    };
};

#endif



#if deal_II_dimension == 2

template <>
void MGDoFHandler<2>::reserve_space () {
  const unsigned int dim = 2;
  
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
				   ////////////////////////////
				   // DESTRUCTION

                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<mg_levels.size(); ++i)
    delete mg_levels[i];
  mg_levels.resize (0);

				   // also delete vector of vertex indices
				   // this calls the destructor which
				   // must free the space
  mg_vertex_dofs.resize (0);
  
  
				   ////////////////////////////
				   // CONSTRUCTION
  
				   // first allocate space for the
				   // lines and quads on each level
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      mg_levels.push_back (new DoFLevel<2>);

      mg_levels.back()->line_dofs = vector<int> (tria->levels[i]->lines.lines.size() *
						 selected_fe->dofs_per_line,
						 -1);
      mg_levels.back()->quad_dofs = vector<int> (tria->levels[i]->quads.quads.size() *
						 selected_fe->dofs_per_quad,
						 -1);
    };

  
				   // now allocate space for the
				   // vertices. To this end, we need
				   // to construct as many objects as
				   // there are vertices and let them
				   // allocate enough space for their
				   // vertex indices on the levels they
				   // live on. We need therefore to
				   // count to how many levels a cell
				   // belongs to, which we do by looping
				   // over all cells and storing the
				   // maximum and minimum level each
				   // vertex we pass by  belongs to
  mg_vertex_dofs.resize (tria->vertices.size());

				   // here again, gcc2.8 fails to
				   // construct the vector properly
				   // using parameters to the constructor
//  vector<unsigned int> min_level (tria->vertices.size(), tria->n_levels());
//  vector<unsigned int> max_level (tria->vertices.size(), 0);
  vector<unsigned int> min_level, max_level;
  min_level.resize (tria->vertices.size(), tria->n_levels());
  max_level.resize (tria->vertices.size(), 0);

  Triangulation<dim>::cell_iterator cell = tria->begin(),
				    endc = tria->end();
  for (; cell!=endc; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	 ++vertex)
      {
	const unsigned int vertex_index = cell->vertex_index(vertex);
	if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
	  min_level[vertex_index] = cell->level();
	if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
	  max_level[vertex_index] = cell->level();
      };

  
				   // now allocate the needed space
  for (unsigned int vertex=0; vertex<tria->vertices.size(); ++vertex)
    {
      Assert (min_level[vertex] < tria->n_levels(),   ExcInternalError());
      Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());

      mg_vertex_dofs[vertex].init (min_level[vertex],
				   max_level[vertex],
				   selected_fe->dofs_per_vertex);
    };
};

#endif


// explicite instantiations
template class MGDoFHandler<deal_II_dimension>;
