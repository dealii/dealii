/* $Id$ */

#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/dof_constraints.h>
#include <grid/tria_levels.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>
#include <lac/dsmatrix.h>
#include <map>
#include <set>
#include <algorithm>




template <int dim>
DoFHandler<dim>::DoFHandler (Triangulation<dim> *tria) :
		tria(tria),
		selected_fe (0),
		used_dofs (0) {};


template <int dim>
DoFHandler<dim>::~DoFHandler () {
  if (selected_fe != 0)
    delete selected_fe;
};






#if deal_II_dimension == 1

template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::begin_raw (const unsigned int level) const {
  return begin_raw_line (level);
};



template <>
DoFHandler<1>::cell_iterator
DoFHandler<1>::begin (const unsigned int level) const {
  return begin_line (level);
};



template <>
DoFHandler<1>::active_cell_iterator
DoFHandler<1>::begin_active (const unsigned int level) const {
  return begin_active_line (level);
};



template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::end () const {
  return end_line ();
};



template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::last_raw () const {
  return last_raw_line ();
};



template <>
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::last_raw (const unsigned int level) const {
  return last_raw_line (level);
};



template <>
DoFHandler<1>::cell_iterator
DoFHandler<1>::last () const {
  return last_line ();
};



template <>
DoFHandler<1>::cell_iterator
DoFHandler<1>::last (const unsigned int level) const {
  return last_line (level);
};



template <>
DoFHandler<1>::active_cell_iterator
DoFHandler<1>::last_active () const {
  return last_active_line ();
};



template <>
DoFHandler<1>::active_cell_iterator
DoFHandler<1>::last_active (const unsigned int level) const {
  return last_active_line (level);
};



template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::begin_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::face_iterator
DoFHandler<1>::begin_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::active_face_iterator
DoFHandler<1>::begin_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::end_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::last_raw_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::raw_face_iterator
DoFHandler<1>::last_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::face_iterator
DoFHandler<1>::last_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::face_iterator
DoFHandler<1>::last_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::active_face_iterator
DoFHandler<1>::last_active_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
DoFDimensionInfo<1>::active_face_iterator
DoFHandler<1>::last_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};




template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::begin_raw_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::quad_iterator
DoFHandler<1>::begin_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::active_quad_iterator
DoFHandler<1>::begin_active_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::end_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::last_raw_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};


template <>
DoFHandler<1>::quad_iterator
DoFHandler<1>::last_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::active_quad_iterator
DoFHandler<1>::last_active_quad (const unsigned int) const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::raw_quad_iterator
DoFHandler<1>::last_raw_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::quad_iterator
DoFHandler<1>::last_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
DoFHandler<1>::active_quad_iterator
DoFHandler<1>::last_active_quad () const {
  Assert (false, ExcNotImplemented());
  return 0;
};

#endif


#if deal_II_dimension == 2

template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::begin_raw (const unsigned int level) const {
  return begin_raw_quad (level);
};



template <>
DoFHandler<2>::cell_iterator
DoFHandler<2>::begin (const unsigned int level) const {
  return begin_quad (level);
};



template <>
DoFHandler<2>::active_cell_iterator
DoFHandler<2>::begin_active (const unsigned int level) const {
  return begin_active_quad (level);
};



template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::end () const {
  return end_quad ();
};



template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::last_raw () const {
  return last_raw_quad ();
};



template <>
DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::last_raw (const unsigned int level) const {
  return last_raw_quad (level);
};



template <>
DoFHandler<2>::cell_iterator
DoFHandler<2>::last () const {
  return last_quad ();
};



template <>
DoFHandler<2>::cell_iterator
DoFHandler<2>::last (const unsigned int level) const {
  return last_quad (level);
};



template <>
DoFHandler<2>::active_cell_iterator
DoFHandler<2>::last_active () const {
  return last_active_quad ();
};



template <>
DoFHandler<2>::active_cell_iterator
DoFHandler<2>::last_active (const unsigned int level) const {
  return last_active_quad (level);
};


template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::begin_raw_face (const unsigned int level) const {
  return begin_raw_line (level);
};



template <>
DoFDimensionInfo<2>::face_iterator
DoFHandler<2>::begin_face (const unsigned int level) const {
  return begin_line (level);
};



template <>
DoFDimensionInfo<2>::active_face_iterator
DoFHandler<2>::begin_active_face (const unsigned int level) const {
  return begin_active_line (level);
};



template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::end_face () const {
  return end_line ();
};



template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::last_raw_face () const {
  return last_raw_line ();
};



template <>
DoFDimensionInfo<2>::raw_face_iterator
DoFHandler<2>::last_raw_face (const unsigned int level) const {
  return last_raw_line (level);
};



template <>
DoFDimensionInfo<2>::face_iterator
DoFHandler<2>::last_face () const {
  return last_line ();
};



template <>
DoFDimensionInfo<2>::face_iterator
DoFHandler<2>::last_face (const unsigned int level) const {
  return last_line (level);
};



template <>
DoFDimensionInfo<2>::active_face_iterator
DoFHandler<2>::last_active_face () const {
  return last_active_line ();
};



template <>
DoFDimensionInfo<2>::active_face_iterator
DoFHandler<2>::last_active_face (const unsigned int level) const {
  return last_active_line (level);
};

#endif




template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::begin_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->begin_raw_line(level)->level(),
			    tria->begin_raw_line(level)->index(),
			    this);
};



template <int dim>
typename DoFHandler<dim>::line_iterator
DoFHandler<dim>::begin_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->begin_line(level)->level(),
			tria->begin_line(level)->index(),
			this);
};



template <int dim>
typename DoFHandler<dim>::active_line_iterator
DoFHandler<dim>::begin_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->begin_active_line(level)->level(),
			       tria->begin_active_line(level)->index(),
			       this);
};



template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::begin_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->begin_raw_quad(level)->level(),
			    tria->begin_raw_quad(level)->index(),
			    this);
};



template <int dim>
typename DoFHandler<dim>::quad_iterator
DoFHandler<dim>::begin_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->begin_quad(level)->level(),
			tria->begin_quad(level)->index(),
			this);
};



template <int dim>
typename DoFHandler<dim>::active_quad_iterator
DoFHandler<dim>::begin_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->begin_active_quad(level)->level(),
			       tria->begin_active_quad(level)->index(),
			       this);
};



template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::end_line () const {
  return raw_line_iterator (tria, -1, -1, this);
};



template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::end_quad () const {
  return raw_quad_iterator (tria, -1, -1, this);
};



template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::last_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->last_raw_line(level)->level(),
			    tria->last_raw_line(level)->index(),
			    this);
};



template <int dim>
typename DoFHandler<dim>::line_iterator
DoFHandler<dim>::last_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->last_line(level)->level(),
			tria->last_line(level)->index(),
			this);
};



template <int dim>
typename DoFHandler<dim>::active_line_iterator
DoFHandler<dim>::last_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->last_active_line(level)->level(),
			       tria->last_active_line(level)->index(),
			       this);
};



template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::last_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->last_raw_quad(level)->level(),
			    tria->last_raw_quad(level)->index(),
			    this);
};




template <int dim>
typename DoFHandler<dim>::quad_iterator
DoFHandler<dim>::last_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->last_quad(level)->level(),
			tria->last_quad(level)->index(),
			this);
};




template <int dim>
typename DoFHandler<dim>::active_quad_iterator
DoFHandler<dim>::last_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->last_active_quad(level)->level(),
			       tria->last_active_quad(level)->index(),
			       this);
};




template <int dim>
typename DoFHandler<dim>::raw_line_iterator
DoFHandler<dim>::last_raw_line () const {
  return last_raw_line (levels.size()-1);
};



template <int dim>
typename DoFHandler<dim>::raw_quad_iterator
DoFHandler<dim>::last_raw_quad () const {
  return last_raw_quad (levels.size()-1);
};



template <int dim>
typename DoFHandler<dim>::line_iterator
DoFHandler<dim>::last_line () const {
  return last_line (levels.size()-1);
};



template <int dim>
typename DoFHandler<dim>::quad_iterator
DoFHandler<dim>::last_quad () const {
  return last_quad (levels.size()-1);
};



template <int dim>
typename DoFHandler<dim>::active_line_iterator
DoFHandler<dim>::last_active_line () const {
  return last_active_line (levels.size()-1);
};




template <int dim>
typename DoFHandler<dim>::active_quad_iterator
DoFHandler<dim>::last_active_quad () const {
  return last_active_quad (levels.size()-1);
};






//------------------------------------------------------------------







template <int dim>
unsigned int DoFHandler<dim>::n_dofs () const {
  return used_dofs;
};



#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::n_boundary_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (false, ExcNotImplemented());
  return 0;
};



template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (false, ExcNotImplemented());
  return 0;
};

#endif





template <int dim>
unsigned int DoFHandler<dim>::n_boundary_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  
  set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  vector<int> dofs_on_face(dofs_per_face);
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      };
  return boundary_dofs.size();
};    



template <int dim>
unsigned int DoFHandler<dim>::n_boundary_dofs (const FunctionMap &boundary_indicators) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());
  
  set<int> boundary_dofs;

  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  vector<int> dofs_on_face(dofs_per_face);
  active_face_iterator face = begin_active_face (),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find(face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  boundary_dofs.insert(dofs_on_face[i]);
      };
  return boundary_dofs.size();
};    



template <int dim>
const Triangulation<dim> & DoFHandler<dim>::get_tria () const {
  return *tria;
};



template <int dim>
const FiniteElementBase<dim> & DoFHandler<dim>::get_selected_fe () const {
  return *selected_fe;
};



template <int dim>
void DoFHandler<dim>::distribute_dofs (const FiniteElementBase<dim> &fe) {
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
  if (selected_fe != 0) delete selected_fe;
  selected_fe = new FiniteElementBase<dim>(fe);
  reserve_space ();

				   // clear user flags because we will
				   // need them
  tria->clear_user_flags ();
  
  unsigned int next_free_dof = 0;   
  active_cell_iterator cell = begin_active(),
		       endc = end();

  for (; cell != endc; ++cell) 
    next_free_dof = distribute_dofs_on_cell (cell, next_free_dof);
  
  used_dofs = next_free_dof;
};



#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::distribute_dofs_on_cell (active_cell_iterator &cell,
						     unsigned int          next_free_dof) {

				   // distribute dofs of vertices
  for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
    {
      cell_iterator neighbor = cell->neighbor(v);

      if (neighbor.state() == valid)
	{
					 // find true neighbor; may be its
					 // a child of #neighbor#
	  while (neighbor->has_children())
	    neighbor = neighbor->child(v==0 ? 1 : 0);

					   // has neighbor already been processed?
	  if (neighbor->user_flag_set())
					   // copy dofs
	    {
	      if (v==0) 
		for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (0, d,
					      neighbor->vertex_dof_index (1, d));
	      else
		for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (1, d,
					      neighbor->vertex_dof_index (0, d));

					       // next neighbor
	      continue;
	    };
	};
            
				       // otherwise: create dofs newly
      for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	cell->set_vertex_dof_index (v, d, next_free_dof++);
    };
  
				   // dofs of line
  for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
    cell->set_dof_index (d, next_free_dof++);

				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif


#if deal_II_dimension == 2

template <>
unsigned int DoFHandler<2>::distribute_dofs_on_cell (active_cell_iterator &cell,
						     unsigned int          next_free_dof) {
  if (selected_fe->dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->vertex_dof_index(vertex, 0) == -1)
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  cell->set_vertex_dof_index (vertex, d, next_free_dof++);
    
  				   // for the four sides
  if (selected_fe->dofs_per_line > 0)
    for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
      {
	line_iterator line = cell->line(side);
	
					 // distribute dofs if necessary:
					 // check whether line dof is already
					 // numbered (check only first dof)
	if (line->dof_index(0) == -1)
					   // if not: distribute dofs
	  for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
	    line->set_dof_index (d, next_free_dof++);	    
      };
  

      				       // dofs of quad
  if (selected_fe->dofs_per_quad > 0)
    for (unsigned int d=0; d<selected_fe->dofs_per_quad; ++d)
      cell->set_dof_index (d, next_free_dof++);

  
				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};

#endif



template <int dim>
void DoFHandler<dim>::renumber_dofs (const RenumberingMethod method,
				     const bool use_constraints,
				     const vector<int> &starting_points) {
				   // make the connection graph
  dSMatrixStruct sparsity (n_dofs(), max_couplings_between_dofs());
  make_sparsity_pattern (sparsity);

  if (use_constraints) 
    {
      ConstraintMatrix constraints;
      make_constraint_matrix (constraints);
      constraints.condense (sparsity);
    };
    
  int n_dofs = sparsity.n_rows();
				   // store the new dof numbers; -1 means
				   // that no new number was chosen yet
				   //
				   // the commented line is what would be the
				   // correct way to do, but gcc2.8 chokes
				   // over that. The other lines are a
				   // workaround
//  vector<int> new_number(sparsity.n_rows(), -1);
  vector<int> new_number;
  new_number.resize (sparsity.n_rows(), -1);
  
				   // store the indices of the dofs renumbered
				   // in the last round. Default to starting
				   // points
  vector<int> last_round_dofs (starting_points);
  
				   // delete disallowed elements
  for (unsigned int i=0; i<last_round_dofs.size(); ++i)
    if ((last_round_dofs[i]<0) || (last_round_dofs[i]>=n_dofs))
      last_round_dofs[i] = -1;
  
  remove_if (last_round_dofs.begin(), last_round_dofs.end(),
	     bind2nd(equal_to<int>(), 0));
  
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
  do_renumbering (new_number);
};


#if deal_II_dimension == 1

template <>
void DoFHandler<1>::do_renumbering (const vector<int> &new_numbers) {
				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // ways it's not only more correct but also
				   // faster
  for (vector<int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
    if (*i != -1)
      *i = new_numbers[*i];

  for (unsigned int level=0; level<levels.size(); ++level) 
    for (vector<int>::iterator i=levels[level]->line_dofs.begin();
	 i!=levels[level]->line_dofs.end(); ++i)
      if (*i != -1)
	*i = new_numbers[*i];
};

#endif


#if deal_II_dimension == 2

template <>
void DoFHandler<2>::do_renumbering (const vector<int> &new_numbers) {
  for (vector<int>::iterator i=vertex_dofs.begin(); i!=vertex_dofs.end(); ++i)
    if (*i != -1)
      *i = new_numbers[*i];

  for (unsigned int level=0; level<levels.size(); ++level) 
    {
      for (vector<int>::iterator i=levels[level]->line_dofs.begin();
	   i!=levels[level]->line_dofs.end(); ++i)
	if (*i != -1)
	  *i = new_numbers[*i];
      for (vector<int>::iterator i=levels[level]->quad_dofs.begin();
	   i!=levels[level]->quad_dofs.end(); ++i)
	if (*i != -1)
	  *i = new_numbers[*i];
    };
};

#endif


#if deal_II_dimension == 1

template <>
void DoFHandler<1>::make_constraint_matrix (ConstraintMatrix &cm) const {
  cm.clear ();
  cm.close ();
};

#endif



#if deal_II_dimension == 2

template <>
void DoFHandler<2>::make_constraint_matrix (ConstraintMatrix &constraints) const {
  const unsigned int dim = 2;
  
  constraints.clear ();

				   // first mark all faces which are subject
				   // to constraints. We do so by looping
				   // over all active cells and checking
				   // whether any of the faces are refined
				   // which can only be from the neighboring
				   // cell because this one is active. In that
				   // case, the face is subject to constraints
  tria->clear_user_flags ();
  Triangulation<dim>::active_cell_iterator cell = tria->begin_active(),
					   endc = tria->end();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->has_children()) 
	cell->face(face)->set_user_flag();
	  
  

  
  line_iterator line = begin_line(),
		endl = end_line();
				   // loop over all lines; only on lines
				   // there can be constraints.
  for (; line != endl; ++line)
				     // if dofs on this line are subject
				     // to constraints
    if (line->user_flag_set() == true)
      {
					 // reserve space to gather
					 // the dof numbers. We could
					 // get them when we need them,
					 // but it seems easier to gather
					 // them only once.
	vector<int> dofs_on_mother;
	vector<int> dofs_on_children;
	dofs_on_mother.reserve (2*selected_fe->dofs_per_vertex+
				selected_fe->dofs_per_line);
	dofs_on_children.reserve (selected_fe->dofs_per_vertex+
				  2*selected_fe->dofs_per_line);

	Assert(2*selected_fe->dofs_per_vertex+selected_fe->dofs_per_line ==
	       selected_fe->constraints().n(),
	       ExcDifferentDimensions(2*selected_fe->dofs_per_vertex+
				      selected_fe->dofs_per_line,
				      selected_fe->constraints().n()));
	Assert(selected_fe->dofs_per_vertex+2*selected_fe->dofs_per_line ==
	       selected_fe->constraints().m(),
	       ExcDifferentDimensions(3*selected_fe->dofs_per_vertex+
				      2*selected_fe->dofs_per_line,
				      selected_fe->constraints().m()));
	
					 // fill the dofs indices. Use same
					 // enumeration scheme as in
					 // #FiniteElement::constraints()#
	for (unsigned int vertex=0; vertex<2; ++vertex)
	  for (unsigned int dof=0; dof!=selected_fe->dofs_per_vertex; ++dof)
	    dofs_on_mother.push_back (line->vertex_dof_index(vertex,dof));
	for (unsigned int dof=0; dof!=selected_fe->dofs_per_line; ++dof)
	  dofs_on_mother.push_back (line->dof_index(dof));

	for (unsigned int dof=0; dof!=selected_fe->dofs_per_vertex; ++dof)
	  dofs_on_children.push_back (line->child(0)->vertex_dof_index(1,dof));
	for (unsigned int child=0; child<2; ++child)
	  for (unsigned int dof=0; dof!=selected_fe->dofs_per_line; ++dof)
	    dofs_on_children.push_back (line->child(child)->dof_index(dof));

					 // for each row in the constraint
					 // matrix for this line:
	for (unsigned int row=0; row!=dofs_on_children.size(); ++row) 
	  {
	    constraints.add_line (dofs_on_children[row]);
	    for (unsigned int i=0; i!=dofs_on_mother.size(); ++i)
	      constraints.add_entry (dofs_on_children[row],
				     dofs_on_mother[i],
				     selected_fe->constraints()(row,i));
	  };
      };

  constraints.close ();
};

#endif



template <int dim>
void DoFHandler<dim>::make_sparsity_pattern (dSMatrixStruct &sparsity) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (sparsity.n_rows() == n_dofs(),
	  ExcDifferentDimensions (sparsity.n_rows(), n_dofs()));
  Assert (sparsity.n_cols() == n_dofs(),
	  ExcDifferentDimensions (sparsity.n_cols(), n_dofs()));

  const unsigned int total_dofs = selected_fe->total_dofs;
  vector<int> dofs_on_this_cell(total_dofs);
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell!=endc; ++cell) 
    {
      cell->get_dof_indices (dofs_on_this_cell);
				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<total_dofs; ++i)
	for (unsigned int j=0; j<total_dofs; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    };
};


#if deal_II_dimension == 1

template <>
void DoFHandler<1>::make_boundary_sparsity_pattern (const vector<int> &,
						    dSMatrixStruct &) const {
    Assert (selected_fe != 0, ExcNoFESelected());
    Assert (false, ExcInternalError());
};



template <>
void DoFHandler<1>::make_boundary_sparsity_pattern (const FunctionMap &,
						    const vector<int> &,
						    dSMatrixStruct &) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (false, ExcInternalError());
};

#endif



template <int dim>
void DoFHandler<dim>::make_boundary_sparsity_pattern (const vector<int> &dof_to_boundary_mapping,
						      dSMatrixStruct &sparsity) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (dof_to_boundary_mapping.size() == n_dofs(), ExcInternalError());
  Assert (sparsity.n_rows() == n_boundary_dofs(),
	  ExcDifferentDimensions (sparsity.n_rows(), n_boundary_dofs()));
  Assert (sparsity.n_cols() == n_boundary_dofs(),
	  ExcDifferentDimensions (sparsity.n_cols(), n_boundary_dofs()));
  Assert (*max_element(dof_to_boundary_mapping.begin(),
		       dof_to_boundary_mapping.end()) == (signed int)sparsity.n_rows()-1,
	  ExcInternalError());

  const unsigned int total_dofs = selected_fe->dofs_per_face;
  vector<int> dofs_on_this_face(total_dofs);
  active_face_iterator face = begin_active_face(),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary())
      {
	face->get_dof_indices (dofs_on_this_face);

					 // make sure all dof indices have a
					 // boundary index
	Assert (*min_element(dofs_on_this_face.begin(),
			     dofs_on_this_face.end()) >=0,
		ExcInternalError());
	
					 // make sparsity pattern for this cell
	for (unsigned int i=0; i<total_dofs; ++i)
	  for (unsigned int j=0; j<total_dofs; ++j) 
	    sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
			  dof_to_boundary_mapping[dofs_on_this_face[j]]);
      };
};



template <int dim>
void DoFHandler<dim>::make_boundary_sparsity_pattern (const FunctionMap &boundary_indicators,
						      const vector<int> &dof_to_boundary_mapping,
						      dSMatrixStruct &sparsity) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (dof_to_boundary_mapping.size() == n_dofs(), ExcInternalError());
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());
  Assert (sparsity.n_rows() == n_boundary_dofs(boundary_indicators),
	  ExcDifferentDimensions (sparsity.n_rows(), n_boundary_dofs(boundary_indicators)));
  Assert (sparsity.n_cols() == n_boundary_dofs(boundary_indicators),
	  ExcDifferentDimensions (sparsity.n_cols(), n_boundary_dofs(boundary_indicators)));
  Assert (*max_element(dof_to_boundary_mapping.begin(),
		       dof_to_boundary_mapping.end()) == (signed int)sparsity.n_rows()-1,
	  ExcInternalError());

  const unsigned int total_dofs = selected_fe->dofs_per_face;
  vector<int> dofs_on_this_face(total_dofs);
  active_face_iterator face = begin_active_face(),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find(face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_this_face);

					 // make sure all dof indices have a
					 // boundary index
	Assert (*min_element(dofs_on_this_face.begin(),
			     dofs_on_this_face.end()) >=0,
		ExcInternalError());
					 // make sparsity pattern for this cell
	for (unsigned int i=0; i<total_dofs; ++i)
	  for (unsigned int j=0; j<total_dofs; ++j)
	    sparsity.add (dof_to_boundary_mapping[dofs_on_this_face[i]],
			  dof_to_boundary_mapping[dofs_on_this_face[j]]);
      };
};





template <int dim>
void DoFHandler<dim>::make_transfer_matrix (const DoFHandler<dim> &transfer_from,
					    dSMatrixStruct        &transfer_pattern) const {
  Assert (tria->n_cells(0) == transfer_from.tria->n_cells(0),
	  ExcGridsDoNotMatch());
  Assert (*selected_fe == *transfer_from.selected_fe,
	  ExcGridsDoNotMatch());
				   // assert for once at the beginning the
				   // the matrices have the correct sizes
#ifdef DEBUG
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      Assert (selected_fe->prolongate(c).m() == selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->prolongate(c).m()));
      Assert (selected_fe->prolongate(c).n() == selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->prolongate(c).n()));
      Assert (selected_fe->restrict(c).m() == selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->restrict(c).m()));
      Assert (selected_fe->restrict(c).n() == selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->restrict(c).n()));
    };
#endif
  
  cell_iterator old_cell = transfer_from.begin(0),
		new_cell = begin(0);
  unsigned int  n_coarse_cells = tria->n_cells(0);

				   // first run: make up sparsity structure
  for (unsigned int j=0; j<n_coarse_cells; ++j, ++old_cell, ++new_cell) 
    transfer_cell (old_cell, new_cell, transfer_pattern);

  transfer_pattern.compress();
};



template <int dim>
void DoFHandler<dim>::make_transfer_matrix (const DoFHandler<dim> &transfer_from,
					    dSMatrix              &transfer_matrix) const {
  cell_iterator old_cell = transfer_from.begin(0),
		new_cell = begin(0);
  unsigned int  n_coarse_cells = tria->n_cells(0);

				   // second run: make up matrix entries
  for (unsigned int j=0; j<n_coarse_cells; ++j, ++old_cell, ++new_cell) 
    transfer_cell (old_cell, new_cell, transfer_matrix);
};



template <int dim>
void DoFHandler<dim>::transfer_cell (const typename DoFHandler<dim>::cell_iterator &old_cell,
				     const typename DoFHandler<dim>::cell_iterator &new_cell,
				     dSMatrixStruct      &transfer_pattern) const {
  if (!new_cell->active() && !old_cell->active())
				     // both cells have children; go deeper
    for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
      transfer_cell (old_cell->child(child), new_cell->child(child),
		     transfer_pattern);
  else
    if (new_cell->active() && old_cell->active())
				       // both cells have no children
      {
	vector<int> old_dofs(selected_fe->total_dofs);
	vector<int> new_dofs(selected_fe->total_dofs);
	old_cell->get_dof_indices (old_dofs);
	new_cell->get_dof_indices (new_dofs);

					 // copy dofs one-by-one
	for (unsigned int j=0; j<old_dofs.size(); ++j)
	  transfer_pattern.add (new_dofs[j], old_dofs[j]);
      }
    else
      if (!new_cell->active() && old_cell->active())
					 // new cell has children, old one has not
	{
	  cell_iterator child[GeometryInfo<dim>::children_per_cell];
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
	    {
	      child[c] = new_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
					   // numbers of old dofs
	  vector<int> old_dof_indices (selected_fe->total_dofs);
	  old_cell->get_dof_indices (old_dof_indices);

	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices(selected_fe->total_dofs);
	      child[c]->get_dof_indices (child_dof_indices);

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->prolongate(c)(k,j) != 0.) 
		    transfer_pattern.add (child_dof_indices[k],
					  old_dof_indices[j]);
	    };
	} else {
					   // old cell has children, new one has not
	  cell_iterator child[GeometryInfo<dim>::children_per_cell];
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
	    {
	      child[c] = old_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
	      					   // numbers of new dofs
	  vector<int> new_dof_indices(selected_fe->total_dofs);
	  new_cell->get_dof_indices(new_dof_indices);
	  
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices (selected_fe->total_dofs);
	      child[c]->get_dof_indices (child_dof_indices);

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->restrict(c)(k,j) != 0.)
		    transfer_pattern.add (new_dof_indices[k],
					  child_dof_indices[j]);
	    };
	};
};



template <int dim>
void DoFHandler<dim>::transfer_cell (const typename DoFHandler<dim>::cell_iterator &old_cell,
				     const typename DoFHandler<dim>::cell_iterator &new_cell,
				     dSMatrix            &transfer_matrix) const {
  if (!new_cell->active() && !old_cell->active())
				     // both cells have children; go deeper
    for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
      transfer_cell (old_cell->child(child), new_cell->child(child),
		     transfer_matrix);
  else
    if (new_cell->active() && old_cell->active())
				       // both cells have no children
      {
	vector<int> old_dofs (selected_fe->total_dofs);
	vector<int> new_dofs (selected_fe->total_dofs);
	old_cell->get_dof_indices (old_dofs);
	new_cell->get_dof_indices (new_dofs);

					 // copy dofs one-by-one
	for (unsigned int j=0; j<old_dofs.size(); ++j)
					   // use the dSMatrix:: as a workaround
					   // for a bug in egcs
	  transfer_matrix.dSMatrix::set (new_dofs[j], old_dofs[j], 1.0);
      }
    else
      if (!new_cell->active() && old_cell->active())
					 // new cell has children, old one has not
	{
	  cell_iterator child[GeometryInfo<dim>::children_per_cell];
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
	    {
	      child[c] = new_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
					   // numbers of old dofs
	  vector<int> old_dof_indices (selected_fe->total_dofs);
	  old_cell->get_dof_indices (old_dof_indices);

	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices (selected_fe->total_dofs);
	      child[c]->get_dof_indices (child_dof_indices);

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->prolongate(c)(k,j) != 0.) 
						     // use the dSMatrix::
						     // as a workaround
						     // for a bug in egcs
		    transfer_matrix.dSMatrix::set (child_dof_indices[k],
						   old_dof_indices[j],
						   selected_fe->prolongate(c)(k,j));
	    };
	} else {
					   // old cell has children, new one has not
	  cell_iterator child[GeometryInfo<dim>::children_per_cell];
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
	    {
	      child[c] = old_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
	      					   // numbers of new dofs
	  vector<int> new_dof_indices (selected_fe->total_dofs);
	  new_cell->get_dof_indices(new_dof_indices);
	  
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices (selected_fe->total_dofs);
	      child[c]->get_dof_indices (child_dof_indices);

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->restrict(c)(k,j) != 0.)
						     // use the dSMatrix:: as
						     // a workaround
						     // for a bug in egcs

		    transfer_matrix.dSMatrix::set (new_dof_indices[k],
						   child_dof_indices[j],
						   selected_fe->restrict(c)(k,j));
	    };
	};
};



#if deal_II_dimension == 1

template <>
unsigned int DoFHandler<1>::max_couplings_between_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  return 3*selected_fe->dofs_per_vertex + 2*selected_fe->dofs_per_line;
};



template <>
unsigned int DoFHandler<1>::max_couplings_between_boundary_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (false, ExcInternalError());
  return 0;
};



template <>
unsigned int DoFHandler<1>::max_transfer_entries (const unsigned int max_level_diff) const {
  Assert (max_level_diff<2, ExcOnlyOnelevelTransferImplemented());
  switch (max_level_diff)
    {
      case 0:
	    return 1;
      case 1:
	    return (2*selected_fe->dofs_per_line+
		    selected_fe->dofs_per_vertex) * 2 + 1;
    };
  return 0;
};

#endif



#if deal_II_dimension == 2

template <>
unsigned int DoFHandler<2>::max_couplings_between_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());

				   // get these numbers by drawing pictures
				   // and counting...
				   // example:
				   //   |     |     |
				   // --x-----x--x--X--
				   //   |     |  |  |
				   //   |     x--x--x
				   //   |     |  |  |
				   // --x--x--*--x--x--
				   //   |  |  |     |
				   //   x--x--x     |
				   //   |  |  |     |
				   // --X--x--x-----x--
				   //   |     |     |
				   // x = vertices connected with center vertex *;
				   //   = total of 19
				   // (the X vertices are connected with * if
				   // the vertices adjacent to X are hanging
				   // nodes)
				   // count lines -> 28 (don't forget to count
				   // mother and children separately!)
  switch (tria->max_adjacent_cells())
    {
      case 4:
	    return (19*selected_fe->dofs_per_vertex +
		    28*selected_fe->dofs_per_line +
		    8*selected_fe->dofs_per_quad);
      case 5:
	    return (21*selected_fe->dofs_per_vertex +
		    31*selected_fe->dofs_per_line +
		    9*selected_fe->dofs_per_quad);
      case 6:
	    return (28*selected_fe->dofs_per_vertex +
		    42*selected_fe->dofs_per_line +
		    12*selected_fe->dofs_per_quad);
      case 7:
	    return (30*selected_fe->dofs_per_vertex +
		    45*selected_fe->dofs_per_line +
		    13*selected_fe->dofs_per_quad);
      case 8:
	    return (37*selected_fe->dofs_per_vertex +
		    56*selected_fe->dofs_per_line +
		    16*selected_fe->dofs_per_quad);
      default:
	    Assert (false, ExcNotImplemented());
	    return 0;
    };
};



template <>
unsigned int DoFHandler<2>::max_couplings_between_boundary_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  return 3*selected_fe->dofs_per_vertex + 2*selected_fe->dofs_per_line;
};



template <>
unsigned int DoFHandler<2>::max_transfer_entries (const unsigned int max_level_diff) const {
  Assert (max_level_diff<2, ExcOnlyOnelevelTransferImplemented());
  switch (max_level_diff)
    {
      case 0:
	    return 1;
      case 1:
	    return (4*selected_fe->dofs_per_quad +
		    12*selected_fe->dofs_per_line+
		    5*selected_fe->dofs_per_vertex) * 4 +1;
    };
  return 0;
};

#endif





template <int dim>
void DoFHandler<dim>::distribute_cell_to_dof_vector (const dVector &cell_data,
						     dVector       &dof_data) const {
  Assert (cell_data.size()==tria->n_active_cells(),
	  ExcWrongSize (cell_data.size(), tria->n_active_cells()));

				   // assign the right size to the output vector
  dof_data.reinit (n_dofs());
				   // count how often we have added a value
				   // in the sum for each dof
  vector<unsigned char> touch_count (n_dofs(), 0);

  active_cell_iterator cell = begin_active(),
		       endc = end();
  unsigned int present_cell = 0;
  const unsigned int total_dofs = selected_fe->total_dofs;
  vector<int> dof_indices (total_dofs);
  
  for (; cell!=endc; ++cell, ++present_cell) 
    {
      cell->get_dof_indices (dof_indices);
      for (unsigned int i=0; i<total_dofs; ++i)
	{
					   // sum up contribution of the
					   // present_cell to this dof
	  dof_data(dof_indices[i]) += cell_data(present_cell);
					   // note that we added another
					   // summand
	  ++touch_count[dof_indices[i]];
	};
    };

				   // compute the mean value on all the
				   // dofs by dividing with the number
				   // of summands.
  for (unsigned int i=0; i<n_dofs(); ++i)
    {
      Assert (touch_count[i]!=0, ExcInternalError());
      dof_data(i) /=  touch_count[i];
    };
};



#if deal_II_dimension == 1

template <>
void DoFHandler<1>::map_dof_to_boundary_indices (vector<int> &) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (false, ExcNotImplemented());
};



template <>
void DoFHandler<1>::map_dof_to_boundary_indices (const FunctionMap &,
						 vector<int> &) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (false, ExcNotImplemented());
};

#endif



template <int dim>
void DoFHandler<dim>::map_dof_to_boundary_indices (vector<int> &mapping) const {
  Assert (selected_fe != 0, ExcNoFESelected());

  mapping.clear ();
  mapping.insert (mapping.end(), n_dofs(), -1);
  
  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  vector<int> dofs_on_face(dofs_per_face);
  int next_boundary_index = 0;
  
  active_face_iterator face = begin_active_face(),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (face->at_boundary()) 
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  if (mapping[dofs_on_face[i]] == -1)
	    mapping[dofs_on_face[i]] = next_boundary_index++;
      };

  Assert (static_cast<unsigned int>(next_boundary_index) == n_boundary_dofs(),
	  ExcInternalError());
};



template <int dim>
void DoFHandler<dim>::map_dof_to_boundary_indices (const FunctionMap &boundary_indicators,
						   vector<int> &mapping) const {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (boundary_indicators.find(255) == boundary_indicators.end(),
	  ExcInvalidBoundaryIndicator());

  mapping.clear ();
  mapping.insert (mapping.end(), n_dofs(), -1);

				   // return if there is nothing to do
  if (boundary_indicators.size() == 0)
    return;
  
  const unsigned int dofs_per_face = selected_fe->dofs_per_face;
  vector<int> dofs_on_face(dofs_per_face);
  int next_boundary_index = 0;
  
  active_face_iterator face = begin_active_face(),
		       endf = end_face();
  for (; face!=endf; ++face)
    if (boundary_indicators.find(face->boundary_indicator()) !=
	boundary_indicators.end())
      {
	face->get_dof_indices (dofs_on_face);
	for (unsigned int i=0; i<dofs_per_face; ++i)
	  if (mapping[dofs_on_face[i]] == -1)
	    mapping[dofs_on_face[i]] = next_boundary_index++;
      };

  Assert (static_cast<unsigned int>(next_boundary_index)
	  == n_boundary_dofs(boundary_indicators),
	  ExcInternalError());
};



#if deal_II_dimension == 1

template <>
void DoFHandler<1>::reserve_space () {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.resize (0);

  vertex_dofs = vector<int>(tria->vertices.size()*
			    selected_fe->dofs_per_vertex,
			    -1);
    
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<1>);

      levels.back()->line_dofs = vector<int>(tria->levels[i]->lines.lines.size() *
					     selected_fe->dofs_per_line,
					     -1);
    };
};

#endif


#if deal_II_dimension == 2

template <>
void DoFHandler<2>::reserve_space () {
  Assert (selected_fe != 0, ExcNoFESelected());
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.resize (0);

  vertex_dofs = vector<int>(tria->vertices.size()*
			    selected_fe->dofs_per_vertex,
			    -1);
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<2>);

      levels.back()->line_dofs = vector<int> (tria->levels[i]->lines.lines.size() *
					      selected_fe->dofs_per_line,
					      -1);
      levels.back()->quad_dofs = vector<int> (tria->levels[i]->quads.quads.size() *
					      selected_fe->dofs_per_quad,
					      -1);
    };
};

#endif



/*-------------- Explicit Instantiations -------------------------------*/
template class DoFHandler<deal_II_dimension>;



