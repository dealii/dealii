/* $Id$ */

#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/dof_constraints.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include "../deal/dsmatrix.h"
#include <multimap.h>
#include <algo.h>


extern TriaRawIterator<1,LineAccessor<1> > _dummy112;
extern TriaRawIterator<2,LineAccessor<2> > _dummy113;
extern TriaIterator<1,LineAccessor<1> > _dummy114;
extern TriaIterator<2,LineAccessor<2> > _dummy115;
extern TriaActiveIterator<1,LineAccessor<1> > _dummy116;
extern TriaActiveIterator<2,LineAccessor<2> > _dummy117;

extern TriaRawIterator<2,QuadAccessor<2> > _dummy118;
extern TriaIterator<2,QuadAccessor<2> > _dummy119;
extern TriaActiveIterator<2,QuadAccessor<2> > _dummy120;

extern TriaRawIterator<1,CellAccessor<1> > _dummy121;
extern TriaRawIterator<2,CellAccessor<2> > _dummy122;
extern TriaIterator<1,CellAccessor<1> > _dummy123;
extern TriaIterator<2,CellAccessor<2> > _dummy124;
extern TriaActiveIterator<1,CellAccessor<1> > _dummy125;
extern TriaActiveIterator<2,CellAccessor<2> > _dummy126;

extern DoFCellAccessor<1>                   _dummy127; //do this to calm down gcc2.7
extern TriaRawIterator<1,DoFCellAccessor<1> >     _dummy128; //wait for gcc2.8
extern TriaIterator<1,DoFCellAccessor<1> >        _dummy129;
extern TriaActiveIterator<1,DoFCellAccessor<1> >  _dummy130;

extern DoFCellAccessor<2>                   _dummy131; //do this to calm down gcc2.7
extern TriaRawIterator<2,DoFLineAccessor<2,LineAccessor<2> > >     _dummy132; //wait for gcc2.8
extern TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >        _dummy133;
extern TriaActiveIterator<2,DoFLineAccessor<2,LineAccessor<2> > >  _dummy134;
extern TriaRawIterator<2,DoFCellAccessor<2> >     _dummy135; //wait for gcc2.8
extern TriaIterator<2,DoFCellAccessor<2> >         _dummy136;
extern TriaActiveIterator<2,DoFCellAccessor<2> >  _dummy137;

extern Triangulation<1> _dummy138;
extern Triangulation<2> _dummy139;



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




DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::begin_raw (const unsigned int level) const {
  return begin_raw_line (level);
};



DoFHandler<1>::cell_iterator
DoFHandler<1>::begin (const unsigned int level) const {
  return begin_line (level);
};



DoFHandler<1>::active_cell_iterator
DoFHandler<1>::begin_active (const unsigned int level) const {
  return begin_active_line (level);
};



DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::end () const {
  return end_line ();
};



DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::last_raw () const {
  return last_raw_line ();
};



DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::last_raw (const unsigned int level) const {
  return last_raw_line (level);
};



DoFHandler<1>::cell_iterator
DoFHandler<1>::last () const {
  return last_line ();
};



DoFHandler<1>::cell_iterator
DoFHandler<1>::last (const unsigned int level) const {
  return last_line (level);
};



DoFHandler<1>::active_cell_iterator
DoFHandler<1>::last_active () const {
  return last_active_line ();
};



DoFHandler<1>::active_cell_iterator
DoFHandler<1>::last_active (const unsigned int level) const {
  return last_active_line (level);
};





DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::begin_raw (const unsigned int level) const {
  return begin_raw_quad (level);
};



DoFHandler<2>::cell_iterator
DoFHandler<2>::begin (const unsigned int level) const {
  return begin_quad (level);
};



DoFHandler<2>::active_cell_iterator
DoFHandler<2>::begin_active (const unsigned int level) const {
  return begin_active_quad (level);
};



DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::end () const {
  return end_quad ();
};



DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::last_raw () const {
  return last_raw_quad ();
};



DoFHandler<2>::raw_cell_iterator
DoFHandler<2>::last_raw (const unsigned int level) const {
  return last_raw_quad (level);
};



DoFHandler<2>::cell_iterator
DoFHandler<2>::last () const {
  return last_quad ();
};



DoFHandler<2>::cell_iterator
DoFHandler<2>::last (const unsigned int level) const {
  return last_quad (level);
};



DoFHandler<2>::active_cell_iterator
DoFHandler<2>::last_active () const {
  return last_active_quad ();
};



DoFHandler<2>::active_cell_iterator
DoFHandler<2>::last_active (const unsigned int level) const {
  return last_active_quad (level);
};




//------------------------------------------------------------------
DoFHandler<1>::raw_cell_iterator
DoFHandler<1>::begin_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->begin_raw_line(level)->level(),
			    tria->begin_raw_line(level)->index(),
			    this);
};



DoFHandler<2>::raw_line_iterator
DoFHandler<2>::begin_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->begin_raw_line(level)->level(),
			    tria->begin_raw_line(level)->index(),
			    this);
};



DoFHandler<1>::line_iterator
DoFHandler<1>::begin_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->begin_line(level)->level(),
			tria->begin_line(level)->index(),
			this);
};



DoFHandler<2>::line_iterator
DoFHandler<2>::begin_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->begin_line(level)->level(),
			tria->begin_line(level)->index(),
			this);
};



DoFHandler<1>::active_line_iterator
DoFHandler<1>::begin_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->begin_active_line(level)->level(),
			       tria->begin_active_line(level)->index(),
			       this);
};



DoFHandler<2>::active_line_iterator
DoFHandler<2>::begin_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->begin_active_line(level)->level(),
			       tria->begin_active_line(level)->index(),
			       this);
};



DoFHandler<2>::raw_quad_iterator
DoFHandler<2>::begin_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->begin_raw_quad(level)->level(),
			    tria->begin_raw_quad(level)->index(),
			    this);
};



DoFHandler<2>::quad_iterator
DoFHandler<2>::begin_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->begin_quad(level)->level(),
			tria->begin_quad(level)->index(),
			this);
};



DoFHandler<2>::active_quad_iterator
DoFHandler<2>::begin_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->begin_active_quad(level)->level(),
			       tria->begin_active_quad(level)->index(),
			       this);
};



DoFHandler<1>::raw_line_iterator
DoFHandler<1>::end_line () const {
  return raw_line_iterator (tria, -1, -1, this);
};



DoFHandler<2>::raw_line_iterator
DoFHandler<2>::end_line () const {
  return raw_line_iterator (tria, -1, -1, this);
};



DoFHandler<2>::raw_quad_iterator
DoFHandler<2>::end_quad () const {
  return raw_quad_iterator (tria, -1, -1, this);
};




DoFHandler<1>::raw_line_iterator
DoFHandler<1>::last_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->last_raw_line(level)->level(),
			    tria->last_raw_line(level)->index(),
			    this);
};



DoFHandler<2>::raw_line_iterator
DoFHandler<2>::last_raw_line (const unsigned int level) const {
  return raw_line_iterator (tria,
			    tria->last_raw_line(level)->level(),
			    tria->last_raw_line(level)->index(),
			    this);
};



DoFHandler<1>::line_iterator
DoFHandler<1>::last_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->last_line(level)->level(),
			tria->last_line(level)->index(),
			this);
};



DoFHandler<2>::line_iterator
DoFHandler<2>::last_line (const unsigned int level) const {
  return line_iterator (tria,
			tria->last_line(level)->level(),
			tria->last_line(level)->index(),
			this);
};



DoFHandler<1>::active_line_iterator
DoFHandler<1>::last_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->last_active_line(level)->level(),
			       tria->last_active_line(level)->index(),
			       this);
};



DoFHandler<2>::active_line_iterator
DoFHandler<2>::last_active_line (const unsigned int level) const {
  return active_line_iterator (tria,
			       tria->last_active_line(level)->level(),
			       tria->last_active_line(level)->index(),
			       this);
};



DoFHandler<2>::raw_quad_iterator
DoFHandler<2>::last_raw_quad (const unsigned int level) const {
  return raw_quad_iterator (tria,
			    tria->last_raw_quad(level)->level(),
			    tria->last_raw_quad(level)->index(),
			    this);
};




DoFHandler<2>::quad_iterator
DoFHandler<2>::last_quad (const unsigned int level) const {
  return quad_iterator (tria,
			tria->last_quad(level)->level(),
			tria->last_quad(level)->index(),
			this);
};




DoFHandler<2>::active_quad_iterator
DoFHandler<2>::last_active_quad (const unsigned int level) const {
  return active_quad_iterator (tria,
			       tria->last_active_quad(level)->level(),
			       tria->last_active_quad(level)->index(),
			       this);
};




DoFHandler<1>::raw_line_iterator
DoFHandler<1>::last_raw_line () const {
  return last_raw_line (levels.size()-1);
};



DoFHandler<2>::raw_line_iterator
DoFHandler<2>::last_raw_line () const {
  return last_raw_line (levels.size()-1);
};



DoFHandler<2>::raw_quad_iterator
DoFHandler<2>::last_raw_quad () const {
  return last_raw_quad (levels.size()-1);
};



DoFHandler<1>::line_iterator
DoFHandler<1>::last_line () const {
  return last_line (levels.size()-1);
};



DoFHandler<2>::line_iterator
DoFHandler<2>::last_line () const {
  return last_line (levels.size()-1);
};



DoFHandler<2>::quad_iterator
DoFHandler<2>::last_quad () const {
  return last_quad (levels.size()-1);
};



DoFHandler<1>::active_line_iterator
DoFHandler<1>::last_active_line () const {
  return last_active_line (levels.size()-1);
};



DoFHandler<2>::active_line_iterator
DoFHandler<2>::last_active_line () const {
  return last_active_line (levels.size()-1);
};



DoFHandler<2>::active_quad_iterator
DoFHandler<2>::last_active_quad () const {
  return last_active_quad (levels.size()-1);
};







template <int dim>
unsigned int DoFHandler<dim>::n_dofs () const {
  return used_dofs;
};



template <int dim>
const Triangulation<dim> & DoFHandler<dim>::get_tria () const {
  return *tria;
};



template <int dim>
const FiniteElement<dim> & DoFHandler<dim>::get_selected_fe () const {
  return *selected_fe;
};



template <int dim>
void DoFHandler<dim>::distribute_dofs (const FiniteElement<dim> &fe) {
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
  reserve_space (fe);
  if (selected_fe != 0) delete selected_fe;
  selected_fe = new FiniteElement<dim>(fe);

				   // clear user flags because we will
				   // need them
  tria->clear_user_flags ();
  
  unsigned int next_free_dof = 0;   
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell != endc; ++cell) 
    next_free_dof = distribute_dofs_on_cell (cell, fe, next_free_dof);
  
  used_dofs = next_free_dof;
};



int DoFHandler<1>::distribute_dofs_on_cell (active_cell_iterator   &cell,
					    const FiniteElement<1> &fe,
					    unsigned int            next_free_dof) {

				   // distribute dofs of vertices
  for (unsigned int v=0; v<2; ++v)
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
		for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (0, d, neighbor->vertex_dof_index (1, d));
	      else
		for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
		  cell->set_vertex_dof_index (1, d, neighbor->vertex_dof_index (0, d));

					       // next neighbor
	      continue;
	    };
	};
            
				       // otherwise: create dofs newly
      for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
	cell->set_vertex_dof_index (v, d, next_free_dof++);
    };
  
				   // dofs of line
  for (unsigned int d=0; d<fe.dofs_per_line; ++d)
    cell->set_dof_index (d, next_free_dof++);

				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};




int DoFHandler<2>::distribute_dofs_on_cell (active_cell_iterator   &cell,
					    const FiniteElement<2> &fe,
					    unsigned int            next_free_dof) {
  if (fe.dofs_per_vertex > 0)
				     // number dofs on vertices
    for (unsigned int vertex=0; vertex<4; ++vertex)
				       // check whether dofs for this
				       // vertex have been distributed
				       // (only check the first dof)
      if (cell->vertex_dof_index(vertex, 0) == -1)
	for (unsigned int d=0; d<fe.dofs_per_vertex; ++d)
	  cell->set_vertex_dof_index (vertex, d, next_free_dof++);
    
  				   // for the four sides
  for (unsigned int side=0; side<4; ++side)
    {
      line_iterator line = cell->line(side);

				       // distribute dofs if necessary
      if ((fe.dofs_per_line > 0) &&
					 // check whether line dof is already
					 // numbered (check only first dof)
	  (line->dof_index(0) == -1))
					 // if not: distribute dofs
	for (unsigned int d=0; d<fe.dofs_per_line; ++d)
	  line->set_dof_index (d, next_free_dof++);	    

				       // note if line is subject to constraints
      if (line->child_index(0) != -1) 
	line->set_user_flag ();
    };
  

      				       // dofs of quad
  if (fe.dofs_per_quad > 0)
    for (unsigned int d=0; d<fe.dofs_per_line; ++d)
      cell->set_dof_index (d, next_free_dof++);

  
				   // note that this cell has been processed
  cell->set_user_flag ();
  
  return next_free_dof;
};




template <int dim>
void DoFHandler<dim>::renumber_dofs (const RenumberingMethod method,
				     bool use_constraints,
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
  vector<int> new_number(sparsity.n_rows(), -1);

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
	  int j;
	  for (j=sparsity.rowstart[row]; j<sparsity.rowstart[row+1]; ++j)
	    if (sparsity.colnums[j] == -1)
	      break;
					   // post: coordination is now
					   // j-rowstart[row]
	  if (j-sparsity.rowstart[row] < (signed int)min_coordination)
	    {
	      min_coordination = j-sparsity.rowstart[row];
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
	for (int j=sparsity.rowstart[last_round_dofs[i]];
	     j<sparsity.rowstart[last_round_dofs[i]+1]; ++j)
	  if (sparsity.colnums[j] == -1)
	    break;
	  else
	    next_round_dofs.push_back (sparsity.colnums[j]);
      
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
      multimap<unsigned int, int, less<unsigned int> > dofs_by_coordination;
      
				       // find coordination number for
				       // each of these dofs
      for (vector<int>::iterator s=next_round_dofs.begin();
	   s!=next_round_dofs.end(); ++s) 
	{
	  unsigned int coordination = 0;
	  for (int j=sparsity.rowstart[*s]; j<sparsity.rowstart[*s+1]; ++j)
	    if (sparsity.colnums[j] == -1)
	      break;
	    else
	      ++coordination;

					   // insert this dof at its
					   // coordination number
	  const pair<const unsigned int, int> new_entry (coordination, *s);
	  dofs_by_coordination.insert (new_entry);
	};
      
				       ////
      multimap<unsigned int, int, less<unsigned int> >::iterator i;
      for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i) 
	new_number[(*i).second] = next_free_number++;

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




void DoFHandler<1>::make_constraint_matrix (ConstraintMatrix &cm) const {
  cm.clear ();
};


  

void DoFHandler<2>::make_constraint_matrix (ConstraintMatrix &constraints) const {
  constraints.clear ();
  
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
	       (unsigned int)selected_fe->constraints().n(),
	       ExcDifferentDimensions(2*selected_fe->dofs_per_vertex+
				      selected_fe->dofs_per_line,
				      selected_fe->constraints().n()));
	Assert(selected_fe->dofs_per_vertex+2*selected_fe->dofs_per_line ==
	       (unsigned int)selected_fe->constraints().m(),
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




void DoFHandler<1>::make_sparsity_pattern (dSMatrixStruct &sparsity) const {
  active_cell_iterator cell = begin_active(),
		       endc = end();

				   // set up an array which dofs are used
				   // on a specific cell
  unsigned int      *dofs_on_this_cell = new unsigned int[selected_fe->total_dofs];
  
  for (; cell!=endc; ++cell) 
    {
      unsigned int dof_number=0;

				       // fill dof indices on vertices
      for (unsigned int vertex=0; vertex<2; ++vertex)
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  dofs_on_this_cell[dof_number++] = cell->vertex_dof_index (vertex,d);

				       // fill dof indices on line
      for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
	dofs_on_this_cell[dof_number++] = cell->dof_index (d);

				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<selected_fe->total_dofs; ++i)
	for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    };

  delete[] dofs_on_this_cell;
};



void DoFHandler<2>::make_sparsity_pattern (dSMatrixStruct &sparsity) const {
  active_cell_iterator cell = begin_active(),
		       endc = end();

				   // set up an array which dofs are used
				   // on a specific cell
  unsigned int      *dofs_on_this_cell = new unsigned int[selected_fe->total_dofs];

  
  for (; cell!=endc; ++cell) 
    {
      int dof_number=0;

				       // fill dof indices on vertices
      for (unsigned int vertex=0; vertex<4; ++vertex)
	for (unsigned int d=0; d<selected_fe->dofs_per_vertex; ++d)
	  dofs_on_this_cell[dof_number++] = cell->vertex_dof_index (vertex,d);

      for (unsigned int line=0; line<4; ++line)
	for (unsigned int d=0; d<selected_fe->dofs_per_line; ++d)
	  dofs_on_this_cell[dof_number++] = cell->line(line)->dof_index (d);
      
				       // fill dof indices on quad
      for (unsigned int d=0; d<selected_fe->dofs_per_quad; ++d)
	dofs_on_this_cell[dof_number++] = cell->dof_index (d);

				       // make sparsity pattern for this cell
      for (unsigned int i=0; i<selected_fe->total_dofs; ++i)
	for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
	  sparsity.add (dofs_on_this_cell[i],
			dofs_on_this_cell[j]);
    };

  delete[] dofs_on_this_cell;
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
  for (unsigned int c=0; c<(1<<dim); ++c)
    {
      Assert ((unsigned int)selected_fe->prolongate(c).m() ==
	      selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->prolongate(c).m()));
      Assert ((unsigned int)selected_fe->prolongate(c).n() ==
	      selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->prolongate(c).n()));
      Assert ((unsigned int)selected_fe->restrict(c).m() ==
	      selected_fe->total_dofs,
	      ExcMatrixHasWrongSize(selected_fe->restrict(c).m()));
      Assert ((unsigned int)selected_fe->restrict(c).n() ==
	      selected_fe->total_dofs,
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
void DoFHandler<dim>::transfer_cell (const TriaIterator<dim,DoFCellAccessor<dim> > &old_cell,
				     const TriaIterator<dim,DoFCellAccessor<dim> > &new_cell,
				     dSMatrixStruct      &transfer_pattern) const {
  if (!new_cell->active() && !old_cell->active())
				     // both cells have children; go deeper
    for (unsigned int child=0; child<(1<<dim); ++child)
      transfer_cell (old_cell->child(child), new_cell->child(child),
		     transfer_pattern);
  else
    if (new_cell->active() && old_cell->active())
				       // both cells have no children
      {
	vector<int> old_dofs, new_dofs;
	old_cell->dof_indices (old_dofs);
	new_cell->dof_indices (new_dofs);
	Assert (old_dofs.size() == selected_fe->total_dofs,
		ExcInternalError ());
	Assert (new_dofs.size() == selected_fe->total_dofs,
		ExcInternalError ());

					 // copy dofs one-by-one
	for (unsigned int j=0; j<old_dofs.size(); ++j)
	  transfer_pattern.add (new_dofs[j], old_dofs[j]);
      }
    else
      if (!new_cell->active() && old_cell->active())
					 // new cell has children, old one has not
	{
	  cell_iterator child[1<<dim];
	  for (unsigned int c=0; c<(1<<dim); ++c) 
	    {
	      child[c] = new_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
					   // numbers of old dofs
	  vector<int> old_dof_indices;
	  old_cell->dof_indices (old_dof_indices);

	  Assert (old_dof_indices.size() == selected_fe->total_dofs,
		  ExcInternalError ());

	  for (unsigned int c=0; c<(1<<dim); ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices;
	      child[c]->dof_indices (child_dof_indices);

	      Assert (child_dof_indices.size() == selected_fe->total_dofs,
		      ExcInternalError ());

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->prolongate(c)(k,j) != 0.) 
		    transfer_pattern.add (child_dof_indices[k],
					  old_dof_indices[j]);
	    };
	} else {
					   // old cell has children, new one has not
	  cell_iterator child[1<<dim];
	  for (unsigned int c=0; c<(1<<dim); ++c) 
	    {
	      child[c] = old_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
	      					   // numbers of new dofs
	  vector<int> new_dof_indices;
	  new_cell->dof_indices(new_dof_indices);

	  Assert (new_dof_indices.size() == selected_fe->total_dofs,
		  ExcInternalError ());
	  
	  for (unsigned int c=0; c<(1<<dim); ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices;
	      child[c]->dof_indices (child_dof_indices);

	      Assert (child_dof_indices.size() == selected_fe->total_dofs,
		      ExcInternalError ());

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->restrict(c)(k,j) != 0.)
		    transfer_pattern.add (new_dof_indices[k],
					  child_dof_indices[j]);
	    };
	};
};



template <int dim>
void DoFHandler<dim>::transfer_cell (const TriaIterator<dim,DoFCellAccessor<dim> > &old_cell,
				     const TriaIterator<dim,DoFCellAccessor<dim> > &new_cell,
				     dSMatrix            &transfer_matrix) const {
  if (!new_cell->active() && !old_cell->active())
				     // both cells have children; go deeper
    for (unsigned int child=0; child<(1<<dim); ++child)
      transfer_cell (old_cell->child(child), new_cell->child(child),
		     transfer_matrix);
  else
    if (new_cell->active() && old_cell->active())
				       // both cells have no children
      {
	vector<int> old_dofs, new_dofs;
	old_cell->dof_indices (old_dofs);
	new_cell->dof_indices (new_dofs);
	Assert (old_dofs.size() == selected_fe->total_dofs,
		ExcInternalError ());
	Assert (new_dofs.size() == selected_fe->total_dofs,
		ExcInternalError ());

					 // copy dofs one-by-one
	for (unsigned int j=0; j<old_dofs.size(); ++j)
	  transfer_matrix.set (new_dofs[j], old_dofs[j], 1.0);
      }
    else
      if (!new_cell->active() && old_cell->active())
					 // new cell has children, old one has not
	{
	  cell_iterator child[1<<dim];
	  for (unsigned int c=0; c<(1<<dim); ++c) 
	    {
	      child[c] = new_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
					   // numbers of old dofs
	  vector<int> old_dof_indices;
	  old_cell->dof_indices (old_dof_indices);

	  Assert (old_dof_indices.size() == selected_fe->total_dofs,
		  ExcInternalError ());

	  for (unsigned int c=0; c<(1<<dim); ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices;
	      child[c]->dof_indices (child_dof_indices);

	      Assert (child_dof_indices.size() == selected_fe->total_dofs,
		      ExcInternalError ());

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->prolongate(c)(k,j) != 0.) 
		    transfer_matrix.set (child_dof_indices[k],
					 old_dof_indices[j],
					 selected_fe->prolongate(c)(k,j));
	    };
	} else {
					   // old cell has children, new one has not
	  cell_iterator child[1<<dim];
	  for (unsigned int c=0; c<(1<<dim); ++c) 
	    {
	      child[c] = old_cell->child(c);
	      Assert (child[c]->active(),
		      ExcOnlyOnelevelTransferImplemented());
	    };
	  
	      					   // numbers of new dofs
	  vector<int> new_dof_indices;
	  new_cell->dof_indices(new_dof_indices);

	  Assert (new_dof_indices.size() == selected_fe->total_dofs,
		  ExcInternalError ());
	  
	  for (unsigned int c=0; c<(1<<dim); ++c)
	    {
					       // numbers of child dofs
	      vector<int> child_dof_indices;
	      child[c]->dof_indices (child_dof_indices);

	      Assert (child_dof_indices.size() == selected_fe->total_dofs,
		      ExcInternalError ());

	      for (unsigned int k=0; k<selected_fe->total_dofs; ++k)
		for (unsigned int j=0; j<selected_fe->total_dofs; ++j)
		  if (selected_fe->restrict(c)(k,j) != 0.)
		    transfer_matrix.set (new_dof_indices[k],
					 child_dof_indices[j],
					 selected_fe->restrict(c)(k,j));
	    };
	};
};




unsigned int DoFHandler<1>::max_couplings_between_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  return 3*selected_fe->dofs_per_vertex + 2*selected_fe->dofs_per_line;
};



unsigned int DoFHandler<2>::max_couplings_between_dofs () const {
  Assert (selected_fe != 0, ExcNoFESelected());
  unsigned int max_adjacent_cells = tria->max_adjacent_cells();

				   // get these numbers by drawing pictures
				   // and counting...
  if (max_adjacent_cells == 4)
    return (13*selected_fe->dofs_per_vertex +
	    20*selected_fe->dofs_per_line +
	    8*selected_fe->dofs_per_quad);
  else
    if (max_adjacent_cells == 5)
      return (15*selected_fe->dofs_per_vertex +
	      23*selected_fe->dofs_per_line +
	      9*selected_fe->dofs_per_quad);
    else
      if (max_adjacent_cells == 6)
					 // are you really sure you
					 // want to use such grids??
	return (19*selected_fe->dofs_per_vertex +
		30*selected_fe->dofs_per_line +
		16*selected_fe->dofs_per_quad);
      else 
	{
	  Assert (false, ExcNotImplemented());
	  return 0;
	};
};



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



void DoFHandler<1>::reserve_space (const FiniteElement<1> &fe) {
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.erase (levels.begin(), levels.end());

  vertex_dofs.erase (vertex_dofs.begin(), vertex_dofs.end());
  vertex_dofs.reserve (tria->vertices.size());
  vertex_dofs.insert (vertex_dofs.end(),
		      tria->vertices.size(),
		      -1);
    
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<1>);

      levels.back()->line_dofs.reserve (tria->levels[i]->lines.lines.size() *
					fe.dofs_per_line);
      levels.back()->line_dofs.insert (levels.back()->line_dofs.end(),
				       (tria->levels[i]->lines.lines.size() *
					fe.dofs_per_line),
				       -1);
    };
};



void DoFHandler<2>::reserve_space (const FiniteElement<2> &fe) {
  Assert (tria->n_levels() > 0, ExcInvalidTriangulation());
  
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.erase (levels.begin(), levels.end());

  vertex_dofs.erase (vertex_dofs.begin(), vertex_dofs.end());
  vertex_dofs.reserve (tria->vertices.size());
  vertex_dofs.insert (vertex_dofs.end(),
		      tria->vertices.size(),
		      -1);
    
  for (unsigned int i=0; i<tria->n_levels(); ++i) 
    {
      levels.push_back (new DoFLevel<2>);

      levels.back()->line_dofs.reserve (tria->levels[i]->lines.lines.size() *
					fe.dofs_per_line);
      levels.back()->line_dofs.insert (levels.back()->line_dofs.end(),
				       (tria->levels[i]->lines.lines.size() *
					fe.dofs_per_line),
				       -1);
      
      levels.back()->quad_dofs.reserve (tria->levels[i]->quads.quads.size() *
					fe.dofs_per_quad);
      levels.back()->quad_dofs.insert (levels.back()->quad_dofs.end(),
				       (tria->levels[i]->quads.quads.size() *
					fe.dofs_per_quad),
				       -1);
    };
};





/*-------------- Explicit Instantiations -------------------------------*/
template class DoFHandler<1>;
template class DoFHandler<2>;



