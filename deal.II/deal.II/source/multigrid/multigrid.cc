// $Id$
// Copyright Guido Kanschat, Universitaet Heidelberg, 1999

#include <grid/tria.h>
#include <grid/mg_dof.h>
#include <grid/mg_dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <numerics/multigrid.h>
#include <lac/vector.h>


MultiGridBase::MultiGridBase(MGTransferBase& transfer,
			     unsigned maxlevel, unsigned minlevel,
			     unsigned pre_smooth, unsigned post_smooth)
		:
		maxlevel(maxlevel), minlevel(minlevel),
		transfer(&transfer),
		n_pre_smooth(pre_smooth), n_post_smooth(post_smooth)
{}


void
MultiGridBase::level_mgstep(unsigned level)
{
  if (level == minlevel)
    {
      coarse_grid_solution(level, d[level], s[level]);
      return;
    }
  
				   // smoothing of the residual by modifying s
  smooth(level, d[level], s[level], n_pre_smooth);

				   // t = d-As
  level_residual(level, t, s[level], d[level]);
  
				   // make t rhs of lower level
  transfer->restrict(level, t, d[level-1]);
				   // do recursion
  level_mgstep(level-1);
				   // do coarse grid correction
  transfer->prolongate(level, s[level], s[level-1]);

				   // smoothing (modify s again)
  post_smooth(level, d[level], s[level], n_post_smooth);
}


void
MultiGridBase::post_smooth(unsigned level,
			   Vector<float>& dst, const Vector<float>& src,
			   unsigned steps)
{
  smooth(level, dst, src, steps);
}



void
MultiGridBase::coarse_grid_solution(unsigned level,
			  Vector<float>& dst,
			  const Vector<float>& src)
{
  smooth(level, dst, src, 10 * (n_pre_smooth + n_post_smooth));
}




// ab hier Wolfgang



/* ------------------------------ MGTransferBase --------------------------- */

MGTransferBase::~MGTransferBase () 
{};



/* ----------------------------- MGTransferPrebuilt ------------------------ */



MGTransferPrebuilt::~MGTransferPrebuilt () 
{};



template <int dim>
void MGTransferPrebuilt::build_matrices (const MGDoFHandler<dim> &mg_dof) 
{
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().total_dofs;
  
				   // reset the size of the array of
				   // matrices
  prolongation_sparsities.clear ();
  prolongation_matrices.clear ();
  prolongation_sparsities.reserve (n_levels-1);
  prolongation_matrices.reserve (n_levels-1);
  

				   // two fields which will store the
				   // indices of the multigrid dofs
				   // for a cell and one of its children
  vector<int> dof_indices_mother (dofs_per_cell);
  vector<int> dof_indices_child (dofs_per_cell);
  
				   // for each level: first build the sparsity
				   // pattern of the matrices and then build the
				   // matrices themselves. note that we only
				   // need to take care of cells on the coarser
				   // level which have children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {
      prolongation_sparsities.push_back (SparseMatrixStruct());
      
      for (typename MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_mother);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().prolongate(child);
	    
		cell->child(child)->get_dof_indices (dof_indices_child);

						 // now tag the entries in the
						 // matrix which will be used
						 // for this pair of mother/child
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			prolongation_sparsities[level].add
			  (dof_indices_child[i],
			   dof_indices_mother[j]);
		      };
	      };
	  };

      prolongation_matrices.push_back (SparseMatrix<float>());
      prolongation_matrices[level].reinit (prolongation_sparsities[level]);


				       // now actually build the matrices
      for (MGDoFHandler<dim>::cell_iterator cell=mg_dof.begin(level);
	   cell != mg_dof.end(level); ++cell)
	if (cell->has_children())
	  {
	    cell->get_mg_dof_indices (dof_indices_mother);

	    for (unsigned int child=0;
		 child<GeometryInfo<dim>::children_per_cell; ++child)
	      {
						 // set an alias to the
						 // prolongation matrix for
						 // this child
		const FullMatrix<double> &prolongation
		  = mg_dof.get_fe().prolongate(child);
	    
		cell->child(child)->get_dof_indices (dof_indices_child);

						 // now set the entries in the
						 // matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (prolongation(i,j) != 0)
		      {
			prolongation_matrices[level].set
			  (dof_indices_child[i],
			   dof_indices_mother[j],
			   prolongation(i,j));
		      };
	      };
	  };
    };
};



void MGTransferPrebuilt::prolongate (const unsigned int   to_level,
				     Vector<float>       &dst,
				     const Vector<float> &src) const 
{
  Assert ((to_level >= 1) && (to_level<prolongation_matrices.size()-1),
	  ExcIndexRange (to_level, 1, prolongation_matrices.size()-1));

  prolongation_matrices[to_level-1].vmult (dst, src);
};



void MGTransferPrebuilt::restrict (const unsigned int   from_level,
				   Vector<float>       &dst,
				   const Vector<float> &src) const 
{
  Assert ((from_level >= 1) && (from_level<prolongation_matrices.size()-1),
	  ExcIndexRange (from_level, 1, prolongation_matrices.size()-1));

  prolongation_matrices[from_level-1].Tvmult (dst, src);
};



// explicit instatations
template
void MGTransferPrebuilt::build_matrices (const MGDoFHandler<deal_II_dimension> &mg_dof);


