//----------------------------  solution_transfer.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solution_transfer.cc  ---------------------------


#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <lac/vector.h>
#include <numerics/solution_transfer.h>


template<int dim, typename number>
SolutionTransfer<dim, number>::SolutionTransfer(const DoFHandler<dim> &dof):
		dof_handler(&dof),
		n_dofs_old(0),
		prepared_for(none)
{}


template<int dim, typename number>
SolutionTransfer<dim, number>::~SolutionTransfer()
{
  clear ();
}


template<int dim, typename number>
void SolutionTransfer<dim, number>::clear ()
{
  if (indices_on_cell.size())
    indices_on_cell.erase(indices_on_cell.begin(), indices_on_cell.end());  
  if (all_pointerstructs.size())
    all_pointerstructs.erase(all_pointerstructs.begin(), all_pointerstructs.end());
  if (dof_values_on_cell.size())
    dof_values_on_cell.erase(dof_values_on_cell.begin(), dof_values_on_cell.end());

  prepared_for=none;
}


template<int dim, typename number>
void SolutionTransfer<dim, number>::prepare_for_pure_refinement()
{ 
  Assert(prepared_for!=pure_refinement, ExcAlreadyPrepForRef());
  Assert(prepared_for!=coarsening_and_refinement, 
	 ExcAlreadyPrepForCoarseAndRef());

  clear();

  const unsigned int n_active_cells = dof_handler->get_tria().n_active_cells();
  const unsigned int dofs_per_cell  = dof_handler->get_fe().dofs_per_cell;
  n_dofs_old=dof_handler->n_dofs();

  indices_on_cell=vector<vector<unsigned int> > (n_active_cells,
						 vector<unsigned int> (dofs_per_cell));

  DoFHandler<dim>::cell_iterator cell = dof_handler->begin(),
				 endc = dof_handler->end();

  for (unsigned int i=0; cell!=endc; ++cell) 
    {
      if (cell->active())
	{
					   // on each cell store the indices
					   // of the dofs. after refining we
					   // get the values on the children
					   // by taking these indices, getting
					   // the respective values out of
					   // the data vectors and prolonging
					   // them to the children
	  cell->get_dof_indices(indices_on_cell[i]);
	  cell->set_user_pointer(&indices_on_cell[i]);

	  ++i;	  
	}
      else
	cell->clear_user_pointer();
    }
  prepared_for=pure_refinement;
}


template<int dim, typename number>  
void
SolutionTransfer<dim, number>::refine_interpolate(const Vector<number> &in,
						  Vector<number>       &out) const
{
  Assert(prepared_for==pure_refinement, ExcNotPrepared());
  Assert(in.size()==n_dofs_old, ExcWrongVectorSize(in.size(),n_dofs_old));
  Assert(out.size()==dof_handler->n_dofs(),
	 ExcWrongVectorSize(out.size(),dof_handler->n_dofs()));

  unsigned int dofs_per_cell=dof_handler->get_fe().dofs_per_cell;  
  Vector<number> local_values(dofs_per_cell);

  DoFHandler<dim>::cell_iterator cell = dof_handler->begin(),
				 endc = dof_handler->end();

  vector<unsigned int> *indexptr;  

  for (; cell!=endc; ++cell) 
    {
      if (cell->user_pointer())
					 // this cell was refined or not
					 // touched at all, so we can get
					 // the new values by just setting
					 // or interpolating to the children,
					 // which is both done by one
					 // function
	{
	  indexptr=static_cast<vector<unsigned int> *>(cell->user_pointer());
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    local_values(i)=in(indexptr->operator[](i));
	  cell->set_dof_values_by_interpolation(local_values, out);
	}
    }
}


template<int dim, typename number>  
void SolutionTransfer<dim, number>::refine_interpolate (Vector<number> &vec) const
{
  Assert(vec.size()==n_dofs_old, ExcWrongVectorSize(vec.size(),n_dofs_old));

  Vector<number> vec_old(vec);
  vec.reinit(dof_handler->n_dofs());
  
  refine_interpolate(vec_old, vec);
}


template<int dim, typename number>
void
SolutionTransfer<dim, number>::
prepare_for_coarsening_and_refinement(const vector<Vector<number> > &all_in)
{
  Assert(prepared_for!=pure_refinement, ExcAlreadyPrepForRef());
  Assert(!prepared_for!=coarsening_and_refinement, 
	 ExcAlreadyPrepForCoarseAndRef());
  
  const unsigned int in_size=all_in.size();
  Assert(in_size!=0, ExcNoInVectorsGiven());

  clear();

  const unsigned int n_active_cells = dof_handler->get_tria().n_active_cells();
  const unsigned int dofs_per_cell  = dof_handler->get_fe().dofs_per_cell;
  n_dofs_old=dof_handler->n_dofs();

  for (unsigned int i=0; i<in_size; ++i)
    {
      Assert(all_in[i].size()==n_dofs_old,
	     ExcWrongVectorSize(all_in[i].size(),n_dofs_old));
    }

				   // first count the number
				   // of cells that'll be coarsened
				   // and that'll stay or be refined
  unsigned int n_cells_to_coarsen=0;
  unsigned int n_cells_to_stay_or_refine=0;
  DoFHandler<dim>::active_cell_iterator act_cell = dof_handler->begin_active(),
					    endc = dof_handler->end();
  for (; act_cell!=endc; ++act_cell) 
    {
      if (act_cell->coarsen_flag_set())
	++n_cells_to_coarsen;
      else
	++n_cells_to_stay_or_refine;
    }
  Assert((n_cells_to_coarsen+n_cells_to_stay_or_refine)==n_active_cells,
	 ExcInternalError());
  Assert(n_cells_to_coarsen%GeometryInfo<dim>::children_per_cell==0,
	 ExcInternalError());
  Assert(n_cells_to_coarsen%GeometryInfo<dim>::children_per_cell==0,
	 ExcTriaPrepCoarseningNotCalledBefore());
  
  const unsigned int n_coarsen_fathers = n_cells_to_coarsen /
					 GeometryInfo<dim>::children_per_cell;

				   // allocate the needed memory
  indices_on_cell    = vector<vector<unsigned int> > (n_cells_to_stay_or_refine,
						      vector<unsigned int> (dofs_per_cell));
  dof_values_on_cell
    = vector<vector<Vector<number> > > (n_coarsen_fathers,
					vector<Vector<number> > (in_size,
								 Vector<number> (dofs_per_cell)));

  all_pointerstructs = vector<Pointerstruct> (n_cells_to_stay_or_refine +
					      n_coarsen_fathers);


				   // we need counters for
				   // the 'to_stay_or_refine' cells 'n_sr',
				   // the 'coarsen_fathers' cells 'n_cf',
				   // and all the cells where a
				   // @p{Pointerstruct} is needed 'n'
  unsigned int n_sr=0, n_cf=0, n=0;
  DoFHandler<dim>::cell_iterator cell = dof_handler->begin();  
  for (; cell!=endc; ++cell) 
    {
      if (cell->active() && !cell->coarsen_flag_set())
	{
					   // cell will not be coarsened,
					   // so we get away by storing the
					   // dof indices and later
					   // interpolating to the children
	  cell->get_dof_indices(indices_on_cell[n_sr]);
	  all_pointerstructs[n].indices_ptr=&indices_on_cell[n_sr];
	  all_pointerstructs[n].dof_values_ptr=0;
	  cell->set_user_pointer(&all_pointerstructs[n]);
	  ++n_sr;
	  ++n;
	}
      else if (cell->has_children() && cell->child(0)->coarsen_flag_set())
	{
					   // note that if one child has the
					   // coaresen flag, then all should
					   // have if Tria::prepare_* has
					   // worked correctly
	  for (unsigned int i=1; i<GeometryInfo<dim>::children_per_cell; ++i)
	    Assert(cell->child(i)->coarsen_flag_set(),
		   ExcTriaPrepCoarseningNotCalledBefore());
	      
	  for (unsigned int j=0; j<in_size; ++j)
	    {
					       // store the data of each of
					       // the input vectors
	      cell->get_interpolated_dof_values(all_in[j],
						dof_values_on_cell[n_cf][j]);
	    }
	  all_pointerstructs[n].indices_ptr=0;
	  all_pointerstructs[n].dof_values_ptr=&dof_values_on_cell[n_cf];
	  cell->set_user_pointer(&all_pointerstructs[n]);    
	  ++n_cf;
	  ++n;
	}
      else
					 // some cell on the lower levels to
					 // which nothing will happen
	cell->clear_user_pointer();
    }
  Assert(n_sr==n_cells_to_stay_or_refine, ExcInternalError());
  Assert(n_cf==n_coarsen_fathers, ExcInternalError());

  prepared_for=coarsening_and_refinement;
}


template<int dim, typename number>
void
SolutionTransfer<dim, number>::prepare_for_coarsening_and_refinement(const Vector<number> &in)
{
  vector<Vector<number> > all_in=vector<Vector<number> >(1, in);
  prepare_for_coarsening_and_refinement(all_in);
}


template<int dim, typename number>
void SolutionTransfer<dim, number>::
interpolate (const vector<Vector<number> > &all_in,
	     vector<Vector<number> >       &all_out) const
{
  Assert(prepared_for==coarsening_and_refinement, ExcNotPrepared());
  for (unsigned int i=0; i<all_in.size(); ++i)
    Assert (all_in[i].size() == n_dofs_old,
	    ExcWrongVectorSize(all_in[i].size(), n_dofs_old));
			      
  unsigned int out_size=all_out.size();

				   // resize the output vector
  if (all_out.size() != all_in.size())
    all_out = vector<Vector<number> >(all_in.size(),
				      Vector<number>(dof_handler->n_dofs()));
  else
    {
      for (unsigned int i=0; i<all_in.size(); ++i)
	if (all_out[i].size() != dof_handler->n_dofs())
	  all_out[i].reinit (dof_handler->n_dofs());
    };


  const unsigned int dofs_per_cell=dof_handler->get_fe().dofs_per_cell;  
  Vector<number> local_values(dofs_per_cell);

  vector<unsigned int>    *indexptr;
  Pointerstruct           *structptr;
  vector<Vector<number> > *valuesptr;
  vector<unsigned int>     dofs(dofs_per_cell);

  DoFHandler<dim>::cell_iterator cell = dof_handler->begin(),
				 endc = dof_handler->end();
  
  for (; cell!=endc; ++cell) 
    {
      if (cell->user_pointer())
	{
	  structptr=static_cast<Pointerstruct *>(cell->user_pointer());
					   // cell stayed or is refined
	  if (indexptr=structptr->indices_ptr)
	    {
	      Assert (structptr->dof_values_ptr == 0,
		      ExcInternalError());

					       // get the values of each of
					       // the input data vectors on
					       // this cell and prolong it
					       // to its children
	      for (unsigned int j=0; j<out_size; ++j)
		{
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    local_values(i)=all_in[j](indexptr->operator[](i));
		  cell->set_dof_values_by_interpolation(local_values,
							all_out[j]);
		}
	    }
					   // children of cell were deleted
	  else if (valuesptr=structptr->dof_values_ptr)
	    {
	      Assert (!cell->has_children(), ExcInternalError());
	      Assert (structptr->indices_ptr == 0,
		      ExcInternalError());

					       // get the local indices
	      cell->get_dof_indices(dofs);

					       // distribute the stored data
					       // to the new vectors
	      for (unsigned int j=0; j<out_size; ++j)
		{
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    all_out[j](dofs[i])=valuesptr->operator[](j)(i);
		}
	    }
					   // undefined status
	  else
	    Assert(false, ExcInternalError());
	}
    }
}


template<int dim, typename number>
void SolutionTransfer<dim, number>::interpolate(const Vector<number> &in,
						Vector<number>       &out) const
{
  vector<Vector<number> > all_in(1);
  all_in[0] = in;
  vector<Vector<number> > all_out(1);
  all_out[0] = out;
  interpolate(all_in,
	      all_out);
  out=all_out[0];
}


template class SolutionTransfer<deal_II_dimension, float>;
template class SolutionTransfer<deal_II_dimension, double>;

/*----------------------------   solutiontransfer.cc     ----------------------*/
