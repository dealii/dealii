/*----------------------------   solutiontransfer.cc     ---------------------*/
/*      $Id$     */
/*           Ralf Hartmann, University of Heidelberg                          */

#include <grid/tria.h>
#include <grid/dof.h>
#include <grid/tria_accessor.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <lac/vector.h>
#include <numerics/solutiontransfer.h>



template<int dim, typename number>
SolutionTransfer<dim, number>::SolutionTransfer(
  const DoFHandler<dim> &dof):
		dof_handler(&dof),
		n_dofs_old(0),
		prepared_for_pure_refinement(0),
		prepared_for_coarsening_and_refinement(0),
		vecs_ptr(0) {}


template<int dim, typename number>
SolutionTransfer<dim, number>::~SolutionTransfer()
{
  reinit();
}


template<int dim, typename number>
void SolutionTransfer<dim, number>::reinit()
{
  if (indices_on_cell.size())
    indices_on_cell.erase(indices_on_cell.begin(), indices_on_cell.end());  
  if (all_pointerstructs.size())
    all_pointerstructs.erase(all_pointerstructs.begin(), all_pointerstructs.end());
  if (dof_values_on_cell.size())
    dof_values_on_cell.erase(dof_values_on_cell.begin(), dof_values_on_cell.end());
  vecs_ptr=0;
  prepared_for_pure_refinement=false;
  prepared_for_coarsening_and_refinement=false;
}



template<int dim, typename number>
void SolutionTransfer<dim, number>::prepare_for_pure_refinement()
{ 
  Assert(!prepared_for_pure_refinement, ExcAlreadyPrepForRef());
  Assert(!prepared_for_coarsening_and_refinement, 
	 ExcAlreadyPrepForRefAndCoarse());

  reinit();

  unsigned int n_cells=dof_handler->get_tria().n_active_cells();
  unsigned int total_dofs=dof_handler->get_fe().total_dofs;
  n_dofs_old=dof_handler->n_dofs();

  indices_on_cell=vector<vector<int> > (
    n_cells, vector<int> (total_dofs));

  DoFHandler<dim>::cell_iterator cell = dof_handler->begin(),
				 endc = dof_handler->end();

  for (unsigned int i=0; cell!=endc; ++cell) 
    {
      if (cell->active())
	{
	  cell->set_user_pointer(&indices_on_cell[i]);
	  cell->get_dof_indices(indices_on_cell[i]);
	  ++i;	  
	}
      else
	cell->clear_user_pointer();
    }
  prepared_for_pure_refinement=true;
}


template<int dim, typename number>  
void SolutionTransfer<dim, number>::refine_interpolate(
  const Vector<number> &in, Vector<number> &out) const
{
  Assert(prepared_for_pure_refinement, ExcNotPrepared());
  Assert(in.size()==n_dofs_old, ExcWrongVectorSize(in.size(),n_dofs_old));
  Assert(out.size()==dof_handler->n_dofs(),
	 ExcWrongVectorSize(out.size(),dof_handler->n_dofs()));

  unsigned int total_dofs=dof_handler->get_fe().total_dofs;  
  Vector<number> local_values(total_dofs);

  DoFHandler<dim>::cell_iterator cell = dof_handler->begin(),
				 endc = dof_handler->end();

  vector<int> *indexptr;  

  for (; cell!=endc; ++cell) 
    {
      if (cell->user_pointer())
	{
	  indexptr=static_cast<vector<int> *>(cell->user_pointer());
	  for (unsigned int i=0; i<total_dofs; ++i)
	    local_values(i)=in(indexptr->operator[](i));
	  cell->set_dof_values_by_interpolation(local_values, out);
	}
    }
}

template<int dim, typename number>  
void SolutionTransfer<dim, number>::refine_interpolate(
  Vector<number> &vec) const
{
  Assert(vec.size()==n_dofs_old, ExcWrongVectorSize(vec.size(),n_dofs_old));

  Vector<number> vec_old(vec);
  vec.reinit(dof_handler->n_dofs());
  
  refine_interpolate(vec_old, vec);
}



template<int dim, typename number>
void SolutionTransfer<dim, number>::prepare_for_coarsening_and_refinement(
  const vector<Vector<number> > &all_in)
{
  Assert(!prepared_for_pure_refinement, ExcAlreadyPrepForRef());
  Assert(!prepared_for_coarsening_and_refinement, 
	 ExcAlreadyPrepForRefAndCoarse());
  unsigned int in_size=all_in.size();
  Assert(in_size!=0, ExcNoInVectorsGiven());

  reinit();

  unsigned int n_cells=dof_handler->get_tria().n_active_cells();
  unsigned int total_dofs=dof_handler->get_fe().total_dofs;
  n_dofs_old=dof_handler->n_dofs();

  for (unsigned int i=0; i<in_size; ++i)
    {
      Assert(all_in[i].size()==n_dofs_old,
	     ExcWrongVectorSize(all_in[i].size(),n_dofs_old));
    }
  vecs_ptr=&all_in;
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
  Assert((n_cells_to_coarsen+n_cells_to_stay_or_refine)==n_cells,
	 ExcInternalError());
  Assert(n_cells_to_coarsen%GeometryInfo<dim>::children_per_cell==0,
	 ExcInternalError());
  Assert(n_cells_to_coarsen%GeometryInfo<dim>::children_per_cell==0,
	 ExcTriaPrepCoarseningNotCalledBefore());
  
  unsigned int n_coarsen_fathers=n_cells_to_coarsen/
				 GeometryInfo<dim>::children_per_cell;

				   // allocate the needed memory
  indices_on_cell=vector<vector<int> > (
    n_cells_to_stay_or_refine, vector<int> (total_dofs));
  dof_values_on_cell=vector<vector<Vector<number> > > (
       n_coarsen_fathers, vector<Vector<number> > (
	 in_size, Vector<number> (total_dofs)));

  all_pointerstructs=vector<Pointerstruct> (
    n_cells_to_stay_or_refine+n_coarsen_fathers);


				   // we need counters for
				   // the 'to_stay_or_refine' cells 'n_sr',
				   // the 'coarsen_fathers' cells 'n_cf',
				   // and all the cells where a
				   // #Pointerstruct# is needed 'n'
  unsigned int n_sr=0, n_cf=0, n=0;
  DoFHandler<dim>::cell_iterator cell = dof_handler->begin();  
  for (; cell!=endc; ++cell) 
    {
      if (cell->active() && !cell->coarsen_flag_set())
	{
	  cell->get_dof_indices(indices_on_cell[n_sr]);
	  all_pointerstructs[n].indices_ptr=&indices_on_cell[n_sr];
	  all_pointerstructs[n].dof_values_ptr=0;
	  cell->set_user_pointer(&all_pointerstructs[n]);
	  ++n_sr, ++n;
	}
      else if (cell->has_children() && cell->child(0)->coarsen_flag_set())
	{
	  for (unsigned int i=1; i<GeometryInfo<dim>::children_per_cell; ++i)
	    Assert(cell->child(i)->coarsen_flag_set(),
		   ExcTriaPrepCoarseningNotCalledBefore());
	      
	  for (unsigned int j=0; j<in_size; ++j)
	    {
	      cell->get_interpolated_dof_values(
		all_in[j], dof_values_on_cell[n_cf][j]);
	    }
	  all_pointerstructs[n].indices_ptr=0;
	  all_pointerstructs[n].dof_values_ptr=&dof_values_on_cell[n_cf];
	  cell->set_user_pointer(&all_pointerstructs[n]);    
	  ++n_cf, ++n;
	}
      else
	cell->clear_user_pointer();
    }
  Assert(n_sr==n_cells_to_stay_or_refine, ExcInternalError());
  Assert(n_cf==n_coarsen_fathers, ExcInternalError());

  prepared_for_coarsening_and_refinement=true;
}

template<int dim, typename number>
void SolutionTransfer<dim, number>::prepare_for_coarsening_and_refinement(
  const Vector<number> &in)
{
  vector<Vector<number> > all_in=vector<Vector<number> >(1, in);
  prepare_for_coarsening_and_refinement(all_in);
}



template<int dim, typename number>
void SolutionTransfer<dim, number>::interpolate (
  vector<Vector<number> > &all_out) const
{
  Assert(prepared_for_coarsening_and_refinement, ExcNotPrepared());
  Assert(vecs_ptr==&all_out, ExcVectorsDifferFromInVectors());
  unsigned int out_size=all_out.size();
				   // local copy of the vectors
  vector<Vector<number> > all_in(all_out);
				   // resize the out_all vectors
  for (unsigned int j=0; j<out_size; ++j)
    all_out[j].reinit(dof_handler->n_dofs());

  unsigned int total_dofs=dof_handler->get_fe().total_dofs;  
  Vector<number> local_values(total_dofs);

  vector<int> *indexptr;
  Pointerstruct *structptr;
  vector<Vector<number> > *valuesptr;
  vector<int> dofs(total_dofs);

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
	      for (unsigned int j=0; j<out_size; ++j)
		{
		  for (unsigned int i=0; i<total_dofs; ++i)
		    local_values(i)=all_in[j](indexptr->operator[](i));
		  cell->set_dof_values_by_interpolation(
		    local_values, all_out[j]);
		}
	    }
					   // children of cell were deleted
	  else if (valuesptr=structptr->dof_values_ptr)
	    {
	      Assert(!cell->has_children(), ExcInternalError());
	      cell->get_dof_indices(dofs);
	      for (unsigned int j=0; j<out_size; ++j)
		{
		  for (unsigned int i=0; i<total_dofs; ++i)
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
void SolutionTransfer<dim, number>::interpolate(Vector<number> &out) const
{
  vector<Vector<number> > all_in=vector<Vector<number> >(1, out);
  interpolate(all_in);
  out=all_in[0];
}


template class SolutionTransfer<2, double>;  

/*----------------------------   solutiontransfer.cc     ----------------------*/
