// $Id$


#include <numerics/time-dependent.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>




template <int dim>
TimeDependent<dim>::~TimeDependent ()
{
  for (typename vector<TimeStepBase<dim>*>::iterator i=timesteps.begin();
       i!=timesteps.end(); ++i)
    {
      (*i)->unsubscribe();
      delete (*i);
    };
  
  timesteps.erase (timesteps.begin(), timesteps.end());
};



template <int dim>
void
TimeDependent<dim>::insert_timestep (TimeStepBase<dim> *new_timestep,
				     const unsigned int position) 
{
  Assert (position<=timesteps.size(),
	  ExcInvalidPosition(position, timesteps.size()));

				   // lock this timestep from deletion
  timesteps.back()->subscribe();

				   // first insert the new time step
				   // into the doubly linked list
				   // of timesteps
  if (position != 0)
    {
      timesteps[position]->set_next_timestep (new_timestep);
      new_timestep->set_previous_timestep (timesteps[position]);
    }
  else
    new_timestep->set_previous_timestep (0);
  
  if (position+1 < timesteps.size())
    {
      timesteps[position+1]->set_previous_timestep (new_timestep);
      new_timestep->set_next_timestep (timesteps[position+1]);
    }
  else
    new_timestep->set_next_timestep (0);

				   // finally enter it into the
				   // array
  timesteps.insert (&timesteps[position], new_timestep);
};



template <int dim>
void
TimeDependent<dim>::add_timestep (TimeStepBase<dim> *new_timestep)
{
  insert_timestep (new_timestep, timesteps.size());
};



template <int dim>
void
TimeDependent<dim>::solve_primal_problem () 
{
  const unsigned int n_timesteps = timesteps.size();

				   // initialize the time steps for
				   // a round of primal problems
  for (unsigned int step=0; step<n_timesteps; ++step)
    timesteps[step]->init_for_primal_problem();
  
  for (unsigned int step=0; step<n_timesteps; ++step)
    {
				       // first thing: wake up as many
				       // timesteps ahead as necessary
      if (step==0)
	{
	  for (unsigned int i=0; i<timestepping_data_primal.look_ahead; ++i)
	    if (i < n_timesteps)
	      timesteps[i]->wake_up();
	}
      else
	if (step+timestepping_data_primal.look_ahead < n_timesteps)
	  timesteps[step+timestepping_data_primal.look_ahead]->wake_up();

				       // actually do the work
      timesteps[step]->solve_primal_problem ();
      
				       // last thing: make those time levels
				       // sleep that are no more needed
      if (step=n_timesteps-1)
	{
	  for (unsigned int i=0; i<timestepping_data_primal.look_back; ++i)
	    if (step >= i)
	      timesteps[step-i]->sleep();
	}
      else
	if (step >= timestepping_data_primal.look_back)
	  timesteps[step-timestepping_data_primal.look_back]->sleep();
    };
};



template <int dim>
void TimeDependent<dim>::start_sweep (const unsigned int s) 
{
  sweep_no = s;

				   // reset the number each
				   // time step has, since some time
				   // steps might have been added since
				   // the last time we visited them
				   //
				   // also set the sweep we will
				   // process in the sequel
  for (unsigned int step=0; step<timesteps.size(); ++step)
    {
      timesteps[step]->set_timestep_no (step);
      timesteps[step]->set_sweep_no (sweep_no);
    };

  for (unsigned int step=0; step<timesteps.size(); ++step)
    timesteps[step]->init_for_sweep ();
};



/* --------------------------------------------------------------------- */


template <int dim>
TimeStepBase<dim>::Flags::Flags (const unsigned int max_refinement_level,
				 const bool delete_and_rebuild_tria):
		max_refinement_level (max_refinement_level),
		delete_and_rebuild_tria (delete_and_rebuild_tria)
{};



template <int dim>
TimeStepBase<dim>::TimeStepBase (const Triangulation<dim> &coarse_grid,
				 const Flags              &flags) :
		tria(0),
		coarse_grid (coarse_grid),
		previous_timestep(0),
		next_timestep (0),
		flags (flags)
{
  coarse_grid.subscribe();
};



template <int dim>
TimeStepBase<dim>::~TimeStepBase () 
{  
  Assert (tria!=0, ExcInternalError());

  coarse_grid.unsubscribe();
};


template <int dim>
void
TimeStepBase<dim>::wake_up () {
  if (flags.delete_and_rebuild_tria)
    restore_grid ();
};



template <int dim>
void
TimeStepBase<dim>::sleep ()
{
  Assert (tria!=0, ExcInternalError());

  if (flags.delete_and_rebuild_tria)
    {
      tria->unsubscribe();
      delete tria;
      tria = 0;
    };
};



template <int dim>
void
TimeStepBase<dim>::init_for_sweep () 
{};



template <int dim>
void
TimeStepBase<dim>::init_for_primal_problem () 
{
  problem_type = primal_problem;
};



template <int dim>
void
TimeStepBase<dim>::init_for_dual_problem () 
{
  problem_type = dual_problem;
};





template <int dim>
void
TimeStepBase<dim>::solve_dual_problem () 
{
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
void TimeStepBase<dim>::save_refine_flags ()
{
  				   // for any of the non-initial grids
				   // store the refinement flags
  refine_flags.push_back (vector<bool>());
  coarsen_flags.push_back (vector<bool>());
  tria->save_refine_flags (refine_flags.back());
  tria->save_coarsen_flags (coarsen_flags.back());
};



template <int dim>
void TimeStepBase<dim>::restore_grid () {
  Assert (tria == 0, ExcGridNotDeleted());
  Assert (refine_flags.size() == coarsen_flags.size(),
	  ExcInternalError());

				   // create a virgin triangulation and
				   // set it to a copy of the coarse grid
  tria = new Triangulation<dim> ();
  tria->subscribe();
  tria->copy_triangulation (coarse_grid);

				   // for each of the previous refinement
				   // sweeps
  for (unsigned int previous_sweep=0; previous_sweep<refine_flags.size();
       ++previous_sweep) 
    {
				       // get flags
      tria->load_refine_flags  (refine_flags[previous_sweep]);
      tria->load_coarsen_flags (coarsen_flags[previous_sweep]);

				       // limit refinement depth if the user
				       // desired so
      if (flags.max_refinement_level != 0)
	{
	  Triangulation<dim>::active_cell_iterator cell, endc;
	  for (cell = tria->begin_active(),
	       endc = tria->end();
	       cell!=endc; ++cell)
	    if (static_cast<unsigned int>(cell->level()) >=
		flags.max_refinement_level)
	      cell->clear_refine_flag();
	};

      tria->execute_coarsening_and_refinement ();
    };
};



template <int dim>
void
TimeStepBase<dim>::set_previous_timestep (const TimeStepBase *previous)
{
  previous_timestep = previous;
};



template <int dim>
void
TimeStepBase<dim>::set_next_timestep (const TimeStepBase *next)
{
  next_timestep     = next;
};


template <int dim>
void
TimeStepBase<dim>::set_timestep_no (const unsigned int step_no)
{
  timestep_no = step_no;
};



template <int dim>
void
TimeStepBase<dim>::set_sweep_no (const unsigned int sweep)
{
  sweep_no = sweep;
};



// explicit instantiations
template class TimeDependent<deal_II_dimension>;
template class TimeStepBase<deal_II_dimension>;


