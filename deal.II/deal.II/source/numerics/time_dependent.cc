// $Id$


#include <numerics/time-dependent.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>



TimeDependent::TimeSteppingData::TimeSteppingData (const unsigned int look_ahead,
						   const unsigned int look_back)
		:
		look_ahead (look_ahead),
		look_back (look_back)
{};




TimeDependent::TimeDependent (const TimeSteppingData &data_primal):
		sweep_no (static_cast<unsigned int>(-1)),
		timestepping_data_primal (data_primal)
{};




TimeDependent::~TimeDependent ()
{
  while (timesteps.size() != 0)
    delete_timestep (0);
};



void
TimeDependent::insert_timestep (TimeStepBase      *new_timestep,
				const unsigned int position) 
{
  Assert (position<=timesteps.size(),
	  ExcInvalidPosition(position, timesteps.size()));

				   // lock this timestep from deletion
  new_timestep->subscribe();

				   // first insert the new time step
				   // into the doubly linked list
				   // of timesteps
  if (position == timesteps.size())
    {
				       // at the end
      new_timestep->set_next_timestep (0);
      if (timesteps.size() > 0)
	{
	  timesteps.back()->set_next_timestep (new_timestep);
	  new_timestep->set_previous_timestep (timesteps.back());
	}
      else
	new_timestep->set_previous_timestep (0);
    }
  else
    if (position == 0)
      {
					 // at the beginning
	new_timestep->set_previous_timestep (0);
	if (timesteps.size() > 0)
	  {
	    timesteps[0]->set_previous_timestep (new_timestep);
	    new_timestep->set_next_timestep (timesteps[0]);
	  }
	else
	  new_timestep->set_next_timestep (0);
      }
    else
      {
					 // inner time step
	timesteps[position-1]->set_next_timestep (new_timestep);
	new_timestep->set_next_timestep (timesteps[position]);
	timesteps[position]->set_previous_timestep (new_timestep);
      };

				   // finally enter it into the
				   // array
  timesteps.insert (&timesteps[position], new_timestep);
};



void
TimeDependent::add_timestep (TimeStepBase *new_timestep)
{
  insert_timestep (new_timestep, timesteps.size());
};



void TimeDependent::delete_timestep (const unsigned int position)
{
  Assert (position<timesteps.size(),
	  ExcInvalidPosition(position, timesteps.size()));

  timesteps[position]->unsubscribe();
  delete timesteps[position];
  timesteps.erase (&timesteps[position]);

				   // reset "next" pointer of previous
				   // time step if possible
				   //
				   // note that if now position==size,
				   // then we deleted the last time step
  if (position != 0)
    timesteps[position-1]->set_next_timestep ((position<timesteps.size()) ?
					      timesteps[position] :
					      0);

				   // same for "previous" pointer of next
				   // time step
  if (position<timesteps.size())
    timesteps[position]->set_previous_timestep ((position!=0) ?
						timesteps[position-1] :
						0);
};



void
TimeDependent::solve_primal_problem () 
{
  const unsigned int n_timesteps = timesteps.size();

				   // initialize the time steps for
				   // a round of primal problems
  for (unsigned int step=0; step<n_timesteps; ++step)
    timesteps[step]->init_for_primal_problem();

				   // wake up the first few time levels
  for (int step=-timestepping_data_primal.look_ahead; step<0; ++step)
    for (int look_ahead=0;
	 look_ahead<=static_cast<int>(timestepping_data_primal.look_ahead); ++look_ahead)
      if (step+look_ahead >= 0)
	timesteps[step+look_ahead]->wake_up(look_ahead);
  
  for (unsigned int step=0; step<n_timesteps; ++step)
    {
				       // first thing: wake up the
				       // timesteps ahead as necessary
      for (unsigned int look_ahead=0;
	   look_ahead<=timestepping_data_primal.look_ahead; ++look_ahead)
	if (step+look_ahead < n_timesteps)
	  timesteps[step+look_ahead]->wake_up(look_ahead);
      
				       // actually do the work
      timesteps[step]->solve_primal_problem ();
      
				       // let the timesteps behind sleep
      for (unsigned int look_back=0;
	   look_back<=timestepping_data_primal.look_back; ++look_back)
	if (step>=look_back)
	  timesteps[step-look_back]->sleep(look_back);
    };

				   // make the last few timesteps sleep
  for (int step=n_timesteps;
       step<static_cast<int>(n_timesteps+timestepping_data_primal.look_back); ++step)
    for (int look_back=0;
	 look_back<=static_cast<int>(timestepping_data_primal.look_back); ++look_back)
      if ((step-look_back>=0) && (step-look_back<static_cast<int>(n_timesteps)))
	timesteps[step-look_back]->sleep(look_back);
};



void TimeDependent::start_sweep (const unsigned int s) 
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



TimeStepBase::TimeStepBase (const double time) :
		previous_timestep(0),
		next_timestep (0),
		sweep_no (static_cast<unsigned int>(-1)),
		timestep_no (static_cast<unsigned int>(-1)),
		time (time)
{};



TimeStepBase::~TimeStepBase () 
{};



void
TimeStepBase::wake_up (const unsigned )
{};



void
TimeStepBase::sleep (const unsigned)
{};



void
TimeStepBase::init_for_sweep () 
{};



void
TimeStepBase::init_for_primal_problem () 
{
  next_action = primal_problem;
};



void
TimeStepBase::init_for_dual_problem () 
{
  next_action = dual_problem;
};




void
TimeStepBase::solve_dual_problem () 
{
  Assert (false, ExcPureVirtualFunctionCalled());
};



double
TimeStepBase::get_backward_timestep () const
{
  Assert (previous_timestep != 0, ExcCantComputeTimestep());
  return time - previous_timestep->time;
};



double
TimeStepBase::get_forward_timestep () const
{
  Assert (next_timestep != 0, ExcCantComputeTimestep());
  return next_timestep->time - time;
};



void
TimeStepBase::set_previous_timestep (const TimeStepBase *previous)
{
  previous_timestep = previous;
};



void
TimeStepBase::set_next_timestep (const TimeStepBase *next)
{
  next_timestep     = next;
};



void
TimeStepBase::set_timestep_no (const unsigned int step_no)
{
  timestep_no = step_no;
};



void
TimeStepBase::set_sweep_no (const unsigned int sweep)
{
  sweep_no = sweep;
};



/* ------------------------------------------------------------------------- */


template <int dim>
TimeStepBase_Tria<dim>::Flags::Flags () {
  Assert (false, ExcInternalError());
};



template <int dim>
TimeStepBase_Tria<dim>::Flags::Flags (const unsigned int max_refinement_level,
				      const bool delete_and_rebuild_tria,
				      const unsigned int wakeup_level_to_build_grid,
				      const unsigned int sleep_level_to_delete_grid):
		max_refinement_level (max_refinement_level),
		delete_and_rebuild_tria (delete_and_rebuild_tria),
		wakeup_level_to_build_grid (wakeup_level_to_build_grid),
		sleep_level_to_delete_grid (sleep_level_to_delete_grid)
{
//   Assert (!delete_and_rebuild_tria || (wakeup_level_to_build_grid>=1),
// 	  ExcInvalidParameter(wakeup_level_to_build_grid));
//   Assert (!delete_and_rebuild_tria || (sleep_level_to_delete_grid>=1),
// 	  ExcInvalidParameter(sleep_level_to_delete_grid));
};




template <int dim>
TimeStepBase_Tria<dim>::TimeStepBase_Tria() :
		TimeStepBase (0),
		tria (0),
		coarse_grid (*reinterpret_cast<Triangulation<dim>*>(0))
{
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
TimeStepBase_Tria<dim>::TimeStepBase_Tria (const double              time,
					   const Triangulation<dim> &coarse_grid,
					   const Flags              &flags) :
		TimeStepBase (time),
		tria(0),
		coarse_grid (coarse_grid),
		flags (flags)
{
  coarse_grid.subscribe();
};



template <int dim>
TimeStepBase_Tria<dim>::~TimeStepBase_Tria () 
{
  if (!flags.delete_and_rebuild_tria)
    {
      tria->unsubscribe ();
      delete tria;
    }
  else
    Assert (tria==0, ExcInternalError());

  coarse_grid.unsubscribe();
};


template <int dim>
void
TimeStepBase_Tria<dim>::wake_up (const unsigned wakeup_level) {
  TimeStepBase::wake_up (wakeup_level);
  
  if (wakeup_level == flags.wakeup_level_to_build_grid)
    if (flags.delete_and_rebuild_tria || !tria)
      restore_grid ();
};



template <int dim>
void
TimeStepBase_Tria<dim>::sleep (const unsigned sleep_level)
{
  if (sleep_level == flags.sleep_level_to_delete_grid)
    {
      Assert (tria!=0, ExcInternalError());
      
      if (flags.delete_and_rebuild_tria)
	{
	  tria->unsubscribe();
	  delete tria;
	  tria = 0;
	};
    };

  TimeStepBase::sleep (sleep_level);
};




template <int dim>
void TimeStepBase_Tria<dim>::save_refine_flags ()
{
  				   // for any of the non-initial grids
				   // store the refinement flags
  refine_flags.push_back (vector<bool>());
  coarsen_flags.push_back (vector<bool>());
  tria->save_refine_flags (refine_flags.back());
  tria->save_coarsen_flags (coarsen_flags.back());
};



template <int dim>
void TimeStepBase_Tria<dim>::restore_grid () {
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


// explicit instantiations
template class TimeStepBase_Tria<deal_II_dimension>;


