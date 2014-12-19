// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/numerics/time_dependent.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parallel.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/lac/vector.h>

#include <functional>
#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN

TimeDependent::TimeSteppingData::TimeSteppingData (const unsigned int look_ahead,
                                                   const unsigned int look_back)
  :
  look_ahead (look_ahead),
  look_back (look_back)
{}


TimeDependent::TimeDependent (const TimeSteppingData &data_primal,
                              const TimeSteppingData &data_dual,
                              const TimeSteppingData &data_postprocess):
  sweep_no (numbers::invalid_unsigned_int),
  timestepping_data_primal (data_primal),
  timestepping_data_dual (data_dual),
  timestepping_data_postprocess (data_postprocess)
{}


TimeDependent::~TimeDependent ()
{
  while (timesteps.size() != 0)
    delete_timestep (0);
}


void
TimeDependent::insert_timestep (const TimeStepBase *position,
                                TimeStepBase       *new_timestep)
{
  Assert ((std::find(timesteps.begin(), timesteps.end(), position) != timesteps.end()) ||
          (position == 0),
          ExcInvalidPosition());
  // first insert the new time step
  // into the doubly linked list
  // of timesteps
  if (position == 0)
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
  else if (position == timesteps[0])
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
      std::vector<SmartPointer<TimeStepBase,TimeDependent> >::iterator insert_position
        = std::find(timesteps.begin(), timesteps.end(), position);

      (*(insert_position-1))->set_next_timestep (new_timestep);
      new_timestep->set_previous_timestep (*(insert_position-1));
      new_timestep->set_next_timestep (*insert_position);
      (*insert_position)->set_previous_timestep (new_timestep);
    };

  // finally enter it into the
  // array
  timesteps.insert ((position == 0 ?
                     timesteps.end() :
                     std::find(timesteps.begin(), timesteps.end(), position)),
                    new_timestep);
}


void
TimeDependent::add_timestep (TimeStepBase *new_timestep)
{
  insert_timestep (0, new_timestep);
}


void TimeDependent::delete_timestep (const unsigned int position)
{
  Assert (position<timesteps.size(),
          ExcInvalidPosition());

  // Remember time step object for
  // later deletion and unlock
  // SmartPointer
  TimeStepBase *t = timesteps[position];
  timesteps[position] = 0;
  // Now delete unsubscribed object
  delete t;

  timesteps.erase (timesteps.begin() + position);

  // reset "next" pointer of previous
  // time step if possible
  //
  // note that if now position==size,
  // then we deleted the last time step
  if (position != 0)
    timesteps[position-1]->set_next_timestep ((position<timesteps.size()) ?
                                              timesteps[position] :
                                              /*null*/SmartPointer<TimeStepBase,TimeDependent>());

  // same for "previous" pointer of next
  // time step
  if (position<timesteps.size())
    timesteps[position]->set_previous_timestep ((position!=0) ?
                                                timesteps[position-1] :
                                                /*null*/SmartPointer<TimeStepBase,TimeDependent>());
}


void
TimeDependent::solve_primal_problem ()
{
  do_loop (std::mem_fun(&TimeStepBase::init_for_primal_problem),
           std::mem_fun(&TimeStepBase::solve_primal_problem),
           timestepping_data_primal,
           forward);
}


void
TimeDependent::solve_dual_problem ()
{
  do_loop (std::mem_fun(&TimeStepBase::init_for_dual_problem),
           std::mem_fun(&TimeStepBase::solve_dual_problem),
           timestepping_data_dual,
           backward);
}


void
TimeDependent::postprocess ()
{
  do_loop (std::mem_fun(&TimeStepBase::init_for_postprocessing),
           std::mem_fun(&TimeStepBase::postprocess_timestep),
           timestepping_data_postprocess,
           forward);
}



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
    timesteps[step]->start_sweep ();
}



void TimeDependent::end_sweep (const unsigned int)
{
  end_sweep ();
}


void TimeDependent::end_sweep ()
{
  void (TimeDependent::*p) (const unsigned int, const unsigned int)
    = &TimeDependent::end_sweep;
  parallel::apply_to_subranges (0U, timesteps.size(),
                                std_cxx11::bind (p, this, std_cxx11::_1, std_cxx11::_2),
                                1);
}



void TimeDependent::end_sweep (const unsigned int begin,
                               const unsigned int end)
{
  for (unsigned int step=begin; step<end; ++step)
    timesteps[step]->end_sweep ();
}



std::size_t
TimeDependent::memory_consumption () const
{
  std::size_t mem = (MemoryConsumption::memory_consumption (timesteps) +
                     MemoryConsumption::memory_consumption (sweep_no) +
                     sizeof(timestepping_data_primal) +
                     sizeof(timestepping_data_dual) +
                     sizeof(timestepping_data_postprocess));
  for (unsigned int i=0; i<timesteps.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*timesteps[i]);

  return mem;
}



/* --------------------------------------------------------------------- */


TimeStepBase::TimeStepBase (const double time) :
  previous_timestep(0),
  next_timestep (0),
  sweep_no (numbers::invalid_unsigned_int),
  timestep_no (numbers::invalid_unsigned_int),
  time (time)
{}



TimeStepBase::~TimeStepBase ()
{}



void
TimeStepBase::wake_up (const unsigned int )
{}



void
TimeStepBase::sleep (const unsigned)
{}



void
TimeStepBase::start_sweep ()
{}



void
TimeStepBase::end_sweep ()
{}



void
TimeStepBase::init_for_primal_problem ()
{
  next_action = primal_problem;
}



void
TimeStepBase::init_for_dual_problem ()
{
  next_action = dual_problem;
}



void
TimeStepBase::init_for_postprocessing ()
{
  next_action = postprocess;
}



void
TimeStepBase::solve_dual_problem ()
{
  Assert (false, ExcPureFunctionCalled());
}



void
TimeStepBase::postprocess_timestep ()
{
  Assert (false, ExcPureFunctionCalled());
}



double
TimeStepBase::get_time () const
{
  return time;
}



unsigned int
TimeStepBase::get_timestep_no () const
{
  return timestep_no;
}



double
TimeStepBase::get_backward_timestep () const
{
  Assert (previous_timestep != 0, ExcCantComputeTimestep());
  return time - previous_timestep->time;
}



double
TimeStepBase::get_forward_timestep () const
{
  Assert (next_timestep != 0, ExcCantComputeTimestep());
  return next_timestep->time - time;
}



void
TimeStepBase::set_previous_timestep (const TimeStepBase *previous)
{
  previous_timestep = previous;
}



void
TimeStepBase::set_next_timestep (const TimeStepBase *next)
{
  next_timestep     = next;
}



void
TimeStepBase::set_timestep_no (const unsigned int step_no)
{
  timestep_no = step_no;
}



void
TimeStepBase::set_sweep_no (const unsigned int sweep)
{
  sweep_no = sweep;
}



std::size_t
TimeStepBase::memory_consumption () const
{
  // only simple data types
  return sizeof(*this);
}



template <int dim>
TimeStepBase_Tria<dim>::TimeStepBase_Tria() :
  TimeStepBase (0),
  tria (0, typeid(*this).name()),
  coarse_grid (0, typeid(*this).name()),
  flags (),
  refinement_flags(0)
{
  Assert (false, ExcPureFunctionCalled());
}



template <int dim>
TimeStepBase_Tria<dim>::TimeStepBase_Tria (const double              time,
                                           const Triangulation<dim> &coarse_grid,
                                           const Flags              &flags,
                                           const RefinementFlags    &refinement_flags) :
  TimeStepBase (time),
  tria(0, typeid(*this).name()),
  coarse_grid (&coarse_grid, typeid(*this).name()),
  flags (flags),
  refinement_flags (refinement_flags)
{}



template <int dim>
TimeStepBase_Tria<dim>::~TimeStepBase_Tria ()
{
  if (!flags.delete_and_rebuild_tria)
    {
      Triangulation<dim> *t = tria;
      tria = 0;
      delete t;
    }
  else
    Assert (tria==0, ExcInternalError());

  coarse_grid = 0;
}



template <int dim>
void
TimeStepBase_Tria<dim>::wake_up (const unsigned int wakeup_level)
{
  TimeStepBase::wake_up (wakeup_level);

  if (wakeup_level == flags.wakeup_level_to_build_grid)
    if (flags.delete_and_rebuild_tria || !tria)
      restore_grid ();
}



template <int dim>
void
TimeStepBase_Tria<dim>::sleep (const unsigned int sleep_level)
{
  if (sleep_level == flags.sleep_level_to_delete_grid)
    {
      Assert (tria!=0, ExcInternalError());

      if (flags.delete_and_rebuild_tria)
        {
          Triangulation<dim> *t = tria;
          tria = 0;
          delete t;
        }
    }

  TimeStepBase::sleep (sleep_level);
}



template <int dim>
void TimeStepBase_Tria<dim>::save_refine_flags ()
{
  // for any of the non-initial grids
  // store the refinement flags
  refine_flags.push_back (std::vector<bool>());
  coarsen_flags.push_back (std::vector<bool>());
  tria->save_refine_flags (refine_flags.back());
  tria->save_coarsen_flags (coarsen_flags.back());
}



template <int dim>
void TimeStepBase_Tria<dim>::restore_grid ()
{
  Assert (tria == 0, ExcGridNotDeleted());
  Assert (refine_flags.size() == coarsen_flags.size(),
          ExcInternalError());

  // create a virgin triangulation and
  // set it to a copy of the coarse grid
  tria = new Triangulation<dim> ();
  tria->copy_triangulation (*coarse_grid);

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
//       if (flags.max_refinement_level != 0)
//      {
//        typename Triangulation<dim>::active_cell_iterator cell, endc;
//        for (cell = tria->begin_active(),
//             endc = tria->end();
//             cell!=endc; ++cell)
//          if (static_cast<unsigned int>(cell->level()) >=
//              flags.max_refinement_level)
//            cell->clear_refine_flag();
//      };

      tria->execute_coarsening_and_refinement ();
    };
}



// have a few helper functions
namespace
{
  template <int dim>
  void
  mirror_refinement_flags (const typename Triangulation<dim>::cell_iterator &new_cell,
                           const typename Triangulation<dim>::cell_iterator &old_cell)
  {
    // mirror the refinement
    // flags from the present time level to
    // the previous if the dual problem was
    // used for the refinement, since the
    // error is computed on a time-space cell
    //
    // we don't mirror the coarsening flags
    // since we want stronger refinement. if
    // this was the wrong decision, the error
    // on the child cells of the previous
    // time slab will indicate coarsening
    // in the next iteration, so this is not
    // so dangerous here.
    //
    // also, we only have to check whether
    // the present cell flagged for
    // refinement and the previous one is on
    // the same level and also active. If it
    // already has children, then there is
    // no problem at all, if it is on a lower
    // level than the present one, then it
    // will be refined below anyway.
    if (new_cell->active())
      {
        if (new_cell->refine_flag_set() && old_cell->active())
          {
            if (old_cell->coarsen_flag_set())
              old_cell->clear_coarsen_flag();

            old_cell->set_refine_flag();
          };

        return;
      };

    if (old_cell->has_children() && new_cell->has_children())
      {
        Assert(old_cell->n_children()==new_cell->n_children(), ExcNotImplemented());
        for (unsigned int c=0; c<new_cell->n_children(); ++c)
          dealii::mirror_refinement_flags<dim> (new_cell->child(c), old_cell->child(c));
      }
  }



  template <int dim>
  bool
  adapt_grid_cells (const typename Triangulation<dim>::cell_iterator &cell1,
                    const typename Triangulation<dim>::cell_iterator &cell2)
  {

    if (cell2->has_children() && cell1->has_children())
      {
        bool grids_changed = false;

        Assert(cell2->n_children()==cell1->n_children(), ExcNotImplemented());
        for (unsigned int c=0; c<cell1->n_children(); ++c)
          grids_changed |= dealii::adapt_grid_cells<dim> (cell1->child(c),
                                                          cell2->child(c));
        return grids_changed;
      };


    if (!cell1->has_children() && !cell2->has_children())
      // none of the two have children, so
      // make sure that not one is flagged
      // for refinement and the other for
      // coarsening
      {
        if (cell1->refine_flag_set() && cell2->coarsen_flag_set())
          {
            cell2->clear_coarsen_flag();
            return true;
          }
        else if (cell1->coarsen_flag_set() && cell2->refine_flag_set())
          {
            cell1->clear_coarsen_flag();
            return true;
          };

        return false;
      };


    if (cell1->has_children() && !cell2->has_children())
      // cell1 has children, cell2 has not
      // -> cell2 needs to be refined if any
      // of cell1's children is flagged
      // for refinement. None of them should
      // be refined further, since then in the
      // last round something must have gone
      // wrong
      //
      // if cell2 was flagged for coarsening,
      // we need to clear that flag in any
      // case. The only exception would be
      // if all children of cell1 were
      // flagged for coarsening, but rules
      // for coarsening are so complicated
      // that we will not attempt to cover
      // them. Rather accept one cell which
      // is not coarsened...
      {
        bool changed_grid = false;
        if (cell2->coarsen_flag_set())
          {
            cell2->clear_coarsen_flag();
            changed_grid = true;
          };

        if (!cell2->refine_flag_set())
          for (unsigned int c=0; c<cell1->n_children(); ++c)
            if (cell1->child(c)->refine_flag_set() ||
                cell1->child(c)->has_children())
              {
                cell2->set_refine_flag();
                changed_grid = true;
                break;
              };
        return changed_grid;
      };

    if (!cell1->has_children() && cell2->has_children())
      // same thing, other way round...
      {
        bool changed_grid = false;
        if (cell1->coarsen_flag_set())
          {
            cell1->clear_coarsen_flag();
            changed_grid = true;
          };

        if (!cell1->refine_flag_set())
          for (unsigned int c=0; c<cell2->n_children(); ++c)
            if (cell2->child(c)->refine_flag_set() ||
                cell2->child(c)->has_children())
              {
                cell1->set_refine_flag();
                changed_grid = true;
                break;
              };
        return changed_grid;
      };

    Assert (false, ExcInternalError());
    return false;
  }



  template <int dim>
  bool
  adapt_grids (Triangulation<dim> &tria1,
               Triangulation<dim> &tria2)
  {
    bool grids_changed = false;

    typename Triangulation<dim>::cell_iterator cell1 = tria1.begin(),
                                               cell2 = tria2.begin();
    typename Triangulation<dim>::cell_iterator endc;
    endc = (tria1.n_levels() == 1 ?
            typename Triangulation<dim>::cell_iterator(tria1.end()) :
            tria1.begin(1));
    for (; cell1!=endc; ++cell1, ++cell2)
      grids_changed |= dealii::adapt_grid_cells<dim> (cell1, cell2);

    return grids_changed;
  }
}


template <int dim>
void TimeStepBase_Tria<dim>::refine_grid (const RefinementData refinement_data)
{
  Vector<float> criteria;
  get_tria_refinement_criteria (criteria);

  // copy the following two values since
  // we may need modified values in the
  // process of this function
  double refinement_threshold = refinement_data.refinement_threshold,
         coarsening_threshold = refinement_data.coarsening_threshold;

  // prepare an array where the criteria
  // are stored in a sorted fashion
  // we need this if cell number correction
  // is switched on.
  // the criteria are sorted in ascending
  // order
  // only fill it when needed
  Vector<float> sorted_criteria;
  // two pointers into this array denoting
  // the position where the two thresholds
  // are assumed
  Vector<float>::const_iterator p_refinement_threshold=0,
                                p_coarsening_threshold=0;


  // if we are to do some cell number
  // correction steps, we have to find out
  // which further cells (beyond
  // refinement_threshold) to refine in case
  // we need more cells, and which cells
  // to not refine in case we need less cells
  // (or even to coarsen, if necessary). to
  // this end, we first define pointers into
  // a sorted array of criteria pointing
  // to the thresholds of refinement or
  // coarsening; moving these pointers amounts
  // to changing the threshold such that the
  // number of cells flagged for refinement
  // or coarsening would be changed by one
  if ((timestep_no != 0) &&
      (sweep_no>=refinement_flags.first_sweep_with_correction) &&
      (refinement_flags.cell_number_correction_steps > 0))
    {
      sorted_criteria = criteria;
      std::sort (sorted_criteria.begin(),
                 sorted_criteria.end());
      p_refinement_threshold = Utilities::lower_bound (sorted_criteria.begin(),
                                                       sorted_criteria.end(),
                                                       static_cast<float>(refinement_threshold));
      p_coarsening_threshold = std::upper_bound (sorted_criteria.begin(),
                                                 sorted_criteria.end(),
                                                 static_cast<float>(coarsening_threshold));
    };


  // actually flag cells the first time
  GridRefinement::refine (*tria, criteria, refinement_threshold);
  GridRefinement::coarsen (*tria, criteria, coarsening_threshold);

  // store this number for the following
  // since its computation is rather
  // expensive and since it doesn't change
  const unsigned int n_active_cells = tria->n_active_cells ();

  // if not on first time level: try to
  // adjust the number of resulting
  // cells to those on the previous
  // time level. Only do the cell number
  // correction for higher sweeps and if
  // there are sufficiently many cells
  // already to avoid "grid stall" i.e.
  // that the grid's evolution is hindered
  // by the correction (this usually
  // happens if there are very few cells,
  // since then the number of cells touched
  // by the correction step may exceed the
  // number of cells which are flagged for
  // refinement; in this case it often
  // happens that the number of cells
  // does not grow between sweeps, which
  // clearly is not the wanted behaviour)
  //
  // however, if we do not do anything, we
  // can get into trouble somewhen later.
  // therefore, we also use the correction
  // step for the first sweep or if the
  // number of cells is between 100 and 300
  // (unlike in the first version of the
  // algorithm), but relax the conditions
  // for the correction to allow deviations
  // which are three times as high than
  // allowed (sweep==1 || cell number<200)
  // or twice as high (sweep==2 ||
  // cell number<300). Also, since
  // refinement never does any harm other
  // than increased work, we allow for
  // arbitrary growth of cell number if
  // the estimated cell number is below
  // 200.
  //
  // repeat this loop several times since
  // the first estimate may not be totally
  // correct
  if ((timestep_no != 0) && (sweep_no>=refinement_flags.first_sweep_with_correction))
    for (unsigned int loop=0;
         loop<refinement_flags.cell_number_correction_steps; ++loop)
      {
        Triangulation<dim> *previous_tria
          = dynamic_cast<const TimeStepBase_Tria<dim>*>(previous_timestep)->tria;

        // do one adaption step if desired
        // (there are more coming below then
        // also)
        if (refinement_flags.adapt_grids)
          dealii::adapt_grids<dim> (*previous_tria, *tria);

        // perform flagging of cells
        // needed to regularize the
        // triangulation
        tria->prepare_coarsening_and_refinement ();
        previous_tria->prepare_coarsening_and_refinement ();


        // now count the number of elements
        // which will result on the previous
        // grid after it will be refined. The
        // number which will really result
        // should be approximately that that we
        // compute here, since we already
        // performed most of the prepare*
        // steps for the previous grid
        //
        // use a double value since for each
        // four cells (in 2D) that we flagged
        // for coarsening we result in one
        // new. but since we loop over flagged
        // cells, we have to subtract 3/4 of
        // a cell for each flagged cell
        Assert(!tria->get_anisotropic_refinement_flag(), ExcNotImplemented());
        Assert(!previous_tria->get_anisotropic_refinement_flag(), ExcNotImplemented());
        double previous_cells = previous_tria->n_active_cells();
        typename Triangulation<dim>::active_cell_iterator cell, endc;
        cell = previous_tria->begin_active();
        endc = previous_tria->end();
        for (; cell!=endc; ++cell)
          if (cell->refine_flag_set())
            previous_cells += (GeometryInfo<dim>::max_children_per_cell-1);
          else if (cell->coarsen_flag_set())
            previous_cells -= (GeometryInfo<dim>::max_children_per_cell-1) /
                              GeometryInfo<dim>::max_children_per_cell;

        // @p{previous_cells} now gives the
        // number of cells which would result
        // from the flags on the previous grid
        // if we refined it now. However, some
        // more flags will be set when we adapt
        // the previous grid with this one
        // after the flags have been set for
        // this time level; on the other hand,
        // we don't account for this, since the
        // number of cells on this time level
        // will be changed afterwards by the
        // same way, when it is adapted to the
        // next time level

        // now estimate the number of cells which
        // will result on this level
        double estimated_cells = n_active_cells;
        cell = tria->begin_active();
        endc = tria->end();
        for (; cell!=endc; ++cell)
          if (cell->refine_flag_set())
            estimated_cells += (GeometryInfo<dim>::max_children_per_cell-1);
          else if (cell->coarsen_flag_set())
            estimated_cells -= (GeometryInfo<dim>::max_children_per_cell-1) /
                               GeometryInfo<dim>::max_children_per_cell;

        // calculate the allowed delta in
        // cell numbers; be more lenient
        // if there are few cells
        double delta_up = refinement_flags.cell_number_corridor_top,
               delta_down = refinement_flags.cell_number_corridor_bottom;

        const std::vector<std::pair<unsigned int,double> > &relaxations
          = (sweep_no >= refinement_flags.correction_relaxations.size() ?
             refinement_flags.correction_relaxations.back() :
             refinement_flags.correction_relaxations[sweep_no]);
        for (unsigned int r=0; r!=relaxations.size(); ++r)
          if (n_active_cells < relaxations[r].first)
            {
              delta_up   *= relaxations[r].second;
              delta_down *= relaxations[r].second;
              break;
            };

        // now, if the number of estimated
        // cells exceeds the number of cells
        // on the old time level by more than
        // delta: cut the top threshold
        //
        // note that for each cell that
        // we unflag we have to diminish the
        // estimated number of cells by
        // @p{children_per_cell}.
        if (estimated_cells > previous_cells*(1.+delta_up))
          {
            // only limit the cell number
            // if there will not be less
            // than some number of cells
            //
            // also note that when using the
            // dual estimator, the initial
            // time level is not refined
            // on its own, so we may not
            // limit the number of the second
            // time level on the basis of
            // the initial one; since for
            // the dual estimator, we
            // mirror the refinement
            // flags, the initial level
            // will be passively refined
            // later on.
            if (estimated_cells>refinement_flags.min_cells_for_correction)
              {
                // number of cells by which the
                // new grid is to be diminished
                double delta_cells = estimated_cells -
                                     previous_cells*(1.+delta_up);

                // if we need to reduce the
                // number of cells, we need
                // to raise the thresholds,
                // i.e. move ahead in the
                // sorted array, since this
                // is sorted in ascending
                // order. do so by removing
                // cells tagged for refinement

                for (unsigned int i=0; i<delta_cells;
                     i += GeometryInfo<dim>::max_children_per_cell-1)
                  if (p_refinement_threshold != sorted_criteria.end())
                    ++p_refinement_threshold;
                  else
                    break;
              }
            else
              // too many cells, but we
              // won't do anything about
              // that
              break;
          }
        else
          // likewise: if the estimated number
          // of cells is less than 90 per cent
          // of those at the previous time level:
          // raise threshold by refining
          // additional cells. if we start to
          // run into the area of cells
          // which are to be coarsened, we
          // raise the limit for these too
          if (estimated_cells < previous_cells*(1.-delta_down))
            {
              // number of cells by which the
              // new grid is to be enlarged
              double delta_cells = previous_cells*(1.-delta_down)-estimated_cells;
              // heuristics: usually, if we
              // add @p{delta_cells} to the
              // present state, we end up
              // with much more than only
              // (1-delta_down)*prev_cells
              // because of the effect of
              // regularization and because
              // of adaption to the
              // following grid. Therefore,
              // if we are not in the last
              // correction loop, we try not
              // to add as many cells as seem
              // necessary at first and hope
              // to get closer to the limit
              // this way. Only in the last
              // loop do we have to take the
              // full number to guarantee the
              // wanted result.
              //
              // The value 0.9 is taken from
              // practice, as the additional
              // number of cells introduced
              // by regularization is
              // approximately 10 per cent
              // of the flagged cells.
              if (loop != refinement_flags.cell_number_correction_steps-1)
                delta_cells *= 0.9;

              // if more cells need to be
              // refined, we need to lower
              // the thresholds, i.e. to
              // move to the beginning
              // of sorted_criteria, which is
              // sorted in ascending order
              for (unsigned int i=0; i<delta_cells;
                   i += (GeometryInfo<dim>::max_children_per_cell-1))
                if (p_refinement_threshold != p_coarsening_threshold)
                  --refinement_threshold;
                else if (p_coarsening_threshold != sorted_criteria.begin())
                  --p_coarsening_threshold, --p_refinement_threshold;
                else
                  break;
            }
          else
            // estimated cell number is ok,
            // stop correction steps
            break;

        if (p_refinement_threshold == sorted_criteria.end())
          {
            Assert (p_coarsening_threshold != p_refinement_threshold,
                    ExcInternalError());
            --p_refinement_threshold;
          };

        coarsening_threshold = *p_coarsening_threshold;
        refinement_threshold = *p_refinement_threshold;

        if (coarsening_threshold>=refinement_threshold)
          coarsening_threshold = 0.999*refinement_threshold;

        // now that we have re-adjusted
        // thresholds: clear all refine and
        // coarsening flags and do it all
        // over again
        cell = tria->begin_active();
        endc  = tria->end();
        for (; cell!=endc; ++cell)
          {
            cell->clear_refine_flag ();
            cell->clear_coarsen_flag ();
          };


        // flag cells finally
        GridRefinement::refine (*tria, criteria, refinement_threshold);
        GridRefinement::coarsen (*tria, criteria, coarsening_threshold);
      };

  // if step number is greater than
  // one: adapt this and the previous
  // grid to each other. Don't do so
  // for the initial grid because
  // it is always taken to be the first
  // grid and needs therefore no
  // treatment of its own.
  if ((timestep_no >= 1) && (refinement_flags.adapt_grids))
    {
      Triangulation<dim> *previous_tria
        = dynamic_cast<const TimeStepBase_Tria<dim>*>(previous_timestep)->tria;


      // if we used the dual estimator, we
      // computed the error information on
      // a time slab, rather than on a level
      // of its own. we then mirror the
      // refinement flags we determined for
      // the present level to the previous
      // one
      //
      // do this mirroring only, if cell number
      // adjustment is on, since otherwise
      // strange things may happen
      if (refinement_flags.mirror_flags_to_previous_grid)
        {
          dealii::adapt_grids<dim> (*previous_tria, *tria);

          typename Triangulation<dim>::cell_iterator old_cell, new_cell, endc;
          old_cell = previous_tria->begin(0);
          new_cell = tria->begin(0);
          endc     = tria->end(0);
          for (; new_cell!=endc; ++new_cell, ++old_cell)
            dealii::mirror_refinement_flags<dim> (new_cell, old_cell);
        };

      tria->prepare_coarsening_and_refinement ();
      previous_tria->prepare_coarsening_and_refinement ();

      // adapt present and previous grids
      // to each other: flag additional
      // cells to avoid the previous grid
      // to have cells refined twice more
      // than the present one and vica versa.
      dealii::adapt_grids<dim> (*previous_tria, *tria);
    };
}



template <int dim>
void TimeStepBase_Tria<dim>::init_for_refinement ()
{
  next_action = grid_refinement;
}



template <int dim>
std::size_t
TimeStepBase_Tria<dim>::memory_consumption () const
{
  return (TimeStepBase::memory_consumption () +
          sizeof(tria) +
          MemoryConsumption::memory_consumption (coarse_grid) +
          sizeof(flags) + sizeof(refinement_flags) +
          MemoryConsumption::memory_consumption (refine_flags) +
          MemoryConsumption::memory_consumption (coarsen_flags));
}



template <int dim>
TimeStepBase_Tria_Flags::Flags<dim>::Flags ()
  :
  delete_and_rebuild_tria (false),
  wakeup_level_to_build_grid (0),
  sleep_level_to_delete_grid (0)
{
  Assert (false, ExcInternalError());
}



template <int dim>
TimeStepBase_Tria_Flags::Flags<dim>::Flags (const bool delete_and_rebuild_tria,
                                            const unsigned int wakeup_level_to_build_grid,
                                            const unsigned int sleep_level_to_delete_grid):
  delete_and_rebuild_tria (delete_and_rebuild_tria),
  wakeup_level_to_build_grid (wakeup_level_to_build_grid),
  sleep_level_to_delete_grid (sleep_level_to_delete_grid)
{
//   Assert (!delete_and_rebuild_tria || (wakeup_level_to_build_grid>=1),
//        ExcInvalidParameter(wakeup_level_to_build_grid));
//   Assert (!delete_and_rebuild_tria || (sleep_level_to_delete_grid>=1),
//        ExcInvalidParameter(sleep_level_to_delete_grid));
}


template <int dim>
typename TimeStepBase_Tria_Flags::RefinementFlags<dim>::CorrectionRelaxations
TimeStepBase_Tria_Flags::RefinementFlags<dim>::default_correction_relaxations
(1,    // one element, denoting the first and all subsequent sweeps
 std::vector<std::pair<unsigned int,double> >(1,    // one element, denoting the upper bound
                                              // for the following
                                              // relaxation
                                              std::make_pair (0U, 0.)));


template <int dim>
TimeStepBase_Tria_Flags::RefinementFlags<dim>::
RefinementFlags (const unsigned int max_refinement_level,
                 const unsigned int first_sweep_with_correction,
                 const unsigned int min_cells_for_correction,
                 const double       cell_number_corridor_top,
                 const double       cell_number_corridor_bottom,
                 const CorrectionRelaxations &correction_relaxations,
                 const unsigned int cell_number_correction_steps,
                 const bool         mirror_flags_to_previous_grid,
                 const bool         adapt_grids) :
  max_refinement_level(max_refinement_level),
  first_sweep_with_correction (first_sweep_with_correction),
  min_cells_for_correction(min_cells_for_correction),
  cell_number_corridor_top(cell_number_corridor_top),
  cell_number_corridor_bottom(cell_number_corridor_bottom),
  correction_relaxations (correction_relaxations.size() != 0 ?
                          correction_relaxations :
                          default_correction_relaxations),
  cell_number_correction_steps(cell_number_correction_steps),
  mirror_flags_to_previous_grid(mirror_flags_to_previous_grid),
  adapt_grids(adapt_grids)
{
  Assert (cell_number_corridor_top>=0, ExcInvalidValue (cell_number_corridor_top));
  Assert (cell_number_corridor_bottom>=0, ExcInvalidValue (cell_number_corridor_bottom));
  Assert (cell_number_corridor_bottom<=1, ExcInvalidValue (cell_number_corridor_bottom));
}


template <int dim>
TimeStepBase_Tria_Flags::RefinementData<dim>::
RefinementData (const double         _refinement_threshold,
                const double         _coarsening_threshold) :
  refinement_threshold(_refinement_threshold),
  // in some rare cases it may happen that
  // both thresholds are the same (e.g. if
  // there are many cells with the same
  // error indicator). That would mean that
  // all cells will be flagged for
  // refinement or coarsening, but some will
  // be flagged for both, namely those for
  // which the indicator equals the
  // thresholds. This is forbidden, however.
  //
  // In some rare cases with very few cells
  // we also could get integer round off
  // errors and get problems with
  // the top and bottom fractions.
  //
  // In these case we arbitrarily reduce the
  // bottom threshold by one permille below
  // the top threshold
  coarsening_threshold((_coarsening_threshold == _refinement_threshold ?
                        _coarsening_threshold :
                        0.999*_coarsening_threshold))
{
  Assert (refinement_threshold >= 0, ExcInvalidValue(refinement_threshold));
  Assert (coarsening_threshold >= 0, ExcInvalidValue(coarsening_threshold));
  // allow both thresholds to be zero,
  // since this is needed in case all indicators
  // are zero
  Assert ((coarsening_threshold < refinement_threshold) ||
          ((coarsening_threshold == 0) && (refinement_threshold == 0)),
          ExcInvalidValue (coarsening_threshold));
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "time_dependent.inst"


DEAL_II_NAMESPACE_CLOSE
