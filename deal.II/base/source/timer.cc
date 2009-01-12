//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/timer.h>
#include <base/exceptions.h>
#include <iostream>

// these includes should probably be properly
// ./configure'd using the AC_HEADER_TIME macro:
#include <sys/time.h>
#include <sys/resource.h>

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#include <mpi.h>
#endif


// on SunOS 4.x, getrusage is stated in the man pages and exists, but
// is not declared in resource.h. declare it ourselves
#ifdef NO_HAVE_GETRUSAGE
extern "C" { 
  int getrusage(int who, struct rusage* ru);
}
#endif

DEAL_II_NAMESPACE_OPEN



Timer::Timer()
                :
                cumulative_time (0.),
		cumulative_wall_time (0.)
{
  start();
}



void Timer::start ()
{
  running    = true;

  struct timeval wall_timer;
  gettimeofday(&wall_timer, NULL);
  start_wall_time = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

  rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  start_time = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;

  rusage usage_children;
  getrusage (RUSAGE_CHILDREN, &usage_children);
  start_time_children = usage_children.ru_utime.tv_sec + 1.e-6 * usage_children.ru_utime.tv_usec;
}



double Timer::stop ()
{
  if (running)
    {
      running = false;
      rusage usage;
      getrusage (RUSAGE_SELF, &usage);
      const double dtime = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
      cumulative_time += dtime - start_time;

      rusage usage_children;
      getrusage (RUSAGE_CHILDREN, &usage_children);
      const double dtime_children =
	usage_children.ru_utime.tv_sec + 1.e-6 * usage_children.ru_utime.tv_usec;
      cumulative_time += dtime_children - start_time_children;

      struct timeval wall_timer;
      gettimeofday(&wall_timer, NULL);
      cumulative_wall_time += wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec 
	- start_wall_time;
    }
  return cumulative_time;
}



double Timer::operator() () const
{
  if (running)
    {
      rusage usage;
      getrusage (RUSAGE_SELF, &usage);
      const double dtime =  usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
      
      rusage usage_children;
      getrusage (RUSAGE_CHILDREN, &usage_children);
      const double dtime_children =
	usage_children.ru_utime.tv_sec + 1.e-6 * usage_children.ru_utime.tv_usec;

      return dtime - start_time + dtime_children - start_time_children + cumulative_time;
    }
  else
    return cumulative_time;
}



double Timer::wall_time () const
{
  if (running)
    {
      struct timeval wall_timer;
      gettimeofday(&wall_timer, NULL);
      return wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start_wall_time + 
	cumulative_wall_time;
    }
  else
    return cumulative_wall_time;
}



void Timer::reset ()
{
  cumulative_time = 0.;
  cumulative_wall_time = 0.;
  running         = false;
}



/* ---------------------------- TimerOutput -------------------------- */

TimerOutput::TimerOutput (std::ostream &stream, 
			  const enum OutputFrequency output_frequency,
			  const enum OutputType output_type)
                          :
                          output_frequency (output_frequency),
			  output_type (output_type),
                          out_stream (stream, true)
{}



TimerOutput::TimerOutput (ConditionalOStream &stream,
			  const enum OutputFrequency output_frequency,
			  const enum OutputType output_type)
                          :
                          output_frequency (output_frequency),
			  output_type (output_type),
                          out_stream (stream)
{}



TimerOutput::~TimerOutput()
{
  while (active_sections.size() > 0)
    exit_section();

  if (output_frequency != every_call)
    print_summary();
}



void 
TimerOutput::enter_section (const std::string &section_name)
{
  Assert (section_name.empty() == false,
	  ExcMessage ("Section string is empty."));

  unsigned int this_section_number = numbers::invalid_unsigned_int;

				   // check whether the requested section
				   // already exists
  for (unsigned int i=0; i<section_names.size(); ++i)
    {
				   // if we found it, we have to have the
				   // index into which we write the
				   // data. ensure that the section was not
				   // active before to avoid any mess
      if (section_names[i].compare (section_name) == 0)
	{
	  this_section_number = i;
	  for (unsigned int j=0; j<active_sections.size(); ++j)
	    {
	      if (active_sections[j] == i)
		{
#ifdef DEBUG
		  std::string exc_text = "Cannot enter the already active section " 
		    + section_name;
		  Assert (active_sections[j] != i, ExcMessage (exc_text.data()));
#endif
		  exit_section(section_name);
		}
	    }
	  break;
	}
    }

				   // in case this is a new section, enlarge
				   // our data storage by one, otherwise
				   // just reset the timer
  if (this_section_number == numbers::invalid_unsigned_int)
    {
      section_names.push_back(section_name);
      section_total_cpu_times.push_back (0.);
      section_total_wall_times.push_back (0.);
      section_timers.push_back (Timer());
      section_n_calls.push_back (0);
      this_section_number = section_names.size() - 1;
    }
  else
    {
      section_timers[this_section_number].reset();
      section_timers[this_section_number].start();
    }

				   // now turn the current section into the
				   // list of active sections and increase
				   // the counter for that section
  active_sections.push_back (this_section_number);
  section_n_calls[this_section_number]++;
}



void 
TimerOutput::exit_section (const std::string &section_name)
{
				// if no string is given, exit the last
				// active section.
  unsigned int active_index_to_delete = numbers::invalid_unsigned_int;
  if (section_name.empty() == true)
    active_index_to_delete = active_sections.size() - 1;
				   // if we got a string, we need to find
				   // the index in the list of active
				   // functions
  else
    {
      for (unsigned int i=0; i<active_sections.size(); ++i)
	if (section_names[active_sections[i]].compare (section_name) == 0)
	  {
	    active_index_to_delete = i;
	    break;
	  }
    }

  Assert (active_index_to_delete != numbers::invalid_unsigned_int,
	  ExcMessage ("Cannot delete a section that was never created."));

  const unsigned int section_to_exit = active_sections[active_index_to_delete];

  Assert (section_to_exit < section_names.size(),
	  ExcInternalError());

				// stop the timer for this section.
  section_timers[section_to_exit].stop();
  section_total_wall_times[section_to_exit] += section_timers[section_to_exit].wall_time();

				// get cpu time. on MPI systems, add the
				// local contributions.
				// 
				// TODO: this should rather be in the Timer
				// class itself, shouldn't it?
  double cpu_time = section_timers[section_to_exit]();
  {

				// On MPI, sum up all the local CPU
				// times.
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    double total_cpu_time = 0.;
    int mpiInitialized;
    MPI_Initialized(&mpiInitialized);

    if( mpiInitialized ) 
      {
	MPI_Allreduce (&cpu_time, &total_cpu_time, 1, MPI_DOUBLE, MPI_SUM, 
		       MPI_COMM_WORLD);
	cpu_time = total_cpu_time;
      }
#endif
    section_total_cpu_times[section_to_exit] += cpu_time;
  }


				   // in case we have to print out
				   // something, do that here...
  if (output_frequency != summary)
    {
      std::string output_time;
      std::ostringstream cpu;
      cpu << cpu_time << "s";
      std::ostringstream wall;
      wall << section_timers[section_to_exit].wall_time() << "s";
      if (output_type == cpu_times)
	output_time = " CPU time: " + cpu.str();
      else if (output_type == wall_times)
	output_time = " wall time: " + wall.str() + ".";
      else
	output_time = ", CPU/wall time: " + cpu.str() + " / " + wall.str() + ".";

      out_stream << section_names[section_to_exit] << output_time
		 << std::endl;
    }

				   // delete the index from the list of
				   // active ones
  std::vector<unsigned int>::iterator position_to_delete 
    = active_sections.begin() + active_index_to_delete;
  active_sections.erase(position_to_delete);
}



void 
TimerOutput::print_summary ()
{
				// in case we want to write CPU times
  if (output_type != wall_times)
    {
      double total_cpu_time;

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      {
	double my_cpu_time = timer_all();
	int mpiInitialized;
	MPI_Initialized(&mpiInitialized);

	if( mpiInitialized ) 
	  {
	    MPI_Allreduce (&my_cpu_time, &total_cpu_time, 1, MPI_DOUBLE, MPI_SUM, 
			   MPI_COMM_WORLD);
	  }
	else
	  total_cpu_time = my_cpu_time;
      }
#else
      total_cpu_time = timer_all();
#endif

				   // check that the sum of all times is
				   // less or equal than the total
				   // time. otherwise, we might have
				   // generated a lot of overhead in this
				   // function.
      double check_time = 0.;
      for (unsigned int i=0; i<section_names.size(); ++i)
	{
	  check_time += section_total_cpu_times[i];
	}
      if (check_time > total_cpu_time)
	{
	  total_cpu_time = check_time;
	  out_stream << std::endl << "Sum of counted times is larger than total time. "
		     << "Timer function may have introduced too much overhead." << std::endl;
	}

				   // now generate a nice table
      out_stream << "\n\n"
		 << "+---------------------------------------------+------------"
		 << "+------------+\n"
		 << "| Total CPU time elapsed since start          |";
      std::cout.width(10);
      std::cout.precision(3);
      out_stream << total_cpu_time << "s |            |\n";
      out_stream << "|                                             |            "
		 << "|            |\n";
      out_stream << "| Section                         | no. calls |";
      std::cout.width(10);
      std::cout.precision(3);
      out_stream << "  CPU time "  << " | % of total |\n";
      out_stream << "+---------------------------------------------+------------"
		 << "+------------+";
      for (unsigned int i=0; i<section_names.size(); ++i)
	{
	  std::string name_out = section_names[i];

				// resize the array so that it is always
				// of the same size
	  unsigned int pos_non_space = name_out.find_first_not_of (" ");
	  name_out.erase(0, pos_non_space);
	  name_out.resize (32, ' ');
	  out_stream << std::endl;
	  out_stream << "| " << name_out;
	  out_stream << "| ";
	  std::cout.width(9);
	  out_stream << section_n_calls[i] << " |";
	  std::cout.width(10);
	  std::cout.precision(3);
	  out_stream << section_total_cpu_times[i] << "s |";
	  std::cout.width(10);
	  std::cout.precision(2);
	  out_stream << section_total_cpu_times[i]/total_cpu_time * 100 << "% |";
	}
      out_stream << std::endl
		 << "+---------------------------------------------+"
		 << "------------+------------+\n"
		 << std::endl;
    }

				// in case we want to write out wallclock times
  if (output_type != cpu_times)
    {
      double total_wall_time = timer_all.wall_time();

 				   // check that the sum of all times is
				   // less or equal than the total
				   // time. otherwise, we might have
				   // generated a lot of overhead in this
				   // function.
      double check_time = 0.;
      for (unsigned int i=0; i<section_names.size(); ++i)
	{
	  check_time += section_total_wall_times[i];
	}
      if (check_time > total_wall_time)
	{
	  total_wall_time = check_time;
	  out_stream << std::endl 
		     << "Sum of counted times is larger than total time. "
		     << "Timer function may have introduced too much overhead." 
		     << std::endl;
	}

				   // now generate a nice table
     out_stream << "\n\n"
		 << "+---------------------------------------------+------------"
		 << "+------------+\n"
		 << "| Total wallclock time elapsed from start     |";
      std::cout.width(10);
      std::cout.precision(3);
      out_stream << total_wall_time << "s |            |\n";
      out_stream << "|                                             |            "
		 << "|            |\n";
      out_stream << "| Section                         | no. calls |";
      std::cout.width(10);
      std::cout.precision(3);
      out_stream << "  CPU time "  << " | % of total |\n";
      out_stream << "+---------------------------------------------+------------"
		 << "+------------+";
      for (unsigned int i=0; i<section_names.size(); ++i)
	{
	  std::string name_out = section_names[i];

				// resize the array so that it is always
				// of the same size
	  unsigned int pos_non_space = name_out.find_first_not_of (" ");
	  name_out.erase(0, pos_non_space);
	  name_out.resize (32, ' ');
	  out_stream << std::endl;
	  out_stream << "| " << name_out;
	  out_stream << "| ";
	  std::cout.width(9);
	  out_stream << section_n_calls[i] << " |";
	  std::cout.width(10);
	  std::cout.precision(3);
	  out_stream << section_total_wall_times[i] << "s |";
	  std::cout.width(10);
	  std::cout.precision(2);
	  out_stream << section_total_wall_times[i]/total_wall_time * 100 << "% |";
	}
      out_stream << std::endl
		 << "+---------------------------------------------+"
		 << "------------+------------+\n"
		 << std::endl;
    }

}


DEAL_II_NAMESPACE_CLOSE
