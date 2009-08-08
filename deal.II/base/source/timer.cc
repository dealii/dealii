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
#include <algorithm>

// these includes should probably be properly
// ./configure'd using the AC_HEADER_TIME macro:
#include <sys/time.h>
#include <sys/resource.h>


// on SunOS 4.x, getrusage is stated in the man pages and exists, but
// is not declared in resource.h. declare it ourselves
#ifdef NO_HAVE_GETRUSAGE
extern "C" { 
  int getrusage(int who, struct rusage* ru);
}
#endif

DEAL_II_NAMESPACE_OPEN



				   // in case we use an MPI compiler, need
				   // to create a communicator just for the
				   // current process
Timer::Timer()
                :
                cumulative_time (0.),
		cumulative_wall_time (0.)
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
		, mpi_communicator (MPI_COMM_SELF)
#endif
{
  start();
}



				   // in case we use an MPI compiler, use
				   // the communicator given from input
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
Timer::Timer(MPI_Comm mpi_communicator)
                :
                cumulative_time (0.),
		cumulative_wall_time (0.),
		mpi_communicator (mpi_communicator)
{
  start();
}
#endif



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

				   // in case of MPI, need to get the time
				   // passed by summing the time over all
				   // processes in the network. works also
				   // in case we just want to have the time
				   // of a single thread, since then the
				   // communicator is MPI_COMM_SELF
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      double local_time = dtime - start_time + dtime_children 
	- start_time_children + cumulative_time;

      int mpiInitialized;
      MPI_Initialized(&mpiInitialized);

      if ( mpiInitialized )
	{
	  double global_time = 0.;
	  MPI_Allreduce (&local_time, &global_time, 1, MPI_DOUBLE, MPI_SUM, 
			 mpi_communicator);
	  return global_time;
	}
      else
	return local_time;
#else
      return dtime - start_time + dtime_children - start_time_children + cumulative_time;
#endif
    }
  else
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    {
      int mpiInitialized;
      MPI_Initialized(&mpiInitialized);

      if ( mpiInitialized )
	{
	  double local_time = cumulative_time;
	  double global_time = 0.;
	  MPI_Allreduce (&local_time, &global_time, 1, MPI_DOUBLE, MPI_SUM, 
			 mpi_communicator);
	  return global_time;
	}
      else
	return cumulative_time;
    }
#else
    return cumulative_time;
#endif
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
                          out_stream (stream, true),
			  output_is_enabled (true)
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
			  , mpi_communicator (MPI_COMM_SELF)
#endif
{}



TimerOutput::TimerOutput (ConditionalOStream &stream,
			  const enum OutputFrequency output_frequency,
			  const enum OutputType output_type)
                          :
                          output_frequency (output_frequency),
			  output_type (output_type),
                          out_stream (stream),
			  output_is_enabled (true)
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
			  , mpi_communicator (MPI_COMM_SELF)
#endif
{}


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI

TimerOutput::TimerOutput (MPI_Comm      mpi_communicator,
			  std::ostream &stream, 
			  const enum OutputFrequency output_frequency,
			  const enum OutputType output_type)
                          :
                          output_frequency (output_frequency),
			  output_type (output_type),
                          out_stream (stream, true),
			  output_is_enabled (true), 
			  mpi_communicator (mpi_communicator)
{}



TimerOutput::TimerOutput (MPI_Comm      mpi_communicator,
			  ConditionalOStream &stream,
			  const enum OutputFrequency output_frequency,
			  const enum OutputType output_type)
                          :
                          output_frequency (output_frequency),
			  output_type (output_type),
                          out_stream (stream),
			  output_is_enabled (true),
			  mpi_communicator (mpi_communicator)
{}

#endif


TimerOutput::~TimerOutput()
{
  while (active_sections.size() > 0)
    exit_section();

  if (output_frequency != every_call && output_is_enabled == true)
    print_summary();
}



void 
TimerOutput::enter_section (const std::string &section_name)
{
  Threads::ThreadMutex::ScopedLock lock (mutex);

  Assert (section_name.empty() == false,
	  ExcMessage ("Section string is empty."));

  Assert (std::find (active_sections.begin(), active_sections.end(),
		     section_name) == active_sections.end(),
	  ExcMessage ("Cannot enter the already active section."));
  
  if (sections.find (section_name) == sections.end())
    {
      sections[section_name].total_cpu_time = 0;
      sections[section_name].total_wall_time = 0;
      sections[section_name].n_calls = 0;
    }

  sections[section_name].timer.reset();
  sections[section_name].timer.start();
  sections[section_name].n_calls++;

  active_sections.push_back (section_name);
}



void 
TimerOutput::exit_section (const std::string &section_name)
{
  Assert (active_sections.size() > 0,
	  ExcMessage("Cannot exit any section because none has been entered!"));

  Threads::ThreadMutex::ScopedLock lock (mutex);

  if (section_name != "")
    {
      Assert (sections.find (section_name) != sections.end(),
	      ExcMessage ("Cannot delete a section that was never created."));
      Assert (std::find (active_sections.begin(), active_sections.end(),
			 section_name) != active_sections.end(),
	      ExcMessage ("Cannot delete a section that has not been entered."));
    }

				   // if no string is given, exit the last
				   // active section.
  const std::string actual_section_name = (section_name == "" ?
					   active_sections.back () :
					   section_name);

  sections[actual_section_name].timer.stop();
  sections[actual_section_name].total_wall_time 
    += sections[actual_section_name].timer.wall_time();

				       // get cpu time. on MPI systems, add
				       // the local contributions. we could
				       // do that also in the Timer class
				       // itself, but we didn't initialize
				       // the Timers here according to that
  double cpu_time = sections[actual_section_name].timer();
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
		       mpi_communicator);
	cpu_time = total_cpu_time;
      }
#endif
    sections[actual_section_name].total_cpu_time += cpu_time;
  }


				   // in case we have to print out
				   // something, do that here...
  if (output_frequency != summary && output_is_enabled == true)
    {
      std::string output_time;
      std::ostringstream cpu;
      cpu << cpu_time << "s";
      std::ostringstream wall;
      wall << sections[actual_section_name].timer.wall_time() << "s";
      if (output_type == cpu_times)
	output_time = ", CPU time: " + cpu.str();
      else if (output_type == wall_times)
	output_time = ", wall time: " + wall.str() + ".";
      else
	output_time = ", CPU/wall time: " + cpu.str() + " / " + wall.str() + ".";

      out_stream << actual_section_name << output_time
		 << std::endl;
    }

				   // delete the index from the list of
				   // active ones
  active_sections.erase (std::find (active_sections.begin(), active_sections.end(),
				    actual_section_name));
}



void 
TimerOutput::print_summary () const
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
	    MPI_Allreduce (&my_cpu_time, &total_cpu_time, 1, MPI_DOUBLE, 
			   MPI_SUM, mpi_communicator);
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
      for (std::map<std::string, Section>::const_iterator
	     i = sections.begin(); i!=sections.end(); ++i)
	check_time += i->second.total_cpu_time;
	
      if (check_time > total_cpu_time)
	total_cpu_time = check_time;
      
				       // generate a nice table
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
      out_stream << "+---------------------------------+-----------+------------"
		 << "+------------+";
      for (std::map<std::string, Section>::const_iterator
	     i = sections.begin(); i!=sections.end(); ++i)
	{
	  std::string name_out = i->first;

					   // resize the array so that it is always
					   // of the same size
	  unsigned int pos_non_space = name_out.find_first_not_of (" ");
	  name_out.erase(0, pos_non_space);
	  name_out.resize (32, ' ');
	  out_stream << std::endl;
	  out_stream << "| " << name_out;
	  out_stream << "| ";
	  std::cout.width(9);
	  out_stream << i->second.n_calls << " |";
	  std::cout.width(10);
	  std::cout.precision(3);
	  out_stream << i->second.total_cpu_time << "s |";
	  std::cout.width(10);
	  std::cout.precision(2);
	  out_stream << i->second.total_cpu_time/total_cpu_time * 100 << "% |";
	}
      out_stream << std::endl
		 << "+---------------------------------+-----------+"
		 << "------------+------------+\n"
		 << std::endl;

      if (check_time > total_cpu_time)
	out_stream << std::endl
		   << "Note: The sum of counted times is larger than the total time.\n"
		   << "(Timer function may have introduced too much overhead, or different\n"
		   << "section timers may have run at the same time.)" << std::endl;
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
      for (std::map<std::string, Section>::const_iterator
	     i = sections.begin(); i!=sections.end(); ++i)
	check_time += i->second.total_wall_time;
      
      if (check_time > total_wall_time)
	total_wall_time = check_time;

				       // now generate a nice table
      out_stream << "\n\n"
		 << "+---------------------------------------------+------------"
		 << "+------------+\n"
		 << "| Total wallclock time elapsed since start    |";
      std::cout.width(10);
      std::cout.precision(3);
      out_stream << total_wall_time << "s |            |\n";
      out_stream << "|                                             |            "
		 << "|            |\n";
      out_stream << "| Section                         | no. calls |";
      std::cout.width(10);
      std::cout.precision(3);
      out_stream << "  wall time | % of total |\n";
      out_stream << "+---------------------------------+-----------+------------"
		 << "+------------+";
      for (std::map<std::string, Section>::const_iterator
	     i = sections.begin(); i!=sections.end(); ++i)
	{
	  std::string name_out = i->first;

				// resize the array so that it is always
				// of the same size
	  unsigned int pos_non_space = name_out.find_first_not_of (" ");
	  name_out.erase(0, pos_non_space);
	  name_out.resize (32, ' ');
	  out_stream << std::endl;
	  out_stream << "| " << name_out;
	  out_stream << "| ";
	  std::cout.width(9);
	  out_stream << i->second.n_calls << " |";
	  std::cout.width(10);
	  std::cout.precision(3);
	  out_stream << i->second.total_wall_time << "s |";
	  std::cout.width(10);
	  std::cout.precision(2);
	  out_stream << i->second.total_wall_time/total_wall_time * 100 << "% |";
	}
      out_stream << std::endl
		 << "+---------------------------------+-----------+"
		 << "------------+------------+\n"
		 << std::endl;

      if (check_time > total_wall_time)
	out_stream << std::endl
		   << "Note: The sum of counted times is larger than the total time.\n"
		   << "(Timer function may have introduced too much overhead, or different\n"
		   << "section timers may have run at the same time.)" << std::endl;
    }
}



void 
TimerOutput::disable_output ()
{
  output_is_enabled = false;
}



void
TimerOutput::enable_output ()
{
  output_is_enabled = true;
}


DEAL_II_NAMESPACE_CLOSE
