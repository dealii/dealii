// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <boost/io/ios_state.hpp>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>

#ifdef DEAL_II_HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
#endif

#ifdef DEAL_II_MSVC
#  include <windows.h>
#endif



DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TimerImplementation
  {
    namespace
    {
      /**
       * Type trait for checking whether or not a type is a
       * std::chrono::duration.
       */
      template <typename T>
      struct is_duration : std::false_type
      {};

      /**
       * Specialization to get the right truth value.
       */
      template <typename Rep, typename Period>
      struct is_duration<std::chrono::duration<Rep, Period>> : std::true_type
      {};

      /**
       * Convert a double precision number with units of seconds into a
       * specified duration type T. Only valid when T is a
       * std::chrono::duration type.
       */
      template <typename T>
      T
      from_seconds(const double time)
      {
        static_assert(is_duration<T>::value,
                      "The template type should be a duration type.");
        return T(std::lround(T::period::den * (time / T::period::num)));
      }

      /**
       * Convert a given duration into a double precision number with units of
       * seconds.
       */
      template <typename Rep, typename Period>
      double
      to_seconds(const std::chrono::duration<Rep, Period> duration)
      {
        return Period::num * double(duration.count()) / Period::den;
      }

      /**
       * Fill a MinMaxAvg struct with default values.
       */
      void
      clear_timing_data(Utilities::MPI::MinMaxAvg &data)
      {
        data.sum       = numbers::signaling_nan<double>();
        data.min       = numbers::signaling_nan<double>();
        data.max       = numbers::signaling_nan<double>();
        data.avg       = numbers::signaling_nan<double>();
        data.min_index = numbers::invalid_unsigned_int;
        data.max_index = numbers::invalid_unsigned_int;
      }
    } // namespace
  }   // namespace TimerImplementation
} // namespace internal



CPUClock::time_point
CPUClock::now() noexcept
{
  double system_cpu_duration = 0.0;
#ifdef DEAL_II_MSVC
  FILETIME   cpuTime, sysTime, createTime, exitTime;
  const auto succeeded = GetProcessTimes(
    GetCurrentProcess(), &createTime, &exitTime, &sysTime, &cpuTime);
  if (succeeded)
    {
      system_cpu_duration =
        (double)(((unsigned long long)cpuTime.dwHighDateTime << 32) |
                 cpuTime.dwLowDateTime) /
        1e7;
    }
    // keep the zero value if GetProcessTimes didn't work
#elif defined(DEAL_II_HAVE_SYS_RESOURCE_H)
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  system_cpu_duration = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
#else
  DEAL_II_WARNING("Unsupported platform. Porting not finished.")
#endif
  return time_point(
    internal::TimerImplementation::from_seconds<duration>(system_cpu_duration));
}



template <typename clock_type_>
Timer::ClockMeasurements<clock_type_>::ClockMeasurements()
  : current_lap_start_time(clock_type::now())
  , accumulated_time(duration_type::zero())
  , last_lap_time(duration_type::zero())
{}



template <typename clock_type_>
void
Timer::ClockMeasurements<clock_type_>::reset()
{
  current_lap_start_time = clock_type::now();
  accumulated_time       = duration_type::zero();
  last_lap_time          = duration_type::zero();
}



Timer::Timer()
  : Timer(MPI_COMM_SELF, /*sync_lap_times=*/false)
{}



Timer::Timer(MPI_Comm mpi_communicator, const bool sync_lap_times_)
  : running(false)
  , mpi_communicator(mpi_communicator)
  , sync_lap_times(sync_lap_times_)
{
  reset();
  start();
}



void
Timer::start()
{
  running = true;
#ifdef DEAL_II_WITH_MPI
  if (sync_lap_times)
    {
      const int ierr = MPI_Barrier(mpi_communicator);
      AssertThrowMPI(ierr);
    }
#endif
  wall_times.current_lap_start_time = wall_clock_type::now();
  cpu_times.current_lap_start_time  = cpu_clock_type::now();
}



double
Timer::stop()
{
  if (running)
    {
      running = false;

      wall_times.last_lap_time =
        wall_clock_type::now() - wall_times.current_lap_start_time;
      cpu_times.last_lap_time =
        cpu_clock_type::now() - cpu_times.current_lap_start_time;

      last_lap_wall_time_data =
        Utilities::MPI::min_max_avg(internal::TimerImplementation::to_seconds(
                                      wall_times.last_lap_time),
                                    mpi_communicator);
      if (sync_lap_times)
        {
          wall_times.last_lap_time =
            internal::TimerImplementation::from_seconds<decltype(
              wall_times)::duration_type>(last_lap_wall_time_data.max);
          cpu_times.last_lap_time =
            internal::TimerImplementation::from_seconds<decltype(
              cpu_times)::duration_type>(
              Utilities::MPI::min_max_avg(
                internal::TimerImplementation::to_seconds(
                  cpu_times.last_lap_time),
                mpi_communicator)
                .max);
        }
      wall_times.accumulated_time += wall_times.last_lap_time;
      cpu_times.accumulated_time += cpu_times.last_lap_time;
      accumulated_wall_time_data =
        Utilities::MPI::min_max_avg(internal::TimerImplementation::to_seconds(
                                      wall_times.accumulated_time),
                                    mpi_communicator);
    }
  return internal::TimerImplementation::to_seconds(cpu_times.accumulated_time);
}



double
Timer::cpu_time() const
{
  if (running)
    {
      const double running_time = internal::TimerImplementation::to_seconds(
        cpu_clock_type::now() - cpu_times.current_lap_start_time +
        cpu_times.accumulated_time);
      return Utilities::MPI::sum(running_time, mpi_communicator);
    }
  else
    {
      return Utilities::MPI::sum(internal::TimerImplementation::to_seconds(
                                   cpu_times.accumulated_time),
                                 mpi_communicator);
    }
}



double
Timer::last_cpu_time() const
{
  return internal::TimerImplementation::to_seconds(cpu_times.last_lap_time);
}



double
Timer::wall_time() const
{
  wall_clock_type::duration current_elapsed_wall_time;
  if (running)
    current_elapsed_wall_time = wall_clock_type::now() -
                                wall_times.current_lap_start_time +
                                wall_times.accumulated_time;
  else
    current_elapsed_wall_time = wall_times.accumulated_time;

  return internal::TimerImplementation::to_seconds(current_elapsed_wall_time);
}



double
Timer::last_wall_time() const
{
  return internal::TimerImplementation::to_seconds(wall_times.last_lap_time);
}



void
Timer::reset()
{
  wall_times.reset();
  cpu_times.reset();
  running = false;
  internal::TimerImplementation::clear_timing_data(last_lap_wall_time_data);
  internal::TimerImplementation::clear_timing_data(accumulated_wall_time_data);
}



/* ---------------------------- TimerOutput -------------------------- */

TimerOutput::TimerOutput(std::ostream &        stream,
                         const OutputFrequency output_frequency,
                         const OutputType      output_type)
  : output_frequency(output_frequency)
  , output_type(output_type)
  , out_stream(stream, true)
  , output_is_enabled(true)
  , mpi_communicator(MPI_COMM_SELF)
{}



TimerOutput::TimerOutput(ConditionalOStream &  stream,
                         const OutputFrequency output_frequency,
                         const OutputType      output_type)
  : output_frequency(output_frequency)
  , output_type(output_type)
  , out_stream(stream)
  , output_is_enabled(true)
  , mpi_communicator(MPI_COMM_SELF)
{}



TimerOutput::TimerOutput(MPI_Comm              mpi_communicator,
                         std::ostream &        stream,
                         const OutputFrequency output_frequency,
                         const OutputType      output_type)
  : output_frequency(output_frequency)
  , output_type(output_type)
  , out_stream(stream, true)
  , output_is_enabled(true)
  , mpi_communicator(mpi_communicator)
{}



TimerOutput::TimerOutput(MPI_Comm              mpi_communicator,
                         ConditionalOStream &  stream,
                         const OutputFrequency output_frequency,
                         const OutputType      output_type)
  : output_frequency(output_frequency)
  , output_type(output_type)
  , out_stream(stream)
  , output_is_enabled(true)
  , mpi_communicator(mpi_communicator)
{}



TimerOutput::~TimerOutput()
{
  auto do_exit = [this]() {
    try
      {
        while (active_sections.size() > 0)
          leave_subsection();
        // don't print unless we leave all subsections
        if ((output_frequency == summary ||
             output_frequency == every_call_and_summary) &&
            output_is_enabled == true)
          print_summary();
      }
    catch (...)
      {}
  };

  // avoid communicating with other processes if there is an uncaught
  // exception
#ifdef DEAL_II_WITH_MPI
#  if __cpp_lib_uncaught_exceptions >= 201411
  // std::uncaught_exception() is deprecated in c++17
  if (std::uncaught_exceptions() > 0 && mpi_communicator != MPI_COMM_SELF)
#  else
  if (std::uncaught_exception() == true && mpi_communicator != MPI_COMM_SELF)
#  endif
    {
      const unsigned int myid =
        Utilities::MPI::this_mpi_process(mpi_communicator);
      if (myid == 0)
        std::cerr
          << "---------------------------------------------------------\n"
          << "TimerOutput objects finalize timed values printed to the\n"
          << "screen by communicating over MPI in their destructors.\n"
          << "Since an exception is currently uncaught, this\n"
          << "synchronization (and subsequent output) will be skipped\n"
          << "to avoid a possible deadlock.\n"
          << "---------------------------------------------------------"
          << std::endl;
    }
  else
    {
      do_exit();
    }
#else
  do_exit();
#endif
}



void
TimerOutput::enter_subsection(const std::string &section_name)
{
  std::lock_guard<std::mutex> lock(mutex);

  Assert(section_name.empty() == false, ExcMessage("Section string is empty."));

  Assert(std::find(active_sections.begin(),
                   active_sections.end(),
                   section_name) == active_sections.end(),
         ExcMessage(std::string("Cannot enter the already active section <") +
                    section_name + ">."));

  if (sections.find(section_name) == sections.end())
    {
      if (mpi_communicator != MPI_COMM_SELF)
        {
          // create a new timer for this section. the second argument
          // will ensure that we have an MPI barrier before starting
          // and stopping a timer, and this ensures that we get the
          // maximum run time for this section over all processors.
          // The mpi_communicator from TimerOutput is passed to the
          // Timer here, so this Timer will collect timing information
          // among all processes inside mpi_communicator.
          sections[section_name].timer = Timer(mpi_communicator, true);
        }


      sections[section_name].total_cpu_time  = 0;
      sections[section_name].total_wall_time = 0;
      sections[section_name].n_calls         = 0;
    }

  sections[section_name].timer.reset();
  sections[section_name].timer.start();
  sections[section_name].n_calls++;

  active_sections.push_back(section_name);
}



void
TimerOutput::leave_subsection(const std::string &section_name)
{
  Assert(!active_sections.empty(),
         ExcMessage("Cannot exit any section because none has been entered!"));

  std::lock_guard<std::mutex> lock(mutex);

  if (!section_name.empty())
    {
      Assert(sections.find(section_name) != sections.end(),
             ExcMessage("Cannot delete a section that was never created."));
      Assert(std::find(active_sections.begin(),
                       active_sections.end(),
                       section_name) != active_sections.end(),
             ExcMessage("Cannot delete a section that has not been entered."));
    }

  // if no string is given, exit the last
  // active section.
  const std::string actual_section_name =
    (section_name.empty() ? active_sections.back() : section_name);

  sections[actual_section_name].timer.stop();
  sections[actual_section_name].total_wall_time +=
    sections[actual_section_name].timer.last_wall_time();

  // Get cpu time. On MPI systems, if constructed with an mpi_communicator
  // like MPI_COMM_WORLD, then the Timer will sum up the CPU time between
  // processors among the provided mpi_communicator. Therefore, no
  // communication is needed here.
  const double cpu_time = sections[actual_section_name].timer.last_cpu_time();
  sections[actual_section_name].total_cpu_time += cpu_time;

  // in case we have to print out something, do that here...
  if ((output_frequency == every_call ||
       output_frequency == every_call_and_summary) &&
      output_is_enabled == true)
    {
      std::string        output_time;
      std::ostringstream cpu;
      cpu << cpu_time << "s";
      std::ostringstream wall;
      wall << sections[actual_section_name].timer.last_wall_time() << "s";
      if (output_type == cpu_times)
        output_time = ", CPU time: " + cpu.str();
      else if (output_type == wall_times)
        output_time = ", wall time: " + wall.str() + ".";
      else
        output_time =
          ", CPU/wall time: " + cpu.str() + " / " + wall.str() + ".";

      out_stream << actual_section_name << output_time << std::endl;
    }

  // delete the index from the list of
  // active ones
  active_sections.erase(std::find(active_sections.begin(),
                                  active_sections.end(),
                                  actual_section_name));
}



std::map<std::string, double>
TimerOutput::get_summary_data(const OutputData kind) const
{
  std::map<std::string, double> output;
  for (const auto &section : sections)
    {
      switch (kind)
        {
          case TimerOutput::OutputData::total_cpu_time:
            output[section.first] = section.second.total_cpu_time;
            break;
          case TimerOutput::OutputData::total_wall_time:
            output[section.first] = section.second.total_wall_time;
            break;
          case TimerOutput::OutputData::n_calls:
            output[section.first] = section.second.n_calls;
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
  return output;
}



void
TimerOutput::print_summary() const
{
  // we are going to change the precision and width of output below. store the
  // old values so the get restored when exiting this function
  const boost::io::ios_base_all_saver restore_stream(out_stream.get_stream());

  // get the maximum width among all sections
  unsigned int max_width = 0;
  for (const auto &i : sections)
    max_width =
      std::max(max_width, static_cast<unsigned int>(i.first.length()));

  // 32 is the default width until | character
  max_width = std::max(max_width + 1, static_cast<unsigned int>(32));
  const std::string extra_dash  = std::string(max_width - 32, '-');
  const std::string extra_space = std::string(max_width - 32, ' ');

  if (output_type != cpu_and_wall_times_grouped)
    {
      // in case we want to write CPU times
      if (output_type != wall_times)
        {
          double total_cpu_time =
            Utilities::MPI::sum(timer_all.cpu_time(), mpi_communicator);

          // check that the sum of all times is less or equal than the total
          // time. otherwise, we might have generated a lot of overhead in this
          // function.
          double check_time = 0.;
          for (const auto &i : sections)
            check_time += i.second.total_cpu_time;

          const double time_gap = check_time - total_cpu_time;
          if (time_gap > 0.0)
            total_cpu_time = check_time;

          // generate a nice table
          out_stream << "\n\n"
                     << "+---------------------------------------------"
                     << extra_dash << "+------------"
                     << "+------------+\n"
                     << "| Total CPU time elapsed since start          "
                     << extra_space << "|";
          out_stream << std::setw(10) << std::setprecision(3) << std::right;
          out_stream << total_cpu_time << "s |            |\n";
          out_stream << "|                                             "
                     << extra_space << "|            "
                     << "|            |\n";
          out_stream << "| Section                         " << extra_space
                     << "| no. calls |";
          out_stream << std::setw(10);
          out_stream << std::setprecision(3);
          out_stream << "  CPU time "
                     << " | % of total |\n";
          out_stream << "+---------------------------------" << extra_dash
                     << "+-----------+------------"
                     << "+------------+";
          for (const auto &i : sections)
            {
              std::string name_out = i.first;

              // resize the array so that it is always of the same size
              unsigned int pos_non_space = name_out.find_first_not_of(' ');
              name_out.erase(0, pos_non_space);
              name_out.resize(max_width, ' ');
              out_stream << std::endl;
              out_stream << "| " << name_out;
              out_stream << "| ";
              out_stream << std::setw(9);
              out_stream << i.second.n_calls << " |";
              out_stream << std::setw(10);
              out_stream << std::setprecision(3);
              out_stream << i.second.total_cpu_time << "s |";
              out_stream << std::setw(10);
              if (total_cpu_time != 0)
                {
                  // if run time was less than 0.1%, just print a zero to avoid
                  // printing silly things such as "2.45e-6%". otherwise print
                  // the actual percentage
                  const double fraction =
                    i.second.total_cpu_time / total_cpu_time;
                  if (fraction > 0.001)
                    {
                      out_stream << std::setprecision(2);
                      out_stream << fraction * 100;
                    }
                  else
                    out_stream << 0.0;

                  out_stream << "% |";
                }
              else
                out_stream << 0.0 << "% |";
            }
          out_stream << std::endl
                     << "+---------------------------------" << extra_dash
                     << "+-----------+"
                     << "------------+------------+\n"
                     << std::endl;

          if (time_gap > 0.0)
            out_stream
              << std::endl
              << "Note: The sum of counted times is " << time_gap
              << " seconds larger than the total time.\n"
              << "(Timer function may have introduced too much overhead, or different\n"
              << "section timers may have run at the same time.)" << std::endl;
        }

      // in case we want to write out wallclock times
      if (output_type != cpu_times)
        {
          double total_wall_time = timer_all.wall_time();

          // now generate a nice table
          out_stream << "\n\n"
                     << "+---------------------------------------------"
                     << extra_dash << "+------------"
                     << "+------------+\n"
                     << "| Total wallclock time elapsed since start    "
                     << extra_space << "|";
          out_stream << std::setw(10) << std::setprecision(3) << std::right;
          out_stream << total_wall_time << "s |            |\n";
          out_stream << "|                                             "
                     << extra_space << "|            "
                     << "|            |\n";
          out_stream << "| Section                         " << extra_space
                     << "| no. calls |";
          out_stream << std::setw(10);
          out_stream << std::setprecision(3);
          out_stream << "  wall time | % of total |\n";
          out_stream << "+---------------------------------" << extra_dash
                     << "+-----------+------------"
                     << "+------------+";
          for (const auto &i : sections)
            {
              std::string name_out = i.first;

              // resize the array so that it is always of the same size
              unsigned int pos_non_space = name_out.find_first_not_of(' ');
              name_out.erase(0, pos_non_space);
              name_out.resize(max_width, ' ');
              out_stream << std::endl;
              out_stream << "| " << name_out;
              out_stream << "| ";
              out_stream << std::setw(9);
              out_stream << i.second.n_calls << " |";
              out_stream << std::setw(10);
              out_stream << std::setprecision(3);
              out_stream << i.second.total_wall_time << "s |";
              out_stream << std::setw(10);

              if (total_wall_time != 0)
                {
                  // if run time was less than 0.1%, just print a zero to avoid
                  // printing silly things such as "2.45e-6%". otherwise print
                  // the actual percentage
                  const double fraction =
                    i.second.total_wall_time / total_wall_time;
                  if (fraction > 0.001)
                    {
                      out_stream << std::setprecision(2);
                      out_stream << fraction * 100;
                    }
                  else
                    out_stream << 0.0;

                  out_stream << "% |";
                }
              else
                out_stream << 0.0 << "% |";
            }
          out_stream << std::endl
                     << "+---------------------------------" << extra_dash
                     << "+-----------+"
                     << "------------+------------+\n"
                     << std::endl;
        }
    }
  else
    // output_type == cpu_and_wall_times_grouped
    {
      const double total_wall_time = timer_all.wall_time();
      double       total_cpu_time =
        Utilities::MPI::sum(timer_all.cpu_time(), mpi_communicator);

      // check that the sum of all times is less or equal than the total time.
      // otherwise, we might have generated a lot of overhead in this function.
      double check_time = 0.;

      for (const auto &i : sections)
        check_time += i.second.total_cpu_time;

      const double time_gap = check_time - total_cpu_time;
      if (time_gap > 0.0)
        total_cpu_time = check_time;

      // generate a nice table
      out_stream << "\n\n+---------------------------------------------"
                 << extra_dash << "+"
                 << "------------+------------+"
                 << "------------+------------+"
                 << "\n"
                 << "| Total CPU/wall time elapsed since start     "
                 << extra_space << "|" << std::setw(10) << std::setprecision(3)
                 << std::right << total_cpu_time << "s |            |"
                 << total_wall_time << "s |            |"
                 << "\n|                                             "
                 << extra_space << "|"
                 << "            |            |"
                 << "            |            |"
                 << "\n| Section                         " << extra_space
                 << "| no. calls |"
                 << "  CPU time  | % of total |"
                 << "  wall time | % of total |"
                 << "\n+---------------------------------" << extra_dash
                 << "+-----------+"
                 << "------------+------------+"
                 << "------------+------------+" << std::endl;

      for (const auto &i : sections)
        {
          std::string name_out = i.first;

          // resize the array so that it is always of the same size
          unsigned int pos_non_space = name_out.find_first_not_of(' ');
          name_out.erase(0, pos_non_space);
          name_out.resize(max_width, ' ');
          out_stream << "| " << name_out << "| ";

          out_stream << std::setw(9);
          out_stream << i.second.n_calls << " |";

          if (output_type != wall_times)
            {
              out_stream << std::setw(10);
              out_stream << std::setprecision(3);
              out_stream << i.second.total_cpu_time << "s |";
              out_stream << std::setw(10);
              if (total_cpu_time != 0)
                {
                  // if run time was less than 0.1%, just print a zero to avoid
                  // printing silly things such as "2.45e-6%". otherwise print
                  // the actual percentage
                  const double fraction =
                    i.second.total_cpu_time / total_cpu_time;
                  if (fraction > 0.001)
                    {
                      out_stream << std::setprecision(2);
                      out_stream << fraction * 100;
                    }
                  else
                    out_stream << 0.0;

                  out_stream << "% |";
                }
              else
                out_stream << 0.0 << "% |";
            }

          if (output_type != cpu_times)
            {
              out_stream << std::setw(10);
              out_stream << std::setprecision(3);
              out_stream << i.second.total_wall_time << "s |";
              out_stream << std::setw(10);

              if (total_wall_time != 0)
                {
                  // if run time was less than 0.1%, just print a zero to avoid
                  // printing silly things such as "2.45e-6%". otherwise print
                  // the actual percentage
                  const double fraction =
                    i.second.total_wall_time / total_wall_time;
                  if (fraction > 0.001)
                    {
                      out_stream << std::setprecision(2);
                      out_stream << fraction * 100;
                    }
                  else
                    out_stream << 0.0;

                  out_stream << "% |";
                }
              else
                out_stream << 0.0 << "% |";
            }
          out_stream << std::endl;
        }

      out_stream << "+---------------------------------" << extra_dash
                 << "+-----------+"
                 << "------------+------------+"
                 << "------------+------------+" << std::endl
                 << std::endl;

      if (output_type != wall_times && time_gap > 0.0)
        out_stream
          << std::endl
          << "Note: The sum of counted times is " << time_gap
          << " seconds larger than the total time.\n"
          << "(Timer function may have introduced too much overhead, or different\n"
          << "section timers may have run at the same time.)" << std::endl;
    }
}



void
TimerOutput::print_wall_time_statistics(const MPI_Comm mpi_comm,
                                        const double   quantile) const
{
  // we are going to change the precision and width of output below. store the
  // old values so the get restored when exiting this function
  const boost::io::ios_base_all_saver restore_stream(out_stream.get_stream());

  AssertDimension(sections.size(),
                  Utilities::MPI::max(sections.size(), mpi_comm));
  Assert(quantile >= 0. && quantile <= 0.5,
         ExcMessage("The quantile must be between 0 and 0.5"));

  // get the maximum width among all sections
  unsigned int max_width = 0;
  for (const auto &i : sections)
    max_width =
      std::max(max_width, static_cast<unsigned int>(i.first.length()));

  // 17 is the default width until | character
  max_width = std::max(max_width + 1, static_cast<unsigned int>(17));
  const std::string extra_dash  = std::string(max_width - 17, '-');
  const std::string extra_space = std::string(max_width - 17, ' ');

  // function to print data in a nice table
  const auto print_statistics = [&](const double given_time) {
    const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(mpi_comm);
    if (n_ranks == 1 || quantile == 0.)
      {
        Utilities::MPI::MinMaxAvg data =
          Utilities::MPI::min_max_avg(given_time, mpi_comm);

        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << data.min << "s ";
        out_stream << std::setw(5) << std::right;
        out_stream << data.min_index << (n_ranks > 99999 ? "" : " ") << "|";
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << data.avg << "s |";
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << data.max << "s ";
        out_stream << std::setw(5) << std::right;
        out_stream << data.max_index << (n_ranks > 99999 ? "" : " ") << "|\n";
      }
    else
      {
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(mpi_comm);
        std::vector<double> receive_data(my_rank == 0 ? n_ranks : 0);
        std::vector<double> result(9);
#ifdef DEAL_II_WITH_MPI
        MPI_Gather(&given_time,
                   1,
                   MPI_DOUBLE,
                   receive_data.data(),
                   1,
                   MPI_DOUBLE,
                   0,
                   mpi_comm);
        if (my_rank == 0)
          {
            // fill the received data in a pair and sort; on the way, also
            // compute the average
            std::vector<std::pair<double, unsigned int>> data_rank;
            data_rank.reserve(n_ranks);
            for (unsigned int i = 0; i < n_ranks; ++i)
              {
                data_rank.emplace_back(receive_data[i], i);
                result[4] += receive_data[i];
              }
            result[4] /= n_ranks;
            std::sort(data_rank.begin(), data_rank.end());

            const unsigned int quantile_index =
              static_cast<unsigned int>(std::round(quantile * n_ranks));
            AssertIndexRange(quantile_index, data_rank.size());
            result[0] = data_rank[0].first;
            result[1] = data_rank[0].second;
            result[2] = data_rank[quantile_index].first;
            result[3] = data_rank[quantile_index].second;
            result[5] = data_rank[n_ranks - 1 - quantile_index].first;
            result[6] = data_rank[n_ranks - 1 - quantile_index].second;
            result[7] = data_rank[n_ranks - 1].first;
            result[8] = data_rank[n_ranks - 1].second;
          }
        MPI_Bcast(result.data(), 9, MPI_DOUBLE, 0, mpi_comm);
#endif
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << result[0] << "s ";
        out_stream << std::setw(5) << std::right;
        out_stream << static_cast<unsigned int>(result[1])
                   << (n_ranks > 99999 ? "" : " ") << "|";
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << result[2] << "s ";
        out_stream << std::setw(5) << std::right;
        out_stream << static_cast<unsigned int>(result[3])
                   << (n_ranks > 99999 ? "" : " ") << "|";
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << result[4] << "s |";
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << result[5] << "s ";
        out_stream << std::setw(5) << std::right;
        out_stream << static_cast<unsigned int>(result[6])
                   << (n_ranks > 99999 ? "" : " ") << "|";
        out_stream << std::setw(10) << std::setprecision(4) << std::right;
        out_stream << result[7] << "s ";
        out_stream << std::setw(5) << std::right;
        out_stream << static_cast<unsigned int>(result[8])
                   << (n_ranks > 99999 ? "" : " ") << "|\n";
      }
  };

  // in case we want to write out wallclock times
  {
    const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(mpi_comm);

    const std::string time_rank_column = "------------------+";
    const std::string time_rank_space  = "                  |";

    // now generate a nice table
    out_stream << "\n"
               << "+------------------------------" << extra_dash << "+"
               << time_rank_column
               << (n_ranks > 1 && quantile > 0. ? time_rank_column : "")
               << "------------+"
               << (n_ranks > 1 && quantile > 0. ? time_rank_column : "")
               << time_rank_column << "\n"
               << "| Total wallclock time elapsed " << extra_space << "|";

    print_statistics(timer_all.wall_time());

    out_stream << "|                              " << extra_space << "|"
               << time_rank_space
               << (n_ranks > 1 && quantile > 0. ? time_rank_space : "")
               << "             "
               << (n_ranks > 1 && quantile > 0. ? time_rank_space : "")
               << time_rank_space << "\n";
    out_stream << "| Section          " << extra_space << "| no. calls "
               << "|   min time  rank |";
    if (n_ranks > 1 && quantile > 0.)
      out_stream << " " << std::setw(5) << std::setprecision(2) << std::right
                 << quantile << "-tile  rank |";
    out_stream << "   avg time |";
    if (n_ranks > 1 && quantile > 0.)
      out_stream << " " << std::setw(5) << std::setprecision(2) << std::right
                 << 1. - quantile << "-tile  rank |";
    out_stream << "   max time  rank |\n";
    out_stream << "+------------------------------" << extra_dash << "+"
               << time_rank_column
               << (n_ranks > 1 && quantile > 0. ? time_rank_column : "")
               << "------------+"
               << (n_ranks > 1 && quantile > 0. ? time_rank_column : "")
               << time_rank_column << "\n";
    for (const auto &i : sections)
      {
        std::string name_out = i.first;

        // resize the array so that it is always of the same size
        unsigned int pos_non_space = name_out.find_first_not_of(' ');
        name_out.erase(0, pos_non_space);
        name_out.resize(max_width, ' ');
        out_stream << "| " << name_out;
        out_stream << "| ";
        out_stream << std::setw(9);
        out_stream << i.second.n_calls << " |";

        print_statistics(i.second.total_wall_time);
      }
    out_stream << "+------------------------------" << extra_dash << "+"
               << time_rank_column
               << (n_ranks > 1 && quantile > 0. ? time_rank_column : "")
               << "------------+"
               << (n_ranks > 1 && quantile > 0. ? time_rank_column : "")
               << time_rank_column << "\n";
  }
}



void
TimerOutput::disable_output()
{
  output_is_enabled = false;
}



void
TimerOutput::enable_output()
{
  output_is_enabled = true;
}

void
TimerOutput::reset()
{
  std::lock_guard<std::mutex> lock(mutex);
  sections.clear();
  active_sections.clear();
  timer_all.restart();
}

TimerOutput::Scope::~Scope()
{
  try
    {
      stop();
    }
  catch (...)
    {}
}


DEAL_II_NAMESPACE_CLOSE
