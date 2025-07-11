// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <Kokkos_Core.hpp>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

#ifdef DEAL_II_WITH_MPI
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <mpi.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif

#ifdef DEAL_II_TRILINOS_WITH_SEACAS
#  include <exodusII.h>
#endif

#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
#  include <execinfo.h>
#endif

#ifdef DEAL_II_HAVE_LIBSTDCXX_DEMANGLER
#  include <cxxabi.h>
#endif

DEAL_II_NAMESPACE_OPEN


namespace deal_II_exceptions
{
  namespace internals
  {
    std::string &
    get_additional_assert_output()
    {
      static std::string additional_assert_output;
      return additional_assert_output;
    }

    bool show_stacktrace = true;

    bool allow_abort_on_exception = true;
  } // namespace internals



  void
  set_additional_assert_output(const std::string &p)
  {
    internals::get_additional_assert_output() = p;
  }



  void
  suppress_stacktrace_in_exceptions()
  {
    internals::show_stacktrace = false;
  }



  void
  disable_abort_on_exception()
  {
    internals::allow_abort_on_exception = false;
  }



  void
  enable_abort_on_exception()
  {
    internals::allow_abort_on_exception = true;
  }
} // namespace deal_II_exceptions



ExceptionBase::ExceptionBase()
  : file("")
  , line(0)
  , function("")
  , cond("")
  , exc("")
  , n_stacktrace_frames(0)
  , what_str("")
{
#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  std::fill(std::begin(raw_stacktrace), std::end(raw_stacktrace), nullptr);
#endif
}



ExceptionBase::ExceptionBase(const ExceptionBase &exc)
  : file(exc.file)
  , line(exc.line)
  , function(exc.function)
  , cond(exc.cond)
  , exc(exc.exc)
  , n_stacktrace_frames(exc.n_stacktrace_frames)
  , what_str("") // don't copy the error message, it gets generated dynamically
                 // by what()
{
#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  // Copy the raw_stacktrace pointers. We don't own them, they just point to the
  // addresses of symbols in the executable's/library's symbol tables -- and as
  // a consequence, it is safe to copy these pointers
  std::copy(std::begin(exc.raw_stacktrace),
            std::end(exc.raw_stacktrace),
            std::begin(raw_stacktrace));
#endif
}



void
ExceptionBase::set_fields(const char *f,
                          const int   l,
                          const char *func,
                          const char *c,
                          const char *e)
{
  file     = f;
  line     = l;
  function = func;
  cond     = c;
  exc      = e;

  // If the system supports this, get a stacktrace how we got here:
  // Note that we defer the symbol lookup done by backtrace_symbols()
  // to when we need it (see what() below). This is for performance
  // reasons, as this requires loading libraries and can take in the
  // order of seconds on some machines.
#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  n_stacktrace_frames = backtrace(raw_stacktrace, 25);
#endif
}



const char *
ExceptionBase::what() const noexcept
{
  // If no error c_string was generated so far, do it now:
  if (what_str.empty())
    generate_message();

  return what_str.c_str();
}



const char *
ExceptionBase::get_exc_name() const
{
  return exc;
}



void
ExceptionBase::print_exc_data(std::ostream &out) const
{
  // print a header for the exception
  out << "An error occurred in line <" << line << "> of file <" << file
      << "> in function" << std::endl
      << "    " << function << std::endl;

  // If the exception stores a string representation of the violated
  // condition, then output it. Not all exceptions do (e.g., when
  // creating an exception inside DEAL_II_NOT_IMPLEMENTED();), so
  // we have to check whether there is anything to print.
  //
  // There are also places where the condition is not very interesting.
  // Specifically, this is the case for places such as
  //   Assert (false, ExcInternalError());
  // Here, the condition is simply 'false'. This is not worth printing,
  // so suppress this case.
  if ((cond != nullptr) && (std::strcmp(cond, "false") != 0))
    out << "The violated condition was: " << std::endl
        << "    " << cond << std::endl;

  // If a string representation of the exception itself is available,
  // consider printing it as well. This is useful if the names of
  // local variables appear in the generation of the error message,
  // because it allows the identification of parts of the error text
  // with what variables may have cause this.
  //
  // On the other hand, this is almost never the case for ExcMessage
  // exceptions which would simply print the same text twice: once for
  // the way the message was composed, and once for the additional
  // information. Furthermore, the former of these two is often spread
  // between numerous "..."-enclosed strings that the preprocessor
  // collates into a single string, making it awkward to read. Consequently,
  // elide this text if the message was generated via an ExcMessage object.
  //
  // There are cases where the exception generation mechanism suppresses
  // the string representation of the exception because it does not add
  // anything -- e.g., DEAL_II_NOT_IMPLEMENTED does this. In those cases,
  // also suppress the output.
  if ((exc != nullptr) &&
      ((cond == nullptr) ||
       (std::strstr(cond, "dealii::ExcMessage") != nullptr)))
    out << "The name and call sequence of the exception was:" << std::endl
        << "    " << exc << std::endl;

  // finally print the additional information the exception provides:
  out << "Additional information: " << std::endl;
}



void
ExceptionBase::print_info(std::ostream &out) const
{
  out << "    (none)" << std::endl;
}



void
ExceptionBase::print_stack_trace(std::ostream &out) const
{
  if (n_stacktrace_frames == 0)
    return;

  if (deal_II_exceptions::internals::show_stacktrace == false)
    return;

  char **stacktrace = nullptr;
#ifdef DEAL_II_HAVE_GLIBC_STACKTRACE
  // We have deferred the symbol lookup to this point to avoid costly
  // runtime penalties due to linkage of external libraries by
  // backtrace_symbols.
  stacktrace = backtrace_symbols(raw_stacktrace, n_stacktrace_frames);
#endif


  // if there is a stackframe stored, print it
  out << std::endl;
  out << "Stacktrace:" << std::endl << "-----------" << std::endl;

  // print the stacktrace. first omit all those frames that have
  // ExceptionBase or deal_II_exceptions in their names, as these
  // correspond to the exception raising mechanism themselves, rather than
  // the place where the exception was triggered
  int frame = 0;
  while ((frame < n_stacktrace_frames) &&
         ((std::string(stacktrace[frame]).find("ExceptionBase") !=
           std::string::npos) ||
          (std::string(stacktrace[frame]).find("deal_II_exceptions") !=
           std::string::npos)))
    ++frame;

  // output the rest
  const unsigned int first_significant_frame = frame;
  for (; frame < n_stacktrace_frames; ++frame)
    {
      out << '#' << frame - first_significant_frame << "  ";

      // the stacktrace frame is actually of the format
      // "filename(functionname+offset) [address]". let's try to get the
      // mangled functionname out:
      std::string        stacktrace_entry(stacktrace[frame]);
      const unsigned int pos_start = stacktrace_entry.find('('),
                         pos_end   = stacktrace_entry.find('+');
      std::string functionname =
        stacktrace_entry.substr(pos_start + 1, pos_end - pos_start - 1);

      std::string demangled_stacktrace_entry =
        stacktrace_entry.substr(0, pos_start);
      demangled_stacktrace_entry += ": ";

      // demangle, and if successful replace old mangled string by
      // unmangled one (skipping address and offset). treat "main"
      // differently, since it is apparently demangled as "unsigned int"
      // for unknown reasons :-) if we can, demangle the function name
#ifdef DEAL_II_HAVE_LIBSTDCXX_DEMANGLER
      int   status;
      char *p =
        abi::__cxa_demangle(functionname.c_str(), nullptr, nullptr, &status);

      if ((status == 0) && (functionname != "main"))
        {
          std::string realname(p);
          // in MT mode, one often gets backtraces spanning several lines
          // because we have so many boost::tuple arguments in the MT
          // calling functions. most of the trailing arguments of these
          // tuples are actually unused boost::tuples::null_type, so we
          // should split them off if they are trailing a template argument
          // list
          while (realname.find(", boost::tuples::null_type>") !=
                 std::string::npos)
            realname.erase(realname.find(", boost::tuples::null_type>"),
                           std::string(", boost::tuples::null_type").size());

          demangled_stacktrace_entry += realname;
        }
      else
        demangled_stacktrace_entry += functionname;

      std::free(p);

#else

      demangled_stacktrace_entry += functionname;
#endif

      // then output what we have
      out << demangled_stacktrace_entry << std::endl;

      // stop if we're in main()
      if (functionname == "main")
        break;
    }

  std::free(stacktrace); // free(nullptr) is allowed
  stacktrace = nullptr;
}



void
ExceptionBase::generate_message() const
{
  // build up a c_string with the error message.
  // Guard this procedure with a try block, we shall not throw at this
  // place...
  try
    {
      std::ostringstream converter;

      converter << std::endl
                << "--------------------------------------------------------"
                << std::endl;

      // Print out general data
      print_exc_data(converter);

      // Print out exception specific data. Put this into another stringstream
      // object for now so that we can break long lines and print them in a
      // more easily read way
      {
        std::ostringstream message;
        print_info(message);

        const auto message_in_lines =
          Utilities::break_text_into_lines(message.str(), 70);

        // Put the message into the stream that will be output.
        for (const auto &line : message_in_lines)
          converter << "    " << line << '\n';
      }


      print_stack_trace(converter);

      if (!deal_II_exceptions::internals::get_additional_assert_output()
             .empty())
        {
          converter
            << "--------------------------------------------------------"
            << std::endl
            << deal_II_exceptions::internals::get_additional_assert_output()
            << std::endl;
        }

      converter << "--------------------------------------------------------"
                << std::endl;

      what_str = converter.str();
    }
  catch (...)
    {
      // On error, resume next. There is nothing better we can do...
      what_str = "ExceptionBase::generate_message () failed";
    }
}



namespace StandardExceptions
{
#ifdef DEAL_II_WITH_MPI
  ExcMPI::ExcMPI(const int error_code)
    : error_code(error_code)
  {}

  void
  ExcMPI::print_info(std::ostream &out) const
  {
    char error_name[MPI_MAX_ERROR_STRING];
    error_name[0]        = '\0';
    int resulting_length = MPI_MAX_ERROR_STRING;

    bool error_name_known = false;
    // workaround for Open MPI 1.6.5 not accepting
    // MPI_ERR_LASTCODE in MPI_Error_class
    if (error_code < MPI_ERR_LASTCODE)
      {
        // get the string name of the error code by first converting it to an
        // error class.
        int error_class  = 0;
        int ierr         = MPI_Error_class(error_code, &error_class);
        error_name_known = (ierr == MPI_SUCCESS);

        // Check the output of the error printing functions. If either MPI
        // function fails we should just print a less descriptive message.
        if (error_name_known)
          {
            ierr = MPI_Error_string(error_class, error_name, &resulting_length);
            error_name_known = error_name_known && (ierr == MPI_SUCCESS);
          }
      }

    out << "deal.II encountered an error while calling an MPI function."
        << std::endl;
    if (error_name_known)
      {
        out << "The description of the error provided by MPI is \""
            << error_name << "\"." << std::endl;
      }
    else
      {
        out
          << "This error code is not equal to any of the standard MPI error codes."
          << std::endl;
      }
    out << "The numerical value of the original error code is " << error_code
        << '.' << std::endl;
  }
#endif // DEAL_II_WITH_MPI



#ifdef DEAL_II_TRILINOS_WITH_SEACAS
  ExcExodusII::ExcExodusII(const int error_code)
    : error_code(error_code)
  {
    // To avoid including yet another header in exceptions.h we assume that
    // EX_NOERR is zero. Check that here:
    static_assert(EX_NOERR == 0,
                  "EX_NOERR is assumed to be zero in all versions of ExodusII");
  }



  void
  ExcExodusII::print_info(std::ostream &out) const
  {
    out << "Error code is " << error_code << '\n';
    out << "String description: " << ex_strerror(error_code) << '\n';
  }
#endif // DEAL_II_TRILINOS_WITH_SEACAS
} // namespace StandardExceptions

namespace deal_II_exceptions
{
  namespace internals
  {
    [[noreturn]] void
    abort(const ExceptionBase &exc) noexcept
    {
      // first print the error
      std::cerr << exc.what() << std::endl;

      // then bail out. if in MPI mode, bring down the entire
      // house by asking the MPI system to do a best-effort
      // operation at also terminating all of the other MPI
      // processes. this is useful because if only one process
      // runs into an assertion, then that may lead to deadlocks
      // if the others don't recognize this, or at the very least
      // delay their termination until they realize that their
      // communication with the job that died times out.
      //
      // Unlike std::abort(), MPI_Abort() unfortunately doesn't break when
      // running inside a debugger like GDB, so only use this strategy if
      // absolutely necessary and inform the user how to use a debugger.
#ifdef DEAL_II_WITH_MPI
      int is_initialized;
      MPI_Initialized(&is_initialized);
      if (is_initialized != 0)
        {
          // do the same as in Utilities::MPI::n_mpi_processes() here,
          // but without error checking to not throw again.
          const int n_proc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
          if (n_proc > 1)
            {
              std::cerr
                << "Calling MPI_Abort now.\n"
                << "To break execution in a GDB session, execute 'break MPI_Abort' before "
                << "running. You can also put the following into your ~/.gdbinit:\n"
                << "  set breakpoint pending on\n"
                << "  break MPI_Abort\n"
                << "  set breakpoint pending auto" << std::endl;

              MPI_Abort(MPI_COMM_WORLD,
                        /* return code = */ 255);
            }
        }
#endif

      // Let's abort the program here. On the host, we need to call
      // std::abort, on devices we need to do something different.
      // Kokkos::abort() does the right thing in all circumstances.
      Kokkos::abort(
        "Abort() was called during dealing with an assertion or exception.");
    }



    void
    do_issue_error_nothrow(const ExceptionBase &exc) noexcept
    {
      if (deal_II_exceptions::internals::allow_abort_on_exception)
        abort(exc);
      else
        {
          // We are not allowed to throw, and not allowed to abort.
          // Just print the exception name to deallog and continue normally:
          deallog << "Exception: " << exc.get_exc_name() << std::endl;
          deallog << exc.what() << std::endl;
        }
    }

  } /*namespace internals*/
} /*namespace deal_II_exceptions*/



DEAL_II_NAMESPACE_CLOSE
