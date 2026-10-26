// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// This is the driver program of the testsuite. It is invoked in two
// different roles:
//
// - As a compile and link launcher (via the RULE_LAUNCH_COMPILE and
//   RULE_LAUNCH_LINK target properties) while building the test
//   executables:
//
//     deal_ii_test_driver launch --log-dir DIR --name NAME [--output FILE]
//                                [--link] -- COMMAND [ARGS...]
//
//   The command is executed and its output is forwarded. If the command
//   fails, its output is stored in DIR/NAME.log and the (possibly stale)
//   output FILE is removed. On success the log file is removed again. A
//   successful link step (--link) removes the whole log directory. The
//   exit code of the command is passed through.
//
// - As the test itself (via ctest):
//
//     deal_ii_test_driver run --test NAME [--expect STAGE] [--log-dir DIR]
//                             [--require FILE]... [--build-only]
//                             [--threads N] [--numdiff EXECUTABLE]
//                             [--comparison FILE] [--run-only]
//                             -- COMMAND [ARGS...]
//
//   The test runs through the stages BUILD, RUN, and DIFF within the
//   current working directory and reports the stage it reached:
//     BUILD  - a build log exists in DIR, or one of the required files
//              (the test executable) is missing,
//     RUN    - the command returned with a nonzero exit code,
//     DIFF   - the normalized output did not match the comparison file
//              (or any of its variants),
//     PASSED - all stages were successful.
//   The program returns with exit code 0 if and only if the reached stage
//   matches the expected stage (default PASSED).
//
// - For debugging purposes, the output normalization can be invoked
//   directly:
//
//     deal_ii_test_driver normalize FILE...


#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#else
#  include <fcntl.h>
#  include <spawn.h>
#  include <sys/wait.h>
#  include <unistd.h>
extern char **environ;
#endif

namespace fs = std::filesystem;

namespace
{
  // ---------------------------------------------------------------------------
  // Small helper functions
  // ---------------------------------------------------------------------------

  [[noreturn]] void
  fatal_error(const std::string &message)
  {
    std::cerr << "deal_ii_test_driver: " << message << std::endl;
    std::exit(2);
  }


  std::string
  read_file(const fs::path &path)
  {
    std::ifstream file(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(file),
                       std::istreambuf_iterator<char>());
  }


  void
  write_file(const fs::path &path, const std::string &content)
  {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    file << content;
    if (!file)
      fatal_error("could not write file \"" + path.string() + "\"");
  }


  void
  remove_file(const fs::path &path)
  {
    std::error_code ec;
    fs::remove(path, ec);
  }


  void
  set_environment(const char *name, const std::string &value)
  {
#ifdef _WIN32
    _putenv_s(name, value.c_str());
#else
    setenv(name, value.c_str(), 1);
#endif
  }


  bool
  starts_with(const std::string &string, const std::string &prefix)
  {
    return string.compare(0, prefix.size(), prefix) == 0;
  }


  bool
  ends_with(const std::string &string, const std::string &suffix)
  {
    return string.size() >= suffix.size() &&
           string.compare(string.size() - suffix.size(),
                          suffix.size(),
                          suffix) == 0;
  }


  // Return the argument following the option at position @p i and advance
  // @p i.
  const std::string &
  next_argument(const std::vector<std::string> &arguments, std::size_t &i)
  {
    if (i + 1 >= arguments.size())
      fatal_error("option " + arguments[i] + " requires an argument");
    return arguments[++i];
  }


  bool
  is_digit(const char c)
  {
    return c >= '0' && c <= '9';
  }


  // Whitespace as understood by perl's "\s":
  bool
  is_space(const char c)
  {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' ||
           c == '\v';
  }


  // ---------------------------------------------------------------------------
  // Process creation
  // ---------------------------------------------------------------------------

  // Run the given command and redirect its stdout and stderr into the file
  // @p output_file. Return the exit code of the command. Similarly to a
  // shell, a command that was terminated by a signal returns 128 plus the
  // signal number.
  int
  run_command(const std::vector<std::string> &command,
              const fs::path                 &output_file)
  {
    if (command.empty())
      fatal_error("no command given");

#ifdef _WIN32
    // Construct a command line that is properly split up again by the
    // usual argument parsing of the C runtime:
    std::string command_line;
    for (const auto &argument : command)
      {
        if (!command_line.empty())
          command_line += ' ';
        if (!argument.empty() &&
            argument.find_first_of(" \t\n\v\"") == std::string::npos)
          {
            command_line += argument;
            continue;
          }
        command_line += '"';
        for (auto it = argument.begin();; ++it)
          {
            unsigned int n_backslashes = 0;
            while (it != argument.end() && *it == '\\')
              {
                ++it;
                ++n_backslashes;
              }
            if (it == argument.end())
              {
                command_line.append(n_backslashes * 2, '\\');
                break;
              }
            else if (*it == '"')
              {
                command_line.append(n_backslashes * 2 + 1, '\\');
                command_line += *it;
              }
            else
              {
                command_line.append(n_backslashes, '\\');
                command_line += *it;
              }
          }
        command_line += '"';
      }

    SECURITY_ATTRIBUTES attributes;
    attributes.nLength              = sizeof(attributes);
    attributes.lpSecurityDescriptor = nullptr;
    attributes.bInheritHandle       = TRUE;

    HANDLE output = CreateFileA(output_file.string().c_str(),
                                GENERIC_WRITE,
                                FILE_SHARE_READ | FILE_SHARE_WRITE,
                                &attributes,
                                CREATE_ALWAYS,
                                FILE_ATTRIBUTE_NORMAL,
                                nullptr);
    if (output == INVALID_HANDLE_VALUE)
      fatal_error("could not open file \"" + output_file.string() + "\"");

    STARTUPINFOA startup_info;
    ZeroMemory(&startup_info, sizeof(startup_info));
    startup_info.cb         = sizeof(startup_info);
    startup_info.dwFlags    = STARTF_USESTDHANDLES;
    startup_info.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
    startup_info.hStdOutput = output;
    startup_info.hStdError  = output;

    PROCESS_INFORMATION process_info;
    ZeroMemory(&process_info, sizeof(process_info));

    if (!CreateProcessA(nullptr,
                        command_line.data(),
                        nullptr,
                        nullptr,
                        TRUE,
                        0,
                        nullptr,
                        nullptr,
                        &startup_info,
                        &process_info))
      {
        CloseHandle(output);
        write_file(output_file,
                   "deal_ii_test_driver: could not execute \"" + command[0] +
                     "\"\n");
        return 127;
      }

    CloseHandle(output);
    WaitForSingleObject(process_info.hProcess, INFINITE);
    DWORD exit_code = 1;
    GetExitCodeProcess(process_info.hProcess, &exit_code);
    CloseHandle(process_info.hProcess);
    CloseHandle(process_info.hThread);
    return static_cast<int>(exit_code);

#else

    posix_spawn_file_actions_t file_actions;
    posix_spawn_file_actions_init(&file_actions);
    posix_spawn_file_actions_addopen(&file_actions,
                                     STDOUT_FILENO,
                                     output_file.c_str(),
                                     O_WRONLY | O_CREAT | O_TRUNC,
                                     0644);
    posix_spawn_file_actions_adddup2(&file_actions,
                                     STDOUT_FILENO,
                                     STDERR_FILENO);

    std::vector<char *> argv;
    for (const auto &argument : command)
      argv.push_back(const_cast<char *>(argument.c_str()));
    argv.push_back(nullptr);

    pid_t     pid;
    const int error =
      posix_spawnp(&pid, argv[0], &file_actions, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&file_actions);

    if (error != 0)
      {
        write_file(output_file,
                   "deal_ii_test_driver: could not execute \"" + command[0] +
                     "\": " + std::strerror(error) + "\n");
        return 127;
      }

    int status = 0;
    while (waitpid(pid, &status, 0) == -1)
      if (errno != EINTR)
        fatal_error("waitpid() failed");

    if (WIFEXITED(status))
      return WEXITSTATUS(status);
    else if (WIFSIGNALED(status))
      return 128 + WTERMSIG(status);
    else
      return 1;
#endif
  }


  // ---------------------------------------------------------------------------
  // Output normalization
  // ---------------------------------------------------------------------------

  // Remove everything starting with the first occurrence of @p pattern.
  void
  truncate_at(std::string &line, const char *pattern)
  {
    const auto position = line.find(pattern);
    if (position != std::string::npos)
      line.erase(position);
  }


  // Return true if the line matches "^\s+<pattern>".
  bool
  matches_indented(const std::string &line, const char *pattern)
  {
    std::size_t i = 0;
    while (i < line.size() && is_space(line[i]))
      ++i;
    return i > 0 && line.compare(i, std::strlen(pattern), pattern) == 0;
  }


  // s/line <\d+> of file <.*\//file </
  void
  normalize_exception_file_name(std::string &line)
  {
    const std::string prefix = "line <";
    const std::string infix  = "> of file <";

    const auto last_slash = line.rfind('/');
    if (last_slash == std::string::npos)
      return;

    for (auto start = line.find(prefix); start != std::string::npos;
         start      = line.find(prefix, start + 1))
      {
        auto i = start + prefix.size();
        if (i >= line.size() || !is_digit(line[i]))
          continue;
        while (i < line.size() && is_digit(line[i]))
          ++i;
        if (line.compare(i, infix.size(), infix) != 0)
          continue;
        i += infix.size();
        // The ".*\/" part matches greedily up to the last slash:
        if (last_slash < i)
          return;
        line.replace(start, last_slash + 1 - start, "file <");
        return;
      }
  }


  // s/-0\.(0+)(?!\d)/0.\1/g
  void
  normalize_negative_zero(std::string &line)
  {
    for (auto start = line.find("-0."); start != std::string::npos;
         start      = line.find("-0.", start))
      {
        auto i = start + 3;
        while (i < line.size() && line[i] == '0')
          ++i;
        if (i > start + 3 && (i == line.size() || !is_digit(line[i])))
          {
            line.erase(start, 1);
            start = i - 1;
          }
        else
          start += 1;
      }
  }


  // s/by the deal.II library on.*// (note that the "." is a wildcard)
  void
  normalize_written_by_library(std::string &line)
  {
    const std::string first  = "by the deal";
    const std::string second = "II library on";
    for (auto start = line.find(first); start != std::string::npos;
         start      = line.find(first, start + 1))
      {
        const auto i = start + first.size();
        if (i < line.size() && line[i] != '\n' &&
            line.compare(i + 1, second.size(), second) == 0)
          {
            line.erase(start);
            return;
          }
      }
  }


  // s/written by deal\.II \d+\.\d+\.\d+(-pre|-rc\d*|)/written by deal.II x.y.z/
  void
  normalize_version_string(std::string &line)
  {
    const std::string prefix = "written by deal.II ";
    for (auto start = line.find(prefix); start != std::string::npos;
         start      = line.find(prefix, start + 1))
      {
        auto i = start + prefix.size();

        bool valid = true;
        for (unsigned int component = 0; component < 3; ++component)
          {
            if (component > 0)
              {
                if (i >= line.size() || line[i] != '.')
                  {
                    valid = false;
                    break;
                  }
                ++i;
              }
            if (i >= line.size() || !is_digit(line[i]))
              {
                valid = false;
                break;
              }
            while (i < line.size() && is_digit(line[i]))
              ++i;
          }
        if (!valid)
          continue;

        if (line.compare(i, 4, "-pre") == 0)
          i += 4;
        else if (line.compare(i, 3, "-rc") == 0)
          {
            i += 3;
            while (i < line.size() && is_digit(line[i]))
              ++i;
          }

        line.replace(start, i - start, "written by deal.II x.y.z");
        return;
      }
  }


  // s/^(\s*)(".*":) \{$/\1\2\n\1\{/
  void
  normalize_json_brace(std::string &line)
  {
    std::size_t indentation = 0;
    while (indentation < line.size() && is_space(line[indentation]))
      ++indentation;

    if (line.size() >= indentation + 5 && line[indentation] == '"' &&
        ends_with(line, "\": {"))
      {
        const std::string whitespace = line.substr(0, indentation);
        line.erase(line.size() - 2);
        line += "\n" + whitespace + "{";
      }
  }


  const char *const x11_authorization_warning =
    "Authorization required, but no authorization protocol specified";


  // Remove insignificant, volatile data from the output of tests. This
  // function reimplements the (formerly used) perl script normalize.pl line
  // by line. The perl script ran with "perl -p", i.e., every line was
  // processed with its trailing newline character.
  std::string
  normalize(const std::string &input)
  {
    std::string output;
    output.reserve(input.size());

    std::size_t position = 0;
    while (position < input.size())
      {
        auto       end         = input.find('\n', position);
        const bool has_newline = (end != std::string::npos);
        if (!has_newline)
          end = input.size();

        std::string line = input.substr(position, end - position);
        position         = has_newline ? end + 1 : end;

        // Convert windows to unix line endings:
        if (!line.empty() && line.back() == '\r')
          line.pop_back();

        // Remove JobID, input file name, and date and time strings:
        truncate_at(line, "JobId");
        truncate_at(line, "# Input file name:");
        truncate_at(line, "%%Creation Date:");
        truncate_at(line, "\"created\"");
        truncate_at(line, "# Time =");
        truncate_at(line, "# Date =");
        if (matches_indented(line, "Time ="))
          line.clear();
        if (matches_indented(line, "Date ="))
          line.clear();
        truncate_at(line, "Time tag:");
        normalize_written_by_library(line);

        // Exceptions:
        normalize_exception_file_name(line);

        // Replace -0.0... (not followed by any other digit) by 0.0...:
        normalize_negative_zero(line);

        // Remove deal.II debug output:
        if (has_newline && starts_with(line, "DEAL") &&
            line.find("::_", 4) != std::string::npos)
          continue;

        // Normalize the version string:
        normalize_version_string(line);

        // Different p4est versions output different text in VTU output.
        // Simply delete these lines:
        if (has_newline && line.find("<PCellData Scalars") != std::string::npos)
          continue;
        if (has_newline)
          {
            const auto i = line.find("<PDataArray type");
            if (i != std::string::npos &&
                (line.find("mpirank", i + 16) != std::string::npos ||
                 line.find("level", i + 16) != std::string::npos))
              continue;
          }

        // Different boost versions output the opening bracket for json
        // output on a new line. Thus always transform
        //     "label": {
        // into
        //     "label":
        //     {
        normalize_json_brace(line);

        // Some of the dependency libraries use the X protocol and print
        // warnings if no authorization is available. Remove the warning
        // and the following line:
        if (line.find(x11_authorization_warning) != std::string::npos)
          {
            if (has_newline)
              {
                const auto next = input.find('\n', position);
                position =
                  (next == std::string::npos) ? input.size() : next + 1;
              }
            continue;
          }

        output += line;
        if (has_newline)
          output += '\n';
      }

    return output;
  }


  // ---------------------------------------------------------------------------
  // The "launch" subcommand
  // ---------------------------------------------------------------------------

  std::string
  sanitize(const std::string &name)
  {
    std::string result = name;
    for (auto &c : result)
      if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || is_digit(c) ||
            c == '.' || c == '-' || c == '_'))
        c = '_';
    return result;
  }


  int
  launch(const std::vector<std::string> &arguments)
  {
    fs::path                 log_dir;
    std::string              name;
    fs::path                 output;
    bool                     link = false;
    std::vector<std::string> command;

    for (std::size_t i = 0; i < arguments.size(); ++i)
      {
        const auto &argument = arguments[i];

        if (argument == "--log-dir")
          log_dir = next_argument(arguments, i);
        else if (argument == "--name")
          name = next_argument(arguments, i);
        else if (argument == "--output")
          output = next_argument(arguments, i);
        else if (argument == "--link")
          link = true;
        else if (argument == "--")
          {
            command.assign(arguments.begin() + i + 1, arguments.end());
            break;
          }
        else
          fatal_error("launch: unknown option \"" + argument + "\"");
      }

    if (log_dir.empty() || name.empty())
      fatal_error("launch: --log-dir and --name are mandatory");

    std::error_code ec;
    fs::create_directories(log_dir, ec);

    const fs::path log_file       = log_dir / (sanitize(name) + ".log");
    const fs::path temporary_file = log_dir / (sanitize(name) + ".log.tmp");

    const int exit_code = run_command(command, temporary_file);

    const std::string command_output = read_file(temporary_file);
    remove_file(temporary_file);
    std::cout << command_output << std::flush;

    if (exit_code == 0)
      {
        if (link)
          fs::remove_all(log_dir, ec);
        else
          remove_file(log_file);
        return 0;
      }

    std::string command_line;
    for (const auto &argument : command)
      command_line += (command_line.empty() ? "" : " ") + argument;

    write_file(log_file,
               "Command: " + command_line + "\n" + "Return code " +
                 std::to_string(exit_code) + "\n" + command_output);

    // Make sure that we do not accidentally pick up a stale output from
    // a previous build:
    if (!output.empty())
      remove_file(output);

    return exit_code;
  }


  // ---------------------------------------------------------------------------
  // The "run" subcommand
  // ---------------------------------------------------------------------------

  enum class Stage
  {
    build,
    run,
    diff,
    passed
  };


  const char *
  stage_name(const Stage stage)
  {
    switch (stage)
      {
        case Stage::build:
          return "BUILD";
        case Stage::run:
          return "RUN";
        case Stage::diff:
          return "DIFF";
        default:
          return "PASSED";
      }
  }


  Stage
  parse_stage(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), [](char c) {
      return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    });
    for (const auto stage :
         {Stage::build, Stage::run, Stage::diff, Stage::passed})
      if (name == stage_name(stage))
        return stage;
    fatal_error("unknown stage \"" + name + "\"");
  }


  // Remove all files in the current directory whose name starts with
  // @p prefix.
  void
  remove_files_with_prefix(const std::string &prefix)
  {
    std::error_code ec;
    for (const auto &entry : fs::directory_iterator(".", ec))
      if (starts_with(entry.path().filename().string(), prefix))
        remove_file(entry.path());
  }


  struct RunOptions
  {
    std::string              test;
    Stage                    expect = Stage::passed;
    fs::path                 log_dir;
    std::vector<fs::path>    required_files;
    bool                     build_only = false;
    std::string              n_threads  = "0";
    std::string              numdiff;
    fs::path                 comparison_file;
    bool                     run_only = false;
    std::vector<std::string> command;
  };


  // The BUILD stage: Check that the build of all required files succeeded.
  bool
  check_build(const RunOptions &options, std::ostream &out)
  {
    const auto &test = options.test;

    std::error_code       ec;
    std::vector<fs::path> log_files;
    if (!options.log_dir.empty())
      for (const auto &entry : fs::directory_iterator(options.log_dir, ec))
        if (ends_with(entry.path().filename().string(), ".log"))
          log_files.push_back(entry.path());
    std::sort(log_files.begin(), log_files.end());

    if (!log_files.empty())
      {
        out << test << ": BUILD failed. ------ Build log:\n";
        for (const auto &log_file : log_files)
          out << read_file(log_file) << "\n";
        return false;
      }

    for (const auto &file : options.required_files)
      if (!fs::exists(file, ec))
        {
          out << test << ": BUILD failed. ------ File " << file.string()
              << " does not exist. The test executables have to be built"
              << " before running the tests, for example by building the"
              << " \"build_tests\" target.\n";
          return false;
        }

    out << test << ": BUILD successful.\n";
    return true;
  }


  // The RUN stage: Run the command and store the normalized output in the
  // file "output".
  bool
  run_test(const RunOptions &options, std::ostream &out)
  {
    const auto &test = options.test;

    //
    // If n_threads is not equal to zero then:
    //  - Export the environment variable DEAL_II_NUM_THREADS. This will
    //    enforce an upper bound of DEAL_II_NUM_THREADS during thread
    //    initialization of the threading pool in deal.II.
    //  - Export TEST_N_THREADS which is internally used in the deal.II
    //    testsuite to explicitly set the number of threads in the header
    //    file tests.h.
    //
    if (options.n_threads != "0")
      {
        set_environment("DEAL_II_NUM_THREADS", options.n_threads);
        set_environment("TEST_N_THREADS", options.n_threads);
      }

    //
    // Limit the OpenMP pool to two threads. Set both variables in case the
    // caller has either one set. OMP_NUM_THREADS permits nesting, while
    // OMP_THREAD_LIMIT limits the total number of threads. In addition we
    // set the OMP_PROC_BIND variable to false to allow free movement of
    // the two worker threads and silence a Kokkos warning.
    //
    set_environment("OMP_NUM_THREADS", "2");
    set_environment("OMP_THREAD_LIMIT", "2");
    set_environment("OMP_PROC_BIND", "false");

    //
    // OpenMPI parameters:
    //  - Allow to oversubscribe the system, meaning that we can issue tests
    //    with more mpi ranks than available cores.
    //  - Ensure that we do not bind mpi tests to specific
    //    cores/processors/sockets. Otherwise we run the risk that multiple
    //    mpi tests are all pinned to the same processor core.
    //
    //    If the test runs exclusively, however, do not override MPI binding
    //    policies. This is important to ensure that performance tests are
    //    scheduled properly.
    //
    // (The first variable of each pair is for OpenMPI 4 and older, the
    // second one for OpenMPI 5.)
    //
    set_environment("OMPI_MCA_rmaps_base_oversubscribe", "1");
    set_environment("PRTE_MCA_rmaps_default_mapping_policy", ":oversubscribe");
    if (std::getenv("TEST_IS_EXCLUSIVE") == nullptr)
      {
        set_environment("OMPI_MCA_hwloc_base_binding_policy", "none");
        set_environment("PRTE_MCA_hwloc_default_binding_policy", "none");
      }

    //
    // Disable Kokkos runtime warnings so that we can oversubscribe threads
    // without Kokkos complaining.
    //
    set_environment("KOKKOS_DISABLE_WARNINGS", "1");

    remove_file("failing_output");
    remove_file("output");
    remove_file("stdout");

    const int exit_code = run_command(options.command, "stdout");

    // Tests either write to the file "output" directly, or to stdout:
    std::error_code ec;
    if (!fs::exists("output", ec))
      fs::rename("stdout", "output", ec);

    if (exit_code != 0)
      {
        fs::rename("output", "failing_output", ec);
        const auto cwd = fs::current_path(ec).string();
        out << test << ": RUN failed. ------ Return code " << exit_code << "\n";
        out << test << ": RUN failed. ------ Result: " << cwd
            << "/failing_output\n";
        out << test << ": RUN failed. ------ Partial output:\n";
        out << read_file("failing_output");
        if (fs::exists("stdout", ec))
          {
            out << "\n";
            out << test
                << ": RUN failed. ------ Additional output on stdout/stderr:\n";
            out << "\n";
            out << read_file("stdout");
          }
        return false;
      }

    write_file("output", normalize(read_file("output")));

    out << test << ": RUN successful.\n";
    return true;
  }


  // The DIFF stage: Compare the file "output" against the comparison file
  // and all its variants.
  bool
  diff_test(const RunOptions &options, std::ostream &out)
  {
    const auto &test = options.test;

    remove_files_with_prefix("failing_diff");
    remove_files_with_prefix("diff");

    //
    // Pick up the main comparison file and all variants. A valid variant
    // name is of the form [...].output.[STRING]. Sort the list the same
    // way as the C locale would do:
    //
    const auto        directory = options.comparison_file.parent_path();
    const std::string base_name = options.comparison_file.filename().string();

    std::vector<std::string> candidates;
    std::error_code          ec;
    for (const auto &entry :
         fs::directory_iterator(directory.empty() ? "." : directory, ec))
      {
        const auto file_name = entry.path().filename().string();
        if (starts_with(file_name, base_name))
          candidates.push_back(file_name);
      }
    std::sort(candidates.begin(), candidates.end());

    for (const auto &candidate : candidates)
      {
        const std::string variant   = candidate.substr(base_name.size());
        const fs::path    diff_file = "diff" + variant;

        //
        // Configure numdiff with
        //   - absolute precision set to 1e-6, everything below is regarded
        //     as being equal to 0
        //   - relative differences of 1e-8
        //   - [space][tab][newline]=,:;<>[](){}^ as separators between
        //     numbers
        //
        const int exit_code = run_command({options.numdiff,
                                           "-a",
                                           "1e-6",
                                           "-r",
                                           "1e-8",
                                           "-s",
                                           " \\t\\r\\n=,:;<>[](){}^",
                                           (directory / candidate).string(),
                                           "output"},
                                          diff_file);

        if (exit_code == 0)
          {
            //
            // Ensure that only a single diff file remains. In case of a
            // successful comparison against a variant, store the fact
            // that we compared against a variant in the diff file:
            //
            remove_files_with_prefix("diff");
            write_file("diff",
                       variant.empty() ?
                         std::string() :
                         test + ": DIFF successful. - Variant: " +
                           (directory / candidate).string() + "\n");
            out << test << ": DIFF successful.\n";
            return true;
          }
      }

    //
    // None of the diffs succeeded, report the diff against the main
    // comparison file:
    //
    for (const auto &entry : fs::directory_iterator(".", ec))
      {
        const auto file_name = entry.path().filename().string();
        if (starts_with(file_name, "diff"))
          fs::rename(entry.path(), "failing_" + file_name, ec);
      }

    const auto cwd = fs::current_path(ec).string();
    out << test
        << ": DIFF failed. ------ Source: " << options.comparison_file.string()
        << "\n";
    out << test << ": DIFF failed. ------ Result: " << cwd << "/output\n";
    out << "Check " << cwd << "/output " << options.comparison_file.string()
        << "\n";
    out << test << ": DIFF failed. ------ Diff:   " << cwd << "/failing_diff\n";
    out << test << ": DIFF failed. ------ First 20 lines of numdiff output:\n";

    std::ifstream failing_diff("failing_diff");
    std::string   line;
    for (unsigned int i = 0; i < 20 && std::getline(failing_diff, line); ++i)
      out << line << "\n";

    return false;
  }


  int
  run(const std::vector<std::string> &arguments)
  {
    RunOptions options;

    for (std::size_t i = 0; i < arguments.size(); ++i)
      {
        const auto &argument = arguments[i];

        if (argument == "--test")
          options.test = next_argument(arguments, i);
        else if (argument == "--expect")
          options.expect = parse_stage(next_argument(arguments, i));
        else if (argument == "--log-dir")
          options.log_dir = next_argument(arguments, i);
        else if (argument == "--require")
          options.required_files.emplace_back(next_argument(arguments, i));
        else if (argument == "--build-only")
          options.build_only = true;
        else if (argument == "--threads")
          options.n_threads = next_argument(arguments, i);
        else if (argument == "--numdiff")
          options.numdiff = next_argument(arguments, i);
        else if (argument == "--comparison")
          options.comparison_file = next_argument(arguments, i);
        else if (argument == "--run-only")
          options.run_only = true;
        else if (argument == "--")
          {
            options.command.assign(arguments.begin() + i + 1, arguments.end());
            break;
          }
        else
          fatal_error("run: unknown option \"" + argument + "\"");
      }

    if (options.test.empty())
      fatal_error("run: --test is mandatory");
    if (!options.build_only && options.command.empty())
      fatal_error("run: no command given");
    if (!options.build_only && !options.run_only &&
        (options.numdiff.empty() || options.comparison_file.empty()))
      fatal_error("run: --numdiff and --comparison are mandatory unless"
                  " --run-only or --build-only is given");

    //
    // Run through all stages and record the stage we reached:
    //

    std::ostringstream log;
    Stage              stage = Stage::build;

    if (check_build(options, log))
      {
        stage = Stage::run;
        if (options.build_only)
          {
            log << options.test << ": RUN skipped.\n";
            log << options.test << ": DIFF skipped.\n";
            stage = Stage::passed;
          }
        else if (run_test(options, log))
          {
            stage = Stage::diff;
            if (options.run_only)
              {
                log << options.test << ": DIFF skipped.\n";
                stage = Stage::passed;
              }
            else if (diff_test(options, log))
              stage = Stage::passed;
          }
      }

    //
    // Print out the test result:
    //

    const std::string test = options.test;
    std::cout << "Test " << test << ": " << stage_name(stage) << "\n";
    std::cout << "===============================   OUTPUT BEGIN  "
                 "===============================\n";

    if (stage == Stage::passed)
      {
        // The diff file is either empty or contains a string explaining
        // which comparison file variant succeeded:
        std::error_code ec;
        if (!options.build_only && !options.run_only && fs::exists("diff", ec))
          std::cout << read_file("diff");
        std::cout << test << ": PASSED.\n";
      }
    else
      {
        std::cout << log.str() << "\n";
        std::cout << test << ": ******    " << stage_name(stage)
                  << " failed    *******\n\n";
      }

    std::cout << "===============================    OUTPUT END   "
                 "===============================\n";

    if (stage != options.expect)
      {
        std::cout << "Expected stage " << stage_name(options.expect)
                  << " - aborting" << std::endl;
        return 1;
      }
    else if (stage != Stage::passed)
      std::cout << "Expected stage " << stage_name(options.expect)
                << " - test considered successful." << std::endl;

    return 0;
  }


  void
  print_usage()
  {
    std::cerr
      << "Usage:\n"
      << "  deal_ii_test_driver launch --log-dir DIR --name NAME\n"
      << "                             [--output FILE] [--link]\n"
      << "                             -- COMMAND [ARGS...]\n"
      << "  deal_ii_test_driver run --test NAME [--expect STAGE]\n"
      << "                          [--log-dir DIR] [--require FILE]...\n"
      << "                          [--build-only] [--threads N]\n"
      << "                          [--numdiff EXECUTABLE]\n"
      << "                          [--comparison FILE] [--run-only]\n"
      << "                          -- COMMAND [ARGS...]\n"
      << "  deal_ii_test_driver normalize FILE...\n";
  }
} // namespace



int
main(int argc, char **argv)
{
  if (argc < 2)
    {
      print_usage();
      return 2;
    }

  const std::string              subcommand = argv[1];
  const std::vector<std::string> arguments(argv + 2, argv + argc);

  if (subcommand == "launch")
    return launch(arguments);
  else if (subcommand == "run")
    return run(arguments);
  else if (subcommand == "normalize")
    {
      for (const auto &file : arguments)
        write_file(file, normalize(read_file(file)));
      return 0;
    }

  print_usage();
  return 2;
}
