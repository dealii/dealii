// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_job_identifier_h
#define dealii_job_identifier_h


#include <deal.II/base/config.h>

#include <string>

DEAL_II_NAMESPACE_OPEN
/**
 * Identification of a program run. <tt>JobIdentifier</tt> determines the
 * start time of a program run and stores it as a program identifier. There
 * exists a library object <tt>dealjobid</tt> of this class. This object can
 * be accessed by all output functions to provide an id for the current job.
 *
 * @ingroup utilities
 */
class JobIdentifier
{
public:
  /**
   * Constructor. Set program identifier to value of <tt>program_id</tt>
   * concatenated with the present time.
   */
  JobIdentifier();

  /**
   * This function returns an identifier for the running program. Currently,
   * the library provides a function returning "JobID".
   *
   * The user may define a replacement of this function in their source code and
   * avoid linking the library version. Unfortunately, this mechanism does not
   * work with shared libraries.
   */
  static const char *
  program_id();

  /**
   * Obtain the base name of the filename passed as argument. That is,
   * if the file is <tt>mypath/file.cc</tt> return just
   * <tt>file</tt>. For example, this function can be called from a
   * user program with argument <tt>__FILE__</tt> to create an
   * identifier for the program being run.
   */
  static std::string
  base_name(const std::string &filename);

  /**
   * Return the value of <tt>id</tt>.
   */
  std::string
  operator()() const;

  /**
   * %Function to identify the presently running program.
   */
  static const JobIdentifier &
  get_dealjobid();

private:
  /**
   * String holding the identifier of the presently running program.
   */
  std::string id;
};

DEAL_II_NAMESPACE_CLOSE

#endif
