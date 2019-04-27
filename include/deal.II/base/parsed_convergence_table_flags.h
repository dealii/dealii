// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
// -------------------------------------------------------------------

#ifndef dealii_base_parsed_convergence_table_flags_h
#define dealii_base_parsed_convergence_table_flags_h

#include <deal.II/base/config.h>

#include <deal.II/base/patterns.h>

DEAL_II_NAMESPACE_OPEN

namespace ParsedConvergenceTableFlags
{
  /**
   * Implemented Error norms.
   */
  enum NormType
  {
    none    = 0x00,
    Linfty  = 0x01,
    L2      = 0x02,
    W1infty = 0x04,
    H1      = 0x08,
    custom  = 0x10,
  };

  /**
   * The ExtraColumns enum.
   */
  enum ExtraColumns
  {
    no_columns = 0x0,
    dofs       = 0x1,
    cells      = 0x2,
    dt         = 0x4
  };
} // namespace ParsedConvergenceTableFlags

namespace Patterns
{
  namespace Tools
  {
    template <>
    struct Convert<ParsedConvergenceTableFlags::NormType>
    {
      /**
       * Return the Correct pattern for NormType.
       */
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return std_cxx14::make_unique<Patterns::Selection>(
          "none|Linfty|L2|W1infty|H1|custom");
      }



      /**
       * Convert a NormType to a string.
       */
      static std::string
      to_string(const ParsedConvergenceTableFlags::NormType & s,
                const std::unique_ptr<Patterns::PatternBase> &p =
                  Convert<ParsedConvergenceTableFlags::NormType>::to_pattern())
      {
        (void)p;
        if (s == ParsedConvergenceTableFlags::Linfty)
          return "Linfty";
        else if (s == ParsedConvergenceTableFlags::L2)
          return "L2";
        else if (s == ParsedConvergenceTableFlags::W1infty)
          return "W1infty";
        else if (s == ParsedConvergenceTableFlags::H1)
          return "H1";
        else if (s == ParsedConvergenceTableFlags::custom)
          return "custom";
        else if (s == ParsedConvergenceTableFlags::none)
          return "none";
        else
          {
            AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
            return "";
          }
      }



      /**
       * Convert a string to a NormType.
       */
      static ParsedConvergenceTableFlags::NormType
      to_value(const std::string &                           norm_string,
               const std::unique_ptr<Patterns::PatternBase> &p =
                 Convert<ParsedConvergenceTableFlags::NormType>::to_pattern())
      {
        (void)p;
        ParsedConvergenceTableFlags::NormType norm =
          ParsedConvergenceTableFlags::none;

        if (norm_string == "Linfty")
          {
            norm = ParsedConvergenceTableFlags::Linfty;
          }
        else if (norm_string == "L2")
          {
            norm = ParsedConvergenceTableFlags::L2;
          }
        else if (norm_string == "W1infty")
          {
            norm = ParsedConvergenceTableFlags::W1infty;
          }
        else if (norm_string == "H1")
          {
            norm = ParsedConvergenceTableFlags::H1;
          }
        else if (norm_string == "custom")
          {
            norm = ParsedConvergenceTableFlags::custom;
          }
        else if (norm_string == "none")
          {
            norm = ParsedConvergenceTableFlags::none;
          }
        else
          {
            AssertThrow(false, ExcMessage("The norm type <" + norm_string + "> is not recognized."));
          }

        return norm;
      }
    };



    template <>
    struct Convert<ParsedConvergenceTableFlags::ExtraColumns>
    {
      /**
       * Return the Correct pattern for ExtraColumns.
       */
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return std_cxx14::make_unique<Patterns::Selection>("none|cells|dofs|dt");
      }



      /**
       * Convert an ExtraColumn object to a string.
       */
      static std::string
      to_string(
        const ParsedConvergenceTableFlags::ExtraColumns &s,
        const std::unique_ptr<Patterns::PatternBase> &   p =
          Convert<ParsedConvergenceTableFlags::ExtraColumns>::to_pattern())
      {
        (void)p;
        if (s == ParsedConvergenceTableFlags::cells)
          return "cells";
        else if (s == ParsedConvergenceTableFlags::dofs)
          return "dofs";
        else if (s == ParsedConvergenceTableFlags::dt)
          return "dt";
        else if (s == ParsedConvergenceTableFlags::no_columns)
          return "none";
        else
          {
            AssertThrow(false, ExcMessage("Didn't recognize a column type."));
            return "";
          }
      }



      /**
       * Convert a string to an ExtraColumn object.
       */
      static ParsedConvergenceTableFlags::ExtraColumns
      to_value(
        const std::string &                           column_string,
        const std::unique_ptr<Patterns::PatternBase> &p =
          Convert<ParsedConvergenceTableFlags::ExtraColumns>::to_pattern())
      {
        (void)p;
        ParsedConvergenceTableFlags::ExtraColumns column =
          ParsedConvergenceTableFlags::no_columns;

        if (column_string == "cells")
          {
            column = ParsedConvergenceTableFlags::cells;
          }
        else if (column_string == "dofs")
          {
            column = ParsedConvergenceTableFlags::dofs;
          }
        else if (column_string == "dt")
          {
            column = ParsedConvergenceTableFlags::dt;
          }
        else if (column_string == "none")
          {
            column = ParsedConvergenceTableFlags::no_columns;
          }
        else
          {
            AssertThrow(false,
                        ExcMessage("The column type <" + column_string + "> is not recognized."));
          }

        return column;
      }
    };
  } // namespace Tools
} // namespace Patterns

DEAL_II_NAMESPACE_CLOSE

#endif
