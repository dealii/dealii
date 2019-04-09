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
// ---------------------------------------------------------------------

#include <deal.II/base/parsed_convergence_table.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  std::vector<std::string>
  get_unique_component_names(const std::vector<std::string> &component_names)
  {
    auto elements = component_names;
    auto last     = std::unique(elements.begin(), elements.end());
    elements.resize(last - elements.begin());
    return elements;
  }



  std::vector<ComponentMask>
  get_unique_component_masks(const std::vector<std::string> &component_names)
  {
    const auto unique_component_names =
      get_unique_component_names(component_names);

    std::vector<ComponentMask> masks;

    std::vector<std::vector<bool>> bools(
      unique_component_names.size(),
      std::vector<bool>(component_names.size(), false));

    unsigned int j = 0;
    for (unsigned int i = 0; i < component_names.size(); ++i)
      {
        if (unique_component_names[j] != component_names[i])
          masks.push_back(ComponentMask(bools[j++]));
        bools[j][i] = true;
      }
    masks.push_back(ComponentMask(bools[j++]));
    AssertDimension(j, unique_component_names.size());
    return masks;
  }
} // namespace

ParsedConvergenceTable::ParsedConvergenceTable(
  const std::vector<std::string> &solution_names,
  const std::vector<std::set<ParsedConvergenceTableFlags::NormType>>
    &list_of_error_norms)
  : ParsedConvergenceTable(solution_names,
                           list_of_error_norms,
                           {ParsedConvergenceTableFlags::cells,
                            ParsedConvergenceTableFlags::dofs},
                           ParsedConvergenceTableFlags::dofs,
                           "",
                           3,
                           true,
                           true)
{}



ParsedConvergenceTable::ParsedConvergenceTable(
  const std::vector<std::string> &component_names,
  const std::vector<std::set<ParsedConvergenceTableFlags::NormType>>
    &list_of_error_norms,
  const std::set<ParsedConvergenceTableFlags::ExtraColumns> &extra_columns,
  const ParsedConvergenceTableFlags::ExtraColumns &          rate_key,
  const std::string &                                        error_file_name,
  const unsigned int &                                       precision,
  const bool &                                               compute_error,
  const bool &                                               output_error)
  : component_names(component_names)
  , unique_component_names(get_unique_component_names(component_names))
  , unique_component_masks(get_unique_component_masks(component_names))
  , norms_per_unique_component(list_of_error_norms)
  , extra_columns(extra_columns)
  , rate_key(rate_key)
  , precision(precision)
  , error_file_name(error_file_name)
  , compute_error(compute_error)
  , output_error(output_error)
{
  AssertDimension(unique_component_names.size(), list_of_error_norms.size());
}



void
ParsedConvergenceTable::add_parameters(ParameterHandler &prm)
{
  prm.add_parameter("Enable output to streams",
                    output_error,
                    "When set to false, printing of the convergence table "
                    "to the stream specified as input to output_table() "
                    "is disabled.");

  prm.add_parameter("Enable computation of the errors",
                    compute_error,
                    "When set to false, no computations are performed.");
  prm.add_parameter("Error precision",
                    precision,
                    "Number of digits to use when printing the error.",
                    Patterns::Integer(0));

  prm.add_parameter("Error file name",
                    error_file_name,
                    "Set this to a filename with extension .txt, .gpl, .org, "
                    "or .tex to enable writing the convergence table to a "
                    "file.");

  prm.add_parameter(
    "List of error norms to compute",
    norms_per_unique_component,
    "Each component is separated by a semicolon, "
    "and each norm by a comma. Implemented norms are Linfty, L2, "
    "W1infty, H1, and custom. If you want to skip a component, use none.");

  prm.add_parameter("Extra columns",
                    extra_columns,
                    "Extra columns to add to the table. Availabl options "
                    "are dofs, cells, and dt.");

  prm.add_parameter("Rate key",
                    rate_key,
                    "Key to use when computing convergence rates. If "
                    "this is set to a column that is not present, or to none, "
                    "then no error rates are computed.");
}



void
ParsedConvergenceTable::output_table(std::ostream &out)
{
  if (compute_error)
    {
      // Add convergence rates if the rate_key is not empty
      if (rate_key != ParsedConvergenceTableFlags::no_columns)
        {
          bool has_key = false;
          for (const auto &col : extra_columns)
            {
              if (rate_key == col)
                has_key = true;

              if (col != ParsedConvergenceTableFlags::no_columns)
                table.omit_column_from_convergence_rate_evaluation(
                  Patterns::Tools::to_string(col));
            }
          if (has_key)
            table.evaluate_all_convergence_rates(
              Patterns::Tools::to_string(rate_key),
              ConvergenceTable::reduction_rate_log2);
        }

      if (output_error)
        table.write_text(out);

      if (error_file_name != "")
        {
          const std::string error_file_format =
            error_file_name.substr(error_file_name.find_last_of(".") + 1);

          std::ofstream table_file(error_file_name);

          if (error_file_format == "tex")
            table.write_tex(table_file);
          else if (error_file_format == "txt")
            table.write_text(table_file);
          else if (error_file_format == "gpl")
            table.write_text(
              table_file, TableHandler::table_with_separate_column_description);
          else if (error_file_format == "org")
            table.write_text(table_file, TableHandler::org_mode_table);
          else
            {
              AssertThrow(false,
                          ExcInternalError(
                            std::string("Unrecognized file format: ") +
                            error_file_format));
            }
          table_file.close();
        }
    }
}
DEAL_II_NAMESPACE_CLOSE
