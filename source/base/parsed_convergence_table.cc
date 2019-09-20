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
          masks.emplace_back(bools[j++]);
        bools[j][i] = true;
      }
    masks.emplace_back(bools[j++]);
    AssertDimension(j, unique_component_names.size());
    return masks;
  }
} // namespace



ParsedConvergenceTable::ParsedConvergenceTable(
  const std::vector<std::string> &                    solution_names,
  const std::vector<std::set<VectorTools::NormType>> &list_of_error_norms)
  : ParsedConvergenceTable(solution_names,
                           list_of_error_norms,
                           2.0,
                           {"cells", "dofs"},
                           "dofs",
                           "reduction_rate_log2",
                           "",
                           3,
                           true)
{}



ParsedConvergenceTable::ParsedConvergenceTable(
  const std::vector<std::string> &                    component_names,
  const std::vector<std::set<VectorTools::NormType>> &list_of_error_norms,
  const double &                                      exponent,
  const std::set<std::string> &                       extra_columns,
  const std::string &                                 rate_key,
  const std::string &                                 rate_mode,
  const std::string &                                 error_file_name,
  const unsigned int &                                precision,
  const bool &                                        compute_error)
  : component_names(component_names)
  , unique_component_names(get_unique_component_names(component_names))
  , unique_component_masks(get_unique_component_masks(component_names))
  , norms_per_unique_component(list_of_error_norms)
  , exponent(exponent)
  , extra_columns(extra_columns)
  , rate_key(rate_key)
  , rate_mode(rate_mode)
  , precision(precision)
  , error_file_name(error_file_name)
  , compute_error(compute_error)
{
  AssertDimension(norms_per_unique_component.size(),
                  unique_component_names.size());
}



void
ParsedConvergenceTable::add_parameters(ParameterHandler &prm)
{
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
    "Each component is separated by a semicolon "
    "and each norm by a comma. See the documentation of VectorTools::NormType "
    "for a list of implemented norms. If you want to skip a component, leave "
    "its entry empty.");


  prm.add_parameter("Exponent for p-norms",
                    exponent,
                    "The exponent to use when computing p-norms.",
                    Patterns::Double(1));

  prm.add_parameter("Extra columns",
                    extra_columns,
                    "Extra columns to add to the table. Available options "
                    "are dofs and cells.",
                    Patterns::List(Patterns::Selection("dofs|cells")));

  prm.add_parameter("Rate key",
                    rate_key,
                    "Key to use when computing convergence rates. If "
                    "this is set to a column that is not present, or to the "
                    "empty string, then no error rates are computed.");

  prm.add_parameter("Rate mode",
                    rate_mode,
                    "What type of error rate to compute. Available options are "
                    "reduction_rate_log2, reduction_rate, and none.",
                    Patterns::Selection(
                      "reduction_rate|reduction_rate_log2|none"));
}



void
ParsedConvergenceTable::prepare_table_for_output()
{
  if (compute_error)
    {
      // Add convergence rates if the rate_key is not empty
      if (rate_key != "")
        {
          bool has_key = false;
          for (const auto &col : extra_columns)
            {
              if (rate_key == col)
                has_key = true;

              if (col != "")
                table.omit_column_from_convergence_rate_evaluation(col);
            }

          for (const auto &extra_col : extra_column_functions)
            if (extra_col.second.second == false)
              {
                if (rate_key == extra_col.first)
                  has_key = true;
                table.omit_column_from_convergence_rate_evaluation(
                  extra_col.first);
              }

          if (has_key)
            {
              if (rate_mode == "reduction_rate_log2")
                table.evaluate_all_convergence_rates(
                  rate_key, ConvergenceTable::reduction_rate_log2);
              else if (rate_mode == "reduction_rate")
                table.evaluate_all_convergence_rates(
                  rate_key, ConvergenceTable::reduction_rate);
              else
                {
                  Assert(rate_mode != "none", ExcInternalError());
                }
            }
          else
            {
              AssertThrow(rate_key != "",
                          ExcMessage(
                            "You specified the key <" + rate_key +
                            "> to compute convergence rates, but that key does "
                            "not exist in the current table."));
            }
        }
    }
}



void
ParsedConvergenceTable::output_table(std::ostream &out)
{
  if (compute_error)
    {
      prepare_table_for_output();
      table.write_text(out);
      output_table();
    }
}



void
ParsedConvergenceTable::output_table()
{
  if (compute_error && error_file_name != "")
    {
      prepare_table_for_output();

      const std::string error_file_format =
        error_file_name.substr(error_file_name.find_last_of('.') + 1);

      std::ofstream table_file(error_file_name);

      if (error_file_format == "tex")
        table.write_tex(table_file);
      else if (error_file_format == "txt")
        table.write_text(table_file);
      else if (error_file_format == "gpl")
        table.write_text(table_file,
                         TableHandler::table_with_separate_column_description);
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



void
ParsedConvergenceTable::add_extra_column(
  const std::string &            column_name,
  const std::function<double()> &custom_function,
  const bool &                   compute_rate)
{
  extra_column_functions[column_name] = {custom_function, compute_rate};
}

DEAL_II_NAMESPACE_CLOSE
