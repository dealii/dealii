// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
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
#include <deal.II/base/table_handler.h>

#include <boost/io/ios_state.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>


DEAL_II_NAMESPACE_OPEN


/*---------------------------------------------------------------------*/

// inline and template functions
namespace internal
{
  double
  TableEntry::get_numeric_value() const
  {
    // we don't quite know the data type in 'value', but
    // it must be one of the ones in the type list of the
    // std::variant. Go through this list and return
    // the value if this happens to be a number
    //
    // first try with int
    try
      {
        return std::get<int>(value);
      }
    catch (...)
      {}


    // ... then with unsigned int...
    try
      {
        return std::get<unsigned int>(value);
      }
    catch (...)
      {}

    // ... then with std::uint64_t...
    try
      {
        return std::get<std::uint64_t>(value);
      }
    catch (...)
      {}

    // ...and finally with double precision:
    try
      {
        return std::get<double>(value);
      }
    catch (...)
      {
        Assert(false,
               ExcMessage("The number stored by this element of the "
                          "table is not a number."));
      }

    return 0;
  }

  void
  TableEntry::cache_string(bool scientific, unsigned int precision) const
  {
    std::ostringstream ss;

    ss << std::setprecision(precision);

    if (scientific)
      ss.setf(std::ios::scientific, std::ios::floatfield);
    else
      ss.setf(std::ios::fixed, std::ios::floatfield);

    if (scientific)
      ss << get_numeric_value();
    else
      std::visit([&ss](const auto &v) { ss << v; }, value);

    cached_value = ss.str();
    if (cached_value.empty())
      cached_value = "\"\"";
  }

  const std::string &
  TableEntry::get_cached_string() const
  {
    return cached_value;
  }

  TableEntry
  TableEntry::get_default_constructed_copy() const
  {
    TableEntry new_entry = *this;
    // Let std::visit figure out which data type is actually stored,
    // and then set the object so stored to a default-constructed
    // one.
    std::visit([](
                 auto &arg) { arg = std::remove_reference_t<decltype(arg)>(); },
               new_entry.value);

    return new_entry;
  }


} // namespace internal

/* ------------------------------------------------ */

TableHandler::Column::Column(const std::string &tex_caption)
  : tex_caption(tex_caption)
  , tex_format("c")
  , precision(4)
  , scientific(false)
  , flag(0)
  , max_length(0)
{}



TableHandler::Column::Column()
  : tex_caption()
  , tex_format("c")
  , precision(4)
  , scientific(false)
  , flag(0)
  , max_length(0)
{}



void
TableHandler::Column::pad_column_below(const unsigned int size)
{
  // we should never have a column that is completely
  // empty and that needs to be padded
  Assert(entries.size() > 0, ExcInternalError());

  // add as many elements as necessary
  while (entries.size() < size)
    {
      entries.push_back(entries.back().get_default_constructed_copy());
      const internal::TableEntry &entry = entries.back();
      entry.cache_string(scientific, precision);
      max_length =
        std::max(max_length,
                 static_cast<unsigned int>(entry.get_cached_string().size()));
    }
}


void
TableHandler::Column::invalidate_cache()
{
  max_length = 0;

  for (const auto &entry : entries)
    {
      entry.cache_string(this->scientific, this->precision);
      max_length =
        std::max(max_length,
                 static_cast<unsigned int>(entry.get_cached_string().size()));
    }
}


/*---------------------------------------------------------------------*/


TableHandler::TableHandler()
  : auto_fill_mode(false)
{}



void
TableHandler::declare_column(const std::string &key)
{
  // see if the column already exists; insert it if not
  Assert(columns.find(key) == columns.end(),
         ExcMessage("You are trying to declare a column with key <" + key +
                    "> but such a column already exists."));

  columns.insert(std::make_pair(key, Column(key)));
  column_order.push_back(key);
}



void
TableHandler::start_new_row()
{
  // figure out the longest current column
  unsigned int max_col_length = 0;
  for (const auto &column : columns)
    max_col_length =
      std::max(max_col_length,
               static_cast<unsigned int>(column.second.entries.size()));


  // then pad all columns to that length with empty strings
  for (auto &column : columns)
    while (column.second.entries.size() < max_col_length)
      {
        column.second.entries.emplace_back("");
        const internal::TableEntry &entry = column.second.entries.back();
        entry.cache_string(column.second.scientific, column.second.precision);
        column.second.max_length =
          std::max(column.second.max_length,
                   static_cast<unsigned int>(entry.get_cached_string().size()));
      }
}



void
TableHandler::set_auto_fill_mode(const bool state)
{
  auto_fill_mode = state;
}


void
TableHandler::add_column_to_supercolumn(const std::string &key,
                                        const std::string &superkey)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));

  if (supercolumns.count(superkey) == 0u)
    {
      std::pair<std::string, std::vector<std::string>> new_column(
        superkey, std::vector<std::string>());
      supercolumns.insert(new_column);
      // replace key in column_order
      // by superkey
      for (auto &column : column_order)
        if (column == key)
          {
            column = superkey;
            break;
          }
    }
  else
    {
      // remove key from column_order
      const auto order_iter =
        std::find(column_order.begin(), column_order.end(), key);
      if (order_iter != column_order.end())
        column_order.erase(order_iter);
    }

  if (supercolumns.count(superkey) != 0u)
    {
      supercolumns[superkey].push_back(key);
      // By default set the
      // tex_supercaption to superkey
      std::pair<std::string, std::string> new_tex_supercaption(superkey,
                                                               superkey);
      tex_supercaptions.insert(new_tex_supercaption);
    }
  else
    DEAL_II_ASSERT_UNREACHABLE();
}



void
TableHandler::set_column_order(const std::vector<std::string> &new_order)
{
  for (const auto &new_column : new_order)
    {
      (void)new_column;
      Assert(supercolumns.count(new_column) || columns.count(new_column),
             ExcColumnOrSuperColumnNotExistent(new_column));
    }

  column_order = new_order;
}


void
TableHandler::set_tex_caption(const std::string &key,
                              const std::string &tex_caption)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].tex_caption = tex_caption;
}



void
TableHandler::set_tex_table_caption(const std::string &table_caption)
{
  tex_table_caption = table_caption;
}



void
TableHandler::set_tex_table_label(const std::string &table_label)
{
  tex_table_label = table_label;
}



void
TableHandler::set_tex_supercaption(const std::string &superkey,
                                   const std::string &tex_supercaption)
{
  Assert(supercolumns.count(superkey), ExcSuperColumnNotExistent(superkey));
  Assert(tex_supercaptions.count(superkey), ExcInternalError());
  tex_supercaptions[superkey] = tex_supercaption;
}



void
TableHandler::set_tex_format(const std::string &key,
                             const std::string &tex_format)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  Assert(tex_format == "l" || tex_format == "c" || tex_format == "r",
         ExcUndefinedTexFormat(tex_format));
  columns[key].tex_format = tex_format;
}



void
TableHandler::set_precision(const std::string &key,
                            const unsigned int precision)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  if (columns[key].precision != precision)
    {
      columns[key].precision = precision;
      columns[key].invalidate_cache();
    }
}


void
TableHandler::set_scientific(const std::string &key, const bool scientific)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  if (columns[key].scientific != scientific)
    {
      columns[key].scientific = scientific;
      columns[key].invalidate_cache();
    }
}


void
TableHandler::write_text(std::ostream &out, const TextOutputFormat format) const
{
  AssertThrow(out.fail() == false, ExcIO());
  boost::io::ios_flags_saver restore_flags(out);

  // first pad the table from below if necessary
  if (auto_fill_mode == true)
    {
      unsigned int max_rows = 0;
      for (std::map<std::string, Column>::const_iterator p = columns.begin();
           p != columns.end();
           ++p)
        max_rows = std::max<unsigned int>(max_rows, p->second.entries.size());

      for (auto &column : columns)
        column.second.pad_column_below(max_rows);
    }

  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  const unsigned int nrows  = n_rows();
  const unsigned int n_cols = sel_columns.size();

  // cache the columns and compute the widths of each column for alignment
  std::vector<const Column *> cols;
  std::vector<unsigned int>   column_widths(n_cols, 0);
  for (unsigned int j = 0; j < n_cols; ++j)
    {
      const std::string                                  &key = sel_columns[j];
      const std::map<std::string, Column>::const_iterator col_iter =
        columns.find(key);
      Assert(col_iter != columns.end(), ExcInternalError());
      cols.push_back(&(col_iter->second));

      column_widths[j] = col_iter->second.max_length;
    }

  switch (format)
    {
      case org_mode_table:
        {
          // write the captions
          out << "| " << std::left;
          for (unsigned int j = 0; j < n_cols; ++j)
            {
              const std::string &key = sel_columns[j];
              column_widths[j] =
                std::max(column_widths[j],
                         static_cast<unsigned int>(key.size()));
              out << std::setw(column_widths[j]);
              out << key << " | ";
            }
          out << '\n';

          // write the body
          for (unsigned int i = 0; i < nrows; ++i)
            {
              out << "| ";
              for (unsigned int j = 0; j < n_cols; ++j)
                {
                  const Column &column = *(cols[j]);

                  out << std::setw(column_widths[j]);
                  out << column.entries[i].get_cached_string();
                  out << " | ";
                }
              out << '\n';
            }

          out << std::flush;
          return;
        }

      case simple_table_with_separate_column_description:
        {
          // write the captions
          for (unsigned int j = 0; j < n_cols; ++j)
            {
              const std::string &key = sel_columns[j];
              out << "# " << j + 1 << ": " << key << '\n';
            }

          // write the body
          for (unsigned int i = 0; i < nrows; ++i)
            {
              for (unsigned int j = 0; j < n_cols; ++j)
                {
                  const Column &column = *(cols[j]);

                  out << column.entries[i].get_cached_string();
                  out << ' ';
                }
              out << '\n';
            }

          out << std::flush;
          return;
        }

      case table_with_separate_column_description:
        {
          // writing the captions for table_with_separate_column_description
          // means that we ignore supercolumns and output the column
          // header for each column. enumerate columns starting with 1
          for (unsigned int j = 0; j < n_cols; ++j)
            {
              const std::string &key = sel_columns[j];
              out << "# " << j + 1 << ": " << key << '\n';
            }
          break;
        }

      case table_with_headers:
        {
          // This format output supercolumn headers and aligns them centered
          // over all the columns that belong to it.
          for (const auto &key : column_order)
            {
              unsigned int width = 0;
              {
                // compute the width of this column or supercolumn
                const std::map<std::string,
                               std::vector<std::string>>::const_iterator
                  super_iter = supercolumns.find(key);
                if (super_iter != supercolumns.end())
                  {
                    const unsigned int n_subcolumns = super_iter->second.size();
                    for (unsigned int k = 0; k < n_subcolumns; ++k)
                      {
                        const std::map<std::string, Column>::const_iterator
                          col_iter = columns.find(super_iter->second[k]);
                        Assert(col_iter != columns.end(), ExcInternalError());

                        width += col_iter->second.max_length;
                      }
                    width += n_subcolumns - 1; // separators between subcolumns
                  }
                else
                  {
                    const std::map<std::string, Column>::const_iterator
                      col_iter = columns.find(key);

                    width = col_iter->second.max_length;
                  }
              }

              // header is longer than the column(s) under it
              if (width < key.size())
                {
                  // make the column or the last column in this
                  // supercolumn wide enough
                  std::string colname;

                  const std::map<std::string,
                                 std::vector<std::string>>::const_iterator
                    super_iter = supercolumns.find(key);
                  if (super_iter != supercolumns.end())
                    colname = super_iter->second.back();
                  else
                    colname = key;

                  // find column and change output width
                  for (unsigned int i = 0; i < n_cols; ++i)
                    {
                      if (sel_columns[i] == colname)
                        {
                          column_widths[i] += key.size() - width;
                          break;
                        }
                    }

                  width = key.size();
                }

              // now write key. try to center it somehow
              const unsigned int front_padding = (width - key.size()) / 2,
                                 rear_padding =
                                   (width - key.size()) - front_padding;
              for (unsigned int i = 0; i < front_padding; ++i)
                out << ' ';
              out << key;
              for (unsigned int i = 0; i < rear_padding; ++i)
                out << ' ';

              out << ' ';
            }
          out << '\n';
          break;
        }

      default:
        DEAL_II_ASSERT_UNREACHABLE();
    }


  // finally output the data itself for
  // table_with_headers or table_with_separate_column_description:
  for (unsigned int i = 0; i < nrows; ++i)
    {
      for (unsigned int j = 0; j < n_cols; ++j)
        {
          const Column &column = *(cols[j]);
          out << std::setw(column_widths[j]);
          out << column.entries[i].get_cached_string();

          // pad after this column
          out << ' ';
        }
      out << '\n';
    }
  out << std::flush;
}


void
TableHandler::write_tex(std::ostream &out, const bool with_header) const
{
  // TODO[TH]: update code similar to
  // write_text() to use the cache
  AssertThrow(out.fail() == false, ExcIO());
  if (with_header)
    out << "\\documentclass[10pt]{report}" << '\n'
        << "\\usepackage{float}" << '\n'
        << '\n'
        << '\n'
        << "\\begin{document}" << '\n';

  out << "\\begin{table}[H]" << '\n'
      << "\\begin{center}" << '\n'
      << "\\begin{tabular}{|";

  // first pad the table from below if necessary
  if (auto_fill_mode == true)
    {
      unsigned int max_rows = 0;
      for (std::map<std::string, Column>::const_iterator p = columns.begin();
           p != columns.end();
           ++p)
        max_rows = std::max<unsigned int>(max_rows, p->second.entries.size());

      for (auto &column : columns)
        column.second.pad_column_below(max_rows);
    }

  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  // write the column formats
  for (const auto &key : column_order)
    {
      // avoid `supercolumns[key]'
      const std::map<std::string, std::vector<std::string>>::const_iterator
        super_iter = supercolumns.find(key);

      if (super_iter != supercolumns.end())
        {
          const unsigned int n_subcolumns = super_iter->second.size();
          for (unsigned int k = 0; k < n_subcolumns; ++k)
            {
              // avoid `columns[supercolumns[key]]'
              const std::map<std::string, Column>::const_iterator col_iter =
                columns.find(super_iter->second[k]);
              Assert(col_iter != columns.end(), ExcInternalError());

              out << col_iter->second.tex_format << "|";
            }
        }
      else
        {
          // avoid `columns[key]';
          const std::map<std::string, Column>::const_iterator col_iter =
            columns.find(key);
          Assert(col_iter != columns.end(), ExcInternalError());
          out << col_iter->second.tex_format << "|";
        }
    }
  out << "} \\hline" << '\n';

  // write the caption line of the table

  for (unsigned int j = 0; j < column_order.size(); ++j)
    {
      std::string key = column_order[j];
      const std::map<std::string, std::vector<std::string>>::const_iterator
        super_iter = supercolumns.find(key);

      if (super_iter != supercolumns.end())
        {
          const unsigned int n_subcolumns = super_iter->second.size();
          // avoid use of `tex_supercaptions[key]'
          std::map<std::string, std::string>::const_iterator
            tex_super_cap_iter = tex_supercaptions.find(key);
          out << '\n'
              << "\\multicolumn{" << n_subcolumns << "}{|c|}{"
              << tex_super_cap_iter->second << "}";
        }
      else
        {
          // col_iter->second=columns[col];
          const std::map<std::string, Column>::const_iterator col_iter =
            columns.find(key);
          Assert(col_iter != columns.end(), ExcInternalError());
          out << col_iter->second.tex_caption;
        }
      if (j < column_order.size() - 1)
        out << " & ";
    }
  out << "\\\\ \\hline" << '\n';

  // write the n rows
  const unsigned int nrows = n_rows();
  for (unsigned int i = 0; i < nrows; ++i)
    {
      const unsigned int n_cols = sel_columns.size();

      for (unsigned int j = 0; j < n_cols; ++j)
        {
          const std::string &key = sel_columns[j];
          // avoid `column[key]'
          const std::map<std::string, Column>::const_iterator col_iter =
            columns.find(key);
          Assert(col_iter != columns.end(), ExcInternalError());

          const Column &column = col_iter->second;

          out << std::setprecision(column.precision);

          if (col_iter->second.scientific)
            out.setf(std::ios::scientific, std::ios::floatfield);
          else
            out.setf(std::ios::fixed, std::ios::floatfield);

          std::visit([&out](const auto &v) { out << v; },
                     column.entries[i].value);

          if (j < n_cols - 1)
            out << " & ";
        }
      out << "\\\\ \\hline" << '\n';
    }

  out << "\\end{tabular}" << '\n' << "\\end{center}" << '\n';
  if (!tex_table_caption.empty())
    out << "\\caption{" << tex_table_caption << "}" << '\n';
  if (!tex_table_label.empty())
    out << "\\label{" << tex_table_label << "}" << '\n';
  out << "\\end{table}" << '\n';
  if (with_header)
    out << "\\end{document}" << '\n';

  // Now flush all of the data we've put into the stream to make it
  // sure it really gets written.
  out << std::flush;
}



void
TableHandler::clear()
{
  columns.clear();
  supercolumns.clear();
  column_order.clear();
  tex_supercaptions.clear();

  tex_table_label.clear();
  tex_table_caption.clear();
}


unsigned int
TableHandler::n_rows() const
{
  if (columns.empty())
    return 0;

  std::map<std::string, Column>::const_iterator col_iter = columns.begin();
  unsigned int n = col_iter->second.entries.size();

  if constexpr (running_in_debug_mode())
    {
      std::string first_name = col_iter->first;
      for (++col_iter; col_iter != columns.end(); ++col_iter)
        Assert(col_iter->second.entries.size() == n,
               ExcWrongNumberOfDataEntries(col_iter->first,
                                           col_iter->second.entries.size(),
                                           first_name,
                                           n));
    }

  return n;
}


void
TableHandler::get_selected_columns(std::vector<std::string> &sel_columns) const
{
  sel_columns.clear();

  for (const auto &key : column_order)
    {
      const std::map<std::string, std::vector<std::string>>::const_iterator
        super_iter = supercolumns.find(key);

      if (super_iter != supercolumns.end())
        {
          // i.e. key is a supercolumn key
          const unsigned int n_subcolumns = super_iter->second.size();
          for (unsigned int k = 0; k < n_subcolumns; ++k)
            {
              const std::string subkey = super_iter->second[k];
              Assert(columns.count(subkey), ExcInternalError());
              sel_columns.push_back(subkey);
            }
        }
      else
        {
          Assert(columns.count(key), ExcInternalError());
          // i.e. key is a column key
          sel_columns.push_back(key);
        }
    }
}


void
TableHandler::clear_current_row()
{
  // Figure out what is the correct (max) length of the columns
  // so that we "shave" one off.
  std::vector<internal::TableEntry>::size_type n = 0;
  for (const auto &column : columns)
    n = std::max(n, column.second.entries.size());

  // shave the top most element
  if (n != 0)
    for (auto &column : columns)
      if (column.second.entries.size() == n)
        column.second.entries.pop_back();
}


DEAL_II_NAMESPACE_CLOSE
