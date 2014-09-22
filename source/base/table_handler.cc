// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/table_handler.h>
#include <deal.II/base/table.h>

#include <sstream>
#include <iostream>
#include <iomanip>


DEAL_II_NAMESPACE_OPEN


/*---------------------------------------------------------------------*/

// inline and template functions
namespace internal
{
  TableEntry::TableEntry ()
  {}


  double TableEntry::get_numeric_value () const
  {
    // we don't quite know the data type in 'value', but
    // it must be one of the ones in the type list of the
    // boost::variant. Go through this list and return
    // the value if this happens to be a number
    //
    // first try with int
    try
      {
        return boost::get<int>(value);
      }
    catch (...)
      {}


    // ... then with unsigned int...
    try
      {
        return boost::get<unsigned int>(value);
      }
    catch (...)
      {}

    // ...and finally with double precision:
    try
      {
        return boost::get<double>(value);
      }
    catch (...)
      {
        Assert (false, ExcMessage ("The number stored by this element of the "
                                   "table is not a number."))
      }

    return 0;
  }

  void TableEntry::cache_string(bool scientific, unsigned int precision) const
  {
    std::ostringstream ss;

    ss << std::setprecision(precision);

    if (scientific)
      ss.setf(std::ios::scientific, std::ios::floatfield);
    else
      ss.setf(std::ios::fixed, std::ios::floatfield);

    ss << value;

    cached_value = ss.str();
    if (cached_value.size()==0)
      cached_value = "\"\"";
  }

  const std::string &TableEntry::get_cached_string() const
  {
    return cached_value;
  }


  namespace Local
  {
    // see which type we can cast to, then use this type to create
    // a default constructed object
    struct GetDefaultValue : public boost::static_visitor<>
    {
      template <typename T>
      void operator()( T &operand ) const
      {
        operand = T();
      }
    };
  }

  TableEntry TableEntry::get_default_constructed_copy () const
  {
    TableEntry new_entry = *this;
    boost::apply_visitor(Local::GetDefaultValue(), new_entry.value);

    return new_entry;
  }


}

/* ------------------------------------------------ */

TableHandler::Column::Column(const std::string &tex_caption)
  :
  tex_caption(tex_caption),
  tex_format("c"),
  precision(4),
  scientific(0),
  flag(0),
  max_length(0)
{}



TableHandler::Column::Column()
  :
  tex_caption(),
  tex_format("c"),
  precision(4),
  scientific(0),
  flag(0),
  max_length(0)
{}



void
TableHandler::Column::pad_column_below (const unsigned int size)
{
  // we should never have a column that is completely
  // empty and that needs to be padded
  Assert (entries.size() > 0, ExcInternalError());

  // add as many elements as necessary
  while (entries.size() < size)
    {
      entries.push_back (entries.back().get_default_constructed_copy());
      internal::TableEntry &entry = entries.back();
      entry.cache_string(scientific, precision);
      max_length = std::max(max_length, static_cast<unsigned int>(entry.get_cached_string().length()));
    }
}


void
TableHandler::Column::invalidate_cache()
{
  max_length = 0;

  for (std::vector<dealii::internal::TableEntry>::iterator it=entries.begin(); it!=entries.end(); ++it)
    {
      it->cache_string(this->scientific, this->precision);
      max_length = std::max(max_length, static_cast<unsigned int>(it->get_cached_string().length()));
    }
}


/*---------------------------------------------------------------------*/


TableHandler::TableHandler()
  :
  auto_fill_mode (false)
{}



void
TableHandler::set_auto_fill_mode (const bool state)
{
  auto_fill_mode = state;
}


void TableHandler::add_column_to_supercolumn (const std::string &key,
                                              const std::string &superkey)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));

  if (!supercolumns.count(superkey))
    {
      std::pair<std::string, std::vector<std::string> >
      new_column(superkey, std::vector<std::string>());
      supercolumns.insert(new_column);
      // replace key in column_order
      // by superkey
      for (unsigned int j=0; j<column_order.size(); ++j)
        if (column_order[j]==key)
          {
            column_order[j]=superkey;
            break;
          }
    }
  else
    {
      // remove key from column_order
      // for erase we need an iterator
      for (std::vector<std::string>::iterator order_iter=column_order.begin();
           order_iter!=column_order.end(); ++order_iter)
        if (*order_iter==key)
          {
            column_order.erase(order_iter);
            break;
          }
    }

  if (supercolumns.count(superkey))
    {
      supercolumns[superkey].push_back(key);
      // By default set the
      // tex_supercaption to superkey
      std::pair<std::string, std::string> new_tex_supercaption(superkey, superkey);
      tex_supercaptions.insert(new_tex_supercaption);
    }
  else
    Assert(false, ExcInternalError());
}



void TableHandler::set_column_order (const std::vector<std::string> &new_order)
{
  for (unsigned int j=0; j<new_order.size(); ++j)
    Assert(supercolumns.count(new_order[j]) || columns.count(new_order[j]),
           ExcColumnOrSuperColumnNotExistent(new_order[j]));

  column_order=new_order;
}


void TableHandler::set_tex_caption (const std::string &key,
                                    const std::string &tex_caption)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].tex_caption=tex_caption;
}



void TableHandler::set_tex_table_caption (const std::string &table_caption)
{
  tex_table_caption=table_caption;
}



void TableHandler::set_tex_table_label (const std::string &table_label)
{
  tex_table_label=table_label;
}



void TableHandler::set_tex_supercaption (const std::string &superkey,
                                         const std::string &tex_supercaption)
{
  Assert(supercolumns.count(superkey), ExcSuperColumnNotExistent(superkey));
  Assert(tex_supercaptions.count(superkey), ExcInternalError());
  tex_supercaptions[superkey]=tex_supercaption;
}



void TableHandler::set_tex_format (const std::string &key,
                                   const std::string &tex_format)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  Assert(tex_format=="l" || tex_format=="c" || tex_format=="r",
         ExcUndefinedTexFormat(tex_format));
  columns[key].tex_format=tex_format;
}



void TableHandler::set_precision (const std::string &key,
                                  const unsigned int precision)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  if (columns[key].precision!=precision)
    {
      columns[key].precision = precision;
      columns[key].invalidate_cache();
    }
}


void TableHandler::set_scientific (const std::string &key,
                                   const bool scientific)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  if (columns[key].scientific!=scientific)
    {
      columns[key].scientific = scientific;
      columns[key].invalidate_cache();
    }
}


void TableHandler::write_text(std::ostream &out,
                              const TextOutputFormat format) const
{
  AssertThrow (out, ExcIO());

  // first pad the table from below if necessary
  if (auto_fill_mode == true)
    {
      unsigned int max_rows = 0;
      for (std::map<std::string, Column>::const_iterator p = columns.begin();
           p != columns.end(); ++p)
        max_rows = std::max<unsigned int>(max_rows, p->second.entries.size());

      for (std::map<std::string, Column>::iterator p = columns.begin();
           p != columns.end(); ++p)
        p->second.pad_column_below (max_rows);
    }

  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  const unsigned int nrows  = n_rows();
  const unsigned int n_cols = sel_columns.size();

  // cache the columns and compute the widths of each column for alignment
  std::vector<const Column *> cols;
  std::vector<unsigned int> column_widths (n_cols, 0);
  for (unsigned int j=0; j<n_cols; ++j)
    {
      std::string key=sel_columns[j];
      const std::map<std::string, Column>::const_iterator
      col_iter=columns.find(key);
      Assert(col_iter!=columns.end(), ExcInternalError());
      cols.push_back(&(col_iter->second));

      column_widths[j] = col_iter->second.max_length;
    }

  switch (format)
    {
    case org_mode_table:
    {
      // write the captions
      out << "| " << std::left;
      for (unsigned int j=0; j<n_cols; ++j)
        {
          const std::string &key = sel_columns[j];
          column_widths[j] = std::max(column_widths[j],
                                      (unsigned int)key.length());
          out << std::setw(column_widths[j]);
          out << key << " | ";
        }
      out << std::endl;

      // write the body
      for (unsigned int i=0; i<nrows; ++i)
        {
          out << "| ";
          for (unsigned int j=0; j<n_cols; ++j)
            {
              const Column &column=*(cols[j]);

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
      for (unsigned int j=0; j<n_cols; ++j)
        {
          const std::string &key = sel_columns[j];
          out << "# " << j+1 << ": " << key << '\n';
        }

      // write the body
      for (unsigned int i=0; i<nrows; ++i)
        {
          for (unsigned int j=0; j<n_cols; ++j)
            {
              const Column &column=*(cols[j]);

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
      for (unsigned int j=0; j<n_cols; ++j)
        {
          std::string key=sel_columns[j];
          out << "# " << j+1 << ": " << key << '\n';
        }
      break;
    }

    case table_with_headers:
    {
      // This format output supercolumn headers and aligns them centered
      // over all the columns that belong to it.
      for (unsigned int j=0; j<column_order.size(); ++j)
        {
          const std::string &key = column_order[j];
          unsigned int width=0;
          {
            // compute the width of this column or supercolumn
            const std::map<std::string, std::vector<std::string> >::const_iterator
            super_iter=supercolumns.find(key);
            if (super_iter!=supercolumns.end())
              {
                const unsigned int n_subcolumns=super_iter->second.size();
                for (unsigned int k=0; k<n_subcolumns; ++k)
                  {
                    const std::map<std::string, Column>::const_iterator
                    col_iter=columns.find(super_iter->second[k]);
                    Assert(col_iter!=columns.end(), ExcInternalError());

                    width += col_iter->second.max_length;
                  }
                width += n_subcolumns - 1; // separators between subcolumns
              }
            else
              {
                const std::map<std::string, Column>::const_iterator
                col_iter=columns.find(key);

                width = col_iter->second.max_length;
              }
          }

          // header is longer than the column(s) under it
          if (width<key.length())
            {
              // make the column or the last column in this
              // supercolumn wide enough
              std::string colname;

              const std::map<std::string, std::vector<std::string> >::const_iterator
              super_iter=supercolumns.find(key);
              if (super_iter!=supercolumns.end())
                colname = super_iter->second.back();
              else
                colname = key;

              // find column and change output width
              for (unsigned int i=0; i<n_cols; ++i)
                {
                  if (sel_columns[i]==colname)
                    {
                      column_widths[i] += key.length() - width;
                      break;
                    }
                }

              width=key.length();
            }

          // now write key. try to center it somehow
          const unsigned int front_padding = (width-key.length())/2,
                             rear_padding  = (width-key.length()) -
                                             front_padding;
          for (unsigned int i=0; i<front_padding; ++i)
            out << ' ';
          out << key;
          for (unsigned int i=0; i<rear_padding; ++i)
            out << ' ';

          out << ' ';
        }
      out << '\n';
      break;
    }

    default:
      Assert (false, ExcInternalError());
    }


  // finally output the data itself for
  // table_with_headers or table_with_separate_column_description:
  for (unsigned int i=0; i<nrows; ++i)
    {
      for (unsigned int j=0; j<n_cols; ++j)
        {
          const Column &column=*(cols[j]);
          out << std::setw(column_widths[j]);
          out << column.entries[i].get_cached_string();

          // pad after this column
          out << ' ';
        }
      out << '\n';
    }
  out << std::flush;
}


void TableHandler::write_tex (std::ostream &out, const bool with_header) const
{
  //TODO[TH]: update code similar to
  //write_text() to use the cache
  AssertThrow (out, ExcIO());
  if (with_header)
    out << "\\documentclass[10pt]{report}" << std::endl
        << "\\usepackage{float}" << std::endl << std::endl << std::endl
        << "\\begin{document}" << std::endl;

  out << "\\begin{table}[H]" << std::endl
      << "\\begin{center}" << std::endl
      << "\\begin{tabular}{|";

  // first pad the table from below if necessary
  if (auto_fill_mode == true)
    {
      unsigned int max_rows = 0;
      for (std::map<std::string, Column>::const_iterator p = columns.begin();
           p != columns.end(); ++p)
        max_rows = std::max<unsigned int>(max_rows, p->second.entries.size());

      for (std::map<std::string, Column>::iterator p = columns.begin();
           p != columns.end(); ++p)
        p->second.pad_column_below (max_rows);
    }

  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  // write the column formats
  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      std::string key=column_order[j];
      // avoid `supercolumns[key]'
      const std::map<std::string, std::vector<std::string> >::const_iterator
      super_iter=supercolumns.find(key);

      if (super_iter!=supercolumns.end())
        {
          const unsigned int n_subcolumns=super_iter->second.size();
          for (unsigned int k=0; k<n_subcolumns; ++k)
            {
              // avoid `columns[supercolumns[key]]'
              const std::map<std::string, Column>::const_iterator
              col_iter=columns.find(super_iter->second[k]);
              Assert(col_iter!=columns.end(), ExcInternalError());

              out << col_iter->second.tex_format << "|";
            }
        }
      else
        {
          // avoid `columns[key]';
          const std::map<std::string, Column>::const_iterator
          col_iter=columns.find(key);
          Assert(col_iter!=columns.end(), ExcInternalError());
          out << col_iter->second.tex_format << "|";
        }
    }
  out << "} \\hline" << std::endl;

  // write the caption line of the table

  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      std::string key=column_order[j];
      const std::map<std::string, std::vector<std::string> >::const_iterator
      super_iter=supercolumns.find(key);

      if (super_iter!=supercolumns.end())
        {
          const unsigned int n_subcolumns=super_iter->second.size();
          // avoid use of `tex_supercaptions[key]'
          std::map<std::string,std::string>::const_iterator
          tex_super_cap_iter=tex_supercaptions.find(key);
          out << std::endl << "\\multicolumn{" << n_subcolumns << "}{|c|}{"
              << tex_super_cap_iter->second << "}";
        }
      else
        {
          // col_iter->second=columns[col];
          const std::map<std::string, Column>::const_iterator
          col_iter=columns.find(key);
          Assert(col_iter!=columns.end(), ExcInternalError());
          out << col_iter->second.tex_caption;
        }
      if (j<column_order.size()-1)
        out << " & ";
    }
  out << "\\\\ \\hline" << std::endl;

  // write the n rows
  const unsigned int nrows=n_rows();
  for (unsigned int i=0; i<nrows; ++i)
    {
      const unsigned int n_cols=sel_columns.size();

      for (unsigned int j=0; j<n_cols; ++j)
        {
          std::string key=sel_columns[j];
          // avoid `column[key]'
          const std::map<std::string, Column>::const_iterator
          col_iter=columns.find(key);
          Assert(col_iter!=columns.end(), ExcInternalError());

          const Column &column=col_iter->second;

          out << std::setprecision(column.precision);

          if (col_iter->second.scientific)
            out.setf(std::ios::scientific, std::ios::floatfield);
          else
            out.setf(std::ios::fixed, std::ios::floatfield);

          out << column.entries[i].value;

          if (j<n_cols-1)
            out << " & ";
        }
      out << "\\\\ \\hline" << std::endl;
    }

  out   << "\\end{tabular}" << std::endl
        << "\\end{center}" << std::endl;
  if (tex_table_caption!="")
    out << "\\caption{"  << tex_table_caption << "}" << std::endl;
  if (tex_table_label!="")
    out << "\\label{"   << tex_table_label << "}" << std::endl;
  out   << "\\end{table}" << std::endl;
  if (with_header)
    out << "\\end{document}" << std::endl;
}


unsigned int TableHandler::n_rows() const
{
  if (columns.size() == 0)
    return 0;

  std::map<std::string, Column>::const_iterator col_iter = columns.begin();
  unsigned int n = col_iter->second.entries.size();
  std::string first_name=col_iter->first;

  for (++col_iter; col_iter!=columns.end(); ++col_iter)
    Assert(col_iter->second.entries.size()==n,
           ExcWrongNumberOfDataEntries(col_iter->first,
                                       col_iter->second.entries.size(),
                                       first_name, n));

  return n;
}


void TableHandler::get_selected_columns(std::vector<std::string> &sel_columns) const
{
  sel_columns.clear();

  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      std::string key=column_order[j];
      const std::map<std::string, std::vector<std::string> >::const_iterator
      super_iter=supercolumns.find(key);

      if (super_iter!=supercolumns.end())
        {
          // i.e. key is a supercolumn key
          const unsigned int n_subcolumns=super_iter->second.size();
          for (unsigned int k=0; k<n_subcolumns; ++k)
            {
              const std::string subkey=super_iter->second[k];
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


DEAL_II_NAMESPACE_CLOSE
