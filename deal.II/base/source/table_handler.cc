//----------------------------  table_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  table_handler.cc  ---------------------------


#include <base/table_handler.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif

#include <iostream>
#include <iomanip>
#include <strings.h>



TableEntryBase::TableEntryBase ()
{}


TableEntryBase::~TableEntryBase ()
{}

/*---------------------------------------------------------------------*/

template <typename number>
TableEntry<number>::TableEntry(const number value):
		val(value)  {}

template <typename number>
TableEntry<number>::~TableEntry()
{}

template <typename number>
number TableEntry<number>::value() const
{
  return val;
}


template <typename number>
void TableEntry<number>::write_tex (std::ostream &out) const
{
  out << val;
}

template <typename number>
void TableEntry<number>::write_text (std::ostream &out) const
{
  out << val;
}


/*---------------------------------------------------------------------*/

TableHandler::Column::Column(const std::string &tex_caption):
		tex_caption(tex_caption),
		tex_format("c"),
		precision(4),
		scientific(0),
		flag(0)  
{}


TableHandler::Column::Column():
		tex_caption(),
		tex_format("c"),
		precision(4),
		scientific(0),
		flag(0)  
{}


TableHandler::Column::~Column()
{
  for (unsigned int i=0; i<entries.size(); ++i)
    delete entries[i];
  entries.clear();
}

/*---------------------------------------------------------------------*/


TableHandler::TableHandler()
{}


template <typename number>
void TableHandler::add_value (const std::string &key, 
			      const number value)
{
  if (!columns.count(key))
    {
      std::pair<std::string, Column> new_column(key, Column(key));
      columns.insert(new_column);
      column_order.push_back(key);
    }

  TableEntry<number> *entry_ptr=new TableEntry<number>(value);
  const std::map<std::string, Column>::iterator col_iter=columns.find(key);
  Assert(col_iter!=columns.end(), ExcInternalError());
  
  col_iter->second.entries.push_back(entry_ptr);
}


void TableHandler::add_column_to_supercolumn (const std::string &key,
					      const std::string &superkey)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));

  if (!supercolumns.count(superkey))
    {
      std::pair<std::string, std::vector<std::string> > new_column(superkey,
								   std::vector<std::string>());
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
  columns[key].precision=precision;
}


void TableHandler::set_scientific (const std::string &key,
				   const bool scientific)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].scientific=scientific;
}


void TableHandler::write_text(std::ostream &out) const
{
  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  const unsigned int nrows  = n_rows();
  const unsigned int n_cols = sel_columns.size();

				   // first compute the widths of each
				   // entry of the table, in order to
				   // have a nicer alignement
  std::vector<std::vector<unsigned int> >
    entry_widths (nrows, std::vector<unsigned int>(n_cols));
  for (unsigned int i=0; i<nrows; ++i)
    for (unsigned int j=0; j<n_cols; ++j)
      {
					 // get key and entry here
	std::string key=sel_columns[j];
	const std::map<std::string, Column>::const_iterator
	  col_iter=columns.find(key);
	Assert(col_iter!=columns.end(), ExcInternalError());
	  
	const Column & column = col_iter->second;

					 // write it into a dummy
					 // stream, just to get its
					 // size upon output
#ifdef HAVE_STD_STRINGSTREAM
	std::ostringstream dummy_out;
#else
	std::ostrstream dummy_out;
#endif
	
	dummy_out << std::setprecision(column.precision);

	if (col_iter->second.scientific)
	  dummy_out.setf (std::ios::scientific, std::ios::floatfield);
	else
	  dummy_out.setf (std::ios::fixed, std::ios::floatfield);
	column.entries[i]->write_text (dummy_out);
#ifndef HAVE_STD_STRINGSTREAM
	dummy_out << std::ends;
#endif
	
					 // get size, note that we are
					 // not interested in the
					 // trailing \0
#ifdef HAVE_STD_STRINGSTREAM
	entry_widths[i][j] = dummy_out.str().length();
#else
	entry_widths[i][j] = strlen(dummy_out.str());
#endif
      };

				   // next compute the width each row
				   // has to have to suit all entries
  std::vector<unsigned int> column_widths (n_cols, 0);
  for (unsigned int i=0; i<nrows; ++i)
    for (unsigned int j=0; j<n_cols; ++j)
      column_widths[j] = std::max(entry_widths[i][j],
				  column_widths[j]);

				   // write the caption line
  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      const std::string & key = column_order[j];
      
				       // if the key of this column is
				       // wider than the widest entry,
				       // then adjust
      if (key.length() > column_widths[j])
	column_widths[j] = key.length();

				       // now write key. try to center
				       // it somehow
      const unsigned int front_padding = (column_widths[j]-key.length())/2,
			  rear_padding = (column_widths[j]-key.length()) -
					 front_padding;
      for (unsigned int i=0; i<front_padding; ++i)
	out << ' ';
      out << key;
      for (unsigned int i=0; i<rear_padding; ++i)
	out << ' ';

				       // finally column break
      out << ' ';
    }
  out << std::endl;
  
  for (unsigned int i=0; i<nrows; ++i)
    {    
      for (unsigned int j=0; j<n_cols; ++j)
	{
	  std::string key=sel_columns[j];
	  const std::map<std::string, Column>::const_iterator
	    col_iter=columns.find(key);
	  Assert(col_iter!=columns.end(), ExcInternalError());
	  
	  const Column &column=col_iter->second;
	  
	  out << std::setprecision(column.precision);

	  if (col_iter->second.scientific)
	    out.setf(std::ios::scientific, std::ios::floatfield);
	  else
	    out.setf(std::ios::fixed, std::ios::floatfield);
	  out << std::setw(column_widths[j]);
	  column.entries[i]->write_text(out);
	  out << " ";
	}
      out << std::endl;
    }
}


void TableHandler::write_tex(std::ofstream &out) const
{
  out << "\\documentclass[10pt]{report}" << std::endl
      << "\\usepackage{float}" << std::endl << std::endl << std::endl
      << "\\begin{document}" << std::endl
      << "\\begin{table}[H]" << std::endl
      << "\\begin{center}" << std::endl
      << "\\begin{tabular}{|";

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
	  
	  column.entries[i]->write_tex(out);
	  
	  if (j<n_cols-1)
	    out << " & ";
	}
      out << "\\\\ \\hline" << std::endl;
    }

  std::string caption="table";
  
  out << "\\end{tabular}" << std::endl
      << "\\end{center}" << std::endl
//      << "\\caption{"  << caption << "}" << std::endl
      << "\\end{table}" << std::endl
      << "\\end{document}" << std::endl;
}


unsigned int TableHandler::n_rows() const
{
  std::map<std::string, Column>::const_iterator col_iter=columns.begin();
  unsigned int n=col_iter->second.entries.size();
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
	  for (unsigned int j=0; j<n_subcolumns; ++j)
	    {
	      const std::string subkey=super_iter->second[j];
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


// explicit instantiations
template class TableEntry<double>;
template class TableEntry<float>;
template class TableEntry<int>;
template class TableEntry<unsigned int>;
template class TableEntry<std::string>;


template
void TableHandler::add_value(const std::string &, const double);

template
void TableHandler::add_value(const std::string &, const int);

template
void TableHandler::add_value(const std::string &, const unsigned int);

template
void TableHandler::add_value(const std::string &, const std::string);
