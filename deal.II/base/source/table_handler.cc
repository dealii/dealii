//----------------------------  table_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  table_handler.cc  ---------------------------


#include <iomanip>


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
void TableEntry<number>::write_tex (ostream &out) const
{
  out << val;
}

template <typename number>
void TableEntry<number>::write_text (ostream &out) const
{
  out << val;
}


/*---------------------------------------------------------------------*/

TableHandler::Column::Column(const string &tex_caption):
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
void TableHandler::add_value (const string &key, 
			      const number value)
{
  if (!columns.count(key))
    {
      pair<string, Column> new_column(key, Column(key));
      columns.insert(new_column);
      column_order.push_back(key);
    }

  TableEntry<number> *entry_ptr=new TableEntry<number>(value);
  const map<string, Column>::iterator col_iter=columns.find(key);
  Assert(col_iter!=columns.end(), ExcInternalError());
  
  col_iter->second.entries.push_back(entry_ptr);
}


void TableHandler::add_column_to_supercolumn (const string &key,
					      const string &superkey)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));

  if (!supercolumns.count(superkey))
    {
      pair<string, vector<string> > new_column(superkey, vector<string>());
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
      for (vector<string>::iterator order_iter=column_order.begin(); 
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
      pair<string, string> new_tex_supercaption(superkey, superkey);
      tex_supercaptions.insert(new_tex_supercaption);
    }
  else
    Assert(false, ExcInternalError());
}


void TableHandler::set_column_order (const vector<string> &new_order)
{
  for (unsigned int j=0; j<new_order.size(); ++j)
    Assert(supercolumns.count(new_order[j]) || columns.count(new_order[j]),
	   ExcColumnOrSuperColumnNotExistent(new_order[j]));

  column_order=new_order;
}


void TableHandler::set_tex_caption (const string &key, 
				    const string &tex_caption)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].tex_caption=tex_caption;
}


void TableHandler::set_tex_supercaption (const string &superkey,
					 const string &tex_supercaption)
{
  Assert(supercolumns.count(superkey), ExcSuperColumnNotExistent(superkey));
  Assert(tex_supercaptions.count(superkey), ExcInternalError());
  tex_supercaptions[superkey]=tex_supercaption;
}


void TableHandler::set_tex_format (const string &key,
				   const string &tex_format)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  Assert(tex_format=="l" || tex_format=="c" || tex_format=="r",
	 ExcUndefinedTexFormat(tex_format));
  columns[key].tex_format=tex_format;
}


void TableHandler::set_precision (const string &key,
				  unsigned int precision)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].precision=precision;
}


void TableHandler::set_scientific (const string &key,
				   bool scientific)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].scientific=scientific;
}


void TableHandler::write_text(ostream &out) const
{
  vector<string> sel_columns;
  get_selected_columns(sel_columns);


// write the caption line
  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      string key=column_order[j];
      const map<string, vector<string> >::const_iterator 
	super_iter=supercolumns.find(key);

      unsigned int max_string_size=6;
      if (super_iter!=supercolumns.end())
	max_string_size+=8*(super_iter->second.size()-1);
	
      if (key.size()>max_string_size)
	key.erase(max_string_size);
      
      out.setf(ios::left);
      out << setw(max_string_size)
//	  << setfill('*') 
	  << key.c_str() << "\t";
    }
  out << endl;

  const unsigned int nrows=n_rows();
  for (unsigned int i=0; i<nrows; ++i)
    {
      const unsigned int n_cols=sel_columns.size();
      
      for (unsigned int j=0; j<n_cols; ++j)
	{
	  string key=sel_columns[j];
					   // avoid `column[key]'
	  const map<string, Column>::const_iterator
	    col_iter=columns.find(key);
	  Assert(col_iter!=columns.end(), ExcInternalError());
	  
	  const Column &column=col_iter->second;
	  
	  out << setprecision(column.precision);

	  if (col_iter->second.scientific)
	    out.setf(ios::scientific, ios::floatfield);
	  else
	    out.setf(ios::fixed, ios::floatfield);
	  column.entries[i]->write_tex(out);
	  
	  if (j<n_cols-1)
	    out << "\t";
	}
      out << endl;
    }
}


void TableHandler::write_tex(ofstream &out) const
{
  out << "\\documentclass[10pt]{report}" << endl
      << "\\usepackage{float}" << endl << endl << endl
      << "\\begin{document}" << endl
      << "\\begin{table}[H]" << endl
      << "\\begin{center}" << endl
      << "\\begin{tabular}{|";

  vector<string> sel_columns;
  get_selected_columns(sel_columns);

				   // write the column formats
  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      string key=column_order[j];
				       // avoid `supercolumns[key]'
      const map<string, vector<string> >::const_iterator 
	super_iter=supercolumns.find(key);

      if (super_iter!=supercolumns.end())
	{
	  const unsigned int n_subcolumns=super_iter->second.size();
	  for (unsigned int k=0; k<n_subcolumns; ++k)
	    {
					       // avoid `columns[supercolumns[key]]'
	      const map<string, Column>::const_iterator
		col_iter=columns.find(super_iter->second[k]);
	      Assert(col_iter!=columns.end(), ExcInternalError());
	      
	      out << col_iter->second.tex_format << "|";
	    }
	}
      else
	{
					   // avoid `columns[key]';
	  const map<string, Column>::const_iterator
	    col_iter=columns.find(key);
	  Assert(col_iter!=columns.end(), ExcInternalError());
	  out << col_iter->second.tex_format << "|";
	}
    }
  out << "} \\hline" << endl;

				   // write the caption line of the table
  
  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      string key=column_order[j];
      const map<string, vector<string> >::const_iterator 
	super_iter=supercolumns.find(key);

      if (super_iter!=supercolumns.end())
	{
	  const unsigned int n_subcolumns=super_iter->second.size();
					   // avoid use of `tex_supercaptions[key]'
	  map<string,string>::const_iterator 
	    tex_super_cap_iter=tex_supercaptions.find(key);
	  out << endl << "\\multicolumn{" << n_subcolumns << "}{|c|}{" 
	      << tex_super_cap_iter->second << "}";
	}
      else
	{
					   // col_iter->second=columns[col];
	  const map<string, Column>::const_iterator
	    col_iter=columns.find(key);
	  Assert(col_iter!=columns.end(), ExcInternalError());
	  out << col_iter->second.tex_caption;
	}
      if (j<column_order.size()-1)
	out << " & ";
    }
  out << "\\\\ \\hline" << endl;

				   // write the n rows
  const unsigned int nrows=n_rows();
  for (unsigned int i=0; i<nrows; ++i)
    {
      const unsigned int n_cols=sel_columns.size();
      
      for (unsigned int j=0; j<n_cols; ++j)
	{
	  string key=sel_columns[j];
					   // avoid `column[key]'
	  const map<string, Column>::const_iterator
	    col_iter=columns.find(key);
	  Assert(col_iter!=columns.end(), ExcInternalError());
	  
	  const Column &column=col_iter->second;
	  
	  out << setprecision(column.precision);

	  if (col_iter->second.scientific)
	    out.setf(ios::scientific, ios::floatfield);
	  else
	    out.setf(ios::fixed, ios::floatfield);
	  
	  column.entries[i]->write_tex(out);
	  
	  if (j<n_cols-1)
	    out << " & ";
	}
      out << "\\\\ \\hline" << endl;
    }

  string caption="table";
  
  out << "\\end{tabular}" << endl
      << "\\end{center}" << endl
//      << "\\caption{"  << caption << "}" << endl
      << "\\end{table}" << endl
      << "\\end{document}" << endl;
}


unsigned int TableHandler::n_rows() const
{
  map<string, Column>::const_iterator col_iter=columns.begin();
  unsigned int n=col_iter->second.entries.size();
  string first_name=col_iter->first;

  for (++col_iter; col_iter!=columns.end(); ++col_iter)
    Assert(col_iter->second.entries.size()==n, 
	   ExcWrongNumberOfDataEntries(col_iter->first, 
				       col_iter->second.entries.size(), 
				       first_name, n));

  return n;
}


void TableHandler::get_selected_columns(vector<string> &sel_columns) const
{
  sel_columns.clear();
  
  for (unsigned int j=0; j<column_order.size(); ++j)
    {
      string key=column_order[j];
      const map<string, vector<string> >::const_iterator 
	super_iter=supercolumns.find(key);

      if (super_iter!=supercolumns.end())
	{
					   // i.e. key is a supercolumn key
	  const unsigned int n_subcolumns=super_iter->second.size();
	  for (unsigned int j=0; j<n_subcolumns; ++j)
	    {
	      const string subkey=super_iter->second[j];
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
  cout << endl;
}


// explicit instantiations
template class TableEntry<double>;
template class TableEntry<float>;
template class TableEntry<int>;
template class TableEntry<unsigned int>;
template class TableEntry<string>;


template
void TableHandler::add_value(const string &, double);

template
void TableHandler::add_value(const string &, int);

template
void TableHandler::add_value(const string &, unsigned int);

template
void TableHandler::add_value(const string &, string);


