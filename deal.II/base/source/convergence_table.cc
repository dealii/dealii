/*----------------------------   convergence_table.cc     ---------------------------*/
/*      $Id$                 */


#include <base/convergence_table.h>
#include <math.h>


ConvergenceTable::ConvergenceTable():
		n_cells_string("n cells"),
		n_dofs_string("n dofs")   
{}


void ConvergenceTable::add_run (unsigned int ncells,
				unsigned int ndofs)
{
  add_value(n_cells_string, ncells);
  add_value(n_dofs_string, ndofs);
}


void ConvergenceTable::evaluate_convergence_rates(const string &key,
						  const RateMode rate_mode)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  
  vector<TableEntryBase *> &entries=columns[key].entries;
  string rate_key=key;

  const unsigned int n=entries.size();
  
  vector<double> values(n);
  for (unsigned int i=0; i<n; ++i)
    {
      if (dynamic_cast<TableEntry<double>*>(entries[i]) != 0)
	values[i]=dynamic_cast<TableEntry<double>*>(entries[i])->value();
      else if (dynamic_cast<TableEntry<float>*>(entries[i]) != 0)
	values[i]=dynamic_cast<TableEntry<float>*>(entries[i])->value();
      else
	Assert(false, ExcWrongValueType());
    }
  
  switch (rate_mode)
    {
      case standard:
	    rate_key+="s";
					     // no value availble for the
					     // first row
	    add_value(rate_key, string("-"));
	    for (unsigned int i=1; i<n; ++i)
	      add_value(rate_key, values[i-1]/values[i]);
	    break;
      case standard_order:
	    rate_key+="so";
					     // no value availble for the
					     // first row
	    add_value(rate_key, string("-"));
	    for (unsigned int i=1; i<n; ++i)
	      add_value(rate_key, log(values[i-1]/values[i])/log(2));
	    break;
      case n_cells:
      case n_dofs:
      case none:
	    break;
      default:
	    ExcNotImplemented();  
    }

  Assert(columns.count(rate_key), ExcInternalError());  
  columns[rate_key].flag=1;
  set_precision(rate_key, 3);

  string superkey="s"+key;
  if (!supercolumns.count(superkey))
    {
      add_column_to_supercolumn(key, superkey);
      set_tex_supercaption(superkey, columns[key].tex_caption);
    }

  add_column_to_supercolumn(rate_key, superkey);
}


void ConvergenceTable::evaluate_convergence_rates(const RateMode rate_mode)
{
				   // make sure that no convergence rates
				   // are evaluated for the n_cells and the
				   // n_dofs column.
  Assert(columns.count(n_cells_string), ExcInternalError());  
  Assert(columns.count(n_dofs_string), ExcInternalError());
  
  columns[n_cells_string].flag=1;
  columns[n_dofs_string].flag=1;
  
  check_n_rows();

  for (map<string, Column>::const_iterator col_iter=columns.begin();
       col_iter!=columns.end(); ++col_iter)
    if (!col_iter->second.flag)
      evaluate_convergence_rates(col_iter->first, rate_mode);
}


/*----------------------------   convergence_table.cc     ---------------------------*/
