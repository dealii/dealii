/* $Id$ */

#include <basic/grid_out.h>


string GridOut::default_suffix (const OutputFormat output_format) 
{
  switch (output_format) 
    {
      case gnuplot: 
	    return ".gnuplot";

      case ucd: 
	    return ".inp";

      case eps: 
	    return ".eps";
	    
      default: 
	    Assert (false, ExcNotImplemented()); 
	    return "";
    };
};
  



GridOut::OutputFormat
GridOut::parse_output_format (const string &format_name)
{
  if (format_name == "ucd")
    return ucd;

  if (format_name == "gnuplot")
    return gnuplot;

  if (format_name == "eps")
    return eps;
  
  AssertThrow (false, ExcInvalidState ());
				   // return something weird
  return OutputFormat(-1);
};



string GridOut::get_output_format_names () {
  return "gnuplot|eps";
};


