//----------------------------  parameter_handler.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  parameter_handler.cc  ---------------------------


//TODO:[WB] (compiler) replace s.c_str() by s when that is possible

#include <base/parameter_handler.h>
#include <base/logstream.h>
#include <base/memory_consumption.h>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <list>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


namespace Patterns
{  

  PatternBase::~PatternBase ()
  {};


  unsigned int
  PatternBase::memory_consumption () const
  {
    if (dynamic_cast<const Integer*>(this) != 0)
      return sizeof(Integer);
    else if (dynamic_cast<const Double*>(this) != 0)
      return sizeof(Double);
    else if (dynamic_cast<const Bool*>(this) != 0)
      return sizeof(Bool);
    else if (dynamic_cast<const Anything*>(this) != 0)
      return sizeof(Anything);
    else
      return sizeof(*this) + 32;
  };
  



  Integer::Integer (const int lower_bound,
		    const int upper_bound) :
		  lower_bound (lower_bound),
		  upper_bound (upper_bound)
  {};



  bool Integer::match (const std::string &test_string) const
  {
#ifdef HAVE_STD_STRINGSTREAM
    std::istringstream str(test_string.c_str());
#else
    std::istrstream str(test_string.c_str());
#endif
    
    int i;
    if (str >> i)
      {
					 // check whether valid bounds
					 // were specified, and if so
					 // enforce their values
	if (lower_bound <= upper_bound)
	  return ((lower_bound <= i) &&
		  (upper_bound >= i));
	else
	  return true;
      };
  
    return false;
  };



  std::string Integer::description () const
  {
				     // check whether valid bounds
				     // were specified, and if so
				     // output their values
    if (lower_bound <= upper_bound)
      {
#ifdef HAVE_STD_STRINGSTREAM
	std::ostringstream description;	
#else
	std::ostrstream description;
#endif
	
	description << "[Integer range "
		    << lower_bound << "..." << upper_bound
		    << " (inclusive)]"
		    << std::ends;
	return description.str();
      }
    else
				       // if no bounds were given, then
				       // return generic string
      return "[Integer]";
  };



  PatternBase *
  Integer::clone () const
  {
    return new Integer(lower_bound, upper_bound);
  };



  Double::Double (const int lower_bound,
		  const int upper_bound) :
		  lower_bound (lower_bound),
		  upper_bound (upper_bound)
  {};



  bool Double::match (const std::string &test_string) const
  {
#ifdef HAVE_STD_STRINGSTREAM
    std::istringstream str(test_string.c_str());
#else
    std::istrstream str(test_string.c_str());
#endif

    double d;
    if (str >> d)
      {
					 // check whether valid bounds
					 // were specified, and if so
					 // enforce their values
	if (lower_bound <= upper_bound)
	  return ((lower_bound <= d) &&
		  (upper_bound >= d));
	else
	  return true;
      };
    return false;
  };



  std::string Double::description () const
  {
				     // check whether valid bounds
				     // were specified, and if so
				     // output their values
    if (lower_bound <= upper_bound)
      {
#ifdef HAVE_STD_STRINGSTREAM
	std::ostringstream description;	
#else
	std::ostrstream description;
#endif

	description << "[Floating point range "
		    << lower_bound << "..." << upper_bound
		    << " (inclusive)]"
		    << std::ends;
	return description.str();
      }
    else
				       // if no bounds were given, then
				       // return generic string
      return "[Double]";
  };


  PatternBase *
  Double::clone () const
  {
    return new Double(lower_bound, upper_bound);
  };



  Selection::Selection (const std::string &seq)
  {
    sequence = seq;

    while (sequence.find(" |") != std::string::npos)
      sequence.replace (sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace (sequence.find("| "), 2, "|");
  };



  bool Selection::match (const std::string &test_string) const
  {
    std::vector<std::string> choices;
    std::string tmp(sequence);
				     // check the different possibilities
    while (tmp.find('|') != std::string::npos) 
      {
	if (test_string == std::string(tmp, 0, tmp.find('|')))
	  return true;
      
	tmp.erase (0, tmp.find('|')+1);
      };
				     // check last choice, not finished by |
    if (test_string == tmp)
      return true;

				     // not found
    return false;
  };



  std::string Selection::description () const
  {
    return sequence;
  };



  PatternBase *
  Selection::clone () const
  {
    return new Selection(sequence);
  };


  unsigned int
  Selection::memory_consumption () const
  {
    return (sizeof(PatternBase) +
	    MemoryConsumption::memory_consumption(sequence));
  };
  

  MultipleSelection::MultipleSelection (const std::string &seq)
  {
    Assert (seq.find (",") == std::string::npos, ExcCommasNotAllowed(seq.find(",")));
  
    sequence = seq;
    while (sequence.find(" |") != std::string::npos)
      sequence.replace (sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace (sequence.find("| "), 2, "|");
  };



  bool MultipleSelection::match (const std::string &test_string_list) const
  {
    std::string tmp = test_string_list;
    std::list<std::string> split_list;

				     // first split the input list
    while (tmp.length() != 0)
      {
	std::string name;
	name = tmp;
      
	if (name.find(",") != std::string::npos)
	  {
	    name.erase (name.find(","), std::string::npos);
	    tmp.erase (0, test_string_list.find(",")+1);
	  }
	else
	  tmp = "";
      
	while ((name.length() != 0) &&
	       (name[0] == ' '))
	  name.erase (0,1);
	while (name[name.length()-1] == ' ')
	  name.erase (name.length()-1, 1);

	split_list.push_back (name);
      };


				     // check the different possibilities
    for (std::list<std::string>::const_iterator test_string = split_list.begin();
	 test_string != split_list.end(); ++test_string) 
      {
	bool string_found = false;
      
	tmp = sequence;
	while (tmp.find('|') != std::string::npos) 
	  {
	    if (*test_string == std::string(tmp, 0, tmp.find('|')))
	      {
						 // string found, quit
						 // loop. don't change
						 // tmp, since we don't
						 // need it anymore.
		string_found = true;
		break;
	      };
	  
	    tmp.erase (0, tmp.find('|')+1);
	  };
					 // check last choice, not finished by |
	if (!string_found)
	  if (*test_string == tmp)
	    string_found = true;

	if (!string_found)
	  return false;
      };

    return true;
  };



  std::string MultipleSelection::description () const
  {
    return sequence;
  };



  PatternBase *
  MultipleSelection::clone () const
  {
    return new MultipleSelection(sequence);
  };


  unsigned int
  MultipleSelection::memory_consumption () const
  {
    return (sizeof(PatternBase) +
	    MemoryConsumption::memory_consumption(sequence));
  };



  Bool::Bool () :
		  Selection ("true|false")
  {};



  PatternBase *
  Bool::clone () const
  {
    return new Bool();
  };



  Anything::Anything ()
  {};



  bool Anything::match (const std::string &) const
  {
    return true;
  };



  std::string Anything::description () const
  {
    return "[Anything]";
  };



  PatternBase *
  Anything::clone () const
  {
    return new Anything();
  };
 
}   // end namespace Patterns



ParameterHandler::ParameterHandler () :
		status(true) {};



ParameterHandler::~ParameterHandler ()
{};



bool ParameterHandler::read_input (std::istream &input)
{
  AssertThrow (input, ExcIO());

  std::string line;
  int lineno=0;
  while (input) 
    {
      ++lineno;
      getline (input, line);
      if (!scan_line (line, lineno)) 
	status = false;
    };

  return status;
};



bool ParameterHandler::read_input (const std::string &filename)
{
  std::ifstream input (filename.c_str());
  if (!input) 
    {
      std::cerr << "ParameterHandler::read_input: could not open file <"
		<< filename << "> for reading." << std::endl
		<< "Trying to make file <"
		<< filename << "> with default values for you." << std::endl;

      std::ofstream output (filename.c_str());
      if (output)
	{
	  print_parameters (output, Text);
	}
      
      return false;
    }
  
  return read_input (input);
}



bool ParameterHandler::read_input_from_string (const char *s)
{
				   // if empty std::string then exit
				   // with success
  if ((s == 0) || ((*s) == 0)) return true;
  
  std::string line;
  std::string input (s);
  int    lineno=0;

				   // if necessary append a newline char
				   // to make all lines equal
  if (input[input.length()-1] != '\n')
    input += '\n';
  
  while (input.size() != 0) 
    {
				       // get one line from Input (=s)
      line.assign (input, 0, input.find('\n'));
				       // delete this part including
				       // the backspace
      input.erase (0, input.find('\n')+1);
      ++lineno;
      
      if (!scan_line (line, lineno)) 
	status = false;
    };

  return status;
};



void ParameterHandler::clear ()
{
  status = true;

  subsection_path.clear ();
  defaults.entries.clear ();
  changed_entries.entries.clear ();

  std::map<std::string, Section*>::iterator p;

  for (p=defaults.subsections.begin(); p!=defaults.subsections.end(); ++p)
    delete p->second;
  for (p=changed_entries.subsections.begin(); p!=changed_entries.subsections.end(); ++p)
    if (p->second) 
      {
	delete p->second;
	Assert (false, ExcInternalError());
      };

  defaults.subsections.clear ();
  changed_entries.subsections.clear ();
};



bool ParameterHandler::declare_entry (const std::string           &entry,
				      const std::string           &default_value,
				      const Patterns::PatternBase &pattern)
{
  Section* p = get_present_defaults_subsection ();

				   // assertions:
				   // entry must not already exist
  Assert (p->entries.find (entry) == p->entries.end(),
	  ExcEntryAlreadyExists (entry));
				   // Default must match Pattern
  Assert (pattern.match (default_value),
	  ExcDefaultDoesNotMatchPattern(default_value,
					pattern.description()));
  
				   // does entry already exist?
  if (p->entries.find (entry) != p->entries.end())
    return false;

  p->entries[entry] = make_pair(default_value, pattern.clone());

				   // check whether default answer matches
				   // the pattern
  if (!pattern.match(default_value))
    return false;

  return true;
};



void ParameterHandler::enter_subsection (const std::string &subsection)
{
  Section* pd = get_present_defaults_subsection ();

				   // does subsection already exist?
  if (pd->subsections.find (subsection) == pd->subsections.end()) 
    {
				       // subsection does not yet exist:
				       // create entry in Defaults and
				       // changed_entries trees
      pd->subsections[subsection] = new Section();

      Section* pc = get_present_changed_subsection ();
      pc->subsections[subsection] = new Section();
    };

				   // finally enter subsection
  subsection_path.push_back (subsection);
};



bool ParameterHandler::leave_subsection ()
{
				   // assert there is a subsection that
				   // we may leave
				   // (use assert since this is a logical
				   // error in a program. When reading input
				   // the scan_line function has to check
				   // whether there is a subsection to be left!)
  Assert (subsection_path.size() != 0, ExcAlreadyAtTopLevel());

  if (subsection_path.size() == 0) 
    return false;

  subsection_path.pop_back ();
  return true;
};



const std::string & ParameterHandler::get (const std::string &entry_string) const
{
  const Section* pd = get_present_defaults_subsection ();
  const Section* pc = get_present_changed_subsection ();

				   // assert that the according entry is already
				   // declared in the defaults tree
  Assert (pd->entries.find (entry_string) != pd->entries.end(),
	  ExcEntryUndeclared(entry_string));
  
  if (pd->entries.find (entry_string) == pd->entries.end()) 
    {
      static const std::string empty_string;
      return empty_string;
    };


// entry exists; now find out whether
				   // it was changed:
  Section::EntryType::const_iterator ptr;
  ptr = pc->entries.find (entry_string);
  if (ptr != pc->entries.end())
    return ptr->second.first;

				   // not changed
  ptr = pd->entries.find (entry_string);
  return ptr->second.first;
};



long int ParameterHandler::get_integer (const std::string &entry_string) const
{
  std::string s = get (entry_string);
  char *endptr;
  long int i = std::strtol (s.c_str(), &endptr, 10);
				   // assert there was no error
  AssertThrow ((s.c_str()!='\0') || (*endptr == '\0'),
	       ExcConversionError(s));

  return i;
};



double ParameterHandler::get_double (const std::string &entry_string) const
{
  std::string s = get (entry_string);
  char *endptr;
  double d = std::strtod (s.c_str(), &endptr);
				   // assert there was no error
  AssertThrow ((s.c_str()!='\0') || (*endptr == '\0'),
	       ExcConversionError(s));

  return d;
};



bool ParameterHandler::get_bool (const std::string &entry_string) const
{
  std::string s = get(entry_string);

  AssertThrow ((s=="true") || (s=="false"), ExcConversionError(s));
  if (s=="true")
    return true;
  else
    return false;
};



std::ostream & ParameterHandler::print_parameters (std::ostream &out,
						   const OutputStyle style)
{
				   // assert that only known formats are
				   // given as "style"
  Assert ((style == Text) || (style == LaTeX), ExcNotImplemented());

  AssertThrow (out, ExcIO());
  
  switch (style) 
    {
      case Text:
	    out << "#Listing of Parameters" << std::endl
		<< "#---------------------" << std::endl;
	    break;
      case LaTeX:
	    out << "\\subsubsection*{Listing of parameters}";
	    out << std::endl << std::endl;
	    out << "\\begin{itemize}"
	        << std::endl;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };
   
				   // dive recursively into the subsections
  print_parameters_section (out, style, 0);

  switch (style) 
    {
      case Text:
	    break;
      case LaTeX:
	    out << "\\end{itemize}" << std::endl;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };
  
  return out;
};


void
ParameterHandler::log_parameters (LogStream &out)
{
  out.push("parameters");
				   // dive recursively into the subsections
  log_parameters_section (out);

  out.pop();
};



void ParameterHandler::print_parameters_section (std::ostream      &out,
						 const OutputStyle  style,
						 const unsigned int indent_level)
{
				   // assert that only known formats are
				   // given as "style"
  Assert ((style == Text) || (style == LaTeX), ExcNotImplemented());

  AssertThrow (out, ExcIO());

  Section *pd = get_present_defaults_subsection ();
  Section *pc = get_present_changed_subsection ();

				   // traverse entry list
  Section::EntryType::const_iterator ptr;

				   // first find out the longest entry name
  unsigned int longest_entry = 0;
  for (ptr = pd->entries.begin(); ptr != pd->entries.end(); ++ptr)
    if (ptr->first.length() > longest_entry)
      longest_entry = ptr->first.length();

				   // print entries one by one
  for (ptr = pd->entries.begin(); ptr != pd->entries.end(); ++ptr)
    {
				       // check whether this entry is listed
				       // in the Changed tree and actually
				       // differs from the default value
      if ((pc->entries.find(ptr->first) != pc->entries.end()) &&
	  (pc->entries[ptr->first].first != pd->entries[ptr->first].first))
	switch (style) 
	  {
	    case Text:
		  out << std::setw(indent_level*2) << ""
		      << "set "
		      << ptr->first
		      << std::setw(longest_entry-ptr->first.length()+3) << " = "
		      << pc->entries[ptr->first].first
		      << "  #"
		      << pd->entries[ptr->first].first
		      << std::endl;
		  break;
	    case LaTeX:
		  out << "\\item {\\bf " << ptr->first << ":} "
		      << pc->entries[ptr->first].first
		      << " ({\\it default:} "
                      << pd->entries[ptr->first].first
		      << ")"
		      << std::endl;
		  break;
	    default:
		  Assert (false, ExcNotImplemented());
	  }
				       // not a changed entry
      else
	switch (style) 
	  {
	    case Text:
		  out << std::setw(indent_level*2) << ""
		      << "set "
		      << ptr->first
		      << std::setw(longest_entry-ptr->first.length()+3) << "= "
		      << ptr->second.first << std::endl;
		  break;
	    case LaTeX:
		  out << "\\item {\\bf " << ptr->first << ":} "
		      << ptr->second.first
		      << std::endl;
		  break;
	    default:
		  Assert (false, ExcNotImplemented());
	  };
    };


				   // now transverse subsections tree
  std::map<std::string, Section*>::const_iterator ptrss;
  for (ptrss = pd->subsections.begin(); ptrss != pd->subsections.end(); ++ptrss)
    {
      switch (style) 
	{
	  case Text:
		out << std::setw(indent_level*2) << ""
		    << "subsection " << ptrss->first << std::endl;
		break;
	  case LaTeX:
		out << std::endl
		    << "\\item {\\bf "
		    << "Subsection " << ptrss->first
		    << "}" << std::endl
		    << "\\begin{itemize}"
		    << std::endl;
		break;
	  default:
		Assert (false, ExcNotImplemented());
	};
      enter_subsection (ptrss->first);
      print_parameters_section (out, style, indent_level+1);
      leave_subsection ();
      switch (style) 
	{
	  case Text:
						 // write end of subsection. one
						 // blank line after each subsection
		out << std::setw(indent_level*2) << ""
		    << "end" << std::endl
		    << std::endl;

						 // if this is a toplevel
						 // subsection, then have two
						 // newlines
		if (indent_level == 0)
		  out << std::endl;
		
		break;
	  case LaTeX:
		out << "\\end{itemize}"
		    << std::endl;
		break;
	  default:
		Assert (false, ExcNotImplemented());
	};
    };
};



void ParameterHandler::log_parameters_section (LogStream &out)
{
  Section *pd = get_present_defaults_subsection ();
  Section *pc = get_present_changed_subsection ();

				   // traverse entry list
  Section::EntryType::const_iterator ptr;

				   // print entries one by one
  for (ptr = pd->entries.begin(); ptr != pd->entries.end(); ++ptr)
    {
      if ((pc->entries.find(ptr->first) != pc->entries.end()) &&
	  (pc->entries[ptr->first].first != pd->entries[ptr->first].first))
					 // check whether this entry is listed
					 // in the Changed tree and actually
					 // differs from the default value
	out << ptr->first << ": "
	    << pc->entries[ptr->first].first << std::endl;
      else
	out << ptr->first << ": "
	    << ptr->second.first << std::endl;
    };


				   // now transverse subsections tree
  std::map<std::string, Section*>::const_iterator ptrss;
  for (ptrss = pd->subsections.begin(); ptrss != pd->subsections.end(); ++ptrss)
    {
      out.push(ptrss->first);
      enter_subsection (ptrss->first);
      log_parameters_section (out);
      leave_subsection ();
      out.pop();
    };
};



bool ParameterHandler::scan_line (std::string        line,
				  const unsigned int lineno)
{
				   // if there is a comment, delete it
  if (line.find('#') != std::string::npos)
    line.erase (line.find("#"), std::string::npos);
				   // replace every whitespace sequence
				   // by " "
  while (line.find('\t') != std::string::npos)
    line.replace (line.find('\t'), 1, " ");
  while (line.find("  ") != std::string::npos)
    line.erase (line.find("  "), 1);
  				   // now every existing whitespace
				   // should be exactly one ' ';
				   // if at end or beginning: delete
  if ((line.length() != 0) && (line[0] == ' '))  line.erase (0, 1);
				   // if line is now empty: leave
  if (line.length() == 0) return true;

  if (line[line.length()-1] == ' ')
    line.erase (line.size()-1, 1);

				   // enter subsection
  if ((line.find ("SUBSECTION ") == 0) ||
      (line.find ("subsection ") == 0))
    {
				       // delete this prefix
      line.erase (0, std::string("subsection").length()+1);

      std::string SecName = line;
      Section* pc = get_present_changed_subsection ();
				       // check whether subsection exists
      if (pc->subsections.find(SecName) == pc->subsections.end()) 
	{
	  std::cerr << "Line " << lineno << ": There is no such subsection to be entered:" << std::endl;
	  for (unsigned int i=0; i<subsection_path.size(); ++i)
	    std::cerr << std::setw(i*2+4) << " "
		      << "subsection " << subsection_path[i] << std::endl;
	  std::cerr << std::setw(subsection_path.size()*2+4) << " "
		    << "subsection " << SecName << std::endl;
	  return false;
	};

				       // subsection exists
      subsection_path.push_back (SecName);
      return true;
    };
  
				   // exit subsection
  if ((line.find ("END") == 0) ||
      (line.find ("end") == 0))
    if (subsection_path.size() == 0) 
      {
	std::cerr << "Line " << lineno
		  << ": There is no subsection to leave here!" << std::endl;
	return false;
      }
    else
      return leave_subsection ();

				   // regular entry
  if ((line.find ("SET ") == 0) ||
      (line.find ("set ") == 0)) 
    {
				       // erase "set" statement and eliminate
				       // spaces around the '='
      line.erase (0, 4);
      if (line.find(" =") != std::string::npos)
	line.replace (line.find(" ="), 2, "=");
      if (line.find("= ") != std::string::npos)
	line.replace (line.find("= "), 2, "=");

				       // extract entry name and value
      std::string entry_name  (line, 0, line.find('='));
      std::string entry_value (line, line.find('=')+1, std::string::npos);

      Section* pd = get_present_defaults_subsection ();

				       // check whether entry was declared
      if (pd->entries.find(entry_name) == pd->entries.end()) 
	{
	  std::cerr << "Line " << lineno
		    << ": No such entry was declared:" << std::endl
		    << "    " << entry_name << std::endl
		    << "    <Present subsection:" << std::endl;
	  for (unsigned int i=0; i<subsection_path.size(); ++i)
	    std::cerr << std::setw(i*2+8) << " "
		      << "subsection " << subsection_path[i] << std::endl;
	  std::cerr << "    >" << std::endl;

	  return false;
	};

				       // if entry was declared:
				       // does it match the regex? if not,
				       // don't enter it into the database
				       // exception: if it contains characters
				       // which specify it as a multiple loop
				       // entry, then ignore content
      if (entry_value.find ('{') == std::string::npos)
	if (!pd->entries[entry_name].second->match(entry_value))
	  {
	    std::cerr << "Line " << lineno << ":" << std::endl
		      << "    The entry value" << std::endl
		      << "        " << entry_value << std::endl
		      << "    for the entry named" << std::endl
		      << "        " << entry_name << std::endl
		      << "    does not match the given pattern" << std::endl
		      << "        " << pd->entries[entry_name].second->description() << std::endl;
	    return false;
	  };
      
      Section* pc = get_present_changed_subsection ();
				       // the following line declares this entry
				       // if not yet existent and overwrites it
				       // otherwise (the pattern is set to a null
				       // pointer, since we don't need the
				       // pattern in the entries section -- only
				       // in the defaults section)
      pc->entries[entry_name] = make_pair(entry_value,
					  static_cast<Patterns::PatternBase*>(0));

      return true;
    };

				   // this line matched nothing known
  std::cerr << "Line " << lineno << ": This line matched nothing known:" << std::endl
	    << "    " << line << std::endl;
  return false;
};



ParameterHandler::Section* ParameterHandler::get_present_defaults_subsection ()
{
  Section* sec = &defaults;
  std::vector<std::string>::const_iterator SecName = subsection_path.begin();
    
  while (SecName != subsection_path.end()) 
    {
      sec = sec->subsections[*SecName];
      ++SecName;
    };

  return sec;
};



const ParameterHandler::Section* ParameterHandler::get_present_defaults_subsection () const
{
  Section* sec = const_cast<Section*>(&defaults); // not nice, but needs to be and
				   // after all: we do not change @p{sec}
  std::vector<std::string>::const_iterator SecName = subsection_path.begin();
    
  while (SecName != subsection_path.end()) 
    {
      sec = sec->subsections[*SecName];
      ++SecName;
    };

  return sec;
};



ParameterHandler::Section* ParameterHandler::get_present_changed_subsection ()
{
  Section* sec = &changed_entries;
  std::vector<std::string>::iterator SecName = subsection_path.begin();
    
  while (SecName != subsection_path.end()) 
    {
      sec = sec->subsections[*SecName];
      ++SecName;
    };

  return sec;
};



const ParameterHandler::Section* ParameterHandler::get_present_changed_subsection () const
{
  Section* sec = const_cast<Section*>(&changed_entries); // same as in get_present_default_s...
  std::vector<std::string>::const_iterator SecName = subsection_path.begin();
    
  while (SecName != subsection_path.end()) 
    {
      sec = sec->subsections[*SecName];
      ++SecName;
    };

  return sec;
};



unsigned int 
ParameterHandler::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (status) +
	  MemoryConsumption::memory_consumption (subsection_path) +
	  MemoryConsumption::memory_consumption (defaults) +
	  MemoryConsumption::memory_consumption (changed_entries));
};



ParameterHandler::Section::~Section () 
{
				   // first release the memory pointed
				   // to by the second component of
				   // the pair, since we became owner
				   // of that memory through the
				   // clone() call
  for (EntryType::iterator q=entries.begin(); q!=entries.end(); ++q)
    delete q->second.second;
				   // then clear entire map
  entries.clear ();

  std::map<std::string, Section*>::iterator p;

  for (p=subsections.begin(); p!=subsections.end(); ++p)
    delete p->second;

  subsections.clear ();
};



unsigned int 
ParameterHandler::Section::memory_consumption () const
{
  unsigned int mem = 0;
  for (EntryType::const_iterator i=entries.begin(); i!=entries.end(); ++i)
    mem += (MemoryConsumption::memory_consumption (i->first) +
	    MemoryConsumption::memory_consumption (i->second.first) +
	    MemoryConsumption::memory_consumption (*i->second.second));
  for (std::map<std::string, Section*>::const_iterator i=subsections.begin(); i!=subsections.end(); ++i)
    mem += (MemoryConsumption::memory_consumption (i->first) +
	    MemoryConsumption::memory_consumption (*(i->second)));
  return mem;
};




MultipleParameterLoop::MultipleParameterLoop() :
		n_branches(0)
{};



MultipleParameterLoop::~MultipleParameterLoop ()
{};



bool MultipleParameterLoop::read_input (std::istream &input)
{
  AssertThrow (input, ExcIO());

  bool x = ParameterHandler::read_input (input);
  if (x)
    init_branches ();
  return x;
};



bool MultipleParameterLoop::read_input (const std::string &filename)
{
  return ParameterHandler::read_input (filename);
				   // don't call init_branches, since this read_input
				   // function calls
				   // MultipleParameterLoop::Readinput(std::istream &, std::ostream &)
				   // which itself calls init_branches.
};



bool MultipleParameterLoop::read_input_from_string (const char *s)
{
  bool x = ParameterHandler::read_input (s);
  init_branches ();
  return x;
};



void MultipleParameterLoop::loop (MultipleParameterLoop::UserClass &uc)
{
  for (unsigned int run_no=0; run_no<n_branches; ++run_no) 
    {
				       // give create_new one-based numbers
      uc.create_new (run_no+1);
      fill_entry_values (run_no);
      uc.run (*this);
    };
};



void MultipleParameterLoop::init_branches ()
{
  multiple_choices.clear ();
  
  ParameterHandler::Section *sec;
				   // first check the defaults entries whether it
				   // contains multiple choices
  sec = &defaults;
  init_branches_section (*sec);

				   // then check changed entries
  sec = &changed_entries;
  init_branches_section (*sec);

				   // split up different values
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    multiple_choices[i].split_different_values ();
  
				   // check whether we have included a multiple
				   // choice entry from the defaults section which
				   // has only one single value after reading the
				   // input
  if (multiple_choices.size() > 0)
    for (std::vector<Entry>::iterator i=multiple_choices.end()-1;
	 i >= multiple_choices.begin(); --i)
      if (i->different_values.size() == 1)
	multiple_choices.erase (i);

				   // finally calculate number of branches
  n_branches = 1;
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    if (multiple_choices[i].type == Entry::variant)
      n_branches *= multiple_choices[i].different_values.size();

				   // check whether array entries have the correct
				   // number of entries
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    if (multiple_choices[i].type == Entry::array)
      if (multiple_choices[i].different_values.size() != n_branches)
	std::cerr << "    The entry value" << std::endl
		  << "        " << multiple_choices[i].entry_value << std::endl
		  << "    for the entry named" << std::endl
		  << "        " << multiple_choices[i].entry_name << std::endl
		  << "    does not have the right number of entries for the " << std::endl
		  << "        " << n_branches << " variant runs that will be performed." << std::endl;


				   // do a first run on filling the values to
				   // check for the conformance with the regexp
				   // (afterwards, this will be lost in the whole
				   // other output)
  for (unsigned int i=0; i<n_branches; ++i) 
    fill_entry_values (i);
};



void MultipleParameterLoop::init_branches_section (const ParameterHandler::Section &sec)
{
				   // check all entries in the present subsection
				   // whether it is a multiple entry
  Section::EntryType::const_iterator e;
  for (e = sec.entries.begin(); e != sec.entries.end(); ++e) 
    if (e->second.first.find('{') != std::string::npos) 
      multiple_choices.push_back (Entry(subsection_path,
					e->first,
					e->second.first));


				   // transverse subsections
  std::map<std::string, Section*>::const_iterator s;
  for (s = sec.subsections.begin(); s != sec.subsections.end(); ++s) 
    {
      enter_subsection (s->first);
      init_branches_section (*s->second);
      leave_subsection ();
    };
};



void MultipleParameterLoop::fill_entry_values (const unsigned int run_no)
{
  int possibilities = 1;
  
  std::vector<Entry>::iterator choice;
  for (choice = multiple_choices.begin();
       choice != multiple_choices.end();
       ++choice) 
    {
				       // temporarily enter the subsection tree
				       // of this multiple entry
      subsection_path.swap (choice->subsection_path);

				       // set entry
      Section* pd = get_present_defaults_subsection ();
      int selection = (run_no/possibilities) % choice->different_values.size();
      std::string entry_value;
      if (choice->type == Entry::variant)
	entry_value = choice->different_values[selection];
      else 
	{
	  if (run_no>=choice->different_values.size()) 
	    {
	      std::cerr << "The given array for entry "
			<< pd->entries[choice->entry_name].first
			<< " does not conatin enough elements! Taking empty string instead." << std::endl;
	      entry_value = "";
	    }
	  else
	    entry_value = choice->different_values[run_no];
	};
      
				       // check conformance with regex
      if (!pd->entries[choice->entry_name].second->match(entry_value))
	{
	  std::cerr << "In run no.  " << run_no+1 << ":" << std::endl
		    << "    The entry value" << std::endl
		    << "        " << entry_value << std::endl
		    << "    for the entry named" << std::endl
		    << "        " << choice->entry_name << std::endl
		    << "    does not match the given pattern" << std::endl
		    << "        " << pd->entries[choice->entry_name].second->description() << std::endl
		    << "    Taking default value" << std::endl
		    << "        " << pd->entries[choice->entry_name].first << std::endl;
	  
					   // select default instead
	  entry_value = pd->entries[choice->entry_name].first;
	};

      Section* pc = get_present_changed_subsection ();
				       // the following line declares this entry
				       // if not yet existent and overwrites it
				       // otherwise (the pattern is set to a null
				       // pointer, since we don't need the
				       // pattern in the entries section -- only
				       // in the defaults section)
      pc->entries[choice->entry_name] = make_pair(entry_value,
						  static_cast<Patterns::PatternBase*>(0));
            
				       // get out of subsection again
      subsection_path.swap (choice->subsection_path);

				       // move ahead if it was a variant entry
      if (choice->type == Entry::variant)
	possibilities *= choice->different_values.size();
    };
  
};




unsigned int 
MultipleParameterLoop::memory_consumption () const
{
  unsigned int mem = ParameterHandler::memory_consumption ();
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    mem += multiple_choices[i].memory_consumption ();
  
  return mem;
};



MultipleParameterLoop::Entry::Entry (const std::vector<std::string> &ssp,
				     const std::string              &Name,
				     const std::string              &Value) :
		subsection_path (ssp), entry_name(Name), entry_value(Value)
{};



void MultipleParameterLoop::Entry::split_different_values ()
{
				   // split string into three parts:
				   // part before the opening "{",
				   // the selection itself, final
				   // part after "}"
  std::string prefix  (entry_value, 0, entry_value.find('{'));
  std::string multiple(entry_value, entry_value.find('{')+1,
		       entry_value.rfind('}')-entry_value.find('{')-1);
  std::string postfix (entry_value, entry_value.rfind('}')+1, std::string::npos);
				   // if array entry {{..}}: delete inner
				   // pair of braces
  if (multiple[0]=='{')
    multiple.erase (0,1);
  if (multiple[multiple.size()-1] == '}')
    multiple.erase (multiple.size()-1, 1);
				   // erase leading and trailing spaces
				   // in multiple
  while (multiple[0] == ' ') multiple.erase (0,1);
  while (multiple[multiple.size()-1] == ' ') multiple.erase (multiple.size()-1,1);
  
				   // delete spaces around '|'
  while (multiple.find(" |") != std::string::npos)
    multiple.replace (multiple.find(" |"), 2, "|");
  while (multiple.find("| ") != std::string::npos)
    multiple.replace (multiple.find("| "), 2, "|");

  while (multiple.find('|') != std::string::npos) 
    {
      different_values.push_back (prefix +
				  std::string(multiple, 0, multiple.find('|'))+
				  postfix);
      multiple.erase (0, multiple.find('|')+1);
    };
				   // make up the last selection ("while" broke
				   // because there was no '|' any more
  different_values.push_back (prefix+multiple+postfix);
				   // finally check whether this was a variant
				   // entry ({...}) or an array ({{...}})
  if ((entry_value.find("{{") != std::string::npos) &&
      (entry_value.find("}}") != std::string::npos))
    type = Entry::array;
  else
    type = Entry::variant;
};


unsigned int 
MultipleParameterLoop::Entry::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (subsection_path) +
	  MemoryConsumption::memory_consumption (entry_name) +
	  MemoryConsumption::memory_consumption (entry_value) +
	  MemoryConsumption::memory_consumption (different_values) +
	  sizeof (type));
};

