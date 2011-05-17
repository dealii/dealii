//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <sstream>
#include <cctype>

#ifdef HAVE_STD_NUMERIC_LIMITS
# include <limits>
#endif

DEAL_II_NAMESPACE_OPEN



//TODO[WB]: various functions here could be simplified by using namespace Utilities

namespace Patterns
{

  namespace
  {
				     /**
				      * Read to the end of the stream and
				      * return whether all there is is
				      * whitespace or whether there are other
				      * characters as well.
				      */
    bool has_only_whitespace (std::istream &in)
    {
      while (in)
	{
	  char c;

					   // skip if we've reached the end of
					   // the line
	  if (!(in >> c))
	    break;

	  if ((c != ' ') && (c != '\t'))
	    return false;
	}
      return true;
    }
  }



  PatternBase* pattern_factory (const std::string& description)
  {
    PatternBase* p;

    p = Integer::create(description);
    if (p != 0)
      return p;

    p = Double::create(description);
    if (p !=0 )
      return p;

    p = Selection::create(description);
    if (p !=0 )
      return p;

    p = List::create(description);
    if (p !=0 )
      return p;

    p = MultipleSelection::create(description);
    if (p !=0 )
      return p;

    p = Bool::create(description);
    if (p!=0 )
      return p;

    p = Anything::create(description);
    if (p !=0 )
      return p;

    p = FileName::create(description);
    if (p !=0 )
      return p;

    p = DirectoryName::create(description);
    if (p!=0 )
      return p;

    Assert(false, ExcNotImplemented());

    return 0;
  }



  PatternBase::~PatternBase ()
  {}


  std::size_t
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
  }



  const int Integer::min_int_value =
#ifdef HAVE_STD_NUMERIC_LIMITS
          std::numeric_limits<int>::min();
#else
          1;
#endif

  const int Integer::max_int_value =
#ifdef HAVE_STD_NUMERIC_LIMITS
          std::numeric_limits<int>::max();
#else
          0;
#endif


  const char* Integer::description_init = "[Integer";

  Integer::Integer (const int lower_bound,
		    const int upper_bound)
                  :
		  lower_bound (lower_bound),
		  upper_bound (upper_bound)
  {}



  bool Integer::match (const std::string &test_string) const
  {
    std::istringstream str(test_string);

    int i;
    if (!(str >> i))
      return false;

    if (!has_only_whitespace (str))
      return false;
				     // check whether valid bounds
				     // were specified, and if so
				     // enforce their values
    if (lower_bound <= upper_bound)
      return ((lower_bound <= i) &&
	      (upper_bound >= i));
    else
      return true;
  }



  std::string Integer::description () const
  {
				     // check whether valid bounds
				     // were specified, and if so
				     // output their values
    if (lower_bound <= upper_bound)
      {
	std::ostringstream description;

	description << description_init
	        <<" range "
		    << lower_bound << "..." << upper_bound
		    << " (inclusive)]";
	return description.str();
      }
    else
				       // if no bounds were given, then
				       // return generic string
      return "[Integer]";
  }



  PatternBase *
  Integer::clone () const
  {
    return new Integer(lower_bound, upper_bound);
  }



  Integer* Integer::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
    {
      std::istringstream is(description);

      if(is.str().size() > strlen(description_init) + 1)
      {
//TODO: verify that description matches the pattern "^\[Integer range \d+\.\.\.\d+\]$"
        int lower_bound, upper_bound;

        is.ignore(strlen(description_init) + strlen(" range "));

        if(!(is >> lower_bound))
          return new Integer();

        is.ignore(strlen("..."));

        if(!(is >> upper_bound))
          return new Integer();

        return new Integer(lower_bound, upper_bound);
      }
      else
        return new Integer();
    }
    else
      return 0;
  }



  const double Double::min_double_value =
#ifdef HAVE_STD_NUMERIC_LIMITS
          -std::numeric_limits<double>::max();
#else
          1;
#endif

  const double Double::max_double_value =
#ifdef HAVE_STD_NUMERIC_LIMITS
          std::numeric_limits<double>::max();
#else
          0;
#endif

  const char* Double::description_init = "[Double";

  Double::Double (const double lower_bound,
		  const double upper_bound)
                  :
		  lower_bound (lower_bound),
		  upper_bound (upper_bound)
  {}



  bool Double::match (const std::string &test_string) const
  {
    std::istringstream str(test_string);

    double d;
    if (!(str >> d))
      return false;

    if (!has_only_whitespace (str))
      return false;
				     // check whether valid bounds
				     // were specified, and if so
				     // enforce their values
    if (lower_bound <= upper_bound)
      return ((lower_bound <= d) &&
	      (upper_bound >= d));
    else
      return true;
  }



  std::string Double::description () const
  {
	std::ostringstream description;

				     // check whether valid bounds
				     // were specified, and if so
				     // output their values
    if (lower_bound <= upper_bound)
      {
	description << description_init
	        << " "
		    << lower_bound << "..." << upper_bound
		    << " (inclusive)]";
	return description.str();
      }
    else
				       // if no bounds were given, then
				       // return generic string
      {
	description << description_init
		    << "]";
	return description.str();
      }
  }


  PatternBase *
  Double::clone () const
  {
    return new Double(lower_bound, upper_bound);
  }



  Double* Double::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
    {
      std::istringstream is(description);

      if(is.str().size() > strlen(description_init) + 1)
      {
        double lower_bound, upper_bound;

        is.ignore(strlen(description_init) + strlen(" range "));

        if(!(is >> lower_bound))
          return new Double();

        is.ignore(strlen("..."));

        if(!(is >> upper_bound))
          return new Double();

        return new Double(lower_bound, upper_bound);
      }
      else
        return new Double();
    }
    else
      return 0;
  }



  const char* Selection::description_init = "[Selection";


  Selection::Selection (const std::string &seq)
  {
    sequence = seq;

    while (sequence.find(" |") != std::string::npos)
      sequence.replace (sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace (sequence.find("| "), 2, "|");
  }



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
  }



  std::string Selection::description () const
  {
    std::ostringstream description;

    description << description_init
                << " "
                << sequence
                << " ]";

    return description.str();
  }



  PatternBase *
  Selection::clone () const
  {
    return new Selection(sequence);
  }


  std::size_t
  Selection::memory_consumption () const
  {
    return (sizeof(PatternBase) +
	    MemoryConsumption::memory_consumption(sequence));
  }



  Selection* Selection::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
    {
      std::string sequence(description);

      sequence.erase(0, std::strlen(description_init) + 1);
      sequence.erase(sequence.length()-2, 2);

      return new Selection(sequence);
    }
    else
      return 0;
  }



  const unsigned int List::max_int_value =
#ifdef HAVE_STD_NUMERIC_LIMITS
          std::numeric_limits<unsigned int>::max();
#else
          numbers::invalid_unsigned_int;
#endif


  const char* List::description_init = "[List";


  List::List (const PatternBase  &p,
              const unsigned int  min_elements,
              const unsigned int  max_elements)
                  :
                  pattern (p.clone()),
                  min_elements (min_elements),
                  max_elements (max_elements)
  {
    Assert (min_elements <= max_elements,
            ExcInvalidRange (min_elements, max_elements));
  }



  List::~List ()
  {
    delete pattern;
    pattern = 0;
  }



  bool List::match (const std::string &test_string_list) const
  {
    std::string tmp = test_string_list;
    std::vector<std::string> split_list;
    split_list.reserve (std::count (tmp.begin(), tmp.end(), ',')+1);

				     // first split the input list
    while (tmp.length() != 0)
      {
        std::string name;
	name = tmp;

	if (name.find(",") != std::string::npos)
	  {
	    name.erase (name.find(","), std::string::npos);
	    tmp.erase (0, tmp.find(",")+1);
	  }
	else
	  tmp = "";

	while ((name.length() != 0) &&
	       (std::isspace (name[0])))
	  name.erase (0,1);

	while (std::isspace (name[name.length()-1]))
	  name.erase (name.length()-1, 1);

	split_list.push_back (name);
      };

    if ((split_list.size() < min_elements) ||
        (split_list.size() > max_elements))
      return false;

				     // check the different possibilities
    for (std::vector<std::string>::const_iterator
           test_string = split_list.begin();
	 test_string != split_list.end(); ++test_string)
      if (pattern->match (*test_string) == false)
        return false;

    return true;
  }



  std::string List::description () const
  {
    std::ostringstream description;

    description << description_init
                << " list of <" << pattern->description() << ">"
                << " of length " << min_elements << "..." << max_elements
                << " (inclusive)"
                << "]";

    return description.str();
  }



  PatternBase *
  List::clone () const
  {
    return new List(*pattern, min_elements, max_elements);
  }


  std::size_t
  List::memory_consumption () const
  {
    return (sizeof(PatternBase) +
	    MemoryConsumption::memory_consumption(*pattern));
  }



  List* List::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
    {
      int min_elements, max_elements;

      std::istringstream is(description);
      is.ignore(strlen(description_init) + strlen(" list of <"));

      std::auto_ptr<char> new_description (new char[is.str().size() + 1]);
      is.getline(&(*new_description), is.str().size(), '>');
      std::string str(&(*new_description));

      std::auto_ptr<PatternBase> base_pattern (pattern_factory(str));

      is.ignore(strlen(" of length "));
      if(!(is >> min_elements))
        return new List(*base_pattern);

      is.ignore(strlen("..."));
      if(!(is >> max_elements))
        return new List(*base_pattern);

      return new List(*base_pattern, min_elements, max_elements);
    }
    else
      return 0;
  }



  const char* MultipleSelection::description_init = "[MultipleSelection";


  MultipleSelection::MultipleSelection (const std::string &seq)
  {
    Assert (seq.find (",") == std::string::npos, ExcCommasNotAllowed(seq.find(",")));

    sequence = seq;
    while (sequence.find(" |") != std::string::npos)
      sequence.replace (sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace (sequence.find("| "), 2, "|");
  }



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
	    tmp.erase (0, tmp.find(",")+1);
	  }
	else
	  tmp = "";

	while ((name.length() != 0) &&
	       (std::isspace (name[0])))
	  name.erase (0,1);
	while (std::isspace (name[name.length()-1]))
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
  }



  std::string MultipleSelection::description () const
  {
    std::ostringstream description;

    description << description_init
                << " "
                << sequence
                << " ]";

    return description.str();
  }



  PatternBase *
  MultipleSelection::clone () const
  {
    return new MultipleSelection(sequence);
  }


  std::size_t
  MultipleSelection::memory_consumption () const
  {
    return (sizeof(PatternBase) +
	    MemoryConsumption::memory_consumption(sequence));
  }



  MultipleSelection* MultipleSelection::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
    {
      std::string sequence(description);

      sequence.erase(0, std::strlen(description_init) + 1);
      sequence.erase(sequence.length()-2, 2);

      return new MultipleSelection(sequence);
    }
    else
      return 0;
  }



  const char* Bool::description_init = "[Bool";


  Bool::Bool ()
                  :
		  Selection ("true|false")
  {}



  std::string Bool::description () const
  {
    std::ostringstream description;

    description << description_init
                << "]";

    return description.str();
  }



  PatternBase *
  Bool::clone () const
  {
    return new Bool();
  }



  Bool* Bool::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
      return new Bool();
    else
      return 0;
  }



  const char* Anything::description_init = "[Anything";


  Anything::Anything ()
  {}



  bool Anything::match (const std::string &) const
  {
    return true;
  }



  std::string Anything::description () const
  {
    std::ostringstream description;

    description << description_init
                << "]";

    return description.str();
  }



  PatternBase *
  Anything::clone () const
  {
    return new Anything();
  }



  Anything* Anything::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
      return new Anything();
    else
      return 0;
  }



  const char* FileName::description_init = "[FileName";


  FileName::FileName (const FileType type)
          : file_type (type)
  {}



  bool FileName::match (const std::string &) const
  {
    return true;
  }



  std::string FileName::description () const
  {
    std::ostringstream description;

    description << description_init;

    if (file_type == input)
      description << " (Type: input)]";
    else
      description << " (Type: output)]";

    return description.str();
  }



  PatternBase *
  FileName::clone () const
  {
    return new FileName(file_type);
  }



  FileName* FileName::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
    {
      std::istringstream is(description);
      std::string file_type;
      FileType type;

      is.ignore(strlen(description_init) + strlen(" (Type:"));

      is >> file_type;

      if (file_type == "input)]")
        type = input;
      else
        type = output;

      return new FileName(type);
    }
    else
      return 0;
  }



  const char* DirectoryName::description_init = "[DirectoryName";


  DirectoryName::DirectoryName ()
  {}



  bool DirectoryName::match (const std::string &) const
  {
    return true;
  }



  std::string DirectoryName::description () const
  {
    std::ostringstream description;

    description << description_init << "]";

    return description.str();
  }



  PatternBase *
  DirectoryName::clone () const
  {
    return new DirectoryName();
  }



  DirectoryName* DirectoryName::create (const std::string& description)
  {
    if(description.compare(0, std::strlen(description_init), description_init) == 0)
      return new DirectoryName();
    else
      return 0;
  }

}   // end namespace Patterns



ParameterHandler::ParameterHandler ()
		:
		entries (new boost::property_tree::ptree())
{}



ParameterHandler::~ParameterHandler ()
{}



std::string
ParameterHandler::mangle (const std::string &s)
{
  std::string u;
  u.reserve (s.size());

  static const std::string allowed_characters
    ("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

				   // for all parts of the string, see
				   // if it is an allowed character or
				   // not
  for (unsigned int i=0; i<s.size(); ++i)
    if (allowed_characters.find (s[i]) != std::string::npos)
      u.push_back (s[i]);
    else
      {
	u.push_back ('_');
	static const char hex[16]
	  = { '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
	u.push_back (hex[static_cast<unsigned char>(s[i])/16]);
	u.push_back (hex[static_cast<unsigned char>(s[i])%16]);
      }

  return u;
}



std::string
ParameterHandler::demangle (const std::string &s)
{
  std::string u;
  u.reserve (s.size());

  for (unsigned int i=0; i<s.size(); ++i)
    if (s[i] != '_')
      u.push_back (s[i]);
    else
      {
	Assert (i+2 < s.size(),
		ExcMessage ("Trying to demangle an invalid string."));

	unsigned char c = 0;
	switch (s[i+1])
	  {
	    case '0':  c = 0 * 16;  break;
	    case '1':  c = 1 * 16;  break;
	    case '2':  c = 2 * 16;  break;
	    case '3':  c = 3 * 16;  break;
	    case '4':  c = 4 * 16;  break;
	    case '5':  c = 5 * 16;  break;
	    case '6':  c = 6 * 16;  break;
	    case '7':  c = 7 * 16;  break;
	    case '8':  c = 8 * 16;  break;
	    case '9':  c = 9 * 16;  break;
	    case 'a':  c = 10 * 16;  break;
	    case 'b':  c = 11 * 16;  break;
	    case 'c':  c = 12 * 16;  break;
	    case 'd':  c = 13 * 16;  break;
	    case 'e':  c = 14 * 16;  break;
	    case 'f':  c = 15 * 16;  break;
	    default:
		  Assert (false, ExcInternalError());
	  }
	switch (s[i+2])
	  {
	    case '0':  c += 0;  break;
	    case '1':  c += 1;  break;
	    case '2':  c += 2;  break;
	    case '3':  c += 3;  break;
	    case '4':  c += 4;  break;
	    case '5':  c += 5;  break;
	    case '6':  c += 6;  break;
	    case '7':  c += 7;  break;
	    case '8':  c += 8;  break;
	    case '9':  c += 9;  break;
	    case 'a':  c += 10;  break;
	    case 'b':  c += 11;  break;
	    case 'c':  c += 12;  break;
	    case 'd':  c += 13;  break;
	    case 'e':  c += 14;  break;
	    case 'f':  c += 15;  break;
	    default:
		  Assert (false, ExcInternalError());
	  }

	u.push_back (static_cast<char>(c));

					 // skip the two characters
	i += 2;
      }

  return u;
}



bool
ParameterHandler::is_parameter_node (const boost::property_tree::ptree &p)
{
  return (p.get_optional<std::string>("value"));
}



std::string
ParameterHandler::get_current_path () const
{
  if (subsection_path.size() > 0)
    {
      std::string p = mangle(subsection_path[0]);
      for (unsigned int i=1; i<subsection_path.size(); ++i)
	{
	  p += path_separator;
	  p += mangle(subsection_path[i]);
	}
      return p;
    }
  else
    return "";
}



std::string
ParameterHandler::get_current_full_path (const std::string &name) const
{
  std::string path = get_current_path ();
  if (path.empty() == false)
    path += path_separator;

  path += mangle(name);

  return path;
}



bool ParameterHandler::read_input (std::istream &input)
{
  AssertThrow (input, ExcIO());

  std::string line;
  int lineno=0;
  bool status = true;

  while (input)
    {
      ++lineno;
      getline (input, line);
      if (!scan_line (line, lineno))
	status = false;
    }

  return status;
}



bool ParameterHandler::read_input (const std::string &filename,
				   const bool optional,
				   const bool write_compact)
{
  PathSearch search("PARAMETERS");

  try
    {
      std::string openname = search.find(filename);
      std::ifstream input (openname.c_str());
      AssertThrow(input, ExcIO());

      return read_input (input);
    }
  catch (const PathSearch::ExcFileNotFound&)
    {
      std::cerr << "ParameterHandler::read_input: could not open file <"
		<< filename << "> for reading." << std::endl;
      if (!optional)
	{
	  std:: cerr << "Trying to make file <"
		     << filename << "> with default values for you." << std::endl;
	  std::ofstream output (filename.c_str());
	  if (output)
	    print_parameters (output, (write_compact ? ShortText : Text));
	}
    }
  return false;
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

  bool status = true;
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
    }

  return status;
}



namespace
{
				   // Recursively go through the 'source' tree
				   // and see if we can find corresponding
				   // entries in the 'destination' tree. If
				   // not, error out (i.e. we have just read
				   // an XML file that has entries that
				   // weren't declared in the ParameterHandler
				   // object); if so, copy the value of these
				   // nodes into the destination object
  bool
  read_xml_recursively (const boost::property_tree::ptree &source,
			const std::string                 &current_path,
			const char                         path_separator,
			const std::vector<std_cxx1x::shared_ptr<const Patterns::PatternBase> > &
			patterns,
			boost::property_tree::ptree       &destination)
  {
    for (boost::property_tree::ptree::const_iterator p = source.begin();
	 p != source.end(); ++p)
      {
					 // a sub-tree must either be a
					 // parameter node or a subsection
	if (p->second.get_optional<std::string>("value"))
	  {
					     // make sure we have a
					     // corresponding entry in the
					     // destination object as well
	    const std::string full_path
	      = (current_path == ""
		 ?
		 p->first
		 :
		 current_path + path_separator + p->first);
	    if (destination.get_optional<std::string> (full_path)
		&&
		destination.get_optional<std::string> (full_path +
						       path_separator +
						       "value"))
	      {
						 // first make sure that the
						 // new entry actually
						 // satisfies its constraints
		const std::string new_value
		  = p->second.get<std::string>("value");
		
		const unsigned int pattern_index
		  = destination.get<unsigned int> (full_path +
						   path_separator +
						   "pattern");
		if (patterns[pattern_index]->match(new_value) == false)
		  {
		    std::cerr << "    The entry value" << std::endl
			      << "        " << new_value << std::endl
			      << "    for the entry named" << std::endl
			      << "        " << full_path << std::endl
			      << "    does not match the given pattern" << std::endl
			      << "        " << patterns[pattern_index]->description()
			      << std::endl;
		    return false;
		  }

						 // set the found parameter in
						 // the destination argument
		destination.put (full_path + path_separator + "value",
				 new_value);
		
						 // this node might have
						 // sub-nodes in addition to
						 // "value", such as
						 // "default_value",
						 // "documentation", etc. we
						 // might at some point in the
						 // future want to make sure
						 // that if they exist that
						 // they match the ones in the
						 // 'destination' tree
	      }
	    else
	      {
		std::cerr << "The entry <" << full_path
			  << "> with value <"
			  << p->second.get<std::string>("value")
			  << "> has not been declared."
			  << std::endl;
		return false;
	      }
	  }
	else
	  {
					     // it must be a subsection
	    const bool result
	      = read_xml_recursively (p->second,
				      (current_path == "" ?
				       p->first :
				       current_path + path_separator + p->first),
				      path_separator,
				      patterns,
				      destination);

					     // see if the recursive read
					     // succeeded. if yes, continue,
					     // otherwise exit now
	    if (result == false)
	      return false;
	  }
      }

    return true;
  }
}



bool ParameterHandler::read_input_from_xml (std::istream &in)
{
				   // read the XML tree assuming that (as we
				   // do in print_parameters(XML) it has only
				   // a single top-level node called
				   // "ParameterHandler"
  boost::property_tree::ptree single_node_tree;
  try
    {
      read_xml (in, single_node_tree);
    }
  catch (...)
    {
      std::cerr << "This input stream appears not to be valid XML"
		<< std::endl;
      return false;
    }

				   // make sure there is a top-level element
				   // called "ParameterHandler"
  if (!single_node_tree.get_optional<std::string>("ParameterHandler"))
    {
      std::cerr << "There is no top-level XML element called \"ParameterHandler\"."
		<< std::endl;
      return false;
    }

				   // ensure that there is only a single
				   // top-level element
  if (std::distance (single_node_tree.begin(), single_node_tree.end()) != 1)
    {
      std::cerr << "The top-level XML element \"ParameterHandler\" is "
		<< "not the only one."
		<< std::endl;
      std::cerr << "(There are "
		<< std::distance (single_node_tree.begin(),
				  single_node_tree.end())
		<< " top-level elements.)"
		<< std::endl;
      return false;
    }

				   // read the child elements recursively
  const boost::property_tree::ptree
    &my_entries = single_node_tree.get_child("ParameterHandler");

  return read_xml_recursively (my_entries, "", path_separator, patterns,
			       *entries);
}



void ParameterHandler::clear ()
{
  entries.reset (new boost::property_tree::ptree());
}



void
ParameterHandler::declare_entry (const std::string           &entry,
                                 const std::string           &default_value,
                                 const Patterns::PatternBase &pattern,
                                 const std::string           &documentation)
{
  entries->put (get_current_full_path(entry) + path_separator + "value",
	       default_value);
  entries->put (get_current_full_path(entry) + path_separator + "default_value",
	       default_value);
  entries->put (get_current_full_path(entry) + path_separator + "documentation",
	       documentation);

				   // clone the pattern and store its
				   // index in the node
  patterns.push_back (std_cxx1x::shared_ptr<const Patterns::PatternBase>
		      (pattern.clone()));
  entries->put (get_current_full_path(entry) + path_separator + "pattern",
	       static_cast<unsigned int>(patterns.size()-1));
				   // also store the description of
				   // the pattern. we do so because we
				   // may wish to export the whole
				   // thing as XML or any other format
				   // so that external tools can work
				   // on the parameter file; in that
				   // case, they will have to be able
				   // to re-create the patterns as far
				   // as possible
  entries->put (get_current_full_path(entry) + path_separator +
	       "pattern_description",
	       patterns.back()->description());
}



void ParameterHandler::enter_subsection (const std::string &subsection)
{
  const std::string current_path = get_current_path ();

				   // if necessary create subsection
  if (!entries->get_child_optional (get_current_full_path(subsection)))
    entries->add_child (get_current_full_path(subsection),
			boost::property_tree::ptree());

				   // then enter it
  subsection_path.push_back (subsection);
}



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
}



std::string
ParameterHandler::get (const std::string &entry_string) const
{
				   // assert that the entry is indeed
				   // declared
  if (boost::optional<std::string> value
      = entries->get_optional<std::string> (get_current_full_path(entry_string) + path_separator + "value"))
    return value.get();
  else
    {
      Assert (false, ExcEntryUndeclared(entry_string));
      return "";
    }
}



long int ParameterHandler::get_integer (const std::string &entry_string) const
{
  std::string s = get (entry_string);
  char *endptr;
  long int i = std::strtol (s.c_str(), &endptr, 10);
				   // assert there was no error
  AssertThrow (*endptr == '\0', ExcConversionError(s));

  return i;
}



double ParameterHandler::get_double (const std::string &entry_string) const
{
  std::string s = get (entry_string);
  char *endptr;
  double d = std::strtod (s.c_str(), &endptr);
				   // assert there was no error
  AssertThrow ((s.c_str()!='\0') || (*endptr == '\0'),
	       ExcConversionError(s));

  return d;
}



bool ParameterHandler::get_bool (const std::string &entry_string) const
{
  std::string s = get(entry_string);

  AssertThrow ((s=="true") || (s=="false") ||
               (s=="yes") || (s=="no"),
               ExcConversionError(s));
  if (s=="true" || s=="yes")
    return true;
  else
    return false;
}



void
ParameterHandler::set (const std::string &entry_string,
                       const std::string &new_value)
{
				   // assert that the entry is indeed
				   // declared
  if (entries->get_optional<std::string>
      (get_current_full_path(entry_string) + path_separator + "value"))
    {
      const unsigned int pattern_index
	= entries->get<unsigned int> (get_current_full_path(entry_string) + path_separator + "pattern");
      AssertThrow (patterns[pattern_index]->match(new_value),
		   ExcValueDoesNotMatchPattern (new_value,
						entries->get<std::string>
						(get_current_full_path(entry_string) +
						 path_separator +
						 "pattern_description")));

      entries->put (get_current_full_path(entry_string) + path_separator + "value",
		   new_value);
    }
  else
    AssertThrow (false, ExcEntryUndeclared(entry_string));
}


void
ParameterHandler::set (const std::string &entry_string,
                       const char        *new_value)
{
                                   // simply forward
  set (entry_string, std::string(new_value));
}


void
ParameterHandler::set (const std::string &entry_string,
                       const double      &new_value)
{
  std::ostringstream s;
  s << std::setprecision(16);
  s << new_value;

                                   // hand this off to the function that
                                   // actually sets the value as a string
  set (entry_string, s.str());
}



void
ParameterHandler::set (const std::string &entry_string,
                       const long int    &new_value)
{
  std::ostringstream s;
  s << new_value;

                                   // hand this off to the function that
                                   // actually sets the value as a string
  set (entry_string, s.str());
}



void
ParameterHandler::set (const std::string &entry_string,
                       const bool        &new_value)
{
  std::ostringstream s;
  s << new_value;

                                   // hand this off to the function that
                                   // actually sets the value as a string
  set (entry_string, s.str());
}



std::ostream &
ParameterHandler::print_parameters (std::ostream     &out,
                                    const OutputStyle style)
{
  AssertThrow (out, ExcIO());

  switch (style)
    {
      case XML:
      {
					 // call the writer
					 // function and exit as
					 // there is nothing
					 // further to do down in
					 // this function
					 //
					 // XML has a requirement that
					 // there can only be one
					 // single top-level entry,
					 // but we may have multiple
					 // entries and sections.  we
					 // work around this by
					 // creating a tree just for
					 // this purpose with the
					 // single top-level node
					 // "ParameterHandler" and
					 // assign the existing tree
					 // under it
	boost::property_tree::ptree single_node_tree;
	single_node_tree.add_child("ParameterHandler",
				   *entries);

	write_xml (out, single_node_tree);
	return out;
      }


      case JSON:
					     // call the writer
					     // function and exit as
					     // there is nothing
					     // further to do down in
					     // this function
	    write_json (out, *entries);
	    return out;

      case Text:
	    out << "# Listing of Parameters" << std::endl
		<< "# ---------------------" << std::endl;
	    break;
      case LaTeX:
	    out << "\\subsubsection*{Listing of parameters}";
	    out << std::endl << std::endl;
	    out << "\\begin{itemize}"
	        << std::endl;
	    break;
      case Description:
	    out << "Listing of Parameters:" << std::endl << std::endl;
	    break;
      case ShortText:
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };

				   // dive recursively into the subsections
  print_parameters_section (out, style, 0);

  switch (style)
    {
      case Text:
      case Description:
      case ShortText:
	    break;
      case LaTeX:
	    out << "\\end{itemize}" << std::endl;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };

  return out;
}



// Print a section in the desired style. The styles are separated into
// several verbosity classes depending on the higher bits.
//
// If bit 7 (128) is set, comments are not printed.
// If bit 6 (64) is set, default values after change are not printed.
void
ParameterHandler::print_parameters_section (std::ostream      &out,
                                            const OutputStyle  style,
                                            const unsigned int indent_level)
{
  AssertThrow (out, ExcIO());

  const boost::property_tree::ptree &current_section
    = entries->get_child (get_current_path());

  switch (style)
    {
      case Text:
      case ShortText:
      {
					 // first find out the longest
					 // entry name to be able to
					 // align the equal signs
					 //
					 // to do this loop over all
					 // nodes of the current tree,
					 // select the parameter nodes
					 // (and discard sub-tree
					 // nodes) and take the
					 // maximum of their lengths
        std::size_t longest_name = 0;
	for (boost::property_tree::ptree::const_iterator
	       p = current_section.begin();
	     p != current_section.end(); ++p)
	  if (is_parameter_node (p->second) == true)
	    longest_name = std::max (longest_name,
				     demangle(p->first).length());

                                         // likewise find the longest
                                         // actual value string to
                                         // make sure we can align the
                                         // default and documentation
                                         // strings
        std::size_t longest_value = 0;
	for (boost::property_tree::ptree::const_iterator
	       p = current_section.begin();
	     p != current_section.end(); ++p)
	  if (is_parameter_node (p->second) == true)
	    longest_value = std::max (longest_value,
				      p->second.get<std::string>("value").length());


					 // print entries one by
					 // one. make sure they are
					 // sorted by using the
					 // appropriate iterators
	bool first_entry = true;
	for (boost::property_tree::ptree::const_assoc_iterator
	       p = current_section.ordered_begin();
	     p != current_section.not_found(); ++p)
	  if (is_parameter_node (p->second) == true)
	    {
	      const std::string value = p->second.get<std::string>("value");

					       // if there is documentation,
					       // then add an empty line (unless
					       // this is the first entry in a
					       // subsection), print the
					       // documentation, and then the
					       // actual entry; break the
					       // documentation into readable
					       // chunks such that the whole
					       // thing is at most 78 characters
					       // wide
	      if ((!(style & 128)) &&
		  !p->second.get<std::string>("documentation").empty())
		{
		  if (first_entry == false)
		    out << std::endl;
		  else
		    first_entry = false;

		  const std::vector<std::string> doc_lines
		    = Utilities::
		    break_text_into_lines (p->second.get<std::string>("documentation"),
					   78 - indent_level*2 - 2);

		  for (unsigned int i=0; i<doc_lines.size(); ++i)
		    out << std::setw(indent_level*2) << ""
			<< "# "
			<< doc_lines[i]
			<< std::endl;
		}



					       // print name and value
					       // of this entry
	      out << std::setw(indent_level*2) << ""
		  << "set "
		  << demangle(p->first)
		  << std::setw(longest_name-demangle(p->first).length()+1) << " "
		  << "= " << value;

					       // finally print the
					       // default value, but
					       // only if it differs
					       // from the actual value
	      if ((!(style & 64)) && value != p->second.get<std::string>("default_value"))
		{
		  out << std::setw(longest_value-value.length()+1) << ' '
		      << "# ";
		  out << "default: " << p->second.get<std::string>("default_value");
		}

	      out << std::endl;
	    }

        break;
      }

      case LaTeX:
      {
					 // print entries one by
					 // one. make sure they are
					 // sorted by using the
					 // appropriate iterators
	for (boost::property_tree::ptree::const_assoc_iterator
	       p = current_section.ordered_begin();
	     p != current_section.not_found(); ++p)
	  if (is_parameter_node (p->second) == true)
	    {
	      const std::string value = p->second.get<std::string>("value");

					       // print name and value
	      out << "\\item {\\bf " << demangle(p->first) << ":} "
		  << value
		  << " (";

					       // if there is a
					       // documenting string,
					       // print it as well
	      if (!p->second.get<std::string>("documentation").empty())
		out <<p->second.get<std::string>("documentation")  << ", ";

					       // finally print default
					       // value
	      out << "{\\it default:} "
		  << p->second.get<std::string>("default_value")
		  << ")"
		  << std::endl;
	    }

        break;
      }

      case Description:
      {
					 // first find out the longest
					 // entry name to be able to
					 // align the equal signs
        std::size_t longest_name = 0;
	for (boost::property_tree::ptree::const_iterator
	       p = current_section.begin();
	     p != current_section.end(); ++p)
	  if (is_parameter_node (p->second) == true)
	    longest_name = std::max (longest_name,
				     demangle(p->first).length());

					 // print entries one by
					 // one. make sure they are
					 // sorted by using the
					 // appropriate iterators
	for (boost::property_tree::ptree::const_assoc_iterator
	       p = current_section.ordered_begin();
	     p != current_section.not_found(); ++p)
	  if (is_parameter_node (p->second) == true)
	    {
	      const std::string value = p->second.get<std::string>("value");

					       // print name and value
	      out << std::setw(indent_level*2) << ""
		  << "set "
		  << demangle(p->first)
		  << std::setw(longest_name-demangle(p->first).length()+1) << " "
		  << " = ";

					       // print possible values:
	      const std::vector<std::string> description_str
		= Utilities::break_text_into_lines (p->second.get<std::string>
						    ("pattern_description"),
						    78 - indent_level*2 - 2, '|');
	      if (description_str.size() > 1)
		{
		  out << std::endl;
		  for (unsigned int i=0; i<description_str.size(); ++i)
		    out << std::setw(indent_level*2+6) << ""
			<< description_str[i] << std::endl;
		}
	      else if (description_str.empty() == false)
		out << "  " << description_str[0] << std::endl;
	      else
		out << std::endl;

					       // if there is a
					       // documenting string,
					       // print it as well
	      if (p->second.get<std::string>("documentation").length() != 0)
		out << std::setw(indent_level*2 + longest_name + 10) << ""
		    << "(" << p->second.get<std::string>("documentation") << ")" << std::endl;
	    }

        break;
      }

      default:
            Assert (false, ExcNotImplemented());
    }


				   // if there was text before and there are
                                   // sections to come, put two newlines
                                   // between the last entry and the first
                                   // subsection
  {
    unsigned int n_parameters = 0;
    unsigned int n_sections   = 0;
    for (boost::property_tree::ptree::const_iterator
	   p = current_section.begin();
	 p != current_section.end(); ++p)
      if (is_parameter_node (p->second) == true)
	++n_parameters;
      else
	++n_sections;

    if ((style != Description)
	&&
	(!(style & 128))
	&&
	(n_parameters != 0)
	&&
	(n_sections != 0))
      out << std::endl << std::endl;
  }

				   // now traverse subsections tree,
				   // in alphabetical order
  for (boost::property_tree::ptree::const_assoc_iterator
	 p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == false)
      {
                                         // first print the subsection header
        switch (style)
          {
            case Text:
	    case Description:
	    case ShortText:
                  out << std::setw(indent_level*2) << ""
                      << "subsection " << demangle(p->first) << std::endl;
                  break;
            case LaTeX:
                  out << std::endl
                      << "\\item {\\bf "
                      << "Subsection " << demangle(p->first)
                      << "}" << std::endl
                      << "\\begin{itemize}"
                      << std::endl;
                  break;
            default:
                  Assert (false, ExcNotImplemented());
          };

                                         // then the contents of the
                                         // subsection
        enter_subsection (demangle(p->first));
        print_parameters_section (out, style, indent_level+1);
        leave_subsection ();
        switch (style)
          {
	    case Text:
						   // write end of
						   // subsection. one
						   // blank line after
						   // each subsection
		  out << std::setw(indent_level*2) << ""
		      << "end" << std::endl
		      << std::endl;

						   // if this is a toplevel
						   // subsection, then have two
						   // newlines
		  if (indent_level == 0)
		    out << std::endl;

		  break;
	    case Description:
		  break;
	    case ShortText:
						   // write end of
						   // subsection.
		  out << std::setw(indent_level*2) << ""
		      << "end" << std::endl;
                  break;
            case LaTeX:
                  out << "\\end{itemize}"
                      << std::endl;
                  break;
            default:
                  Assert (false, ExcNotImplemented());
          }
      }
}



void
ParameterHandler::log_parameters (LogStream &out)
{
  out.push("parameters");
				   // dive recursively into the
				   // subsections
  log_parameters_section (out);

  out.pop();
}



void
ParameterHandler::log_parameters_section (LogStream &out)
{
  const boost::property_tree::ptree &current_section
    = entries->get_child (get_current_path());

				   // print entries one by
				   // one. make sure they are
				   // sorted by using the
				   // appropriate iterators
  for (boost::property_tree::ptree::const_assoc_iterator
	 p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == true)
      out << demangle(p->first) << ": "
          << p->second.get<std::string>("value") << std::endl;

				   // now transverse subsections tree
				   // now traverse subsections tree,
				   // in alphabetical order
  for (boost::property_tree::ptree::const_assoc_iterator
	 p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == false)
      {
	out.push (demangle(p->first));
        enter_subsection (demangle(p->first));
        log_parameters_section (out);
        leave_subsection ();
	out.pop ();
      }
}



bool
ParameterHandler::scan_line (std::string        line,
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
  if ((line.length() != 0) && (std::isspace (line[0])))  line.erase (0, 1);
				   // if line is now empty: leave
  if (line.length() == 0) return true;

  if (std::isspace (line[line.length()-1]))
    line.erase (line.size()-1, 1);

				   // enter subsection
  if ((line.find ("SUBSECTION ") == 0) ||
      (line.find ("subsection ") == 0))
    {
				       // delete this prefix
      line.erase (0, std::string("subsection").length()+1);

      const std::string subsection = line;

				       // check whether subsection exists
      if (!entries->get_child_optional (get_current_full_path(subsection)))
	{
	  std::cerr << "Line " << lineno
                    << ": There is no such subsection to be entered:" << get_current_full_path(subsection) <<std::endl;
	  for (unsigned int i=0; i<subsection_path.size(); ++i)
	    std::cerr << std::setw(i*2+4) << " "
		      << "subsection " << subsection_path[i] << std::endl;
	  std::cerr << std::setw(subsection_path.size()*2+4) << " "
		    << "subsection " << subsection << std::endl;
	  return false;
	}

				       // subsection exists
      subsection_path.push_back (subsection);
      return true;
    }

				   // exit subsection
  if ((line.find ("END") == 0) ||
      (line.find ("end") == 0))
    {
      if (subsection_path.size() == 0)
	{
	  std::cerr << "Line " << lineno
		    << ": There is no subsection to leave here!" << std::endl;
	  return false;
	}
      else
	return leave_subsection ();
    }

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

      const std::string current_path = get_current_path ();

				       // assert that the entry is indeed
				       // declared
      if (entries->get_optional<std::string> (get_current_full_path(entry_name) + path_separator + "value"))
	{
					   // if entry was declared:
					   // does it match the regex? if not,
					   // don't enter it into the database
					   // exception: if it contains characters
					   // which specify it as a multiple loop
					   // entry, then ignore content
	  if (entry_value.find ('{') == std::string::npos)
	    {
	      const unsigned int pattern_index
		= entries->get<unsigned int> (get_current_full_path(entry_name) + path_separator + "pattern");
	      if (!patterns[pattern_index]->match(entry_value))
		{
		  std::cerr << "Line " << lineno << ":" << std::endl
			    << "    The entry value" << std::endl
			    << "        " << entry_value << std::endl
			    << "    for the entry named" << std::endl
			    << "        " << entry_name << std::endl
			    << "    does not match the given pattern" << std::endl
			    << "        " << patterns[pattern_index]->description()
			    << std::endl;
		  return false;
		}
	    }

	  entries->put (get_current_full_path(entry_name) + path_separator + "value",
		       entry_value);
	  return true;
	}
      else
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
	}
    }

				   // this line matched nothing known
  std::cerr << "Line " << lineno
            << ": This line matched nothing known ('set' or 'subsection' missing!?):" << std::endl
	    << "    " << line << std::endl;
  return false;
}



std::size_t
ParameterHandler::memory_consumption () const
{
//TODO: add to this an estimate of the memory in the property_tree
  return (MemoryConsumption::memory_consumption (subsection_path));
}



bool
ParameterHandler::operator == (const ParameterHandler & prm2)  const
{
  if(patterns.size() != prm2.patterns.size())
    return false;

  for(unsigned int j=0; j<patterns.size(); ++j)
    if (patterns[j]->description() != prm2.patterns[j]->description())
      return false;

				   // instead of walking through all
				   // the nodes of the two trees
				   // entries and prm2.entries and
				   // comparing them for equality,
				   // simply dump the content of the
				   // entire structure into a string
				   // and compare those for equality
  std::ostringstream o1, o2;
  write_json (o1, *entries);
  write_json (o2, *prm2.entries);
  return (o1.str() == o2.str());
}




MultipleParameterLoop::UserClass::~UserClass ()
{}



MultipleParameterLoop::MultipleParameterLoop()
                :
		n_branches(0)
{}



MultipleParameterLoop::~MultipleParameterLoop ()
{}



bool MultipleParameterLoop::read_input (std::istream &input)
{
  AssertThrow (input, ExcIO());

  bool x = ParameterHandler::read_input (input);
  if (x)
    init_branches ();
  return x;
}



bool MultipleParameterLoop::read_input (const std::string &filename,
					bool optional,
					bool write_compact)
{
  return ParameterHandler::read_input (filename, optional, write_compact);
				   // don't call init_branches, since this read_input
				   // function calls
				   // MultipleParameterLoop::Readinput(std::istream &, std::ostream &)
				   // which itself calls init_branches.
}



bool MultipleParameterLoop::read_input_from_string (const char *s)
{
  bool x = ParameterHandler::read_input (s);
  init_branches ();
  return x;
}



void MultipleParameterLoop::loop (MultipleParameterLoop::UserClass &uc)
{
  for (unsigned int run_no=0; run_no<n_branches; ++run_no)
    {
				       // give create_new one-based numbers
      uc.create_new (run_no+1);
      fill_entry_values (run_no);
      uc.run (*this);
    };
}



void MultipleParameterLoop::init_branches ()
{
  multiple_choices.clear ();
  init_branches_current_section ();

				   // split up different values
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    multiple_choices[i].split_different_values ();

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
		  << "        " << n_branches << " variant runs that will be performed."
                  << std::endl;


				   // do a first run on filling the values to
				   // check for the conformance with the regexp
				   // (later on, this will be lost in the whole
				   // other output)
  for (unsigned int i=0; i<n_branches; ++i)
    fill_entry_values (i);
}



void MultipleParameterLoop::init_branches_current_section ()
{
  const boost::property_tree::ptree &current_section
    = entries->get_child (get_current_path());

				   // check all entries in the present
				   // subsection whether they are
				   // multiple entries
				   //
				   // we loop over entries in sorted
				   // order to guarantee backward
				   // compatibility to an earlier
				   // implementation
  for (boost::property_tree::ptree::const_assoc_iterator
	 p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == true)
      {
	const std::string value = p->second.get<std::string>("value");
	if (value.find('{') != std::string::npos)
	  multiple_choices.push_back (Entry(subsection_path,
					    demangle(p->first),
					    value));
      }

				   // then loop over all subsections
  for (boost::property_tree::ptree::const_iterator
	 p = current_section.begin();
       p != current_section.end(); ++p)
    if (is_parameter_node (p->second) == false)
      {
	enter_subsection (demangle(p->first));
	init_branches_current_section ();
	leave_subsection ();
      }
}




void MultipleParameterLoop::fill_entry_values (const unsigned int run_no)
{
  unsigned int possibilities = 1;

  std::vector<Entry>::iterator choice;
  for (choice = multiple_choices.begin();
       choice != multiple_choices.end();
       ++choice)
    {
      const unsigned int selection
	= (run_no/possibilities) % choice->different_values.size();
      std::string entry_value;
      if (choice->type == Entry::variant)
	entry_value = choice->different_values[selection];
      else
	{
	  if (run_no>=choice->different_values.size())
	    {
	      std::cerr << "The given array for entry <"
			<< choice->entry_name
			<< "> does not contain enough elements! Taking empty string instead."
                        << std::endl;
	      entry_value = "";
	    }
	  else
	    entry_value = choice->different_values[run_no];
	}

				       // temporarily enter the
				       // subsection tree of this
				       // multiple entry, set the
				       // value, and get out
				       // again. the set() operation
				       // also tests for the
				       // correctness of the value
				       // with regard to the pattern
      subsection_path.swap (choice->subsection_path);
      set (choice->entry_name, entry_value);
      subsection_path.swap (choice->subsection_path);

				       // move ahead if it was a variant entry
      if (choice->type == Entry::variant)
	possibilities *= choice->different_values.size();
    }
}




std::size_t
MultipleParameterLoop::memory_consumption () const
{
  std::size_t mem = ParameterHandler::memory_consumption ();
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    mem += multiple_choices[i].memory_consumption ();

  return mem;
}



MultipleParameterLoop::Entry::Entry (const std::vector<std::string> &ssp,
				     const std::string              &Name,
				     const std::string              &Value)
                :
		subsection_path (ssp), entry_name(Name), entry_value(Value), type (Entry::array)
{}



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
  while (std::isspace (multiple[0])) multiple.erase (0,1);
  while (std::isspace (multiple[multiple.size()-1])) multiple.erase (multiple.size()-1,1);

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
}


std::size_t
MultipleParameterLoop::Entry::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (subsection_path) +
	  MemoryConsumption::memory_consumption (entry_name) +
	  MemoryConsumption::memory_consumption (entry_value) +
	  MemoryConsumption::memory_consumption (different_values) +
	  sizeof (type));
}

DEAL_II_NAMESPACE_CLOSE
