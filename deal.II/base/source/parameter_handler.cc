/* $Id$ */

#include <basic/parameter_handler.h>
#include <fstream.h>
#include <iomanip.h>

extern "C" {
#include <stdlib.h>
}



const String ParameterHandler::RegularExpressions::WhiteSpace = "[ \n\t]+",
	     ParameterHandler::RegularExpressions::Integer    = "\\(-?[0-9]+[ \n\t]*\\)",
	     ParameterHandler::RegularExpressions::Double     = "\\("
								    "-?"
								    "\\("
									"\\([0-9]+\\.[0-9]*\\)"
									"\\|"
									"\\([0-9]+\\)"
									"\\|"
									"\\(\\.[0-9]+\\)"
									"\\)"
								    "\\([eE][---+]?[0-9]+\\)?"
								    "[ \n\t]*"
								    "\\)",
	     ParameterHandler::RegularExpressions::AlphaNum   = "\\([0-9A-Za-z]+[ \n\t]*\\)";



ParameterHandler::ParameterHandler () :
		status(true) {};


bool ParameterHandler::read_input (istream &input) {
  String line;
  int lineno=0;
  while (input) 
    {
      ++lineno;
				       // get one line; limit of 10000
				       // chars is arbitrary but should
				       // be enough
      char c[10000];
      input.getline ((char *)&c, 10000);
      line = c;
      if (!scan_line (line, lineno)) 
	status = false;
    };

  return status;
};


bool ParameterHandler::read_input (const String &filename) {
  ifstream input (filename);
  if (!input) 
    {
      cerr << "ParameterHandler::read_input: could not open file <"
	   << filename << ">. Aborting." << endl;
      return false;
    };
  
  return read_input (input);
};



bool ParameterHandler::read_input_from_string (const char *s) {
				   // if empty string then exit
				   // with success
  if ((s == 0) || ((*s) == 0)) return true;
  
  String line;
  String input (s);
  int    lineno=0;

				   // if necessary append a newline char
				   // to make all lines equal
  if (input[input.length()-1] != '\n')
    input += '\n';
  
  while (input) 
    {
				       // get one line from Input (=s)
      line = input.before ('\n');
				       // delete this part including
				       // the backspace
      input.del (line);
      input.del (input[0]);
      ++lineno;
      
      if (!scan_line (line, lineno)) 
	status = false;
    };

  return status;
};




void ParameterHandler::clear () {
  status = true;

  subsection_path.erase       (subsection_path.begin(), subsection_path.end());

  defaults.entries.erase     (defaults.entries.begin(), defaults.entries.end());
  changed_entries.entries.erase (changed_entries.entries.begin(),
				 changed_entries.entries.end());

  map<String, Section*, less<String> >::iterator p;

  for (p=defaults.subsections.begin(); p!=defaults.subsections.end(); ++p)
    delete (*p).second;
  for (p=changed_entries.subsections.begin(); p!=changed_entries.subsections.end(); ++p)
    delete (*p).second;

  defaults.subsections.erase (defaults.subsections.begin(), defaults.subsections.end());
  changed_entries.subsections.erase (changed_entries.subsections.begin(),
				     changed_entries.subsections.end());
};



bool ParameterHandler::declare_entry    (const String &entry,
					 const String &default_value,
					 const String &pattern) {
  Section* p = get_present_defaults_subsection ();

				   // assertions:
				   // entry must not already exist
  Assert (p->entries.find (entry) == p->entries.end(),
	  ExcEntryAlreadyExists (entry));
				   // Default must match Pattern
  Assert (default_value.matches (Regex(pattern)),
	  ExcDefaultDoesNotMatchRegex(default_value, pattern));
  
				   // does entry already exist?
  if (p->entries.find (entry) != p->entries.end())
    return false;

  p->entries[entry] = make_pair(default_value, pattern);

				   // check whether default answer matches
				   // the pattern
  if (!default_value.matches (Regex(pattern)))
    return false;

  return true;
};
  
  


void ParameterHandler::enter_subsection (const String &subsection) {
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

    

bool ParameterHandler::leave_subsection () {
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



const String & ParameterHandler::get (const String &entry_string) const {
  const Section* pd = get_present_defaults_subsection ();
  const Section* pc = get_present_changed_subsection ();

				   // assert that the according entry is already
				   // declared in the defaults tree
  Assert (pd->entries.find (entry_string) != pd->entries.end(),
	  ExcEntryUndeclared(entry_string));
  
  if (pd->entries.find (entry_string) == pd->entries.end()) 
    {
      static const String empty_string;
      return empty_string;
    };


				   // entry exists; now find out whether
				   // it was changed:
  map<String, pair<String,String>, less<String> >::const_iterator ptr;
  ptr = pc->entries.find (entry_string);
  if (ptr != pc->entries.end())
    return (*ptr).second.first;

				   // not changed
  ptr = pd->entries.find (entry_string);
  return (*ptr).second.first;
};



long int ParameterHandler::get_integer (const String &entry_string) const {
  String s = get (entry_string);
  char *endptr;
  long int i = strtol (s.chars(), &endptr, 10);
				   // assert there was no error
  Assert ((s.chars()!='\0') || (*endptr == '\0'),
	  ExcConversionError(s));

  return i;
};



double ParameterHandler::get_double (const String &entry_string) const {
  String s = get (entry_string);
  char *endptr;
  double d = strtod (s.chars(), &endptr);
				   // assert there was no error
  Assert ((s.chars()!='\0') || (*endptr == '\0'),
	  ExcConversionError(s));

  return d;
};




ostream & ParameterHandler::print_parameters (ostream &out, OutputStyle style) {
				   // assert that only known formats are
				   // given as "style"
  Assert ((style == Text) || (style == LaTeX), ExcNotImplemented());

  switch (style) 
    {
      case Text:
	    out << "Listing of Parameters" << endl
		<< "---------------------" << endl;
	    break;
      case LaTeX:
	    out << "\\subsubsection*{Listing of parameters}";
	    out << endl << endl;
	    out << "\\begin{itemize}"
	        << endl;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };
   
				   // dive recursively into the subsections
  print_parameters_section (out, style, 1);

  switch (style) 
    {
      case Text:
	    break;
      case LaTeX:
	    out << "\\end{itemize}" << endl;
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };
  
  return out;
};



void ParameterHandler::print_parameters_section (ostream &out,
						 const OutputStyle style,
						 const unsigned int indent_level) {
				   // assert that only known formats are
				   // given as "style"
  Assert ((style == Text) || (style == LaTeX), ExcNotImplemented());
  
  Section *pd = get_present_defaults_subsection ();
  Section *pc = get_present_changed_subsection ();

				   // traverse entry list
  map<String, pair<String,String>, less<String> >::const_iterator ptr;

				   // first find out the longest entry name
  unsigned int longest_entry = 0;
  for (ptr = pd->entries.begin(); ptr != pd->entries.end(); ++ptr)
    if ((*ptr).first.length() > longest_entry)
      longest_entry = (*ptr).first.length();

				   // print entries one by one
  for (ptr = pd->entries.begin(); ptr != pd->entries.end(); ++ptr)
    {
				       // check whether this entry is listed
				       // in the Changed tree
      if (pc->entries.find((*ptr).first) != pc->entries.end()) 
	switch (style) 
	  {
	    case Text:
		  out << setw(indent_level*2) << ""
		      << (*ptr).first
		      << setw(longest_entry-(*ptr).first.length()+3) << " = "
		      << pc->entries[(*ptr).first].first
		      << "  <"
		      << pd->entries[(*ptr).first].first
		      << ">"
		      << endl;
		  break;
	    case LaTeX:
		  out << "\\item {\\bf " << (*ptr).first << ":} "
		      << pc->entries[(*ptr).first].first
		      << " ({\\it default:} "
                      << pd->entries[(*ptr).first].first
		      << ")"
		      << endl;
		  break;
	    default:
		  Assert (false, ExcNotImplemented());
	  }
				       // not a changed entry
      else
	switch (style) 
	  {
	    case Text:
		  out << setw(indent_level*2) << ""
		      << (*ptr).first
		      << setw(longest_entry-(*ptr).first.length()+2) << "= "
		      << (*ptr).second.first << endl;
		  break;
	    case LaTeX:
		  out << "\\item {\\bf " << (*ptr).first << ":} "
		      << (*ptr).second.first
		      << endl;
		  break;
	    default:
		  Assert (false, ExcNotImplemented());
	  };
    };
  

				   // now transverse subsections tree
  map<String, Section*, less<String> >::const_iterator ptrss;
  for (ptrss = pd->subsections.begin(); ptrss != pd->subsections.end(); ++ptrss)
    {
      switch (style) 
	{
	  case Text:
		out << setw(indent_level*2) << ""
		    << "subsection " << (*ptrss).first << endl;
		break;
	  case LaTeX:
		out << endl
		    << "\\item {\\bf "
		    << "Subsection " << (*ptrss).first
		    << "}" << endl
		    << "\\begin{itemize}"
		    << endl;
		break;
	  default:
		Assert (false, ExcNotImplemented());
	};
      enter_subsection ((*ptrss).first);
      print_parameters_section (out, style, indent_level+1);
      leave_subsection ();
      switch (style) 
	{
	  case Text:
		break;
	  case LaTeX:
		out << "\\end{itemize}"
		    << endl;
		break;
	  default:
		Assert (false, ExcNotImplemented());
	};
    };
};


	  
  

bool ParameterHandler::scan_line (String line, const unsigned int lineno) {
				   // if there is a comment, delete it
  if (line.contains ('#'))
    line = line.before ("#");
				   // replace every whitespace sequence
				   // by ' '
  line.gsub (Regex(RegularExpressions::WhiteSpace), ' ');
				   // now every existing whitespace
				   // should be exactly on ' ';
				   // if at end or beginning: delete
  if ((line.length() != 0) && (line[0] == ' '))  line.del (' ');
				   // if line is now empty: leave
  if (line.length() == 0) return true;

  if (line[line.length()-1] == ' ')  line.at (' ', -1) = "";

  static Regex re_enter_subsection ("\\(SUBSECTION\\|subsection\\)"
				    "[ \n\t]+"
				    "\\(.*\\)");
  static Regex re_entry           ("\\(SET\\|set\\)"   // "set" command
				   "[ \n\t]+"
				   "\\(.*[^ \t\n]\\)"   // 2: entry name
				   "[ \n\t]*"
				   "="
				   "[ \n\t]*"
				   "\\(.*\\)");        // 3: entry value
  static Regex re_exit_subsection  ("\\(END\\|end\\)");

				   // enter subsection
  if (line.matches (re_enter_subsection)) 
    {
				       // first find out name of subsection
      int start, len;
      re_enter_subsection.match (line, line.length());
      re_enter_subsection.match_info (start, len, 2);
      
      String SecName = line.at (start, len);
      Section* pc = get_present_changed_subsection ();
				       // check whether subsection exists
      if (pc->subsections.find(SecName) == pc->subsections.end()) 
	{
	  cerr << "Line " << lineno << ": There is no such subsection to be entered:" << endl;
	  for (unsigned int i=0; i<subsection_path.size(); ++i)
	    cerr << setw(i*2+4) << " "
		 << "subsection " << subsection_path[i] << endl;
	  cerr << setw(subsection_path.size()*2+4) << " "
	       << "subsection " << SecName << endl;
	  return false;
	};

				       // subsection exists
      subsection_path.push_back (SecName);
      return true;
    };
  

				   // exit subsection
  if (line.matches (re_exit_subsection)) 
    if (subsection_path.size() == 0) 
      {
	cerr << "Line " << lineno << ": There is no subsection to leave here!" << endl;
	return false;
      }
    else
      return leave_subsection ();


				   // regular entry
  if (line.matches (re_entry)) 
    {
      int start, end;
      re_entry.match (line, line.length());

				       // extract entry name
      re_entry.match_info (start, end, 2);
      String entry_name = line.at (start, end);

      				       // extract entry value
      re_entry.match_info (start, end, 3);
      String entry_value = line.at (start, end);

      Section* pd = get_present_defaults_subsection ();

				       // check whether entry was declared
      if (pd->entries.find(entry_name) == pd->entries.end()) 
	{
	  cerr << "Line " << lineno
	       << ": No such entry was declared:" << endl
	       << "    " << entry_name << endl
	       << "    <Present subsection:" << endl;
	  for (unsigned int i=0; i<subsection_path.size(); ++i)
	    cerr << setw(i*2+8) << " "
		 << "subsection " << subsection_path[i] << endl;
	  cerr << "    >" << endl;

	  return false;
	};

				       // if entry was declared:
				       // does it match the regex? if not,
				       // don't enter it into the database
				       // exception: if it contains characters
				       // which specify it as a multiple loop
				       // entry, then ignore content
      if (!entry_value.contains ('{'))
	if (!entry_value.matches (Regex(pd->entries[entry_name].second)))
	  {
	    cerr << "Line " << lineno << ":" << endl
		 << "    The entry value" << endl
		 << "        " << entry_value << endl
		 << "    for the entry named" << endl
		 << "        " << entry_name << endl
		 << "    does not match the given regular expression" << endl
		 << "        " << pd->entries[entry_name].second << endl;
	    return false;
	  };
      
      Section* pc = get_present_changed_subsection ();
				       // the following line declares this entry
				       // if not yet existent and overwrites it
				       // otherwise
      pc->entries[entry_name] = make_pair(entry_value, String(""));

      return true;
    };

				   // this line matched nothing known
  cerr << "Line " << lineno << ": This line matched nothing known:" << endl
       << "    " << line << endl;
  return false;
};



ParameterHandler::Section* ParameterHandler::get_present_defaults_subsection () {
  Section* sec = &defaults;
  vector<String>::const_iterator SecName;
  SecName = subsection_path.begin();
    
  while (SecName != subsection_path.end()) 
    {
      sec = sec->subsections[*SecName];
      ++SecName;
    };

  return sec;
};



const ParameterHandler::Section* ParameterHandler::get_present_defaults_subsection () const {
				   // simply call the non-const version
				   // of this function
  typedef Section* (ParameterHandler::*xptr) ();
  xptr x = &get_present_defaults_subsection;

  return x();
};



ParameterHandler::Section* ParameterHandler::get_present_changed_subsection () {
  Section* sec = &changed_entries;
  vector<String>::iterator SecName;
  SecName = subsection_path.begin();
    
  while (SecName != subsection_path.end()) 
    {
      sec = sec->subsections[*SecName];
      ++SecName;
    };

  return sec;
};



const ParameterHandler::Section* ParameterHandler::get_present_changed_subsection () const {
				   // simply call the non-const version
				   // of this function
  typedef Section* (ParameterHandler::*xptr) ();
  xptr x = &get_present_changed_subsection;

  return x();
};



ParameterHandler::Section::~Section () {
  entries.erase (entries.begin(), entries.end());

  map<String, Section*, less<String> >::iterator p;

  for (p=subsections.begin(); p!=subsections.end(); ++p)
    delete (*p).second;


  subsections.erase (subsections.begin(), subsections.end());
};





MultipleParameterLoop::MultipleParameterLoop() :
		n_branches(0) {};


bool MultipleParameterLoop::read_input (istream &input) {
  bool x = ParameterHandler::read_input (input);
  if (x) init_branches ();
  return x;
};


bool MultipleParameterLoop::read_input (const String &filename) {
				   // I don't know why it is necessary to
				   // declare this function: simply not
				   // overloading it in MultipleParameterLoop
				   // should suffice, but then the compiler can't
				   // find the inherited version and complains that
				   // it can't find a function MultipleParameterLoop::
				   // read_input (String, ostream) instead of trying the
				   // base class (maybe wait for gcc2.8)
  return ParameterHandler::read_input (filename);
				   // don't call init_branches, since this read_input
				   // function calls
				   // MultipleParameterLoop::Readinput(istream &, ostream &)
				   // which itself calls init_branches.
};


bool MultipleParameterLoop::read_input_from_string (const char *s) {
  bool x = ParameterHandler::read_input (s);
  init_branches ();
  return x;
};





void MultipleParameterLoop::loop (MultipleParameterLoop::UserClass &uc) {
  for (int run_no=0; run_no<n_branches; ++run_no) 
    {
				       // give create_new one-based numbers
      uc.create_new (run_no+1);
      fill_entry_values (run_no);
      uc.run (*this);
    };
};


void MultipleParameterLoop::init_branches () {
  if (multiple_choices.size() != 0)
    multiple_choices.erase (multiple_choices.begin(), multiple_choices.end());
  
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
    for (vector<Entry>::iterator i=multiple_choices.end()-1;
	 i >= multiple_choices.begin(); --i)
      if ((*i).different_values.size() == 1)
	multiple_choices.erase (i);

				   // finally calculate number of branches
  n_branches = 1;
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    if (multiple_choices[i].type == variant)
      n_branches *= multiple_choices[i].different_values.size();

				   // check whether array entries have the correct
				   // number of entries
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    if (multiple_choices[i].type == array)
      if (multiple_choices[i].different_values.size() != (unsigned int)n_branches)
	cerr << "    The entry value" << endl
	     << "        " << multiple_choices[i].entry_value << endl
	     << "    for the entry named" << endl
	     << "        " << multiple_choices[i].entry_name << endl
	     << "    does not have the right number of entries for the " << endl
	     << "        " << n_branches << " variant runs that will be performed." << endl;

  
  
				   // do a first run on filling the values to
				   // check for the conformance with the regexp
				   // (afterwards, this will be lost in the whole
				   // other output)
  for (int i=0; i<n_branches; ++i) 
    fill_entry_values (i);
};



void MultipleParameterLoop::init_branches_section (const ParameterHandler::Section &sec) {
				   // check all entries in the present subsection
				   // whether it is a multiple entry
  map<String, pair<String,String>, less<String> >::const_iterator e;
  for (e = sec.entries.begin(); e != sec.entries.end(); ++e) 
    if ((*e).second.first.contains('{')) 
      multiple_choices.push_back (Entry(subsection_path,
					(*e).first,
					(*e).second.first));
	        

				   // transverse subsections
  map<String, Section*, less<String> >::const_iterator s;
  for (s = sec.subsections.begin(); s != sec.subsections.end(); ++s) 
    {
      enter_subsection ((*s).first);
      init_branches_section (*(*s).second);
      leave_subsection ();
    };
};


void MultipleParameterLoop::fill_entry_values (const unsigned int run_no) {
  int possibilities = 1;
  
  vector<Entry>::iterator choice;
  for (choice = multiple_choices.begin();
       choice != multiple_choices.end();
       ++choice) 
    {
				       // temporarily enter the subsection tree
				       // of this multiple entry
      subsection_path.swap ((*choice).subsection_path);

				       // set entry
      Section* pd = get_present_defaults_subsection ();
      int selection = (run_no/possibilities) % (*choice).different_values.size();
      String entry_value;
      if ((*choice).type == variant)
	entry_value = (*choice).different_values[selection];
      else 
	{
	  if (run_no>=(*choice).different_values.size()) 
	    {
	      cerr << "The given array for entry "
		<< pd->entries[(*choice).entry_name].first
		<< " does not conatin enough elements! Taking empty string instead." << endl;
	      entry_value = "";
	    }
	  else
	    entry_value = (*choice).different_values[run_no];
	};
      
				       // check conformance with regex
      if (!entry_value.matches (Regex(pd->entries[(*choice).entry_name].second)))
	{
	  cerr << "In run no.  " << run_no+1 << ":" << endl
	       << "    The entry value" << endl
	       << "        " << entry_value << endl
	       << "    for the entry named" << endl
	       << "        " << (*choice).entry_name << endl
	       << "    does not match the given regular expression" << endl
	       << "        " << pd->entries[(*choice).entry_name].second << endl
	       << "    Taking default value" << endl
	       << "        " << pd->entries[(*choice).entry_name].first << endl;
	  
					   // select default instead
	  entry_value = pd->entries[(*choice).entry_name].first;
	};

      Section* pc = get_present_changed_subsection ();
      pc->entries[(*choice).entry_name] = make_pair(entry_value, String(""));
            
				       // get out of subsection again
      subsection_path.swap ((*choice).subsection_path);

				       // move ahead if it was a variant entry
      if ((*choice).type == variant)
	possibilities *= (*choice).different_values.size();
    };
  
};

  

MultipleParameterLoop::Entry::Entry (const vector<String> &ssp,
				     const String& Name,
				     const String& Value) :
    subsection_path (ssp), entry_name(Name), entry_value(Value) {};



void MultipleParameterLoop::Entry::split_different_values () {
  Regex re_entry ("\\([^{]*\\){{? *\\([^}]*[^ }]\\) *}}?\\(.*\\)");
  int start, len;

				   // split string into three parts:
				   // part before the opening "{",
				   // the selection itself, final
				   // part after "}"
  String prefix, multiple, postfix;
  re_entry.match (entry_value, entry_value.length());

  re_entry.match_info (start, len, 1);
  prefix   = entry_value.at (start, len);

  re_entry.match_info (start, len, 2);
  multiple = entry_value.at (start, len);

  re_entry.match_info (start, len, 3);
  postfix  = entry_value.at (start, len);


				   // delete spaces around '|'
  multiple.gsub (" |", '|');
  multiple.gsub ("| ", '|');
  
  while (multiple.contains ('|')) 
    {
      different_values.push_back (prefix+multiple.before('|')+postfix);
      multiple.through ('|') = "";
    };
				   // make up the last selection ("while" broke
				   // because there was no '|' any more
  different_values.push_back (prefix+multiple+postfix);

				 // finally check whether this was a variant
				 // entry ({...}) or an array ({{...}})
  if (entry_value.contains(Regex(".*{{.*}}.*")))
    type = array;
  else
    type = variant;
};


