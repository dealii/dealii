/* $Id$ */
/* Author: Luca Heltai, Wolfgang Bangerth, 2005 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2005 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


                                 // As usual, we start with include
                                 // files. This program is content with really
                                 // few of these -- we only need two files
                                 // from the library (one for input and output
                                 // of graphical data, one for parameter
                                 // handling), and a few C++ standard headers:
#include <base/data_out_base.h>
#include <base/parameter_handler.h>

#include <list>
#include <iostream>
#include <fstream>


                                 // Before we start with the actual program,
                                 // let us declare a few global variables that
                                 // will be used to hold the parameters this
                                 // program is going to use. Usually, global
                                 // variables are frowned upon for a good
                                 // reason, but since we have such a short
                                 // program here that does only a single
                                 // thing, we may stray from our usual line
                                 // and make these variables global, rather
                                 // than passing them around to all functions
                                 // or encapsulating them into a class.
                                 //
                                 // The variables we have are: first, an
                                 // object that will hold parameters of
                                 // operation, such as output format (unless
                                 // given on the command line); second, the
                                 // names of input and output files; third,
                                 // the format in which the output is to be
                                 // written; and fourth the name of a
                                 // parameter file:
ParameterHandler         prm;
std::vector<std::string> input_file_names;
std::string              output_file; 
std::string              output_format;
std::string              parameter_file;


                                 // All the stuff this program does can be
                                 // done from here on. As described in the
                                 // introduction, what we have to do is
                                 // declare what values the parameter file can
                                 // have, parse the command line, read the
                                 // input files, then write the output. We
                                 // will do this in this order of operation,
                                 // but before that let us declare a function
                                 // that prints a message about how this
                                 // program is to be used; the function first
                                 // prints a general message, and then goes on
                                 // to list the parameters that are allowed in
                                 // the parameter file (the
                                 // ``ParameterHandler'' class has a function
                                 // to do exactly this; see the results
                                 // section for what it prints):
void
print_usage_message ()
{
  static const char* message
    =
    "\n"
    "Converter from deal.II intermediate format to other graphics formats.\n"
    "\n"
    "Usage:\n"
    "    ./step-19 [-p parameter_file] list_of_input_files \n"
    "              [-x output_format] [-o output_file]\n"
    "\n"
    "Parameter sequences in brackets can be omitted a parameter file is\n"
    "specified on the command line and if it provides values for these\n"
    "missing parameters.\n"
    "\n"
    "The parameter file has the following format and allows the following\n"
    "values (you can cut and paste this and use it for your own parameter\n"
    "file):\n"
    "\n";
  std::cout << message;

  prm.print_parameters (std::cout, ParameterHandler::Text);
}


                                 // The second function is used to declare the
                                 // parameters this program accepts from the
                                 // input file. While we don't actually take
                                 // many parameters from the input file except
                                 // for, possibly, the output file name and
                                 // format, we nevertheless want to show how
                                 // to work with parameter files.
                                 //
                                 // In short, the ``ParameterHandler'' class
                                 // works as follows: one declares the entries
                                 // of parameters that can be given in input
                                 // files together, and later on one can read
                                 // an input file in which these parameters
                                 // are set to their values. If a parameter is
                                 // not listed in the input file, the default
                                 // value specified in the declaration of that
                                 // parameter is used. After that, the program
                                 // can query the values assigned to certain
                                 // parameters from the ``ParameterHandler''
                                 // object.
                                 //
                                 // Declaring parameters can be done using the
                                 // ``ParameterHandler::declare_entry''
                                 // function. It's arguments are the name of a
                                 // parameter, a default value (given as a
                                 // string, even if the parameter is numeric
                                 // in nature, and thirdly an object that
                                 // describes constraints on values that may
                                 // be passed to this parameter. In the
                                 // example below, we use an object of type
                                 // ``Patterns::Anything'' to denote that
                                 // there are no constraints on file names
                                 // (this is, of course, not true -- the
                                 // operating system does have constraints,
                                 // but from an application standpoint, almost
                                 // all names are valid). In other cases, one
                                 // may, for example, use
                                 // ``Patterns::Integer'' to make sure that
                                 // only parameters are accepted that can be
                                 // interpreted as integer values (it is also
                                 // possible to specify bounds for integer
                                 // values, and all values outside this range
                                 // are rejected), ``Patterns::Double'' for
                                 // floating point values, classes that make
                                 // sure that the given parameter value is a
                                 // comma separated list of things, etc. Take
                                 // a look at the ``Patterns'' namespace to
                                 // see what is possible.
                                 //
                                 // The fourth argument to ``declare_entry''
                                 // is a help string that can be printed to
                                 // document what this parameter is meant to
                                 // be used for and other information you may
                                 // consider important when declaring this
                                 // parameter. The default value of this
                                 // fourth argument is the empty string.
void declare_parameters ()
{
  prm.declare_entry ("Output file", "",
                     Patterns::Anything(),
                     "The name of the output file to be generated");

  DataOutInterface<1>::declare_parameters (prm);
}

  

void
parse_command_line (const int                argc,
                    char             *const* argv)
{
  if (argc < 2)
    {
      print_usage_message ();
      exit (0);
    }

  std::list<std::string> args;
  for (int i=1; i<argc; ++i)
    args.push_back (argv[i]);
  
  while (args.size())
    {
      if (args.front() == std::string("-p"))
        {
          if (args.size() == 1)
            {
              std::cerr << "Error: flag '-p' must be followed by the "
                        << "name of a parameter file."
                        << std::endl;
              print_usage_message ();
              exit (1);
            }
          args.pop_front ();
          parameter_file = args.front ();
          args.pop_front ();
        }
      else if (args.front() == std::string("-x"))
        {
          if (args.size() == 1)
            {
              std::cerr << "Error: flag '-x' must be followed by the "
                        << "name of an output format."
                        << std::endl;
              print_usage_message ();
              exit (1);
            }
          args.pop_front ();
          output_format = args.front();
          args.pop_front ();
        }
      else if (args.front() == std::string("-o"))
        {
          if (args.size() == 1)
            {
              std::cerr << "Error: flag '-o' must be followed by the "
                        << "name of an output file."
                        << std::endl;
              print_usage_message ();
              exit (1);
            }
          args.pop_front ();
          output_file = args.front();
          args.pop_front ();
        }
      else
        {
          input_file_names.push_back (args.front());
          args.pop_front ();
        }
    }

  if (input_file_names.size() == 0)
    {
      std::cerr << "Error: No input file specified." << std::endl;
      print_usage_message ();
      exit (1);
    }
}



template <int dim, int spacedim>
void do_convert ()
{
  std::ifstream intermediate_stream (input_file_names[0].c_str());
  AssertThrow (intermediate_stream, ExcIO());
  
  DataOutReader<dim,spacedim> intermediate_data;  
  intermediate_data.read (intermediate_stream);

				   // for all following input files
  for (unsigned int i=1; i<input_file_names.size(); ++i)
    {
      std::ifstream additional_stream (input_file_names[i].c_str());
      AssertThrow (additional_stream, ExcIO());

      DataOutReader<dim,spacedim> additional_data;  
      additional_data.read (additional_stream);
      intermediate_data.merge (additional_data);
    } 
  
  std::ofstream output_stream (output_file.c_str());
  AssertThrow (output_stream, ExcIO());

  const DataOutBase::OutputFormat format
    = intermediate_data.parse_output_format (output_format);
  
  intermediate_data.write(output_stream, format); 
}



void convert ()
{
  AssertThrow (input_file_names.size() > 0,
               ExcMessage ("No input files specified."));
  
  std::ifstream input(input_file_names[0].c_str());
  AssertThrow (input, ExcIO());
  
  const std::pair<unsigned int, unsigned int>
    dimensions = DataOutBase::determine_intermediate_format_dimensions (input);

  switch (dimensions.first)
    {
      case 1:
            switch (dimensions.second)
              {
                case 1:
                      do_convert <1,1> ();
                      return;
                      
                case 2:
                      do_convert <1,2> ();
                      return;
              }
            AssertThrow (false, ExcNotImplemented());
            
      case 2:
            switch (dimensions.second)
              {
                case 2:
                      do_convert <2,2> ();
                      return;

                case 3:
                      do_convert <2,3> ();
                      return;
              }
            AssertThrow (false, ExcNotImplemented());

      case 3:
            switch (dimensions.second)
              {
                case 3:
                      do_convert <3,3> ();
                      return;
              }
            AssertThrow (false, ExcNotImplemented());
    }
  
  AssertThrow (false, ExcNotImplemented());
}



int main (int argc, char ** argv)
{
  try
    {
      declare_parameters ();
      parse_command_line (argc, argv);

      convert ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
 
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
 
  return 0;
}                                
