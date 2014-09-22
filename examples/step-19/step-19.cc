/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2005 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Luca Heltai, Wolfgang Bangerth, 2005
 */


// @sect4{Preliminaries}

// As usual, we start with include files. This program is content with really
// few of these -- we only need two files from the library (one for input and
// output of graphical data, one for parameter handling), and a few C++
// standard headers:
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_handler.h>

#include <list>
#include <iostream>
#include <fstream>

// As mentioned in the first few tutorial programs, all names in deal.II are
// declared in a namespace <code>dealii</code>. To make using these function
// and class names simpler, we import the entire content of that namespace
// into the global scope. As done for all previous programs already, we'll
// also place everything we do here into a namespace of its own:
namespace Step19
{
  using namespace dealii;

  // Before we start with the actual program, let us declare a few global
  // variables that will be used to hold the parameters this program is going
  // to use. Usually, global variables are frowned upon for a good reason, but
  // since we have such a short program here that does only a single thing, we
  // may stray from our usual line and make these variables global, rather
  // than passing them around to all functions or encapsulating them into a
  // class.
  //
  // The variables we have are: first, an object that will hold parameters of
  // operation, such as output format (unless given on the command line);
  // second, the names of input and output files; and third, the format in
  // which the output is to be written:
  ParameterHandler         prm;
  std::vector<std::string> input_file_names;
  std::string              output_file;
  std::string              output_format;


  // All the stuff this program does can be done from here on. As described in
  // the introduction, what we have to do is declare what values the parameter
  // file can have, parse the command line, read the input files, then write
  // the output. We will do this in this order of operation, but before that
  // let us declare a function that prints a message about how this program is
  // to be used; the function first prints a general message, and then goes on
  // to list the parameters that are allowed in the parameter file (the
  // <code>ParameterHandler</code> class has a function to do exactly this;
  // see the results section for what it prints):
  void
  print_usage_message ()
  {
    static const char *message
      =
        "\n"
        "Converter from deal.II intermediate format to other graphics formats.\n"
        "\n"
        "Usage:\n"
        "    ./step-19 [-p parameter_file] list_of_input_files \n"
        "              [-x output_format] [-o output_file]\n"
        "\n"
        "Parameter sequences in brackets can be omitted if a parameter file is\n"
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


  // @sect4{Declaring parameters for the input file}

  // The second function is used to declare the parameters this program
  // accepts from the input file. While we don't actually take many parameters
  // from the input file except for, possibly, the output file name and
  // format, we nevertheless want to show how to work with parameter files.
  //
  // In short, the <code>ParameterHandler</code> class works as follows: one
  // declares the entries of parameters that can be given in input files
  // together, and later on one can read an input file in which these
  // parameters are set to their values. If a parameter is not listed in the
  // input file, the default value specified in the declaration of that
  // parameter is used. After that, the program can query the values assigned
  // to certain parameters from the <code>ParameterHandler</code> object.
  //
  // Declaring parameters can be done using the
  // <code>ParameterHandler::declare_entry</code> function. It's arguments are
  // the name of a parameter, a default value (given as a string, even if the
  // parameter is numeric in nature, and thirdly an object that describes
  // constraints on values that may be passed to this parameter. In the
  // example below, we use an object of type <code>Patterns::Anything</code>
  // to denote that there are no constraints on file names (this is, of
  // course, not true -- the operating system does have constraints, but from
  // an application standpoint, almost all names are valid). In other cases,
  // one may, for example, use <code>Patterns::Integer</code> to make sure
  // that only parameters are accepted that can be interpreted as integer
  // values (it is also possible to specify bounds for integer values, and all
  // values outside this range are rejected), <code>Patterns::Double</code>
  // for floating point values, classes that make sure that the given
  // parameter value is a comma separated list of things, etc. Take a look at
  // the <code>Patterns</code> namespace to see what is possible.
  //
  // The fourth argument to <code>declare_entry</code> is a help string that
  // can be printed to document what this parameter is meant to be used for
  // and other information you may consider important when declaring this
  // parameter. The default value of this fourth argument is the empty string.
  //
  // I always wanted to have an example program describing the
  // <code>ParameterHandler</code> class, because it is so particularly
  // useful. It would have been useful in a number of previous example
  // programs (for example, in order to let the tolerance for linear solvers,
  // or the number of refinement steps be determined by a run-time parameter,
  // rather than hard-coding them into the program), but it turned out that
  // trying to explain this class there would have overloaded them with things
  // that would have distracted from the main purpose. However, while writing
  // this program, I realized that there aren't all that many parameters this
  // program can usefully ask for, or better, it turned out: declaring and
  // querying these parameters was already done centralized in one place of
  // the library, namely the <code>DataOutInterface</code> class that handles
  // exactly this -- managing parameters for input and output.
  //
  // So the second function call in this function is to let the
  // <code>DataOutInterface</code> declare a good number of parameters that
  // control everything from the output format to what kind of output should
  // be generated if output is written in a specific graphical format. For
  // example, when writing data in encapsulated postscript (EPS) format, the
  // result is just a 2d projection, not data that can be viewed and rotated
  // with a viewer. Therefore, one has to choose the viewing angle and a
  // number of other options up front, when output is generated, rather than
  // playing around with them later on. The call to
  // <code>DataOutInterface::declare_parameters</code> declares entries that
  // allow to specify them in the parameter input file during run-time. If the
  // parameter file does not contain entries for them, defaults are taken.
  //
  // As a final note: <code>DataOutInterface</code> is a template, because it
  // is usually used to write output for a specific space dimension. However,
  // this program is supposed to be used for all dimensions at the same time,
  // so we don't know at compile time what the right dimension is when
  // specifying the template parameter. Fortunately, declaring parameters is
  // something that is space dimension independent, so we can just pick one
  // arbitrarily. We pick <code>1</code>, but it could have been any other
  // number as well.
  void declare_parameters ()
  {
    prm.declare_entry ("Output file", "",
                       Patterns::Anything(),
                       "The name of the output file to be generated");

    DataOutInterface<1>::declare_parameters (prm);

    // Since everything that this program can usefully request in terms of
    // input parameters is already handled by now, let us nevertheless show
    // how to use input parameters in other circumstances. First, parameters
    // are like files in a directory tree: they can be in the top-level
    // directory, but you can also group them into subdirectories to make it
    // easier to find them or to be able to use the same parameter name in
    // different contexts.
    //
    // Let us first declare a dummy parameter in the top-level section; we
    // assume that it will denote the number of iterations, and that useful
    // numbers of iterations that a user should be able to specify are in the
    // range 1...1000, with a default value of 42:
    prm.declare_entry ("Dummy iterations", "42",
                       Patterns::Integer (1,1000),
                       "A dummy parameter asking for an integer");

    // Next, let us declare a sub-section (the equivalent to a
    // subdirectory). When entered, all following parameter declarations will
    // be within this subsection. To also visually group these declarations
    // with the subsection name, I like to use curly braces to force my editor
    // to indent everything that goes into this sub-section by one level of
    // indentation. In this sub-section, we shall have two entries, one that
    // takes a Boolean parameter and one that takes a selection list of
    // values, separated by the '|' character:
    prm.enter_subsection ("Dummy subsection");
    {
      prm.declare_entry ("Dummy generate output", "true",
                         Patterns::Bool(),
                         "A dummy parameter that can be fed with either "
                         "'true' or 'false'");
      prm.declare_entry ("Dummy color of output", "red",
                         Patterns::Selection("red|black|blue"),
                         "A dummy parameter that shows how one can define a "
                         "parameter that can be assigned values from a finite "
                         "set of values");
    }
    prm.leave_subsection ();
    // After this, we have left the subsection again. You should have gotten
    // the idea by now how one can nest subsections to separate
    // parameters. There are a number of other possible patterns describing
    // possible values of parameters; in all cases, if you try to pass a
    // parameter to the program that does not match the expectations of the
    // pattern, it will reject the parameter file and ask you to fix it. After
    // all, it does not make much sense if you had an entry that contained the
    // entry "red" for the parameter "Generate output".
  }


  // @sect4{Parsing the command line}

  // Our next task is to see what information has been provided on the command
  // line. First, we need to be sure that there is at least one parameter: an
  // input file. The format and the output file can be specified in the
  // parameter file, but the list of input files can't, so at least one
  // parameter needs to be there. Together with the name of the program (the
  // zeroth parameter), <code>argc</code> must therefore be at least 2. If
  // this is not the case, we print an error message and exit:
  void
  parse_command_line (const int     argc,
                      char *const *argv)
  {
    if (argc < 2)
      {
        print_usage_message ();
        exit (1);
      }

    // Next, collect all parameters in a list that will be somewhat simpler to
    // handle than the <code>argc</code>/<code>argv</code> mechanism. We omit
    // the name of the executable at the zeroth index:
    std::list<std::string> args;
    for (int i=1; i<argc; ++i)
      args.push_back (argv[i]);

    // Then process all these parameters. If the parameter is <code>-p</code>,
    // then there must be a parameter file following (which we should then
    // read), in case of <code>-x</code> it is the name of an output
    // format. Finally, for <code>-o</code> it is the name of the output
    // file. In all cases, once we've treated a parameter, we remove it from
    // the list of parameters:
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
            const std::string parameter_file = args.front ();
            args.pop_front ();

            // Now read the input file:
            prm.read_input (parameter_file);

            // Both the output file name as well as the format can be
            // specified on the command line. We have therefore given them
            // global variables that hold their values, but they can also be
            // set in the parameter file. We therefore need to extract them
            // from the parameter file here, because they may be overridden by
            // later command line parameters:
            if (output_file == "")
              output_file = prm.get ("Output file");

            if (output_format == "")
              output_format = prm.get ("Output format");

            // Finally, let us note that if we were interested in the values
            // of the parameters declared above in the dummy subsection, we
            // would write something like this to extract the value of the
            // Boolean flag (the <code>prm.get</code> function returns the
            // value of a parameter as a string, whereas the
            // <code>prm.get_X</code> functions return a value already
            // converted to a different type):
            prm.enter_subsection ("Dummy subsection");
            {
              prm.get_bool ("Dummy generate output");
            }
            prm.leave_subsection ();
            // We would assign the result to a variable, of course, but don't
            // here in order not to generate an unused variable that the
            // compiler might warn about.
            //
            // Alas, let's move on to handling of output formats:
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

        // Otherwise, this is not a parameter that starts with a known minus
        // sequence, and we should consider it to be the name of an input
        // file. Let us therefore add this file to the list of input files:
        else
          {
            input_file_names.push_back (args.front());
            args.pop_front ();
          }
      }

    // Next check a few things and create errors if the checks fail. Firstly,
    // there must be at least one input file
    if (input_file_names.size() == 0)
      {
        std::cerr << "Error: No input file specified." << std::endl;
        print_usage_message ();
        exit (1);
      }
  }


  // @sect4{Generating output}

  // Now that we have all the information, we need to read all the input
  // files, merge them, and generate a single output file. This, after all,
  // was the motivation, borne from the necessity encountered in the step-18
  // tutorial program, to write this program in the first place.
  //
  // So what we do first is to declare an object into which we will merge the
  // data from all the input file, and read in the first file through a
  // stream. Note that every time we open a file, we use the
  // <code>AssertThrow</code> macro to check whether the file is really
  // readable -- if it isn't then this will trigger an exception and
  // corresponding output will be generated from the exception handler in
  // <code>main()</code>:
  template <int dim, int spacedim>
  void do_convert ()
  {
    DataOutReader<dim,spacedim> merged_data;

    {
      std::ifstream input (input_file_names[0].c_str());
      AssertThrow (input, ExcIO());

      merged_data.read (input);
    }

    // For all the other input files, we read their data into an intermediate
    // object, and then merge that into the first object declared above:
    for (unsigned int i=1; i<input_file_names.size(); ++i)
      {
        std::ifstream input (input_file_names[i].c_str());
        AssertThrow (input, ExcIO());

        DataOutReader<dim,spacedim> additional_data;
        additional_data.read (input);
        merged_data.merge (additional_data);
      }

    // Once we have this, let us open an output stream, and parse what we got
    // as the name of the output format into an identifier. Fortunately, the
    // <code>DataOutBase</code> class has a function that does this parsing
    // for us, i.e. it knows about all the presently supported output formats
    // and makes sure that they can be specified in the parameter file or on
    // the command line. Note that this ensures that if the library acquires
    // the ability to output in other output formats, this program will be
    // able to make use of this ability without having to be changed!
    std::ofstream output_stream (output_file.c_str());
    AssertThrow (output_stream, ExcIO());

    const DataOutBase::OutputFormat format
      = DataOutBase::parse_output_format (output_format);

    // Finally, write the merged data to the output:
    merged_data.write(output_stream, format);
  }


  // @sect4{Dispatching output generation}

  // The function above takes template parameters relating to the space
  // dimension of the output, and the dimension of the objects to be
  // output. (For example, when outputting whole cells, these two dimensions
  // are the same, but the intermediate files may contain only data pertaining
  // to the faces of cells, in which case the first parameter will be one less
  // than the space dimension.)
  //
  // The problem is: at compile time, we of course don't know the dimensions
  // used in the input files. We have to plan for all cases, therefore. This
  // is a little clumsy, since we need to specify the dimensions statically at
  // compile time, even though we will only know about them at run time.
  //
  // So here is what we do: from the first input file, we determine (using a
  // function in <code>DataOutBase</code> that exists for this purpose) these
  // dimensions. We then have a series of switches that dispatch, statically,
  // to the <code>do_convert</code> functions with different template
  // arguments. Not pretty, but works. Apart from this, the function does
  // nothing -- except making sure that it covered the dimensions for which it
  // was called, using the <code>AssertThrow</code> macro at places in the
  // code that shouldn't be reached:
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
}



// @sect4{main()}

// Finally, the main program. There is not much more to do than to make sure
// parameters are declared, the command line is parsed (which includes reading
// parameter files), and finally making sure the input files are read and
// output is generated. Everything else just has to do with handling
// exceptions and making sure that appropriate output is generated if one is
// thrown.
int main (int argc, char **argv)
{
  try
    {
      using namespace Step19;

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
