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

// Base libraries
#include <base/data_out_base.h>

// C++ libraries
#include <iostream>
#include <fstream>




template <int dim, int spacedim>
void do_convert (const std::vector<std::string> &input_files,
                 const std::string              &output_file,
                 const std::string              &/*parameter_file*/,
                 const std::string              &output_format_extension)
{
  std::ifstream intermediate_stream (input_files[0].c_str());
  AssertThrow (intermediate_stream, ExcIO());
  
  DataOutReader<dim,spacedim> intermediate_data;  
  intermediate_data.read (intermediate_stream);

				   // for all following input files
  for (unsigned int i=1; i<input_files.size(); ++i)
    {
      std::ifstream additional_stream (input_files[i].c_str());
      AssertThrow (additional_stream, ExcIO());

      DataOutReader<dim,spacedim> additional_data;  
      additional_data.read (additional_stream);
      intermediate_data.merge (additional_data);
    } 
  
  std::ofstream output_stream (output_file.c_str());
  AssertThrow (output_stream, ExcIO());

  const typename DataOutInterface<dim,spacedim>::OutputFormat
    output_format
    = intermediate_data.parse_output_format (output_format_extension);
  
  intermediate_data.write(output_stream, output_format); 
}


void convert (const std::vector<std::string> &input_files,
              const std::string              &output_file,
              const std::string              &parameter_file,
              const std::string              &output_format_extension)
{
  std::ifstream input(input_files[0].c_str());
  AssertThrow (input, ExcIO());
  
  const std::pair<unsigned int, unsigned int>
    dimensions = DataOutBase::determine_intermediate_format_dimensions (input);

  switch (dimensions.first)
    {
      case 1:
            switch (dimensions.second)
              {
                case 1:
                      do_convert <1,1> (input_files,
                                        output_file,
                                        parameter_file,
                                        output_format_extension);
                      return;
                      
                case 2:
                      do_convert <1,2> (input_files,
                                        output_file,
                                        parameter_file,
                                        output_format_extension);
                      return;
              }
            AssertThrow (false, ExcNotImplemented());
            
      case 2:
            switch (dimensions.second)
              {
                case 2:
                      do_convert <2,2> (input_files,
                                        output_file,
                                        parameter_file,
                                        output_format_extension);
                      return;

                case 3:
                      do_convert <2,3> (input_files,
                                        output_file,
                                        parameter_file,
                                        output_format_extension);
                      return;
              }
            AssertThrow (false, ExcNotImplemented());

      case 3:
            switch (dimensions.second)
              {
                case 3:
                      do_convert <3,3> (input_files,
                                        output_file,
                                        parameter_file,
                                        output_format_extension);
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
      if (argc == 1)
        {
          std::cout << std::endl
		    << "Converter from deal.II intermediate format to "
                    << "other graphics formats."
                    << std::endl << std::endl
                    << "Usage: ./step-19 [-p parameter_file] "
                    << "list_of_input_files -x output_format "
                    << "output_file"
                    << std::endl
		    << std::endl;
          exit (1);
        }
      
      std::string              parameter_file;
      std::string              output_format_extension;
      std::vector<std::string> file_names;
      
      for (int i=1; i<argc; ++i)
        {
          if (std::string(argv[i]) == std::string("-p"))
            {
              if (i+1 == argc)
                {
                  std::cerr << "Error: flag '-p' must be followed by the "
                            << "name of a parameter file."
                            << std::endl;
                  exit (1);
                }
              parameter_file = argv[i+1];
              ++i;
            }
          else if (std::string(argv[i]) == std::string("-x"))
            {
              if (i+1 == argc)
                {
                  std::cerr << "Error: flag '-x' must be followed by the "
                            << "name of an output format."
                            << std::endl;
                  exit (1);
                }
              output_format_extension = argv[i+1];
              ++i;
            }
          else
            file_names.push_back (argv[i]);
        }

      if (file_names.size() < 2)
        {
          std::cerr << "Error: No input and/or output files specified."
                    << std::endl;
          exit (1);
        }

      const std::vector<std::string> input_files (file_names.begin(),
                                                  file_names.end()-1);
      const std::string              output_file (file_names.back());

      if (output_format_extension.size() == 0)
        {
          output_format_extension
            = std::string(output_file.begin() + output_file.rfind ('.') + 1,
                          output_file.end());

          if (output_format_extension.size() == 0)
            {
              std::cerr << "Error: could not determine output file format "
                        << "from name of last file name."
                        << std::endl;
              exit (1);
            }
        }

      convert (input_files,
               output_file,
               parameter_file,
               output_format_extension);
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
