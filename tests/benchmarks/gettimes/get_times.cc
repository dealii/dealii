// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{

  ifstream input;
  stringstream ss;
  const int IGNORE_LINES = 3; //Number of initial lines to ignore
  const string DELIM = "*"; //String to divide test data

  input.open("temp.txt");

  vector<string> names;
  vector<double> times;

  //If no revision number, retrieve column names
  if (argc <= 1)
    {
      string curr_line;

      //Ignore first n lines + delimeter- they don't hold any useful data
      for (int i = 0; i < IGNORE_LINES + 1; i++)
        {
          getline(input, curr_line);
        }

      while (!input.eof())
        {

          getline(input, curr_line);
          if (curr_line == "")
            continue;
          if (curr_line == DELIM)
            break;

          ////cout << "curr line: " << curr_line << endl;
          //Looking for the string after the '|' and ' '
          int first_char = 0;
          if (curr_line[0] == '|')
            first_char++;
          if (curr_line[1] == ' ')
            first_char++;


          ////cout << "first char at pos: " << first_char << endl;

          //Find end of string
          int num_chars = 0;
          while (curr_line[first_char + num_chars] != '|')
            num_chars++;

          num_chars--;

          while (curr_line[first_char + num_chars] == ' ')
            num_chars--;

          names.push_back(curr_line.substr(first_char, num_chars+1));
        }
    }

  else //Else, extract execution time from each line
    {

      int time_index = 0;

      while (!input.eof())
        {

          string curr_line = "";

          //Read in line
          getline(input,curr_line);

          //Check for delimeter
          if (curr_line == DELIM)
            {
              ////cout << "Delimeter detected" << endl;
              time_index = 0;
              //Ignore first n lines- they don't hold any useful data
              for (int i = 0; i < IGNORE_LINES; i++)
                {
                  string dummy;
                  getline(input, dummy);
                  ////cout << "Skipping: " << dummy << endl;
                }
              continue;
            }

          ////cout << "Reading: " << curr_line << endl;

          if (curr_line == "" && DELIM != "")
            continue;


          //Looking for a number that ends with 's'
          int last_s = curr_line.rfind('s');
          //In case 's' is not used in the future
          assert(isdigit(curr_line[last_s - 1])); // Test: s preceded by number

          //Find time string and convert to double
          int num_start = last_s - 1;
          while (curr_line[num_start] != ' ')
            num_start--;
          string timestr = curr_line.substr(num_start+1, last_s - num_start);
          double time = (double)atof(timestr.substr(0,timestr.size()-1).c_str());

          assert(times.size() >= time_index);
          // first addition of times to vector; each loop, times.size() should be one less than time_index
          if (times.size() == time_index)
            times.push_back(time);
          else
            times[time_index] = (time < times[time_index]) ? time : times[time_index]; //else, determines minimum time and stores it
          time_index++;

        }
    }


  //Output individual names
  if (argc <= 1)
    {
      for (int i = 0; i < names.size(); i++)
        {
          cout << names[i];
          if (i < names.size() - 1)
            cout << endl;
        }
    }

  //Output individual times
  if (argc > 1)
    {
      cout << argv[1];
      for (int i = 0; i < times.size(); i++)
        cout << " " << times[i];
    }

  cout << endl;


  return 0;
}
