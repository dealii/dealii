//Jordan Cazamias
//9-19-2012

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[]) {
	
	ifstream input;
	stringstream ss;
	const int IGNORE_LINES = 3; //Number of initial lines to ignore

	input.open("temp.txt");

	////cout << "Revision " << argv[1] << endl;
	
	//Ignore first n lines- they don't hold any useful data
	for(int i = 0; i < IGNORE_LINES; i++)
	{
		string dummy;
		getline(input, dummy);
	}

	vector<string> names;
	vector<double> times;
	
	if(argc <= 1) {
		while(!input.eof()) {
			//If no revision number, retrieve column names 
			string curr_line;
			getline(input, curr_line);
			if(curr_line == "")
				continue;
		
			////cout << "curr line: " << curr_line << endl;
			//Looking for the string after the '|' and ' '
			int first_char = 0;
			if(curr_line[0] == '|')
				first_char++;
			if (curr_line[1] == ' ')
				first_char++;

			////cout << "first char at pos: " << first_char << endl;
	
			//Test case: char at position first_char is a letter
			assert(isalpha(curr_line[first_char]));
			
			//Find position of first space that occurs after the word
			int num_chars = 0;
			while(!isspace(curr_line[first_char + num_chars]))
				num_chars++;
	
			names.push_back(curr_line.substr(first_char, num_chars));
		}
	}


	//Extract and execution time from each line
	while(!input.eof())
	{
		string curr_line;
		getline(input, curr_line);
		////cout << curr_line << endl;
		if(curr_line == "")
			continue;
		
		//Looking for a number that ends with 's'
		int last_s = curr_line.rfind('s');
		//In case 's' is not used in the future
		assert(isdigit(curr_line[last_s - 1])); // Test: s preceded by number
		
		//Find time string and convert to double
		int num_start = last_s - 1;
		while(curr_line[num_start] != ' ')//isdigit(curr_line[num_start]) || curr_line[num_start] == '.')
			num_start--;
		string timestr = curr_line.substr(num_start+1, last_s - num_start);
		////cout << timestr << endl;
		double time = (double)atof(timestr.substr(0,timestr.size()-1).c_str());
		////cout << time << endl;
		times.push_back(time);
	}


	//Output individual names
	if(argc <= 1) {
		for(int i = 0; i < names.size(); i++) {
			cout << names[i];
			if(i < names.size() - 1)
				cout << endl;
		}
	}

	//Output individual times
	if(argc > 1) {
		cout << argv[1]; 
		for(int i = 0; i < times.size(); i++)
			cout << " " << times[i];
	}

	cout << endl;

	
	return 0;
}
