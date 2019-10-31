#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cctype> 
#include <sstream>
#include <regex>
using namespace std;

 
/*
 * It will iterate through all the lines in file and
 * put them in given vector
 */
bool getFileContent(std::string fileName, std::vector<std::string> & vecOfStrs)
{
 
	// Open the File
	std::ifstream in(fileName.c_str());
 
	// Check if object is valid
	if(!in)
	{
		std::cerr << "Cannot open the File : "<<fileName<<std::endl;
		return false;
	}
 
	std::string str;	
	// Read the next line from File untill it reaches the end.
	while (std::getline(in, str))
	{
		// Line contains string of length > 0 then save it in vector
		if(str.size() > 0)
			vecOfStrs.push_back(str);
	}
	//Close The File
	in.close();
	return true;
}
 
vector<string> process_line(string input)
{
	// remove all spaces in the line
	input.erase(std::remove(input.begin(),input.end(),' '),input.end());
	
	vector<string> line_array;
	std::stringstream ss(input);
	string str;
	while (getline(ss, str, ',')) 
	{
		line_array.push_back(str);
	}

	return line_array;
}
 
			
			
vector<int> node_blocks;

int main()
{
	std::vector<std::string> vecOfStr;
 
	// Get the contents of file in a vector
	bool result = getFileContent("1ElementTest_Dload_OP.inp", vecOfStr);
 	
	// Process keywords
	if(result)
	{
		int line_counter = 0;
		for(auto & raw_line : vecOfStr)
		{
			
			auto line_array = process_line(raw_line); // there is move semantics going on here				
			
			
			for ( auto word: line_array ) {
				
				//cout << "word = " << word << endl;
				//if ( word.find("*NODE") != string::npos && word.length() == 5 ) {
				//	node_blocks.push_back( line_counter );
				//} 
				
				// using regular expression
				if (std::regex_search(word,std::regex("(\\*NODE)\\b"))) {
					cout << "got node\n";
					node_blocks.push_back( line_counter );
				}
				
				if (std::regex_search(word,std::regex("(\\*NODEOUTPUT)\\b"))) {
					cout << "got node output\n";
					node_blocks.push_back( line_counter );
				}
				
			} // end loop
				
				
			line_counter += 1;
		
		} // end loop over lines
		
		
		// =====================
		// process node blocks
		// =====================
		for ( auto nblock: node_blocks )
		{
			bool BreakFlag = false;
			
			//cout << "block = " << nblock << endl;
			// header
			auto line_array = process_line(vecOfStr[nblock]);			
			
			for ( auto i = nblock+1; ; i++ )
			{
				auto line_array = process_line(vecOfStr[i]);
				
				for ( auto word: line_array ) 
				{
					if ( word.find("*") != string::npos ) { 
						BreakFlag = true;
						break;
					}
					cout << vecOfStr[i] << endl;
					
						
						
				}

				if ( BreakFlag ) break;
			} // end loop over node block
			
			
		}
		
			
	} // end if file open OK
	
	return 0;
}
	