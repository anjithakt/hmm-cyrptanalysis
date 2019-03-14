#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <map>
#include <fstream>
#include <sstream>


using namespace std;

#include <boost/algorithm/string.hpp>    


string get_file_contents(char* filename) {
	cout << "Opening " << filename << "\n";
	ifstream in(filename, ios::in);
	std::string str((std::istreambuf_iterator<char>(in)),
                 std::istreambuf_iterator<char>());
	return str;
}

bool is_not_valid_state(char s) {
	return (s < 'a' || s >'z') && s != ' ';
}


string sanitize_observation(string observation) {
	boost::algorithm::to_lower(observation);
	observation.erase(std::remove_if(observation.begin(), observation.end(), &is_not_valid_state), observation.end());
	return observation;
}

double** array_alloc(int I, int J) {
	double** arr = new double*[I];
	for(int i =0; i <I; i++) {
		arr[i] = new double[J];
	}
	return arr;
}


int char_to_index(char i) {
	if (i == ' ') 
		return 0;
	else return 1 + i - 'a';
}




double** get_digraph_matrix(char *filename) {
	string line;
	double** digraph=array_alloc(27,27);
	for(int i=0; i< 27; ++i) {
		for(int j=0; j<27; ++j) {
			digraph[i][j] = 5;
		}
	}

	ifstream in(filename, ios::in);
	while(getline(in, line)) {
		line = sanitize_observation(line);
		for(int i=1; i < line.length(); ++i) {
			digraph[char_to_index(line[i-1])][char_to_index(line[i])]++;
		}
	}

	for(int i=0; i< 27; ++i) {
		double sum =0;
		for(int j=0; j<27; ++j) {
			sum += digraph[i][j];
		}
		for(int j=0; j< 27; ++j) {
			digraph[i][j] = digraph[i][j] / sum;
		}
	}
	return digraph;
}

void print_matrix(double** &arr, int n, int m) {
	cout << "\n";
	for(int i=0; i < n ; i++) {
		for (int j =0; j <m; j++) {
			cout << arr[i][j] <<" ";
		}
		cout << "\n";
	}
}

int main(int argc, char *argv[]) {
	srand(0);
	double** digraph_matrix = get_digraph_matrix(argv[1]);
	print_matrix(digraph_matrix, 27,27);

}