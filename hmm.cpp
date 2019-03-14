#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <map>
#include <fstream>
#include <sstream>


#include <boost/algorithm/string.hpp>    


using namespace std;


set<char> get_states(string observation) {
	set<char> state_set(observation.begin(), observation.end());
	return state_set;
}



double** array_alloc(int I, int J) {
	double** arr = new double*[I];
	for(int i =0; i <I; i++) {
		arr[i] = new double[J];
	}
	return arr;
}


double*** array_alloc(int I, int J, int K) {
	double*** arr = new double**[I];
	for(int i =0; i <I; i++) {
		arr[i] = array_alloc(J, K);
	}
	return arr;
}

void array_dealloc(double** &array, int N, int M) {
	for(int i =0; i <M ; i++) {
		delete array[i];
	}
	delete array;
	array = NULL;
}


double approx_value(double value) {
	return value + ( value * .001  * ((rand() % 1000) - 500) / 500);
}


double** row_stochastic(int N, int M) {
	double base_value=1.0/M;

	double** arr = array_alloc(N,M);

	for(int i =0; i < N; i++) {
		double sum_so_far =0;
		for(int j =0; j < M-1; j++) {
			arr[i][j] = approx_value(base_value);
			sum_so_far += arr[i][j];
		}
		arr[i][M-1] = 1 - sum_so_far;
	}

	return arr;
}


bool is_not_valid_state(char s) {
	return (s < 'a' || s >'z') && s != ' ';
}

string sanitize_observation(string observation) {
	boost::algorithm::to_lower(observation);
	observation.erase(std::remove_if(observation.begin(), observation.end(), &is_not_valid_state), observation.end());
	return observation;
}



map<char, int> get_state_map(set<char> state_set) {
	map<char, int> states;
	int index=0;
	for(set<char>::iterator it=state_set.begin(); it!= state_set.end(); ++it) {
		states[*it]=index++;
	}
	return states;
}


//***************** Debug utils **************

void print_matrix(double** &arr, int n, int m) {
	cout << "\n";
	for(int i=0; i < n ; i++) {
		for (int j =0; j <m; j++) {
			cout << arr[i][j] <<" ";
		}
		cout << "\n";
	}
}

void print_map(map<char, int> &myMap) {
	for(map<char, int>::const_iterator it = myMap.begin(); it != myMap.end(); ++it) {
		cout << it -> first <<" -> " << it -> second <<"\n";
	}
	cout << "\n";
}

string get_file_contents(char* filename) {
	cout << "Opening " << filename << "\n";
	ifstream in(filename, ios::in);
	std::string str((std::istreambuf_iterator<char>(in)),
                 std::istreambuf_iterator<char>());
	return str;
}
// ***************** END Debug Utils *****


class HmmState {

public:
	int N;
	int M;
	int T;
	string observation;

	double** A;
	double** B;
	double* pi;
	double* c;

	double** alpha;
	double** beta;
	double** gamma;
	double*** digamma;

	int MAX_ITERS= 100000;
	int MIN_ITERS =100;

	double EPISON = 0.0000001;

	map<char, int> states;

	// public:
	void init(int n, string &obs) {
		N = n;
		observation = sanitize_observation(obs);
		T= observation.size();
		
		set<char> state_set = get_states(observation);
		M = state_set.size();
		states = get_state_map(state_set);

		print_map(states);

		A = row_stochastic(N, N);
		B = row_stochastic(N, M);
		double** pih = row_stochastic(1, N);
		pi=pih[0];
		//delete pih;

		c = new double[T];
		alpha = array_alloc(T, N);
		beta = array_alloc(T, N);
		gamma = array_alloc(T, N);
		digamma = array_alloc(T, N, N);
	}

	void forward_pass() {
		c[0] =0;
		for(int i =0; i < N; ++i) {
			alpha[0][i] = pi[i] * B[i][states[observation[0]]];
			c[0]+=alpha[0][i];
		}
		c[0] = 1.0/ c[0];
		for(int i =0; i < N; ++i) {
			alpha[0][i] *= c[0];
		}

		for(int t=1; t < T; ++t) {
			c[t] =0;
			for(int i=0; i<N; ++i) {
				alpha[t][i] = 0;
				for(int j=0; j<N; ++j) {
					alpha[t][i] += alpha[t-1][j] * A[j][i];
				}
				alpha[t][i] *= B[i][states[observation[t]]];
				c[t] += alpha[t][i];
			}
			c[t] = 1.0/c[t];
			for(int i=0; i<N; ++i) {
				alpha[t][i] *= c[t];
			}
		}
	}


	void backward_pass() {
		for(int i=0; i<N; ++i) {
			beta[T-1][i] = c[T-1];
		}
		for(int t = T-2; t>=0 ; --t) {
			for(int i =0; i < N; ++i) {
				beta[t][i] =0;
				for(int j=0; j<N ; ++j) {
					beta[t][i] += A[i][j] * B[j][states[observation[t+1]]] * beta[t+1][j];
				}
				beta[t][i] *= c[t];
			}
		}
	}

	void compute_gamma() {
		for(int t=0; t<T-1; ++t){
			for(int i =0; i< N; i++) {
				gamma[t][i] = 0;
				for(int j=0; j <N ; j++) {
					digamma[t][i][j] = alpha[t][i] * A[i][j] * B[j][states[observation[t+1]]] * beta[t+1][j];
					gamma[t][i] += digamma[t][i][j];
				}
			}
		}
		for(int i=0; i < N; i++) {
			gamma[T-1][i] = alpha[T-1][i];
		}
	}


	void re_estimate() {
		for (int i =0; i < N; ++i) {
			pi[i] = gamma[0][i];
		}

		for(int i=0; i < N; ++i) {
			double denom =0;
			for(int t=0; t<T-1 ; ++t) {
				denom += gamma[t][i];
			}
			for(int j=0; j <N; ++j) {
				double numer = 0;
				for(int t=0; t < T-1; ++t) {
					numer += digamma[t][i][j];
				}
				A[i][j] = numer / denom;
			}
		}

		for(int i=0; i <N; ++i) {
			double denom =0;
			for(int t =0; t< T; ++t) {
				denom += gamma[t][i];
			}
			for(int j=0; j< M; ++j) {
				double numer =0;
				for(int t=0; t< T; t++) {
					if(states[observation[t]] == j) {
						numer = numer + gamma[t][i];
					}
				}
				B[i][j] = numer/ denom;
			}
		}
	}


	double compute_log_prob() {
		double log_prob =0;
		for(int i=0; i <T; ++i) {
			log_prob += log(c[i]);
		}
		log_prob *=-1;

		return log_prob;
	}

	void converge() {
		double oldLogProb = std::numeric_limits<double>::min();

		int iters =0;
		double log_prob = oldLogProb;
		do{
			oldLogProb = log_prob;
			forward_pass();
			backward_pass();
			compute_gamma();
			re_estimate();
			log_prob = compute_log_prob();
			iters++;
			cout.precision(17);
			cout << iters << " " << "log_prob: " << fixed << log_prob << "\n";
		}while( (iters < MIN_ITERS || log_prob > oldLogProb + EPISON ) && iters < MAX_ITERS);
	}

	void print() {


		cout << "Pi is \n";
		for(int i =0; i< N; i++) {
			cout << pi[i] <<" ";
		}
		cout << "\n";

		cout << "Matrix A  \n"; 
		print_matrix(A, N, N);

		cout << "Matrix B  \n"; 
		print_matrix(B, N, M);
		
		/*cout << "Matrix alpha  \n"; 
		print_matrix(alpha, T, N);
		
		cout << "Matrix beta  \n"; 
		print_matrix(beta, T, N);

		cout << "C is \n";
		for(int i =0; i< T; i++) {
			cout << c[i] <<" ";
		}
		cout << "\n";

		cout <<"Matrix gamma \n";
		print_matrix(gamma, T, N);*/
	}
};



int main(int argc, char *argv[]) {
	srand(1);

	string content = get_file_contents(argv[1]);
	HmmState hmm = HmmState();

	hmm.init(2, content);
	hmm.converge();
	hmm.print();
}