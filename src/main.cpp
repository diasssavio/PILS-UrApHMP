//============================================================================
// Name        : main.cpp
// Author      : Sávio S. Dias
// Version     : 1.0
// Copyright   :
// Description : A GRASP/ILS-based heuristic for the UrApHMP
//============================================================================

#include <iostream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <ctime>

#include "../include/FWChrono.h"
#include "../include/mt19937ar.h"
#include "../include/grasp.h"
#include "../include/ils.h"
#include "../include/solution.h"
#include "../include/UrApHMP.h"

using namespace std;

template<typename T>
T string_to(const string& s){
	istringstream i(s);
	T x;
	if (!(i >> x)) return 0;
	return x;
}

int main(int argc, char* args[]) {
	FWChrono timer;
	timer.start();
	// time_t seed = time(NULL);
	unsigned long seed = string_to<unsigned long>(args[3]);
	// srand(seed);
	init_genrand(seed);

	int n;

	scanf("%d", &n);
	double X = 3.0, alpha_1 = 0.75, delta = 2.0;
	uraphmp instance(n, X, alpha_1, delta);
	int p = string_to<int>(args[1]);
	int r = string_to<int>(args[2]);

	// Reading file information
	vector< vector< double> > aux(n);
	for(int i = 0; i < n; i++){
		vector< double > aux2(n);
		for(int j = 0; j < n; j++)
			scanf("%lf", &aux2[j]);
		aux[i] = aux2;
	}
	instance.set_traffics(aux);

	for(int i = 0; i < n; i++){
		vector< double > aux2(n);
		for(int j = 0; j < n; j++)
			scanf("%lf", &aux2[j]);
		aux[i] = aux2;
	}
	instance.set_distances(aux);
	aux.clear();


	int max_iterations = 0.1 * n;
	int max_r = max_iterations / 2;
	double alpha_2 = 0.2;

	// grasp ILS(instance, max_iterations, p, r, alpha_2, timer);
	// solution result = ILS.execute();

	// ils ILS(instance, max_iterations, p, r, timer);
	ils ILS(instance, max_iterations, max_r, alpha_2, p, r, timer);
	solution result = ILS.execute();

	timer.stop();
	// printf("TOTAL EXECUTION TIME: %.2lf", timer.getStopTime());
	printf("%.2lf,", timer.getStopTime());
	// result.show_data();

	// printf("\nIT_LOG:\n");
	vector< pair< double, unsigned> > it_log = ILS.get_it_log();
	vector< double > times = ILS.get_times();
	double min_time = 0.0;
	for(unsigned i = 0; i < it_log.size(); i++){
		// printf("#%d:\t%.2lf\t%.2lf\n", it_log[i].second, it_log[i].first, times[i]);

		if(it_log[i].first == result.get_total_cost() && min_time == 0.0)
			min_time = times[i];
	}

	// vector< int > path = GRASP.get_path();
	// printf("\nPATH RELINKING USAGE:\n");
	// for(unsigned i = 0; i < path.size(); i++)
	// 	printf("%d\t", path[i]);

	printf("%.2lf,%.2lf\n", min_time, result.get_total_cost());
	// cout << "," << seed << endl;

	return 0;
}
