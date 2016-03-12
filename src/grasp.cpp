/*
 * grasp.cpp
 *
 *  Created on: Apr 30, 2015
 *      Author: Sávio S. Dias
 */

#include "../include/grasp.h"

grasp::grasp( uraphmp& instance, size_t max_iterations, int p, int r, double alpha, FWChrono& timer ) : max_iterations(max_iterations), p(p), r(r), alpha(alpha) {
	this->set_instance(instance);
	this->timer = timer;
}

grasp::~grasp() { }

void grasp::set_max_iterations( size_t max_iterations ){ this->max_iterations = max_iterations; }

void grasp::set_instance( uraphmp& instance ){ this->instance = instance; }

void grasp::set_best( const solution& sol ){ this->best = sol; }

size_t grasp::get_max_iterations(){ return this->max_iterations; }

uraphmp& grasp::get_instance(){ return this->instance; }

solution& grasp::get_best(){ return this->best; }

vector< pair< double, unsigned > >& grasp::get_it_log(){ return this->it_log; }

vector< double >& grasp::get_times(){ return this->times; }

vector< int >& grasp::get_path(){ return this->path; }

bool grasp::my_comparison( pair< double, int > p1, pair< double, int > p2 ){ return (p1.first < p2.first); }

solution grasp::greedy_randomized_construction(){
	vector< vector< double > > traffics = instance.get_traffics();
	vector< vector< double > > distances = instance.get_distances();

	solution sol(instance, p, r);

	// Step 1 - Locating the hubs - alloc_hubs
	set< unsigned > hubs;
	for(int p = 0; p < sol.get_p(); p++){
		// Calculating the g(h) for each node unselected as hub - The Candidate List
		vector< pair< double, unsigned > > g;
		for(unsigned h = 0; h < mesh.size(); h++){
			// Adaptive part of the procedure
			if(hubs.find(mesh[h]) != hubs.end()) continue;

			// Creating the e(j,h)
			vector< double > e;
			for(int j = 0; j < instance.get_n(); j++){
				// Adaptive part of the procedure
				if(hubs.find(j) != hubs.end()) continue;

				double aux = 0.0;
				for(int i = 0; i < instance.get_n(); i++)
					aux += traffics[j][i];
				aux *= distances[j][ mesh[h] ];
				e.push_back(aux);
			}

			// Sorting the elements of e(j,h)
			sort(e.begin(), e.end());

			// Choosing the k = |n/p| elements to compose g(h)
			int k = instance.get_n()/sol.get_p();
			double sum = 0.0;
			for(int j = 0; j < k; j++)
				sum += e[j];
			g.push_back(make_pair(sum, mesh[h]));
		}

		// Creating the RCL based on g(h) -- GRASP
		double g_min = (*min_element(g.begin(), g.end(), my_comparison)).first;
		double g_max = (*max_element(g.begin(), g.end(), my_comparison)).first;
		vector< unsigned > RCL;
		for(unsigned i = 0; i < g.size(); i++)
			if(g[i].first <= (g_min + this->alpha*(g_max - g_min)))
				RCL.push_back(g[i].second);

		// Selecting randomly the hub from RCL
		hubs.insert(RCL[ genrand_int32() % RCL.size() ]);

//		pair< double, int > g_min = *min_element(g.begin(), g.end(), my_comparison);
//		hubs.push_back(g_min.second);
	}
	sol.set_alloc_hubs(hubs);

	// Step 2 - Assignment of the r hubs - assigned_hubs
	sol.assign_hubs();

	// Step 3 - Route traffics
	sol.route_traffics();

	return sol;
}

solution grasp::local_search_n1(solution& p_sol){
	// Evaluating the neighbors
	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		solution improved = neighborhood1(partial);
		// Checking the best elements in the Neighborhood 1
		if(improved.get_total_cost() < p_sol.get_total_cost())
			return improved;
		else is_improved = false;
	}

	return improved;
}

solution grasp::local_search_rn1( solution& p_sol ){
	// Checking the best elements in the Neighborhood 1
	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		improved = r_neighborhood1(partial);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
		else is_improved = false;
	}

	return partial;
}

solution grasp::local_search_c2n1( solution& p_sol ){
	// 1.2.Test the use w/ RN1

	// Checking the best elements in the Adjacent points
	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		improved = closest2_n1(partial);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
		else is_improved = false;
	}

	return partial;
}

solution grasp::local_search_na( solution& p_sol ){
	// Checking the best solution in the Neighborhood A

	solution partial = p_sol;
	solution improved;
	bool is_improved = true;
	while(is_improved){
		improved = neighborhood_a(partial);
//		improved = neighbors[ neighbors.size() - 1 ];
		if(improved.get_total_cost() < partial.get_total_cost())
			partial = improved;
		else is_improved = false;
	}

	return improved;
}

solution grasp::neighborhood1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of the worst evaluated
	 * Hub in the p_sol (partial solution) given.
	 */

	vector< double > temp = p_sol.get_hubs_cost();
	vector< double >::iterator temp_it = max_element(temp.begin(), temp.end());
	unsigned h = temp_it - temp.begin(); // hub to be exchanged

	// Generating Neighborhood
	vector< solution > neighbors;
	for(unsigned i = 0; i < mesh.size(); i++){
		if(p_sol.is_hub(mesh[i])) continue;
		set< unsigned > hubs(p_sol.get_alloc_hubs());
		set< unsigned >::iterator it = hubs.begin();
		advance(it, h);
		hubs.erase(it);
		hubs.insert(mesh[i]);
//		hubs[h] = mesh[i];

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
//		s1.set_assigned_hubs(p_sol.get_assigned_hubs());
//		s1.assign_partial_hubs(h, i); // Changing assigned_hubs in one position
		s1.route_traffics();
		neighbors.push_back(s1);
	}

	return *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
}

solution grasp::r_neighborhood1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of a random Hub
	 * in p_sol (partial solution) given.
	 */

	// Generating Neighborhood
	vector< solution > neighbors;
	for(unsigned i = 0; i < mesh.size(); i++){
		if(p_sol.is_hub(mesh[i])) continue;
		int h = genrand_int32() % p; // hub to be exchanged
		set< unsigned > hubs(p_sol.get_alloc_hubs());
		set< unsigned >::iterator it = hubs.begin();
		advance(it, h);
		hubs.erase(it);
		hubs.insert(mesh[i]);
//		hubs[h] = mesh[i];

		solution s1(instance, p, r);
		s1.set_alloc_hubs(hubs);
		s1.assign_hubs();
		s1.route_traffics();
		neighbors.push_back(s1);
//		if(s1.get_cost() < p_sol.get_cost()) break;
	}

	return *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
}

solution grasp::closest2_n1( solution& p_sol ){
	/**
	 * Generate a Neighborhood based on a one-point exchange of the two closest points of Hubs
	 * in p_sol (partial solution) given.
	 */
	//  1.1.Check the generation of previously calculated solutions as neighbors
	//		There's many intersection of neighbors

	vector< vector< double > > distances = instance.get_distances();

	// Generating Neighborhood
	vector< solution > neighbors; // Analyzing the 2 nearest points for each hub
	set< unsigned > alloc_hubs = p_sol.get_alloc_hubs();
	set< unsigned >::iterator i;
	unsigned count = 0;
	for(i = alloc_hubs.begin(); i != alloc_hubs.end(); i++, count++){
		// Finding the first two non-hub nodes
		unsigned j = 0;
		while(p_sol.is_hub(mesh[j]))
			j++;
		unsigned k = j + 1;
		while(p_sol.is_hub(mesh[k]))
			k++;
		pair< unsigned, double > min1, min2;
		if(distances[*i][mesh[k]] < distances[*i][mesh[j]]){
			min1 = make_pair(mesh[k], distances[*i][mesh[k]]);
			min2 = make_pair(mesh[j], distances[*i][mesh[j]]);
		} else {
			min2 = make_pair(mesh[k], distances[*i][mesh[k]]);
			min1 = make_pair(mesh[j], distances[*i][mesh[j]]);
		}
		// Getting the two Euclidian closest non-hub nodes
		for(unsigned j = k + 1; j < mesh.size(); j++){
			if(*i == mesh[j] || p_sol.is_hub(mesh[j])) continue; // main diagonal & hub nodes avoidance
			if(distances[*i][mesh[j]] < min1.second){
				min2 = min1;
				min1 = make_pair(mesh[j], distances[*i][mesh[j]]);
			}
		}

		// Making the solution objects
		set< unsigned > hubs1(alloc_hubs);
		set< unsigned > hubs2(alloc_hubs);
		set< unsigned >::iterator it = hubs1.begin();
		advance(it, count);
		hubs1.erase(it);
		hubs1.insert(min1.first);
		it = hubs2.begin();
		advance(it, count);
		hubs2.erase(it);
		hubs2.insert(min2.first);

		solution s1(instance, p, r);
		solution s2(instance, p, r);
		s1.set_alloc_hubs(hubs1);
		s1.assign_hubs();
		s1.route_traffics();
		s2.set_alloc_hubs(hubs2);
		s2.assign_hubs();
		s2.route_traffics();

		neighbors.push_back(s1);
		neighbors.push_back(s2);
	}

	return *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
}

solution grasp::neighborhood_a( solution& p_sol ){
	// Generating Neighborhood
	//      3.Implement the partial evaluation of the solution for this neighborhood
	//		The issue here is that the solution quality will drop

	vector< solution > neighbors;
	set< unsigned > hubs(p_sol.get_alloc_hubs());
	for(int i = 0; i < instance.get_n(); i++){
		vector< unsigned > assigned_hubs(p_sol.get_assigned_hubs(i));

		// Finding the symmetric difference
		vector< int > to_assign;
		for(set< unsigned >::iterator j = hubs.begin(); j != hubs.end(); j++)
			if(find(assigned_hubs.begin(), assigned_hubs.end(), *j) == assigned_hubs.end())
				to_assign.push_back(*j);

		// Generating the neighbors
		for(int j = 0; j < this->r; j++)
			for(unsigned k = 0; k < to_assign.size(); k++){
				solution s1(instance, p, r);
				s1.set_alloc_hubs(hubs);
				s1.set_assigned_hubs(p_sol.get_assigned_hubs());
				s1.set_assigned_hub(i, j, to_assign[k]);
				s1.set_cost(p_sol.get_cost());
				s1.set_f_chosen(p_sol.get_f_chosen());
				s1.set_s_chosen(p_sol.get_s_chosen());
				s1.route_partial_traffics(i);
				neighbors.push_back(s1);
//				if(s1.get_cost() < p_sol.get_cost())
//					return neighbors;
			}
	}

	return *min_element(neighbors.begin(), neighbors.end(), solution::my_sol_comparison);
}

void grasp::preprocessing(){
	vector< vector< double > > traffics = instance.get_traffics();
	vector< vector< double > > distances = instance.get_distances();
	int n = instance.get_n();

	vector< pair< double, int > > aux;
	for(int h = 0; h < n; h++){
		// Creating the e(j,h)
		vector< double > e;
		for(int j = 0; j < n; j++){
			if(h == j) continue;

			double aux = 0.0;
			for(int i = 0; i < n; i++)
				aux += traffics[j][i];
			aux *= distances[j][h];
			e.push_back(aux);
		}

		// Sorting the elements of e(j,h)
		sort(e.begin(), e.end());

		// Chosing the k = |n/p| elements to compose g(h)
		int k = n/p;
		double sum = 0.0;
		for(int j = 0; j < k; j++)
			sum += e[j];
		aux.push_back(make_pair(sum, h));
	}

	// Eliminating the x percent most expensive
	vector< int > mesh;
	int x = (0.1 * n);
	sort(aux.begin(), aux.end(), my_comparison);
	for(unsigned i = 0; i < aux.size() - x; i++)
		mesh.push_back(aux[i].second);

	this->mesh = mesh;
}

solution grasp::path_relinking(solution& origin, solution& destination){
	vector< unsigned > hubs_dif;
	set< unsigned > hubs_o = origin.get_alloc_hubs();
	for(set< unsigned >::iterator i = hubs_o.begin(); i != hubs_o.end(); i++){
		if( destination.is_hub(*i) ) continue;
		hubs_dif.push_back(*i);
	}

	if(hubs_dif.size() > 1){
		int k = 0;
		set< unsigned > hubs_d(destination.get_alloc_hubs());
		for(set< unsigned >::iterator i = hubs_d.begin(); i != hubs_d.end(); i++){
			if( origin.is_hub(*i) ) continue;

			solution s1(instance, p, r);
			hubs_d.erase(i);
			hubs_d.insert(hubs_dif[k]);
//			hubs_d[i] = hubs_dif[k];
			k++;

			s1.set_alloc_hubs(hubs_d);
			s1.assign_hubs();
			s1.route_traffics();
			solution improved = local_search_rn1(s1);

			if( improved.get_cost() < destination.get_cost() ) return improved;
		}
	}

	return destination;
}

solution& grasp::execute(){
	// Pre-processing
	preprocessing();

	// processing
	vector< double > partial_improment;
	//	int t_max = 0.1 * max_iterations;
	//	double mean = 0.0;
	unsigned i = 1, k = 0;
	bool first = true;
	while(i < max_iterations){
		solution initial = greedy_randomized_construction();
//		solution improved = local_search_c2n1(initial);
		solution improved = local_search_rn1(initial);

		// Acceptance criterion
		if(!first){
			if(improved.get_total_cost() < best.get_total_cost()) {
				set_best(improved);
				i = 1;
			}else{
				improved = path_relinking(improved, best);
				if(improved.get_total_cost() < best.get_total_cost()){
					set_best(improved);
					i = 1;
					path.push_back(k);
				}else i++;
			}/*else i++;
			else{ // Testing the LS_a only when LS_h doesn't improve the solution
				improved = local_search_na(improved);
				if(improved.get_total_cost() < best.get_total_cost()){
					set_best(improved);
					i = 1;
				}else i++;
			}*/
		}else{
			best = improved;
			first = false;
		}

		it_log.push_back(make_pair(best.get_total_cost(), k++));
		times.push_back(((double) timer.getMilliSpan() / 1000));
//		printf("#%d:\t%.2lf\t%.2lf\n", k++, best.get_total_cost(), ((double) timer.getMilliSpan() / 1000));
	}


	// Post-processing
	//	set_best(local_search_rn1(best));

	return best;
}
