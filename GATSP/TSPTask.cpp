#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm> 
#include <random>
#include <time.h>
#include <assert.h>

using namespace std;

struct Chromosome : public vector<int>
{
	double fitness = 0;
	double totalDistance = 0;

	bool operator<(const Chromosome& chromo) const
	{
		if (totalDistance == chromo.totalDistance) return (fitness <= chromo.fitness) ? 1 : 0;
		return (totalDistance < chromo.totalDistance) ? 1 : 0;
	}
	bool operator==(const Chromosome& chromo) const
	{
		return (totalDistance == chromo.totalDistance) ? 1 : 0;
	}

	Chromosome() = default;
	Chromosome(const Chromosome&) = default;
	Chromosome(Chromosome&&) = default;
	Chromosome& operator=(const Chromosome&) = default;
	Chromosome& operator=(Chromosome&&) = default;
	Chromosome(int size) : vector<int>(size) {}
	~Chromosome() = default;
};

typedef vector<vector<double>> Distances;
typedef vector<Chromosome> Population;

// help functions
string city_name(int i);
void print_population(const Population& popul);
void print_distance_matrix(const vector<string>& cities, const Distances& dist);
void initialize_all(const string& filename,
					unsigned int& city_num,
					unsigned int& generation_num,
					Distances& dist,
					Population& popul,
					vector<string>& cities
					);

// assist function
int find_best(const Population& popul);
int find_worst(const Population& popul);

// algorithm functions
void calc_fitness(Population& popul);
Population newGeneration(const Population& old_popul, const Distances& dist);
Chromosome create_chromo(int i, const Distances& dist);
Chromosome select_one(Population& popul);
Chromosome crossover(const Chromosome& chromo_a, const Chromosome& chromo_b, const Distances& dist);
Chromosome shuffle_crossover(const Chromosome& chromo_a, const Distances& dist);
Chromosome mutate(const Chromosome& chromo, const Distances& dist);


// probabilities
double const MUTATE = 0.3;
double const CROSSOVER = 0.6;
double const SHUFFLE = 0.1;


void main()
{
	srand(time(NULL));

	unsigned int city_num, generation_num;
	double best_dist = INFINITY;
	
	Distances dist;
	Population popul;
	vector<string> cities;

	// file structure:
	// line 1: number of cities
	// line 2: names of cities
	// line 3: distance matrix
	initialize_all("Distance_Matrix.txt", city_num, generation_num, dist, popul, cities);

	print_distance_matrix(cities, dist);

	// executing genetic algorithm
	for (int k = 0; k < generation_num; ++k)
	{
		// tracking the best path
		if (popul[find_best(popul)].totalDistance < best_dist)
		{
			best_dist = popul[find_best(popul)].totalDistance;
			cout << "\nCurrent best distance: " << best_dist << " km" << endl;
		}

		// creating new generation based on the previous
		Population temp_popul = newGeneration(popul, dist);
		popul = move(temp_popul);

		// filling with fitness value ( better fitness - better chance to be chosen )

		calc_fitness(popul);

		// optional output of all chromosomes in generation (slowing down the execution)
		/*print_population(popul);*/
	}

	cout << "\nBest distance: " << best_dist << " km" << endl;
	cout << "Path: ";
	for (size_t j = 0; j <  city_num + 1; ++j)
	{
		cout << city_name(popul[find_best(popul)].at(j));
		if (j != city_num) cout << "->";
	}
	cout << endl;
	system("pause");
}

void initialize_all(const string& filename,
	unsigned int& city_num,
	unsigned int& generation_num,
	Distances& dist, Population& popul, vector<string>& cities)
{
	cout << "Input the number of generations: ";
	cin >> generation_num;
	cout << endl;

	size_t i, j;

	ifstream fin(filename);
	assert(fin.is_open());

	fin >> city_num;

	cities.resize(city_num);
	dist.resize(city_num);
	popul.resize(city_num);

	for (i = 0; i < city_num; ++i) fin >> cities[i];

	for (i = 0; i < city_num; ++i)
	{
		dist[i].resize(city_num);
		for (j = 0; j < city_num; ++j)
		{
			fin >> dist[i][j];
			if (dist[i][j] == 0) dist[i][j] = INFINITY;
		}
	}

	for (i = 0; i < city_num; ++i) popul[i] = create_chromo(i, dist);
	calc_fitness(popul);

	fin.close();
}

void print_distance_matrix(const vector<string>& cities, const Distances& dist)
{
	size_t i, j , size = cities.size();

	cout.width(13);
	cout << ' ';
	for (i = 0; i < size; ++i)
	{
		cout.width(13);
		cout << left << cities[i];
	}
	cout << endl;
	for (i = 0; i < size; ++i)
	{
		cout.width(13);
		cout << left << cities[i];
		for (j = 0; j < size; ++j)
		{
			cout.width(13);
			cout << dist[i][j];
		}

		cout << endl;
	}
}

void print_population(const Population& popul)
{
	int size = popul.size();

	cout << "\n----------------------------------------------------------" << endl;

	for (size_t i = 0; i < size; ++i)
	{
		cout << "Path " << i + 1 << ": ";
		for (size_t j = 0; j < size + 1; ++j)
		{
			cout << city_name(popul[i].at(j));
			if (j != size) cout << "->";
		}
		cout << "\nTotal distance: " << popul[i].totalDistance << " km" << endl;
	}
}

Chromosome create_chromo(int city, const Distances& dist)
{
	int size = dist.size(), i_city = city, min = 0;
	Chromosome chromo;
	
	// pushing start location to the path
	chromo.push_back(i_city);

	// initialize first generation using greedy algorthim
	// choosing locally optimal ( shortest ) path
	for (size_t i = 1; i < size; ++i)
	{
		min = i_city;

		for (size_t j = 0; j < size; ++j)
		{
			if (dist[i_city][min] > dist[i_city][j] && find(chromo.begin(), chromo.end(), j) == chromo.end())
			{
				min = j;
			}
		}

		chromo.totalDistance += dist[i_city][min];
		i_city = min;
		chromo.push_back(i_city);
	}

	//calculating the totall distance for comparing pathes
	chromo.totalDistance += dist[i_city][city];

	// closing the circle with the start location
	chromo.push_back(city);

	return move(chromo);
}

string city_name(int i)
{
	string city;
	switch (i)
	{
	case 0: city = "Lviv";         break;
	case 1: city = "Donetsk";      break;
	case 2: city = "Kiev";         break;
	case 3: city = "Odessa";       break;
	case 4: city = "Dnepr";        break;
	case 5: city = "Kharkiv";      break;
	case 6: city = "Zaporizhzhia"; break;
	case 7: city = "Krivoy-Rog";   break;
	case 8: city = "Ternopil";     break;
	case 9: city = "Vinnitsa";     break;
	default: city = "Unknown";	   break;
	}
	return city;
}

void calc_fitness(Population& popul)
{
	vector<double> fitness_rank;
	double max_probability = 0;
	size_t i;
	int fit_best = find_best(popul);

	// summing up the reverse value for distances,
	// so that the shorter is path - the bigger is value
	for (i = 0; i < popul.size(); ++i)
	{
		fitness_rank.push_back((1.0) / popul[i].totalDistance);
		max_probability += fitness_rank[i];
	}

	// getting the percents of each value for every path (shorter paths have bigger probability percent)
	for (i = 0; i < popul.size(); ++i)
	{
		popul[i].fitness = fitness_rank[i] / max_probability;
	}
}

int find_best(const Population& popul)
{
	int best = 0;
	for (size_t i = 1; i < popul.size(); ++i)
	{
		if ((popul[i]) < (popul[best])) best = i;
	}

	return best;
}

int find_worst(const Population& popul)
{
	int worst = 0;
	for (size_t i = 1; i < popul.size(); ++i)
	{
		if ((popul[worst]) >(popul[i]))
			worst = i;
	}

	return worst;
}

Chromosome select_one(const Population& popul)
{
	double choice = static_cast<double>((rand() % 10000)) / 10000.0;

	int index = 0;

	// every path has its specific percent area
	// so we are moving area by area until we reach
	while (choice > 0)
	{
		choice -= popul[index].fitness;
		index++;
	}

	if(--index < 0) index = 0;

	return popul[index];
}

Chromosome crossover(const Chromosome& chromo_a, const Chromosome& chromo_b, const Distances& dist)
{
	// this crossover is swaping randomly chosen frist part of the chromo_a
	// and filling the tail by values from chormo_b that doesn't to be appeared in chromo_a part
	Chromosome newChromo;
	int i,
		size = chromo_a.size(),
		part = rand() % size;

	for (i = 0; i < part; ++i)
	{
		newChromo.push_back(chromo_a.at(i));
	}
	for (i = 0; i < size; ++i)
	{
		if (find(newChromo.begin(), newChromo.end(), chromo_b.at(i)) == newChromo.end())
			newChromo.push_back(chromo_b.at(i));
	}
	newChromo.push_back(newChromo.at(0));
	
	// calculating distance for the path
	for (i = 0; i < size - 2; ++i)
	{
		int i_1 = newChromo.at(i),
			i_2 = newChromo.at(i + 1);
		newChromo.totalDistance += dist[i_1][i_2];
	}

	return move(newChromo);
}

Chromosome mutate(const Chromosome& chromo, const Distances& dist)
{
	// in this mutate function we are swaping the starting point with the randomly chosen value
	int i, j,
		size = chromo.size(),
		shift = rand() % (size - 1);

	Chromosome newChromo;

	for (i = 0; i < size; ++i) newChromo.push_back(chromo.at(i));
	swap(newChromo.at(0), newChromo.at(shift));
	newChromo.at(size - 1) = newChromo.at(0);

	// calculating the distance
	for (i = 0; i < size - 2; ++i)
	{
		int i_1 = newChromo.at(i),
			i_2 = newChromo.at(i + 1);
		newChromo.totalDistance += dist[i_1][i_2];
	}

	return move(newChromo);
}

Chromosome shuffle_crossover(const Chromosome& chromo_a, const Distances& dist)
{
	// using value from chomo_a as a indexe in chomo_b
	// to get the value from chomo_b that will be a index of value to swap
	// with the starting value
	Chromosome newChromo;
	int i, i_1, i_2, size = chromo_a.size();
	for (i = 0; i < chromo_a.size(); ++i)
	{
		i_1 = chromo_a[i];
		i_2 = chromo_a[i_1];
		newChromo.push_back(chromo_a[i_2]);
	}

	// calculating distances
	for (i = 0; i < size - 2; ++i)
	{
		i_1 = newChromo.at(i);
		i_2 = newChromo.at(i + 1);
		newChromo.totalDistance += dist[i_1][i_2];
	}

	return move(newChromo);
}

Population newGeneration(const Population& old_popul, const Distances& dist)
{
	Population newGener;

	// always adding the best one
	newGener.push_back(old_popul[find_best(old_popul)]);

	// fill the new population based on the randomly chosen chromoes and actions
	for (size_t i = 1; i < old_popul.size(); ++i)
	{
		Chromosome first_selected = select_one(old_popul),
				   second_selected = select_one(old_popul);

		double choice = (double)((rand() % 100)) / 100;

		if (first_selected == second_selected && SHUFFLE - choice >= 0)
		{
			first_selected = shuffle_crossover(first_selected, dist);
		}
		else if (MUTATE - choice >= 0)
		{
			first_selected = mutate(first_selected, dist);
		}
		else if (CROSSOVER - choice >= 0)
		{
			first_selected = crossover(first_selected, second_selected, dist);
		}

		newGener.push_back(first_selected);
	}

	return move(newGener);
}