#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm> 
#include <random>
#include <time.h>

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

string city_name(int i);
void calc_fitness(Population& popul);
void print_population(const Population& popul);
int find_best(const Population& popul);
int find_worst(const Population& popul);
Population newGeneration(const Population& old_popul, const Distances& dist);
Chromosome create_chromo(int i, const Distances& dist);
Chromosome select_one(Population& popul);
Chromosome crossover(const Chromosome& chromo_a, const Chromosome& chromo_b, const Distances& dist);
Chromosome mutate(const Chromosome& chromo, const Distances& dist);


double const MUTATE = 0.3;
double const CROSSOVER = 0.6;
double const SHUFFLE = 0.1;


void main()
{
	srand(time(NULL));

	unsigned int city_num, const generation_num = 1000;
	double best_dist = INFINITY;
	size_t i, j;
	Distances dist;
	Population popul;
	vector<string> cities;

	// file structure:
	// line 1: number of cities
	// line 2: names of cities
	// line 3: distance matrix
	ifstream fin("Distance_Matrix.txt");
	if (!fin.is_open()) return;
	else
	{
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
	}

	for (i = 0; i < city_num; ++i) popul[i] = create_chromo(i, dist);
	calc_fitness(popul);

	cout.width(13);
	cout << ' ';
	for (i = 0; i < city_num; ++i)
	{
		cout.width(13);
		cout << left << cities[i];
	}
	cout << endl;
	for (i = 0; i < city_num; ++i)
	{
		cout.width(13);
		cout << left << cities[i];
		for (j = 0; j < city_num; ++j)
		{
			cout.width(13);
			cout << dist[i][j];
		}

		cout << endl;
	}

	for (int k = 0; k < generation_num; ++k)
	{
		if (popul[find_best(popul)].totalDistance < best_dist)
		{
			best_dist = popul[find_best(popul)].totalDistance;
			cout << "\nCurrent best distance: " << best_dist << " km" << endl;
		}
		Population temp_popul = newGeneration(popul, dist);
		popul = move(temp_popul);
		calc_fitness(popul);
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

	chromo.push_back(i_city);

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

	chromo.totalDistance += dist[i_city][city];
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

	for (i = 0; i < popul.size(); ++i)
	{
		fitness_rank.push_back((1.0) / popul[i].totalDistance);
		max_probability += fitness_rank[i];
	}

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
	int i, j,
		size = chromo.size(),
		shift = rand() % (size - 1);

	Chromosome newChromo;

	for (i = 0; i < size; ++i) newChromo.push_back(chromo.at(i));
	swap(newChromo.at(0), newChromo.at(shift));
	newChromo.at(size - 1) = newChromo.at(0);

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
	Chromosome newChromo;
	int i, i_1, i_2, size = chromo_a.size();
	for (i = 0; i < chromo_a.size(); ++i)
	{
		i_1 = chromo_a[i];
		i_2 = chromo_a[i_1];
		newChromo.push_back(chromo_a[i_2]);
	}
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

	newGener.push_back(old_popul[find_best(old_popul)]);

	for (size_t i = 1; i < old_popul.size(); ++i)
	{
		Chromosome first_selected = select_one(old_popul),
				   second_selected = select_one(old_popul);

		double choice = (double)((rand() % 100)) / 100;

		if (first_selected == second_selected &&
			SHUFFLE - choice >= 0)
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