#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <algorithm>
#include <functional>
#include <set>
#include <fstream>
using namespace std;

struct gene {
    double x_1, x_2;
};

struct chromosome {
    gene person;
    double fitness_value;
    double selection_chance;

};

struct generation {

    generation(int iterations, int total_population){
        population = vector<vector<chromosome>>(iterations, vector<chromosome>(total_population));
    };
    vector<vector<chromosome>> population;
};

void calc_fitness(chromosome & x){
    x.fitness_value = 21.5 + x.person.x_1 * sin(4* M_PI *x.person.x_1) + x.person.x_2 * sin(20* M_PI *x.person.x_2);
}

void sort_population(vector<chromosome> & population_input) {
    sort(population_input.begin(), population_input.end(), [](chromosome &lhs, chromosome rhs) {
        return lhs.fitness_value < rhs.fitness_value;
    });
}

void generate_selection_chance(vector<chromosome> & population_input, const double s) {
    int total_population = population_input.size();
    for (int rank = 0; rank < population_input.size(); ++rank) {
        population_input[rank].selection_chance = (2. - s)/(total_population) + 2. * rank * (s - 1.)/(total_population * (total_population - 1.));
    }
};

//The function below is on lecture 3, slide 8

int parent_selection_rank_based(vector<chromosome> & population_input, function<double(void)> rand01) {
    double cumulative_probability = 0;
    vector<double> cumulative_probabilities(population_input.size());
    for(int i = 0; i < population_input.size(); ++i) {
        cumulative_probability += population_input[i].selection_chance;
        cumulative_probabilities[i] = cumulative_probability;
    }
    double r = rand01();        //r is a random number from 0 to 1
    int position = 0;
    while(position < population_input.size() && r > cumulative_probabilities[position]) {
        //while you are not exceeding population size and you check the first cumulative probability that is more than r
        position++;
    }
    //cout << "position is " << position << endl;
    return position;    //return the first position that's greater than r
}


// The function below is implementing tournament selection for the population
// It is located on Lecture 3, slide 11

int parent_selection_tournament(vector<chromosome> & population_input, const int k, function<double(void)> rand01) {
    int winner = -1;
    double best_fitness = -numeric_limits<double>::max_digits10; // set best_fitness minus infinity
    set<int> indices;
    int count = 0;
    // we want k number of unique indices
    while(count < k) {
        int random_index = round(rand01() * population_input.size());
        if(indices.find(random_index) == indices.end()) { // find whether random_index already exists in indices.
            indices.insert(random_index);
            ++count; // success, count up by + 1.
        } else {
            // random_index already exists in indices. Thus, don't insert it and retry another random_index.
            continue;
        }
    }

    //For loop below is intended to take the fittest out of the randomly selected k number of parents
    vector<int> indices_vec(indices.begin(), indices.end());
    for(int i = 0; i < k; ++i) {
        int index = indices_vec[i];
        double fitness = population_input[index].fitness_value;
        if(fitness > best_fitness) {
            best_fitness = fitness;
            winner = index;
        }
    }
    return winner;    //return the first position that's greater than r
}


void crossover_and_mutation_function(generation & gen, int this_generation, const double crossover_probability, const double mutation_probability, const int k, function<double(void)> rand01, int selector) {
    int position1,position2 = 0;
    int number_of_crossovers = round(gen.population[this_generation].size() * crossover_probability);
    for(int i = 0; i < number_of_crossovers; ++i) {

        //Figuring out index of parent 1
        if(selector == 0) {
            position1 = parent_selection_tournament(gen.population[this_generation], k, rand01);
        } else{
            generate_selection_chance(gen.population[this_generation],2);
            position1 = parent_selection_rank_based(gen.population[this_generation], rand01);
        }
        //Taking parent 1
        chromosome X1 = gen.population[this_generation][position1];

        //Figuring out index of parent 2
        if(selector == 0) {
            position2 = parent_selection_tournament(gen.population[this_generation], k, rand01);
        } else{
            generate_selection_chance(gen.population[this_generation],2);
            position2 = parent_selection_rank_based(gen.population[this_generation], rand01);
        }
        //Taking parent 2
        chromosome X2 = gen.population[this_generation][position2];

        //Creating a child
        chromosome X_new;

        //go to lecture 2, slide 25
        //For crossover

        double crossover_x1_rand = rand01();
        X_new.person.x_1 = X1.person.x_1 * crossover_x1_rand + (1 - crossover_x1_rand) * X2.person.x_1;
        //performing intermediate crossover between the x_1 of both parents
        double crossover_x2_rand = rand01();
        X_new.person.x_2 = X1.person.x_2 * crossover_x2_rand + (1 - crossover_x2_rand) * X2.person.x_2;
        //performing intermediate crossover between the x_2 of both parents




        double check_mutation = rand01();
        //go to lecture 2, slide 18
        //For mutation

        if(check_mutation < mutation_probability) {
            double mutate1 = rand01();
            double mutate2 = rand01();
            X_new.person.x_1 + (mutate1 - mutate2) * 12;
            if(X_new.person.x_1 > 12) {
                X_new.person.x_1 = X_new.person.x_1 - 12;
            }
            else if(X_new.person.x_1 < -12){
                X_new.person.x_1 = X_new.person.x_1 + 12;
            }
            X_new.person.x_2 + (mutate1 - mutate2) * 6;
            if(X_new.person.x_2 > 6) {
                X_new.person.x_2 = X_new.person.x_2 - 6;
            }
            else if(X_new.person.x_2 < -6){
                X_new.person.x_2 = X_new.person.x_2 + 6;
            }
        }
        //putting them in the next generation which is gen.population[current generation + 1][individual within that generation]
        calc_fitness(X_new);
        //populating the new generation with the crossovers
        gen.population[this_generation + 1][i] = X_new;
    }
    // Put the best genes from previous generation into the new generation.
    //The for loop decrements from lowest up until the number of crossovers
    for(int i = gen.population[this_generation].size() - 1; i >= number_of_crossovers; --i) {
        gen.population[this_generation + 1][i] = gen.population[this_generation][i];
    }
}


double calculate_mean(vector<chromosome> & population_input, int current_generation){
    double sum = 0;
    for(int index = 0; index < population_input.size(); ++index){
        sum += population_input[index].fitness_value;
    }
    double mean = sum/population_input.size();
    return mean;
}


double calculate_standard_deviation(vector<chromosome> & population_input, double mean, int current_generation){
    double standard_deviation = 0.0;
    double subtraction = 0;
    double squared = 0;
    double sum = 0.0;
    for(int index = 0; index < population_input.size(); ++index){
        subtraction = (population_input[index].fitness_value - mean);
        squared = pow(subtraction,2.0);
        standard_deviation += squared;
    }

    standard_deviation = sqrt(standard_deviation/population_input.size());
    return standard_deviation;
}


int main() {
    int nb_of_generations = 500;
    int total_population = 10000;
    int selector;
//    cout << "Insert population size: " << endl;
//    cin >> total_population;
//    cout << "Insert how many generations: " << endl;
//    cin >> nb_of_generations;
    cout << "Select which version to Implement: " << endl << "0: First version where tournament based selection was used " << endl
    << "1: Second version where rank based selection was used" << endl;
    cin >> selector;


    ofstream data_file_mean("mean.csv");
    ofstream data_file_standard_deviation("standardDeviation.csv");
    //data_file_mean.open("data.csv");
    std::mt19937_64 rng; // random generator
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    //Reference for the line of code above to generate random numbers
    // https://sort-care.github.io/C++-Random-Number/

    uniform_real_distribution<double> distr01(0.,1.);   // generate a random number between 0 and 1
    uniform_real_distribution<double> distr_x1(-12.,12.);   //for the rand function, creating constraints (x1)
    uniform_real_distribution<double> distr_x2(-6.,6.);     //for the rand function, creating constraints (x2)
    //Reference for the line above to create a random number between the constraints described in the assignment
    // http://www.cplusplus.com/reference/random/uniform_real_distribution/



    double mean = 0;
    double standard_deviation = 0;



    generation gen(nb_of_generations,total_population);

    //populating initial generation (generation 0)
    for(int i = 0; i < total_population; ++i) {
        chromosome c;                      //create a person and then in the next few lines, provide its x1,x2, and fitness
        c.person.x_1 = distr_x1(rng);     //obtain x1 for a single individual
        c.person.x_2 = distr_x2(rng);     //obtain x2 for a single individual
        calc_fitness(c);        //this is for one individual which then you place inside a vector in the next line
        gen.population[0][i] = c;   //put chromosome c in generation 0 at chromosome i
       // cout << "x1 = " << c.person.x_1 << " x2 = " << c.person.x_2 << " fitness = " << c.fitness_value << endl;
    }



    //sort_population is essentially sorting the vector from 0 to 1000 (lowest fitness to highest fitness
    sort_population(gen.population[0]); //ranking population which is the vector within the generation 0


    mean = calculate_mean(gen.population[0] ,0);

    data_file_mean << "Generation" << "," << "Mean" << "\n" << 0 << "," << mean << "\n";

    standard_deviation = calculate_standard_deviation(gen.population[0],mean,0);

    data_file_standard_deviation << "Generation" << "," << "Standard Deviation" << "\n" << 0 << "," << standard_deviation << "\n";

    // generate_selection_chance(gen.population[0],2);     //you calculate the chance based off of ranking
    crossover_and_mutation_function(gen, 0, 0.8, 0.1, 4, [&]() { return distr01(rng); },selector);

    //ranked population 1
    sort_population(gen.population[1]);

    mean = calculate_mean(gen.population[1], 1);
    data_file_mean << 1 << "," << mean << "\n";


    standard_deviation = calculate_standard_deviation(gen.population[1],mean,1);
    data_file_standard_deviation << 1 << "," << standard_deviation << "\n";


    //took the best fit in generation 1 which was ranked previously
    //Consider removing
    double best_fitness = gen.population[1][total_population - 1].fitness_value;


    //The block of code below considers generation 2 and onwards. Generation 1 was initialized in the beginning
    //and generation 1 was dealt with separately.
    //Now (on the block below) we will call crossover_and_mutation for all upcoming generations




    int count_best_fitness = 0;
    for(int i = 1; i < nb_of_generations - 1; i++) {

        crossover_and_mutation_function(gen, i, 0.8, 0.1, 4, [&]() { return distr01(rng); },selector);
        sort_population(gen.population[i + 1]);

        mean = calculate_mean(gen.population[i], i);
        data_file_mean << i << "," << mean << "\n";
        standard_deviation = calculate_standard_deviation(gen.population[i],mean,i);
        data_file_standard_deviation << i << "," << standard_deviation << "\n";

        //Consider removing
        double fitness = gen.population[i + 1][total_population - 1].fitness_value;

    }
    cout << "x1 = " << gen.population[nb_of_generations - 1][total_population - 1].person.x_1 << " x2 = " << gen.population[nb_of_generations - 1][total_population - 1].person.x_2 << " fitness = " << gen.population[nb_of_generations - 1][total_population - 1].fitness_value << endl;
   // data_file_mean.close();

    return 0;
}
