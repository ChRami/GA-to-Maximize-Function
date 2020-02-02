#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
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

void rank_population(vector<chromosome> & population_input) {
    sort(population_input.begin(), population_input.end(), [](chromosome &lhs, chromosome rhs) {
        return lhs.fitness_value < rhs.fitness_value;
    });
}

void generate_selection_chance(vector<chromosome> & population_input, const double s) {
    int total_population = population_input.size();
    for(int rank = 0; rank < population_input.size(); ++rank){
        population_input[rank].selection_chance = 2. * s/(total_population) + 2. * rank * (s - 1.)/(total_population * (total_population - 1.));
    }
};

int parent_selection(vector<chromosome> & population_input){
    double cumulative_probability = 0;
    vector<double> cumulative_probabilities(population_input.size());
    for(int i = 0; i < population_input.size(); ++i) {
        cumulative_probability += population_input[i].selection_chance;
        cumulative_probabilities[i] = cumulative_probability;
    }
    random_device rand_dev;
    mt19937 generator(rand_dev());      //random number generator (reference later)
    uniform_real_distribution<double> distr(0, 1);
    double r = distr(generator);        //r is a random number from 0 to 1
    int position = 0;
    while(position < population_input.size() && r > cumulative_probabilities[position]) {
        //while you are not exceeding population size and you check the first cumulative probability that is more than r
        position++;
    }
    return position;    //return the first position that's greater than r
}

void crossover_and_mutation(generation & population_input, int this_generation, const double crossover_fraction, const double mutation_chance) {
    random_device rand_dev;
    mt19937 generator(rand_dev());      //random number generator (reference later)
    uniform_real_distribution<double> distr(0., 1.);
    int number_of_crossovers = population_input.size() * crossover_fraction;
    for(int i = 0; i < number_of_crossovers; ++i) {
        int position1 = parent_selection(population_input[this_generation]);
        chromosome X1 = population_input[this_generation][position1];
        int position2 = parent_selection(population_input[this_generation]);
        chromosome X2 = population_input[this_generation][position2];
        chromosome X_new;
        double random_number = distr(generator);

        //go to lecture 2, slide 25

        X_new.person.x_1 = X1.person.x_1 * random_number + (1 - random_number) * X2.person.x_1;
        //performing intermediate crossover between the x_1 of both parents
        random_number = distr(generator);
        X_new.person.x_2 = X1.person.x_2 * random_number + (1 - random_number) * X2.person.x_2;
        //performing intermediate crossover between the x_2 of both parents
        double check_mutation = distr(generator);
        //go to lecture 2, slide 18
        if(check_mutation < mutation_chance) {
            double mutate1 = distr(generator);
            double mutate2 = distr(generator);
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
        //putting them in the next generation which is population_input[current generation + 1][individual within that generation]
        calc_fitness(X_new);
        population_input[this_generation + 1][i] = X_new;
    }
    int elitism = (1 - crossover_fraction) * population_input[this_generation].size();
    for(int i = population_input[this_generation].size(); i > number_of_crossovers; --i) {
        population_input[this_generation + 1][number_of_crossovers + 1 + i]
    }

}

int main() {
    random_device rand_dev;
    mt19937 generator(rand_dev());      //random number generator (reference later)
    uniform_real_distribution<double> distR_x1(-12.,12.);   //for the rand function, creating constraints (x1)
    uniform_real_distribution<double> distR_x2(-6.,6.);     //for the rand function, creating constraints (x2)
    generation gen(10,1000);
    for(int i = 0; i < 1000; ++i) {
        chromosome c;                      //create a person and then in the next few lines, provide its x1,x2, and fitness
        c.person.x_1 = distR_x1(generator);     //obtain x1 for a single individual
        c.person.x_2 = distR_x2(generator);     //obtain x2 for a single individual
        calc_fitness(c);        //this is for one individual which then you place inside a vector in the next line
        gen.population[0][i] = c;   //put chromosome c in generation 0 at chromosome i
       // cout << "x1 = " << c.person.x_1 << " x2 = " << c.person.x_2 << " fitness = " << c.fitness_value << endl;
    }
    //rank_population is essentially sorting the vector from 0 to 1000 (lowest fitness to highest fitness
    rank_population(gen.population[0]); //ranking population which is the vector within the generation 0
    generate_selection_chance(gen.population[0],2);     //you calculate the chance based off of ranking
    for(int i = 0; i < gen.population[0].size(); ++i){
        cout << "x1 = " << gen.population[0][i].person.x_1 << " x2 = " << gen.population[0][i].person.x_2 << " fitness = " << gen.population[0][i].fitness_value <<
        " Chance = "<< gen.population[0][i].selection_chance <<endl;
    }

    cout << "testing" << endl;

    return 0;
}
