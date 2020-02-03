#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <algorithm>
#include <functional>
#include <set>
using namespace std;

struct gene {
    double x_1, x_2;
};

struct chromosome {
    gene person;
    double fitness_value;
    //double selection_chance;

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

//void generate_selection_chance(vector<chromosome> & population_input, const double s) {
//    int total_population = population_input.size();
//    for (int rank = 0; rank < population_input.size(); ++rank) {
//        population_input[rank].selection_chance = (2. - s)/(total_population) + 2. * rank * (s - 1.)/(total_population * (total_population - 1.));
//    }
//};

//int parent_selection_rank_based(vector<chromosome> & population_input, function<double(void)> rand01) {
//    double cumulative_probability = 0;
//    vector<double> cumulative_probabilities(population_input.size());
//    for(int i = 0; i < population_input.size(); ++i) {
//        cumulative_probability += population_input[i].selection_chance;
//        cumulative_probabilities[i] = cumulative_probability;
//    }
//    double r = rand01();        //r is a random number from 0 to 1
//    int position = 0;
//    while(position < population_input.size() && r > cumulative_probabilities[position]) {
//        //while you are not exceeding population size and you check the first cumulative probability that is more than r
//        position++;
//    }
//    cout << "position is " << position << endl;
//    return position;    //return the first position that's greater than r
//}

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

//
//void crossover_and_mutation_rank_based(generation & gen, int this_generation, const double crossover_fraction, const double mutation_chance, function<double(void)> rand_x1, function<double(void)> rand_x2, function<double(void)> rand01) {
//    int number_of_crossovers = round(gen.population[this_generation].size() * crossover_fraction);
//    for(int i = 0; i < number_of_crossovers; ++i) {
//        int position1 = parent_selection_rank_based(gen.population[this_generation], rand01);
//        chromosome X1 = gen.population[this_generation][position1];
//        int position2 = parent_selection_rank_based(gen.population[this_generation], rand01);
//        chromosome X2 = gen.population[this_generation][position2];
//        chromosome X_new;
//        double random_number_x1 = rand_x1();
//
//        //go to lecture 2, slide 25
//        double crossover_x1_rand = rand01();
//        X_new.person.x_1 = X1.person.x_1 * crossover_x1_rand + (1 - crossover_x1_rand) * X2.person.x_1;
//        //performing intermediate crossover between the x_1 of both parents
//        double crossover_x2_rand = rand01();
//        X_new.person.x_2 = X1.person.x_2 * crossover_x2_rand + (1 - crossover_x2_rand) * X2.person.x_2;
//        //performing intermediate crossover between the x_2 of both parents
//        double check_mutation = rand01();
//        //go to lecture 2, slide 18
//        if(check_mutation < mutation_chance) {
//            double mutate1 = rand01();
//            double mutate2 = rand01();
//            X_new.person.x_1 + (mutate1 - mutate2) * 12;
//            if(X_new.person.x_1 > 12) {
//                X_new.person.x_1 = X_new.person.x_1 - 12;
//            }
//            else if(X_new.person.x_1 < -12){
//                X_new.person.x_1 = X_new.person.x_1 + 12;
//            }
//            X_new.person.x_2 + (mutate1 - mutate2) * 6;
//            if(X_new.person.x_2 > 6) {
//                X_new.person.x_2 = X_new.person.x_2 - 6;
//            }
//            else if(X_new.person.x_2 < -6){
//                X_new.person.x_2 = X_new.person.x_2 + 6;
//            }
//        }
//        //putting them in the next generation which is gen.population[current generation + 1][individual within that generation]
//        calc_fitness(X_new);
//        //populating the new generation with the crossovers
//        gen.population[this_generation + 1][i] = X_new;
//    }
//    // Put the best genes from previous generation into the new generation.
//    //The for loop decrements from lowest up until the number of crossovers
//    for(int i = gen.population[this_generation].size() - 1; i >= number_of_crossovers; --i) {
//        gen.population[this_generation + 1][i] = gen.population[this_generation][i];
//    }
//}

void crossover_and_mutation_tournament(generation & gen, int this_generation, const double crossover_fraction, const double mutation_chance, const int k, function<double(void)> rand_x1, function<double(void)> rand_x2, function<double(void)> rand01) {
    int number_of_crossovers = round(gen.population[this_generation].size() * crossover_fraction);
    for(int i = 0; i < number_of_crossovers; ++i) {
        int position1 = parent_selection_tournament(gen.population[this_generation], k, rand01);
        chromosome X1 = gen.population[this_generation][position1];
        int position2 = parent_selection_tournament(gen.population[this_generation], k, rand01);
        chromosome X2 = gen.population[this_generation][position2];
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

        if(check_mutation < mutation_chance) {
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

int main() {
    std::mt19937_64 rng; // random generator
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // rng.seed(666); // use black magic to speed up code
    uniform_real_distribution<double> distr01(0.,1.);   // generate a random number between 0 and 1
    uniform_real_distribution<double> distr_x1(-12.,12.);   //for the rand function, creating constraints (x1)
    uniform_real_distribution<double> distr_x2(-6.,6.);     //for the rand function, creating constraints (x2)



    int nb_of_generations = 1000;
    int total_population = 1000;



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

    //rank_population is essentially sorting the vector from 0 to 1000 (lowest fitness to highest fitness
    rank_population(gen.population[0]); //ranking population which is the vector within the generation 0
    // generate_selection_chance(gen.population[0],2);     //you calculate the chance based off of ranking
    crossover_and_mutation_tournament(gen, 0, 0.8, 0.1, 4,
            [&](){return distr_x1(rng);}, [&](){return distr_x2(rng);}, [&](){return distr01(rng);});
//    for(int i = 0; i < gen.population[1].size(); ++i){
//        cout << "x1 = " << gen.population[1][i].person.x_1 << " x2 = " << gen.population[1][i].person.x_2 << " fitness = " << gen.population[1][i].fitness_value << endl;
//    }


    //ranked population 1
    rank_population(gen.population[1]);

    //took the best fit in generation 1 which was ranked previously
    double best_fitness = gen.population[1][total_population - 1].fitness_value;


    int count_best_fitness = 0;
    for(int i = 1; i < nb_of_generations - 1; i++) {
        // rank_population(gen.population[i]); //ranking population which is the vector within the generation 0
        // generate_selection_chance(gen.population[i],2);     //you calculate the chance based off of ranking
        crossover_and_mutation_tournament(gen, i, 0.8, 0.1, 4,
                               [&](){return distr_x1(rng);}, [&](){return distr_x2(rng);}, [&](){return distr01(rng);});
        rank_population(gen.population[i + 1]);
        double fitness = gen.population[i + 1][total_population - 1].fitness_value;

        //checking if the best one changed since the previous generation
        //checking to see if it didn't changed by subtracting them and seeing if it is 0 (or in this case below a threshold)

        if(abs(best_fitness - fitness) < 1e-6) { // best_fitness == fitness
            count_best_fitness++;
        } else {
            best_fitness = fitness;
            count_best_fitness = 0;
        }

        //check to see if it changed after 100 times, if it didn't then stop
        if(count_best_fitness == 100) {
            //i, at this point, is the generation where it didn't change anymore (we are still inside the main for loop)
            nb_of_generations = i;
            cout << "found best after " << i << " generations" << endl;

            //break here because you found it

            break;
        }
    }
    cout << "x1 = " << gen.population[nb_of_generations - 1][total_population - 1].person.x_1 << " x2 = " << gen.population[nb_of_generations - 1][total_population - 1].person.x_2 << " fitness = " << gen.population[nb_of_generations - 1][total_population - 1].fitness_value << endl;

    return 0;
}
