#include <iostream>
#include <list>
#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>
#include <chrono>
#include <ctime>
#include<random>
#include<algorithm>
#include<cstdio>
using namespace std;
#define PI 3.141592 

class random_variables{
    int length;
    list<float> val;
};
class gaussian_variable{
    public:
    float mean;
    float std;
    gaussian_variable(float m=0, float sigma=1){
        mean=m;
        std = sigma;
    }
};

class initial_conditions{
    public:
        int T;
        float delta;
        int N_delta;
        float K;
        float B;
        float S_0;
    initial_conditions(int t=0, float delt=0, int N_delt=0, float K_t=0, float B_t=0, float S_ =0){
        T=t; delta = delt; N_delta=N_delt; K=K_t; B=B_t; S_0=S_;
    }
};

class monte_carlo_output{
    public:
        float mean_estim;
        float std_estim;
        float inf_born;
        float sup_born;
        float err ;
    monte_carlo_output(float mean=0, float std=0, float inf=0, float sup=0, float er=0){
        mean_estim = mean; std_estim=std; inf_born = inf; sup_born =sup; err = er;
    }
};
float gaussian_output(float T);
float abramowitz_stegun_int(float);
float vanilla_eur_option( gaussian_variable ,float ,float ,float );
monte_carlo_output monte_carlo_option_down_put(initial_conditions ,gaussian_variable, int);
float mini(list<float> );
monte_carlo_output anti_monte_carlo_option_down_put(initial_conditions ,gaussian_variable , int );
float h(initial_conditions , float , gaussian_variable );
monte_carlo_output monte_carlo_vanilla_options(initial_conditions , gaussian_variable , int );
float h_anti(initial_conditions , float , gaussian_variable );
monte_carlo_output control_monte_carlo_down_in(initial_conditions, gaussian_variable, int);
