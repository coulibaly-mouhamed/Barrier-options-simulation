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

float gaussian_output(float T);
float abramowitz_stegun_int(float);
float vanilla_eur_option( gaussian_variable ,float ,float ,float );