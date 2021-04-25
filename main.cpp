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
#include"options.hpp"
#include<cstdio>

using namespace std;
#define PI 3.141592 

int main(){
    srand (1999);
    //const float sigma = 0.15;
   
    float T=2;
    float K =1;
    float S_0=1;
    ofstream fichier("Q_14.txt", ios::out| ios::trunc);
    if(fichier){
        float sigma = 0.15;
        gaussian_variable W(0.015,sigma);
        initial_conditions cond(2,1./52,104,1,0.7,S_0);
        for( int N=100; N<=1E7 ; N *=2){
            monte_carlo_output output = control_monte_carlo_down_in( cond, W, N);
            fichier<<N<<'\t'<<output.mean_estim<<'\t'<<output.inf_born<<'\t'<<output.sup_born<<'\t'<<output. err <<endl;
            //cout<<N<<endl;
        }
        
       /*
        for (int i=0; i<= 1000; i++){
            float t = (float) i/1000;
            float B = (1-t)*0.5+ t;
            float sigma =0.15;
            gaussian_variable W(0.015,sigma);
            initial_conditions cond(2,1./52,104,1,B,S_0);
            monte_carlo_output output = monte_carlo_option_down_put( cond, W, 10000);
            fichier<<B<<'\t'<<output.mean_estim<<'\t'<<output.inf_born<<'\t'<<output.sup_born<<'\t'<<output. err <<endl;
            cout<<i<<endl;

        } 
        //fichier<<vanilla_eur_option( W, T, K, S_0)<<endl;
        */

    }

    fichier.close();
}