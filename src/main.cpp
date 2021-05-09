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
    ofstream fichier("Q_14_2.txt", ios::out| ios::trunc);
    if(fichier){
        /*
        float sigma = 0.15;
        gaussian_variable W(0.015,sigma);
        initial_conditions cond(2,1./52,52,1,0.7,S_0);
        //#pragma omp parallel for
        for( int N=100; N<=1E6 ; N*=2){
            monte_carlo_output output = monte_carlo_option_down_put( cond, W, N);
            fichier<<N<<'\t'<<output.mean_estim<<'\t'<<output.err<<endl;
            //cout<<N<<endl;
        };
        */
        //cout<<vanilla_eur_option( W, 2,1,1)<<endl;
        
        float sigma = 0.15;
        gaussian_variable W(0.015,sigma);
        initial_conditions cond(2,1./52,52,1,0.7,S_0);
        //#pragma omp parallel for
        for( int N=100; N<=1E6 ; N*=2){
            monte_carlo_output output = monte_carlo_option_down_put( cond, W, N);
            monte_carlo_output output2 = control_monte_carlo_down_in( cond, W, N);
            fichier<<N<<'\t'<<output.mean_estim<<'\t'<<output2.mean_estim<<'\t'<<output.err<<'\t'<<output2. err <<endl;
            //cout<<N<<endl;
        };
        
        
        /*
        float sigma =0.15;
        gaussian_variable W(0.015,sigma);
        //#pragma omp parallel for
        for(int i=0; i<=200; i++){
            float t =(float) i/200;
            float sigma = t*0.8;
            gaussian_variable W(0.015,sigma);
            initial_conditions cond(2,1./52,52,1,0.7,S_0);
            monte_carlo_output output = anti_monte_carlo_option_down_put( cond, W, 30000);
            fichier<< sigma <<'\t'<<output.mean_estim<<'\t'<<output.err<<endl;
        }
        */
        /*
        list<int> inv_delta= list<int> (6);
        list<int>::iterator it= inv_delta.begin();
        *it = 250;
        it++;
        *it = 52;
        it++;
        *it = 12;
        it++;
        *it = 4;
        it++;
        *it = 1;
        it++;
        *it = 3;
        list<int>::iterator iter= inv_delta.begin();
        for(; iter != inv_delta.end(); iter++){
            float delta = 1./(*iter) ;
            int N_delta = *iter;
            float sigma =0.15;
            gaussian_variable W(0.015,sigma);
            initial_conditions cond(2,delta,N_delta,1,0.7,S_0);
            monte_carlo_output output = anti_monte_carlo_option_down_put( cond, W, 200000);
            monte_carlo_output output2 = monte_carlo_option_down_out_proba_non_sortie( cond, W, 200000);
            fichier<<delta<<'\t'<<output.mean_estim<<'\t'<<output2.mean_estim<<'\t'<<N_delta<<'\t'<<output.err<<'\t'<<output2.err<<endl;
            cout<<N_delta<<endl;
        } 
        */
        //fichier<<vanilla_eur_option( W, T, K, S_0)<<endl;
        

    }

    fichier.close();
}