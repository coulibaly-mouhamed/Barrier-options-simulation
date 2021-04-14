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
    const float sigma = 0.15;
    gaussian_variable W(0.015,0.15);
    float T=2;
    float K =1;
    float S_0=1;
    
    ofstream fichier("european_vanilla_options.txt", ios::out| ios::trunc);
    if(fichier){
        for( int N=30; N<=30000000000 ; N =N*2){
            float mean_estim= 0.;
            float std_estim = 0.;
            for(int i=1; i<=N;i++){
                float c =gaussian_output(T);
                float S_T = S_0*exp((W.mean-pow(sigma,2)/2)*T + sigma*c);
                if (K-S_T>0){mean_estim += (K-S_T);};
                std_estim += pow(fmax(0,K-S_T)*exp(-W.mean*T),2);
            }
            mean_estim = exp(-W.mean*T)*(mean_estim/(N));
            std_estim = (std_estim/N -pow(mean_estim,2)*exp(-2*W.mean*T))*(N/(N-1));
            std_estim = sqrt(std_estim);
            float inf_born = mean_estim - 1.645*std_estim/sqrt(N);
            float sup_born = mean_estim + 1.645*std_estim/sqrt(N);
            float err = 1.96*std_estim/sqrt(N);
            fichier<<N<<'\t'<<mean_estim<<'\t'<<inf_born<<'\t'<<sup_born<<'\t'<< err <<endl;
        }
        fichier<<vanilla_eur_option( W, T, K, S_0)<<endl;
    
    }

    fichier.close();
}