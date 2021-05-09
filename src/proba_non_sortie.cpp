#include"options.hpp"
#include <iostream>
#include <list>
#include <cstdlib>
#include <string>
#include <fstream>
#include <math.h>
#include <chrono>
#include <ctime>
#include<algorithm>
#include<random>
#include<cstdio>
using namespace std;

#define PI 3.141592 

float ev_non_sortie_delta(initial_conditions cond,float X_T,gaussian_variable W){
    float r = W.mean;
    float sigma = W.std;
    int T = cond.T;
    float inter1 = cond.S_0*exp((r-pow(sigma,2)/2)*T+sigma*sqrt(T)*X_T); //Calcul S_T
    list<float> S_list = list<float> (cond.N_delta);
    list<float>::iterator iter= S_list.begin();
    int i =1;
    // First we create a  list of W_T_i output 
    for (iter= S_list.begin() ; iter != S_list.end(); iter++){
        if (iter==S_list.end()){
            *iter = inter1;
        }
        else{
            float X_i = gaussian_output(1);
            float t = i*cond.delta;
            float c = cond.S_0*exp((r-pow(sigma,2)/2)*t+sigma*sqrt(t)*X_i);
            *iter = c;
            i ++;
        }
    };
    float S_min = mini( S_list);
    if (S_min>= cond.B) {
        return 1.;
    }
    return 0.;

}


monte_carlo_output monte_carlo_proba_non_sortie_delta(initial_conditions cond,gaussian_variable W, int trajec){
    float mean_estim =0.;
    float std_estim = 0.;
    float inf_born=0.;
    float sup_born =0.;
    float err =0.;
    //#pragma omp parallel for
    for(int i=1; i<= trajec; i++){
        float X_T = gaussian_output(1);
        float inter = ev_non_sortie_delta( cond, X_T, W);
        mean_estim += inter;
        std_estim += pow(inter,2);
    }
    mean_estim = mean_estim/(trajec);
    std_estim = std_estim/trajec -pow(mean_estim,2);
    std_estim = (trajec/(trajec-1)) * std_estim;
    std_estim = sqrt(std_estim);
    err = 1.645*std_estim/sqrt(trajec);
    inf_born = mean_estim - err;
    sup_born = mean_estim + err;
    monte_carlo_output output(mean_estim, std_estim, inf_born , sup_born, err);
    return output;

}

float proba_non_sortie(initial_conditions cond,float X_T,gaussian_variable W){
    float r = W.mean;
    float sigma = W.std;
    int T = cond.T;
    float inter1 = cond.S_0*exp((r-pow(sigma,2)/2)*T+sigma*sqrt(T)*X_T);
    list<float> S_list = list<float> (cond.N_delta);
    list<float>::iterator iter= S_list.begin();
    int i =1;
    // First we create a  list of W_T_i output 
    for (iter= S_list.begin() ; iter != S_list.end(); iter++){
        if (iter==S_list.end()){
            *iter =  inter1;
        }
        else{
            float X_i = gaussian_output(1);
            float t = i*cond.delta;
            float c = cond.S_0*exp((r-pow(sigma,2)/2)*t+sigma*sqrt(t)*X_i);
            *iter = c;
            i ++;
        }
    };
    float S_min = mini( S_list);
    if (S_min>=cond.B){
        float prod_proba =1.;
        list<float>::iterator it= S_list.begin();
        float previous = *it;
        it ++;
        for(;it != S_list.end(); it ++){
            if (*it > cond.B && previous > cond.B){
                //cout<<(exp(-2*(*it-cond.B)*(previous- cond.B)*cond.delta))<<endl;
                prod_proba = prod_proba* (exp(-2*(*it-cond.B)*(previous- cond.B)*cond.delta));
            }
            previous = *it;
        }
        //cout<<exp(-r*T)*(cond.K-inter1)*prod_proba<<endl;
        return prod_proba;
    }
    else {
        return 0.;
    }
}

monte_carlo_output monte_carlo_proba_non_sortie(initial_conditions cond,gaussian_variable W, int trajec){
    float mean_estim =0.;
    float std_estim = 0.;
    float inf_born=0.;
    float sup_born =0.;
    float err =0.;
    //#pragma omp parallel for
    for(int i=0; i<= trajec; i++){
        float X_T = gaussian_output(1);
        float inter=proba_non_sortie( cond, X_T, W);
        mean_estim += inter;
        std_estim += pow(inter,2);
    }
    mean_estim = mean_estim/(trajec);
    std_estim = std_estim/trajec -pow(mean_estim,2);
    std_estim = (trajec/(trajec-1)) * std_estim;
    std_estim = sqrt(std_estim);
    err = 1.645*std_estim/sqrt(trajec);
    inf_born = mean_estim - err;
    sup_born = mean_estim + err;
    monte_carlo_output output(mean_estim, std_estim, inf_born , sup_born, err);
    return output;
}