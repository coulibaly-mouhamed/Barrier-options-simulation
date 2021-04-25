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
///////////////////////////////////////////////////////////////////////////////////////////////////
/// We first compute two independent random variables which are  uniformly distributed/////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

float gaussian_output(float T){
    // We use Box-Muller Technique to compute random gaussian variable
    float epsilon =1E-4;
    float U =((float)(rand()+1.0))/(RAND_MAX + 1.0) ;
    float V = ((float)(rand()+1.0))/(RAND_MAX+1.0) ;
    float X = sqrt(-2*log(U))*cos(2*PI*V);
    //float Y = sqrt(-2*log(U))*sin(2*PI*V);
    return X*sqrt(T); //On retourne une gaussienne centrée réduite faut multiplier par T aprés
}

//////////////////////////////////////////////////////////////
//We compute the Abramowitz and Stegun approximation//////////
//////////////////////////////////////////////////////////////

float abramowitz_stegun_int( float x){
    float b_0 =0.2316419;
    float b_1 = 0.319381530;
    float b_2 = -0.356563782;
    float  b_3 = 1.781477937 ;
    float b_4 = -1.821255978;
    float  b_5 = 1.330274429;
    float t = 1./(1+b_0*x);
    float poly_member = b_1*pow(t,1)+b_2*pow(t,2)+b_3*pow(t,3)+b_4*pow(t,4)+b_5*pow(t,5);
    float right_member = (1./sqrt(2*PI))*exp(-1./2*pow(x,2))*poly_member;
    return 1-right_member;   
}

//////////////////////////////////////////
/// Fonctions Auxiliaires ///////////////
/////////////////////////////////////////

float mini(list<float> S){
    list<float>::iterator it = S.begin();
    float min = *it;
    for(;it != S.end(); it++){
        if (*it < min){min = *it;}
    }
    return min;
}

//////////////////////////////
// Vanilla european option////
float vanilla_eur_option(gaussian_variable W,float T,float K,float S_0){
    /*
    Compute the price of a vanilla european option 
        Input: T the time when we want to apply the call
               sigma: the standar deviation of the borwnian path 
               r : the mean of the \mu in the Black Scholes model 
        Output : 
            The price P_euro 
    */
    float r = W.mean;
    float sigma  = W.std;
    float theta_1 = (1./sigma)*(log(S_0/K)+(r+pow(sigma,2)/2)*T)*(1/sqrt(T));
    float theta_2 = (1./sigma)*(log(K/S_0)-(r-pow(sigma,2)/2)*T)*(1/sqrt(T));
    float P_1 = S_0*(1- abramowitz_stegun_int(theta_1));
    float P_2 = K*exp(-r*T)*abramowitz_stegun_int(theta_2);
    return P_2-P_1;
}
///////////////////////////////////////////////////
//////Monte-Carlo method for Vanilla options//////
//////////////////////////////////////////////////

monte_carlo_output monte_carlo_vanilla_options(initial_conditions cond, gaussian_variable W, int N){
      int T = cond.T;
      float sigma = W.std;
      float S_0 = cond.S_0;
      float K = cond.K;
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
      float err = 1.645*std_estim/sqrt(N);
      float inf_born = mean_estim - err;
      float sup_born = mean_estim + err;
      monte_carlo_output output(mean_estim, std_estim, inf_born , sup_born, err);
      return output;
    
}

////////////////////////////////////////////////////
/// Monte-Carlo method for  Down and out options////
////////////////////////////////////////////////////

//Fonction auxiliaire h
float h(initial_conditions cond, float X, gaussian_variable W){
   float r = W.mean;
   float sigma = W.std;
   int T = cond.T;
   float inter1 = cond.S_0*exp((r-pow(sigma,2)/2)*T+sigma*sqrt(T)*X);
   if (cond.K<inter1){
       return 0.;
   } 
   else {
       // Simulation de l'indicatrice du minimum
        list<float> S_list = list<float> (cond.N_delta);
        list<float>::iterator iter= S_list.begin();
        int i =1;
        // First we create a  list of W_T_i output 
        for (iter= S_list.begin() ; iter != S_list.end(); iter++){
            float t = i*cond.delta;
            float c = cond.S_0*exp((r-pow(sigma,2)/2)*t+sigma*sqrt(t)*X);
            *iter = c;
            i ++;
        };
        float S_min = mini( S_list);
        if (S_min>=cond.B){
            return exp(-r*T)*(cond.K-inter1);
        }
        else {
            return 0.;
        }
   }
}

// Fonction auxiliaire h_anti qui sert à calculer P_DI
float h_anti(initial_conditions cond, float X, gaussian_variable W){
   float r = W.mean;
   float sigma = W.std;
   int T = cond.T;
   float inter1 = cond.S_0*exp((r-pow(sigma,2)/2)*T+sigma*sqrt(T)*X);
   if (cond.K<inter1){
       return 0.;
   } 
   else {
       // Simulation de l'indicatrice du minimum
        list<float> S_list = list<float> (cond.N_delta);
        list<float>::iterator iter= S_list.begin();
        int i =1;
        // First we create a  list of W_T_i output 
        for (iter= S_list.begin() ; iter != S_list.end(); iter++){
            float t = i*cond.delta;
            float c = cond.S_0*exp((r-pow(sigma,2)/2)*t+sigma*sqrt(t)*X);
            *iter = c;
            i ++;
        };
        float S_min = mini( S_list);
        if (S_min<=cond.B){
            return exp(-r*T)*(cond.K-inter1);
        }
        else {
            return 0.;
        }
   }
}

monte_carlo_output monte_carlo_option_down_put(initial_conditions cond,gaussian_variable W, int trajec){
    float r = W.mean;
    float sigma = W.std;
    float mean_estim =0;
    float std_estim =0;
    float inf_born =0;
    float sup_born =0;
    float err =0;
    for(int j=1; j<= trajec ;j++){
        float X_i = gaussian_output(1);
        float inter = h(cond,X_i,W);
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

// Reducion de la variance par variables antithétiques 


monte_carlo_output anti_monte_carlo_option_down_put(initial_conditions cond,gaussian_variable W, int trajec){
    float mean_estim =0.;
    float std_estim =0.;
    float err =0.;
    float inf_born=0.;
    float sup_born =0.;
    for (int i=1; i<= trajec/2; i++){
        float X_i = gaussian_output(1); // Gaussienne centrée réduite 
        float inter = h(cond, X_i,W)+h(cond, -X_i, W);
        mean_estim += inter;
        std_estim  += pow(inter,2);
    }
    mean_estim = 1./trajec* mean_estim;
    std_estim = 1./trajec* std_estim;
    std_estim = std_estim - pow(mean_estim,2);
    std_estim = (trajec/(trajec-1) )* std_estim;
    err = 1.645*sqrt(std_estim)/sqrt(trajec);
    inf_born = mean_estim - err;
    sup_born = mean_estim +err;
    monte_carlo_output output(mean_estim, std_estim, inf_born , sup_born, err);
    return output;
}

// Reduction de la variance par variable de controle 
monte_carlo_output control_monte_carlo_down_in(initial_conditions cond, gaussian_variable W, int trajec){
    float r = W.mean;
    float sigma = W.std;
    float mean_estim =0;
    float std_estim =0;
    float inf_born =0;
    float sup_born =0;
    float err =0;
    float price =vanilla_eur_option( W, cond.T, cond.K, cond.S_0);
    #pragma omp parallel for
    for(int j=1; j<= trajec ;j++){
        float X_i = gaussian_output(1);
        float inter = h_anti(cond,X_i,W);
        mean_estim += inter;
        std_estim += pow(inter,2);
    }
    float inter = mean_estim/trajec;
    mean_estim = price - mean_estim/(trajec);
    std_estim = pow(price,2)-2*price*inter +std_estim/trajec -pow(mean_estim,2);
    //cout<<std_estim<<endl;
    std_estim = (trajec/(trajec-1)) * std_estim;
    std_estim = sqrt(std_estim);
    err = 1.645*std_estim/sqrt(trajec);
    inf_born = mean_estim - err;
    sup_born = mean_estim + err;
    monte_carlo_output output(mean_estim, std_estim, inf_born , sup_born, err);
    return output;
}