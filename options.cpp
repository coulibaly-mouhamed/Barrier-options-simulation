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