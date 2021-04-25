
import numpy as np  
#import matplotlib
import matplotlib.pyplot as plt

#lecture du fichier 
def traitements( filename ):
    file = open(filename)
    output = file.readlines()
    file.close()
    it=[]
    mean_est=[]
    inf_born =[]
    sup_born =[]
    err =[]
    for char in output:
        L = []
        n =""
        i=0
        while(char[i] != '\n'):
            if (char[i]=='\t'):
                n = float (n)
                #print(n)
                L.append(n)
                n = ""
            else:    
                n += char[i]
            i+=1
        n = float(n)
        L.append(n)
        #print(L)
        if (len(L)>1):
            it.append(L[0])
            mean_est.append(L[1])
            inf_born.append(L[2])
            sup_born.append(L[3])
            err.append(L[4])
        else:
            true_value = L[0]
    it = np.array(it)
    mean_est = np.array(mean_est)
    inf_born = np.array(inf_born)
    sup_born = np.array(sup_born)
    err  = np.array(err)
    return it, mean_est,err

it,mean_est,err = traitements("Q_14.txt" )
it2, mean_est2, err2 = traitements("Q_9_2.txt")

def display_monte_carlo_comp(it,mean_est,err,it2,mean_est2,err2):
    plt.figure()
    plt.scatter(np.log(it)/np.log(10),mean_est)
    plt.scatter(np.log(it2)/np.log(10),mean_est2)
    plt.errorbar(np.log(it)/np.log(10),mean_est,yerr=err,fmt = 'none', capsize = 5, ecolor = 'blue', zorder = 1,label='Sans réduction de la variance ')
    plt.errorbar(np.log(it2)/np.log(10),mean_est2,yerr=err2,fmt = 'none', capsize = 5, ecolor = 'orange', zorder = 1,label= 'Avec réduction de la Variance')
    plt.xlabel('log(N)')
    plt.ylabel('Estimation de la moyenne par une méthode de Monte-Carlo' )
    plt.legend()
    plt.show()

def display_monte_carlo(it,mean_est,err):
    plt.figure()
    plt.scatter(np.log(it)/np.log(10),mean_est)
    plt.errorbar(np.log(it)/np.log(10),mean_est,yerr=err,fmt = 'none', capsize = 5, ecolor = 'blue', zorder = 1,label='Sans réduction de la variance ')
    plt.xlabel('log(N)')
    plt.ylabel('Estimation de la moyenne par une méthode de Monte-Carlo' )
    plt.show()

def display_monte_carlo_B(it,mean_est):
    plt.figure()
    plt.scatter(it,mean_est,c='blue',label='S_0=1')
    plt.errorbar(it,mean_est,yerr=err)
    plt.xlabel('B')
    plt.ylabel('Estimation de la moyenne par une méthode de Monte-Carlo')
    plt.legend()
    plt.show()
    
def display_monte_carlo_B_2(it,mean_est,it2,mean_est2):
    plt.figure()
    plt.scatter(it,mean_est,c='blue',label='S_0=1')
    plt.scatter(it2,mean_est2,c='green',label='S_0=0.8')
    #plt.errorbar(it,mean_est,yerr=err)
    plt.xlabel('Sigma')
    plt.ylabel('Estimation de la moyenne par une méthode de Monte-Carlo')
    plt.legend()
    plt.show()

display_monte_carlo(it,mean_est,err)
#display_monte_carlo_B(it,mean_est)
#display_monte_carlo_B_2(it,mean_est,it2,mean_est2)
""" plt.figure()
plt.scatter(np.log(it),err)
plt.scatter(np.log(it2),err2)
plt.show() """
#T = np.linspace(0,N,N+1)
#plt.plot(T,res)
#plt.hist(res,40,density=True)              

