
import numpy as np  
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

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

#it,mean_est,err = traitements("Q_4.txt" )
#it2, mean_est2, err2 = traitements("Q_7_normal.txt")

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
   # plt.errorbar(it,mean_est,yerr=err)
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

def display_q_4(filename):
    file = open(filename)
    output = file.readlines()
    file.close()
    N=[]
    est=[]
    err=[]
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
            N.append(L[0])
            est.append(L[1])
            err.append(L[2])
        else:
            true_value = L[0]
    est1 = np.array(est)
    N = np.array(N)
    err = np.array(err)
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(np.log(N)/np.log(10),est,c='blue',label=r"$P^{euro}$")
    plt.errorbar(np.log(N)/np.log(10),est,yerr=err,fmt = 'none',zorder = 1,capsize = 5)
    plt.axhline(y=0.0692719,c='red')
    plt.xlabel(r'$\log_{10}(N)$')
    plt.ylabel(r'Estimation du prix')
    plt.legend()
    plt.show()  
    
def display_q_7(filename):
    file = open(filename)
    output = file.readlines()
    file.close()
    N=[]
    est1=[]
    est2 =[]
    err1=[]
    err2=[]
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
            N.append(L[0])
            est1.append(L[1])
            est2.append(L[2])
            err1.append(L[3])
            err2.append(L[4])
            
        else:
            true_value = L[0]
    est1 = np.array(est1)
    N = np.array(N)
    est2 = np.array(est2)
    err1 = np.array(err1)
    err2 = np.array(err2)
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(np.log(N)/np.log(10),est1,label=r"Méthode Classique")
    plt.scatter(np.log(N)/np.log(10),est2,label=r"Avec variable de contrôle")
    plt.errorbar(np.log(N)/np.log(10),est1,yerr=err1,fmt = 'none',zorder = 1,capsize = 5)
    plt.errorbar(np.log(N)/np.log(10),est2,yerr=err2,fmt = 'none',zorder = 1,capsize = 5)
    plt.xlabel(r'$\log_{10}(N)$')
    plt.ylabel(r'Estimation du prix')
    plt.legend()
    plt.show()  
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(np.log(N)/np.log(10),err1,label=r"Méthode Classique")
    plt.scatter(np.log(N)/np.log(10),err2,label=r"Avec variable de contrôle")
    plt.xlabel(r'$\log_{10}(N)$')
    plt.ylabel(r'Erreur')
    plt.legend()
    plt.show()  

def display_q_8(filename):
    file = open(filename)
    output = file.readlines()
    file.close()
    B=[]
    est=[]
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
            B.append(L[0])
            est.append(L[1])
            err.append(L[2])
            
        else:
            true_value = L[0]
    B = np.array(B)
    est = np.array(est)
    err = np.array(err)
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(B,est,c='blue',label=r"$P^{DO \Delta}$ $S_0=0.8$")
    #plt.errorbar(delta,est2,yerr=err2,fmt = 'none')
    plt.xlabel(r'$\sigma$')
    plt.ylabel(r'Estimation du prix')
    plt.legend()
    plt.show()
    
def display_q_11(filename):
    file = open(filename)
    output = file.readlines()
    file.close()
    tau=[]
    est1=[]
    est2 =[]
    delta =[]
    err1=[]
    err2=[]
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
            tau.append(L[0])
            est1.append(L[1])
            est2.append(L[2])
            delta.append(L[3])
            err1.append(L[4])
            err2.append(L[5])
            
        else:
            true_value = L[0]
    tau = np.array(tau)
    est1 = np.array(est1)
    est2 = np.array(est2)
    delta  = np.array(delta)
    err1 = np.array(err1)
    err2 = np.array(err2)
    
    
    
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(delta,est2,c='blue',label=r"$P^{DO \Delta}$")
    plt.scatter(delta,est1,c='red',label=r"$P^{DO}$")
    plt.errorbar(delta,est1,yerr=err1,fmt = 'none')
    plt.errorbar(delta,est2,yerr=err2,fmt = 'none')
    plt.xlabel(r'$1/\Delta$')
    plt.ylabel(r'Estimation de la probabilté de non sortie')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(delta, np.abs(est2-est1))
    plt.xlabel(r'$1/\Delta$')
    plt.ylabel(r'Eccart')
    plt.show()
    
def display_q_12(filename):
    file = open(filename)
    output = file.readlines()
    file.close()
    tau=[]
    est1=[]
    est2 =[]
    delta =[]
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
            tau.append(L[0])
            est1.append(L[1])
            est2.append(L[2])
            delta.append(L[3])
            
        else:
            true_value = L[0]
    tau = np.array(tau)
    est1 = np.array(est1)
    est2 = np.array(est2)
    delta  = np.array(delta)
    
    
    
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.scatter(delta,est2,c='blue',label=r"$\zeta^{DO \Delta}$")
    plt.scatter(delta,est1,c='red',label=r"$\zeta^{DO}$")
    plt.xlabel(r'$1/\Delta$')
    plt.ylabel(r'Estimation de la probabilté de non sortie')
    plt.legend()
    plt.show()
    

display_q_7('Q_14_2.txt')
#display_monte_carlo(it,mean_est,err)
#display_monte_carlo_B(it,mean_est)
#display_monte_carlo_comp(it2,mean_est2,err2,it,mean_est,err)
""" plt.figure()
plt.scatter(np.log(it),err)
plt.scatter(np.log(it2),err2)
plt.show() """
#T = np.linspace(0,N,N+1)
#plt.plot(T,res)
#plt.hist(res,40,density=True)              

