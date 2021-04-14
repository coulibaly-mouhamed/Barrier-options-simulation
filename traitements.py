
import numpy as np  
import matplotlib
import matplotlib.pyplot as plt

#lecture du fichier 
file = open("gaussian_test.txt")
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
plt.scatter(np.log(it),mean_est)
plt.errorbar(np.log(it),mean_est,yerr=err,fmt = 'none', capsize = 10, ecolor = 'red', zorder = 1)
plt.xlabel('log(N)')
plt.ylabel('Estimation de la moyenne par une méthode de Monte-Carlo' )
#res =np.array(L)
##print(np.sqrt(2))
N=30000
#T = np.linspace(0,N,N+1)
#plt.plot(T,res)
#plt.hist(res,40,density=True)              
plt.show()
