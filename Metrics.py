###################################################################################################
############################################ IMPORTS ##############################################
###################################################################################################

from definitions import *
from matplotlib import pyplot as plt
import numpy as np
from TOV import XI, THETA, MU

###################################################################################################
############################################ METRICS ##############################################
###################################################################################################

data = open("metric.dat","w+")
data.write("\n\nxi\tlamda(xi)\tnu(xi)\n")
lamda = [] 
nu = []
count = 0
for i in XI:
    lamda.append(np.log(1/(1-2*MU[count]/i))/2)
    nu.append(np.log((a_0/(a_0+1))**(2*(n+1))*a_0/(a_0+2*XI[-1]*(n-1)*dT(XI[-1], THETA[-1], MU[-1]))*((a_0+1)/(a_0+THETA[count])**(2*(n+1))))/2)
    data.write(str(i)+"\t"+str(lamda[count])+"\t"+str(nu[count])+"\n")
    count += 1

# Plots
plt.plot(XI, lamda)
plt.xlabel("xi")
plt.ylabel("lamda")
plt.title("Lamda metric - Distance")
name = "Plots/Lambda_n_{}_K_{}.pdf".format(n,K)
plt.savefig(name)
plt.clf()
plt.plot(XI, nu)
plt.xlabel("xi")
plt.ylabel("nu")
plt.title("Nu metric - Distance")
name = "Plots/Nu_n_{}_K_{}.pdf".format(n,K)
plt.savefig(name)
data.close()
