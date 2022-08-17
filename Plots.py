###################################################################################################
############################################ IMPORTS ##############################################
###################################################################################################

from matplotlib import pyplot as plt
from TOV import R, M, P
from definitions import n,K

###################################################################################################
############################################# PLOTS ###############################################
###################################################################################################

plt.plot(R,M)
plt.xlabel("R(km)")
plt.ylabel("M(Solar)")
plt.title("Mass - Distance")
name = "Plots/M_n_{}_K_{}.pdf".format(n,K)
plt.savefig(name)
plt.clf()
plt.plot(R,P)
plt.xlabel("R(km)")
plt.ylabel("Pressure (dyn/cm^2)")
plt.title("Pressure - Distance")
name = "Plots/Pressure_n_{}_K_{}.pdf".format(n,K)
plt.savefig(name)
