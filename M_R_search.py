###################################################################################################
############################################ IMPORTS ##############################################
###################################################################################################

from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from matplotlib.collections import LineCollection
import numpy as np
import csv

###################################################################################################
#################################### CONSTANTS INITIALIZATION #####################################
###################################################################################################

h = 0.0001
m_3 = 4*np.pi/3

###################################################################################################
##################################### DIFFERENTIAL EQUATIONS ######################################
###################################################################################################

def Theta(xi):
    """Gives Taylor expansion for dimensionless pressure.

    :xi: Dimensionless distance.
    :returns: T where p = p_0*T(xi)^(n+1), p0 = K*rho_0^gamma, where gamma = n + 1/n

    """
    return 1 + theta_2*(xi**2) + theta_4*(xi**4)

def Mu(xi):
    """Gives Taylor expansion for dimensionless mass.

    :xi: Dimensionless distance.
    :returns: Mu where Mu = sqrt(epsilon_0)*m, m: mass of the star.

    """
    return m_3*(xi**3) + m_5*(xi**5)

def dM(xi, t):
    """Differential equation for dimensionless mass

    :xi: Dimensionless distance
    :t: Dimensionless pressure
    :returns: dM/dxi = 4*Pi*xi^2*(T(xi)^n)

    """
    return 4*np.pi*xi**2*t**n

def dT(xi, t, mu):
    """Differential equation for dimensionless pressure

    :xi: Dimensionless distance
    :t: Dimensionless pressure
    :mu: Dimensionless mass
    :returns: dT/dxi = -[(a_0+T(xi))/(n+1)]*[(M(xi)+4Pi*xi^3(T(xi)^(n+1))/a_0)/(xi(xi-2M(xi)))]

    """
    return -1.0*(a_0+t)/(n+1)*(mu+4*np.pi*xi**3*t**(n+1)/a_0)/(xi*(xi-2*mu))

###################################################################################################
######################################## UNIT CONVERSION ##########################################
###################################################################################################

def M_to_solar(M):
    """Conversion of Mass units from dimensionless in G = c = 1 to Solar Masses

    :M: Dimensionless mass
    :returns: m: Mass in Solar Masses

    """
    return M/np.sqrt(epsilon_0)*1.3466*10**28*5*10**(-34)
def R_to_km(xi):
    """Conversion of Distance from dimensionless in G = c = 1 to km

    :xi: Dimensionless Distance
    :returns: R: Distance in km

    """
    return xi/np.sqrt(epsilon_0)*10**(-5)

###################################################################################################
############################################ MAIN FILE ############################################
###################################################################################################

def main():
# Double loop for Plots
    for K in np.arange(3*10**8, 4*10**8, 2*10**7):
        print("\n\nK={}".format(K))
        for n in np.arange(1.40, 1.52, 0.01):
            XI = []
            MU = []
            print("\n\nn={}".format(n))
            for a_0 in np.arange(0.01, 0.04, 0.0001):
                xi = 0.00001
                theta_2 = -2*np.pi*((1+a_0)*(3+a_0))/(3*a_0*(n+1))
                m_5 = 4*np.pi*n*theta_2/5
                theta_4 = -theta_2*(m_3+4*np.pi/a_0)/(2*(n+1))
                theta = Theta(xi)
                M = Mu(xi)
                epsilon_0 = (a_0*K)**(-1.0*n)

# TOV Solutions
                while theta>=0.0001:
                    M_k1 = dM(xi, theta)
                    theta_k1 = dT(xi, theta, M)
                    M_k2 = dM(xi+h, theta+theta_k1*h)
                    theta_k2 = dT(xi+h, theta+theta_k1*h, M+M_k1*h)
                    M = M + (M_k1+M_k2)*h/2
                    theta = theta + (theta_k1+theta_k2)*h/2
                    xi = xi + h
                print("R={}\tM={}\ta_0={}".format(R_to_km(xi), M_to_solar(M), a_0))

# Keep total Mass and Radius in XI and MU respectively
                XI.append(R_to_km(xi))
                MU.append(M_to_solar(M))

# Finding the Max Value of MU
            i = 0
            while True:
                if MU[i] > MU[i+1]:
                    break
                i = i + 1

            # Finding the value of R for which M[R] = 1.4
            # if MU[i] > 1.4:
                # j = 0
                # while True:
                    # if MU[j] > 1.4:
                        # break
                    # j = j + 1
                # if 1.4-MU[j-1] < MU[j] - 1.4:
                    # j = j-1

# Plots
            leg = "n = {}".format(n)
            plt.plot(XI,MU,label=leg)
            ax = plt.axes()
            title = "Classic M R Plot in proper units"
            plt.title(title)
            plt.legend(loc="lower right")
            plt.xlabel("R (km)")
            plt.ylabel("M (Msolar)")
            ax.xaxis.set_major_locator(ticker.FixedLocator([XI[i]]))
            ax.yaxis.set_major_locator(ticker.FixedLocator([1.4,MU[i]]))
            l1 = [(XI[i],0),(XI[i],MU[i])]
            l2 = [(0,MU[i]),(XI[i],MU[i])]
            lc1 = LineCollection([l1,l2])
            plt.gca().add_collection(lc1)
            # if MU[i] > 1.4:
                # ax.xaxis.set_major_locator(ticker.FixedLocator([XI[i],XI[j]]))
                # l3 = [(XI[j],0),(XI[j],1.4)]
                # l4 = [(0,1.4),(XI[j],1.4)]
                # lc2 = LineCollection([l3,l4])
                # plt.gca().add_collection(lc2)
            # else:
                # ax.xaxis.set_major_locator(ticker.FixedLocator([XI[i]]))
            name = "Plots/n_{}_K_{}_proper_units.pdf".format(n, K)
            plt.savefig(name)
            plt.clf()
            print("\n\n")
if __name__=='__main__':
    main()
