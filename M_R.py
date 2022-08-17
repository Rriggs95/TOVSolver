###################################################################################################
############################################ IMPORTS ##############################################
###################################################################################################

from matplotlib import pyplot as plt
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

# Differential Equations
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
############################################ MAIN FILE ############################################
###################################################################################################

def main():
# Double loop for Plots
    for n in np.arange(1.5, 1.8, 0.01):
        XI = []
        MU = []
        print("n={}\n".format(n))
        for a_0 in np.arange(0.00001, 0.0006, 0.00001):
            xi = 0.00001
            theta_2 = -2*np.pi*((1+a_0)*(3+a_0))/(3*a_0*(n+1))
            m_5 = 4*np.pi*n*theta_2/5
            theta_4 = -theta_2*(m_3+4*np.pi/a_0)/(2*(n+1))
            theta = Theta(xi)
            M = Mu(xi)

# TOV Solutions
            while theta>=0.0001:
                M_k1 = dM(xi, theta)
                theta_k1 = dT(xi, theta, M)
                M_k2 = dM(xi+h, theta+theta_k1*h)
                theta_k2 = dT(xi+h, theta+theta_k1*h, M+M_k1*h)
                M = M + (M_k1+M_k2)*h/2
                theta = theta + (theta_k1+theta_k2)*h/2
                xi = xi + h
            print("xi={}\tM={}\ta_0={}".format(xi, M, a_0))

# Keep total Mass and Radius in XI and MU respectively
            XI.append(xi)
            MU.append(M)

# Plots
        lbl = "n = {}".format(n)
        plt.plot(XI,MU,label=lbl)
        print("\n\n")
        title = "Classic_M_R_Plot"
        plt.title(title)
        plt.legend(loc="lower right")
        plt.xlabel("xi")
        plt.ylabel("M")
        name = "Plots/n_{}.pdf".format(n)
        plt.savefig(name)
        plt.clf()

if __name__ == '__main__':
    main()
