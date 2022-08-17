###################################################################################################
############################################ IMPORTS ##############################################
###################################################################################################

import numpy as np

###################################################################################################
#################################### CONSTANTS INITIALIZATION #####################################
###################################################################################################

K = 3.0*10**(22)
n = 4/7
a_0 = 0.02
m_3 = 4*np.pi/3
theta_2 = -2*np.pi*((1+a_0)*(3+a_0))/(3*a_0*(n+1))
m_5 = 4*np.pi*n*theta_2/5
theta_4 = -theta_2*(m_3+4*np.pi/a_0)/(2*(n+1))

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

def MUnits(M):
    """Conversion of Mass units from dimensionless in G = c = 1 to Solar Masses

    :M: Dimensionless mass
    :returns: m: Mass in Solar Masses

    """
    return M/np.sqrt((a_0*K)**(-1.0*n))*1.3466*10**28*5*10**(-34)

def RUnits(xi):
    """Conversion of Distance from dimensionless in G = c = 1 to km

    :xi: Dimensionless Distance
    :returns: R: Distance in km

    """
    return xi/np.sqrt((a_0*K)**(-1.0*n))*10**(-5)

def PrUnits(t):
    """Converion of Pressure from dimensionless in G = c = 1 to cgs

    :t: Dimensionless Pressure
    :returns: P: Pressure in cgs

    """
    return (a_0*K)**(-1.0*(n+1))*t**(n+1)/(8.263*10**(-50))
