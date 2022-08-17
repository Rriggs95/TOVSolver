###################################################################################################
############################################ IMPORTS ##############################################
###################################################################################################

import numpy as np
from definitions import *

###################################################################################################
############################################## DATA ###############################################
###################################################################################################

data = open("TOV.dat","w+")
data.write("xi: Dimensionless distance\ntheta: Dimensionless pressure T\nM: Dimensionless mass M\n\n\nxi\ttheta\tM\n")

###################################################################################################
##################################### NUMERICAL INTEGRATIONS ######################################
###################################################################################################

h = 0.00001
xi = 0.00001
theta = Theta(xi)
m = Mu(xi)
XI = []
THETA = []
MU = []
while theta>=0.0001:
    XI.append(xi)
    MU.append(m)
    THETA.append(theta)
    data.write(str(xi)+"\t"+str(theta)+"\t"+str(m)+"\n")
    m_k1 = dM(xi, theta)
    theta_k1 = dT(xi, theta, m)
    m_k2 = dM(xi+h, theta+theta_k1*h)
    theta_k2 = dT(xi+h, theta+theta_k1*h, m+m_k1*h)
    m = m + (m_k1+m_k2)*h/2
    theta = theta + (theta_k1+theta_k2)*h/2
    xi = xi + h
data.close()

R = []
M = []
P = []
for i in range(0, len(XI)-1, 1):
    R.append(RUnits(XI[i]))
    M.append(MUnits(MU[i]))
    P.append(PrUnits(THETA[i]))
