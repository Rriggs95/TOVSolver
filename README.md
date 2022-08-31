# TOVSolver
<p>
TOV Solver is the Python code used for my undergraduate thesis for the University of Thessaly, were we explored solutions for the <a href="https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation">Tolman–Oppenheimer–Volkoff equation</a> using Polytropic Equations of State.</p>
<p>
The files above have the following function:
<ul>
<li>definitions.py: defines all the functions used for the solving of the Equations</li>
<li>TOV.py: solves the equations numerically using the second order Runge-Kutta and Dimensionless Units for distance (G = c = 1). Returns distance in km, mass in Solar Masses and Pressure in cgs</li>
<li>Plots.py: plots the Mass and Pressure as a function of distance</li>
<li>Metrics.py: Caluclates and plots the &lambda; and &nu; metrics</li>
<li>M_R_search.py: Searches for through the different values of K and n, plots the Mass and Pressure and saves the plot files for further review.</li>
<li>TOV.ipynb: All of the above added in a Jupyter notebook, used for testing code and minor tweaks
</p>
