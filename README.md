# Tuberculosis
Simulation of the dynamics of tuberculosis lesions in mice lungs

A total of 5 MATLAB scripts are included.
4 of them are functions and 1 of them is a script for plotting.

OrganLung.m
	Main function that requires the initial parameters as inputs.
	It simulates the dynamics of tuberculosis including all hypotheses.
	It returns the different quantities to be plotted.

OrganLung2.m
	Main function that requires the initial parameters as inputs.
	It simulates the dynamics of tuberculosis removing only reinfection hypothesis.
	It returns the different quantities to be plotted.

OrganLung3.m
	Main function that requires the initial parameters as inputs.
	It simulates the dynamics of tuberculosis removing only coalescence hypothesis.
	It returns the different quantities to be plotted.

MC.m
	Auxiliary Monte Carlo function that allows to calculate areas.

OrganLung_Plot.m
	Script where 8 simulations are performed by calling the 4 aforementioned functions.
	Running this script takes of the order of one hour (time depends on the initial parameters).
	The results are the 5 different figures that are included in the article.
