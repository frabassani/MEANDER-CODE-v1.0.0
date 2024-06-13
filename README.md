# MEANDER-CODE-v1.0.0
First release of the Matlab MEANDER CODE for river meandering dynamics.

This Matlab code was edited by Carlo Camporeale and Francesca Bassani (contacts information below).
It implements the two meandering models:
(1) HIPS, developed by Hasegawa (1977) and Ikeda, Parker and Sawai (1981) and
(2) ZS, developed by Zolezzi and Seminara (2001).

It combines the codes for meandering dynamics of Camporeale et al. (2005) and the [Meander Centerline Migration Model](https://github.com/FluidMechanicsUNIPD/Meander-Centerline-Migration-Model) of Bogoni et al. (2017), which were both implemented in Fortran. 
Particularly, the bed characterization, the implementation of the ZS model and all related sub-functions are re-adapted from Bogoni et al. (2017), while all the other features (e.g., the numerical algorithm for the temporal evolution of the river planform, the local curvature, the shifting of the points according to the excess bank longitudinal velocities, the algorithm for cutoff search, etc.) are derived from Camporeale et al. (2005).
In addition, the code is parallelized and extended to the case of varying flows.

# CODE STRUCTURE
Define the settings to run the code as follows:
1. “Camporeale_Bassani_main”: set here (i) the input parameters (e.g., the bankfull discharge defining the river geometry, the mean sediment size, the coefficient of erodibility); (ii) the time step for meandering evolution (e.g., days); (iii) the option defining the meandering model (HIPS or ZS) and whether the forcing discharge is constant or variable; and (v) the number to be associated with the name of the simulation.  
2. “input_data”: set here (i) the bedload sediment transport formula; (ii) the option to extract the variable discharge from a Compound Poisson Process or from a Gamma-distributed noise; (iii) the option for savings; (iv) the simulation end time and time lapse for savings; (v) when to stop the simulation (short- or long-term).
The main file to be launched is: “Camporeale_Bassani_main”.

Once the previous two steps have been executed, the code starts the iteration over time to simulate the evolution of the river centerline. This is done in the script “meander”, where at each time step the computation of the excess bank velocity according to the HIPS or ZS model, the corresponding shifting of river points and the search for cutoffs are performed.
Note that savings can be adapted at lines 51-80 of the function Itcurv1. 



# Additional information

CONTACTS

For more information about this code, you can reach out to:

carlo.camporeale@polito.it

francesca.bassani@epfl.ch 




REFERENCES 

Bogoni, M., Putti M., and Lanzoni S. (2017). Modeling meandermorphodynamics over self-formed heterogeneous floodplains. Water Resources Research, 53, 5137–5157.

Camporeale, C., Perona, P., Porporato, A., and Ridolfi, L. (2005). On the long-term behavior of meandering rivers. Water Resources Research, 41 (12).

Hasegawa, K. (1977). Computer simulation of the gradual migration of meandering channels. Proceedings of the Hokkaido Branch, Japan Society of Civil Engineering, 197–202.

Ikeda, S., Parker, G., and Sawai, K. (1981). Bend theory of river meanders. Part 1. Linear development. Journal of Fluid Mechanics, 112, 363–377.

Zolezzi, G., & Seminara, G. (2001). Downstream and upstream influence in river meandering. Part 1. General theory and application to overdeepening. Journal of Fluid Mechanics, 438 (13), 183–211.
