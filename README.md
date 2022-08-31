# TilmanModelSim
Simulation of the Tilman model for 3 species and 2 resources, implemented in the R programming language.

The code used for the simulations is shown, comprising several functions:

- hazparam(): Function that generates the values of the necessary parameters to solve the differential equations according to the criteria described in the Supplementary material.
      
Receives: 

nsp (int): number of species to use in the simulation; for the work described in the appendix, 3 species were used.
      
Returns:

list with the values of the parameters vr: vector of birth rates, vD: vector with mortality rates, mK: matrix with the values of k_ij for the 3 species with the 2 resources, mY: matrix with the values of Y_ij of the 3 species with the 2 resources; Fr: vector with the supply rates of the two resources; mRa: matrix with the values of R* for each pair of resource species.
 
 
- plotTilman(): Function that determines the ZNGI, the consumption vectors, the supply line for sampling and the coexistence regions.
       
Receives:
           
vr (vector): vector with birth rates.

vD (vector): vector with death rates.

mK (array): matrix with the k_ij values of the 3 species with the 2 resources.

mY (array): matrix with the Y_ij values of the 3 species with the 2 resources.
       
Returns:
           
list with xord: supply line coordinates; xC1 and yC1: coordinates of the consumption vectors of resource 1; xC2 and yC2: coordinates of the consumption vectors of resource 1.
 
 
 
 
 
 
- muestra(): Function that performs the sampling in the supply line.
        
Receives:
           
xord (): supply line
           
intens (int): number of points to sample
        
        
         
         
      
- minsoft(): Function to minimize the minimum for solving differential equations.





- convpar(): Function that names the parameters generated in hazparam().
        
Receives:
           
par (list): list of parameters generated in hazparam()
        
Returns:
           
list of named parameters.           
     
     
     
- Tilman(): Function that establishes the differential equations of the model for its resolution.
         
Receives:
           
t (seq): time sequence for which output is desired; the first time value is the start time.
           
state (list): initial state values of resource concentrations and species densities.
           
parameters (list): list with generated parameter values.
         
Returns:
         
object of class deSolve with an array containing the values of the species densities and resource concentrations at the specified output times.
   
   
   

- hazdat(): Function that generates the parameters of an artificial community, determines the supply region for the sampling of supply points (to establish the initial states of the resources) and solves the differential equations to obtain the dynamics of the species and resources.    
          
Returns:
          
list with par (list): parameter values generated in hazparam(); npar (list): list with the sample data generated in plotTilman(); y salt (array): array containing the values of the species densities and resource concentrations at the specified output times.
           
           
           
           
- Ra_est(): Function that estimates the values of R* for each species-resource pair according to the criteria described in the methods with different threshold values of the maximum density.       
           
Receives: 
           
comp (list): output of hazdat()
            
umbrales (vector): vector with threshold values for estimation.
                      
Returns:
           
Ras_mat (array): matrix with the R* values of each species-resource pair estimated for each threshold value.
            
            
            
       
- explore(): Function that generates the dynamics for the three species and the two resources using the previous functions to generate an artificial community, estimates the R* values of each species-resource pair, calculates the error associated with the estimates of each threshold, and calculates the correlation coefficient between the simulated R* values and the estimated R* values.
                  
Returns:
            
list with: dens_max: maximum densities of the species; real_Ra (matrix): simulated R* values of each species for each resource; est_RAS (matrix): estimated values of R* with each threshold; error_RAS (matrix): error of each estimation; and cor_coefs (vector): correlation coefficients between the R*, both simulated and estimated.
              
              
          
