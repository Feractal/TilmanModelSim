# TilmanModelSim
Simulación del modelo de Tilman para 3 especies y 2 recursos, implementado en el lenguaje de programación R.

Se muestra el código utilizado para las simulaciones, que comprende varias funciones:


- hazparam(): Función que genera los valores de los parámetros necesarios para la resolución de las ecuaciones diferenciales de acuerdo a los criterios descritos en el apéndice A.
      
Recibe: 

nsp (int): número de especies a utilizar en la simulación; para el trabajo descrito en el apéndice, se utilizaron 3 especies.
      
Regresa:

lista con los valores de los parámetros vr: vector de tasas de natalidad, vD: vector con tasas de mortalidad, mK: matriz con los valores de k_ij para  las 3 especies con los 2 recursos, mY: matriz con los valores de Y_ij de las 3 especies con los 2 recursos; Fr: vector con las tasas de suministro de              los dos recursos; mRa: matriz con los valores de R* para cada par de especie recurso.
 
 
- plotTilman(): Función que determina las ZNGI, los vectores de consumo, la línea de abastecimiento para el muestro y las regiones de coexistencia.
       
Recibe:
           
vr (vector): vector con tasas de natalidad

vD (vector): vector con tasas de mortalidad

mK (array): matriz con los valores de k_ij de las 3 especies con los 2 recursos

mY (array): matriz con los valores de Y_ij de las 3 especies con los 2 recursos
       
       
Regresa:
           
lista con xord: coordenadas de la línea de abastecimiento; xC1 y yC1: coordenadas de los vectores de consumo del recurso 1; xC2 y yC2: coordenadas de                los vectores de consumo del recurso 1
 
 
- muestra(): Función que hace el muestreo en la línea de abastecimiento-
        
Recibe:
           
xord (): línea de abastecimiento
           
intens (int): número de puntos a muestrear
           
           
      
- minsoft(): Función para minimizar el mínimo para la resolución de las ecuaciones diferenciales.



- convpar(): Función que nombra los parámetros generados en hazparam().
        
Recibe:
           
par (list): lista de parámetros generados en hazparam()
        
Regresa:
           
lista de parámetros nombrados
           
           
- Tilman(): Función que establece las ecuaciones diferenciales del modelo para su resolución.
         
Recibe:
           
t (seq): secuencia de tiempo para la que se desea salida; el primer valor de tiempos es el tiempo inicial
           
state (list): valores de los estados iniciales de las concentraciones de recursos y densidades de especies
           
parameters (list): lista con los valores de los parámetros generados
         
Regresa:
         
objeto de clase deSolve con una matriz que contiene los valores de las densidades de especies y concentraciones de recursos en los tiempos de salida especificados
             

- hazdat(): Función que genera los parámetros de una comunidad artificial, determina la región de abastecimiento para el muestreo de puntos de suministro (para    establecer los estados iniciales de los recursos) y resuelve las ecuaciones diferenciales para obtener la dinámica de las especies y recursos.
          
          
Regresa:
          
lista con par (list): valores de los parámetros generados en hazparam(); npar (list): lista con los datos del muestreo generados en plotTilman(); y sal (array): matriz que contiene los valores de las densidades de especies y concentraciones de recursos en los tiempos de salida especificados
           
           
- Ra_est(): Función que estima los valores de R* para cada par de especie-recurso de acuerdo a los criterios descritos en los métodos con diferentes valores     umbrales de la densidad máxima.
           
           
Recibe: 
           
comp (list): salida de hazdat()
            
umbrales (vector): vector con los valores umbrales para la estimación
           
           
Regresa:
           
Ras_mat (array): matriz con los valores de R* de cada par de especie-recurso estimados para cada valor umbral
             
       
- explore(): Función que genera la dinámica para las tres especies y los dos recursos utilizando las funciones anteriores para generar una comunidad artificial,    estima los valores de R* de cada par de especie-recurso, calcula el error asociado a las estimaciones de cada umbral, y calcula el coeficiente de correlación entre los valores de R* simulados y los valores de R* estimados.  
            
            
Regresa:
            
lista con: dens_max: densidades maximas de las especies; real_Ra (matrix): valores de R* simulados de cada especie para cada recurso; est_RAS (matrix): valores estimados de R* con cada umbral; error_RAS (matrix): error de cada estimación; y cor_coefs (vector): coeficientes de correlación entre los R*, tanto los simulados como los estimados.
              
              
          
