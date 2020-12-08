library(deSolve)
library(easyNCDF)

#FUNCIONES

hazparam = function(nsp){
  #Genera parametros de acuerdo a los criterios descritos en el apendice
  #Recibe:
  ### nsp: numero de especies
  #Regresa:
  ### Lista con los valores de los parametros
  D = runif(1)*2                             #Tasa de mortalidad (es igual para todas las especies; D_i)
  r = D + runif(nsp)*2 + 1                   #Tasa de crecimiento de las tres especies (deben ser mayores que D; r_i)
  Ra1 = sort(runif(nsp)*2)                   #Valor de R* del recurso 1 para cada especie (R*_i1)
  Ra2 = sort(runif(nsp)*2,decreasing = T)    #Valor de R* del recurso 2 para cada especie (R*_i2)
  K1 = Ra1*(r-D)/D                           #Constante de saturacion media de cada especie limitada por recurso 1 (k_i1)
  K2 = Ra2*(r-D)/D                           #Constante de saturacion media de cada especie limitada por recurso 2 (k_i2)
  m = sort(tan(runif(nsp)*(pi/2)))           #Pendiente de vector de consumo para calculo de Y_ij
  Y2 = runif(nsp)*3                          #Cantidad de recurso 2 necesario para producir un individuo de cada especie (Y_i1)
  Y1 = ((K2+m)*Y2)/(K1+1)                    #Cantidad de recurso 1 necesario para producir un individuo de cada especie (Y_i2)
  Fr = runif(1)*2                            #Tasa de suministro de recurso (es la misma para ambos recursos)
  list(vr = r, vD = rep(D,nsp), vF = rep(Fr,2), mK = cbind(K1,K2), mY = cbind(Y1,Y2), mRa = cbind(Ra1,Ra2))
}


plotTilman = function(vr,vD,mK,mY){
  #Genera graficos con las ZNGI, vectores de consumo y la region de abastecimiento
  #Recibe:
  ### vr: vector con las tasas de crecimiento de las especies (r_i)
  ### vD: vector con las tasas de mortalidad de las especies (D_i)
  ### mK: matriz con las constantes de saturacion media de cada especie para cada recurso (k_ij)
  ### mY: matriz con las cantidades de recurso necesarias para producir un individuo de cada especie (Y_ij)
  
  nsp = dim(mK)[1]
  mr = cbind(vr, vr)
  mD = cbind(vD, vD)
  Ra = mD*mK/(mr-mD)
  mr = mr[order(Ra[,1]),]
  mD = mD[order(Ra[,1]),]
  mY = mY[order(Ra[,1]),]
  mK = mK[order(Ra[,1]),]
  Ra = Ra[order(Ra[,1]),]
  
  #Distancias entre los puntos de interseccion de las ZNGI y
  eles = sqrt(Ra[,1]^2 + Ra[,2]^2)
  inter = sqrt(Ra[2:nsp,1]^2 + Ra[1:(nsp-1),2]^2)
  distmax = c(eles, inter)[which.max(c(eles, inter))]
  
  #Linea con pendiente -1 del 1.5 de la distancia de la interseccion mas lejana al origen
  xab = sqrt((distmax[1]*1.5)^2/2)
  yab = xab
  xord = xab + yab
  b1 = -1
  
  plot(-1,-1, xlim = c(0, xord), ylim = c(0, xord), xlab = "R1", ylab = "R2")
  for(i in 1:nsp){
    lines(c(Ra[i,1], Ra[i,1]), c(Ra[i,2], 1000000), col=i)
    lines(c(Ra[i,1], 1000000), c(Ra[i,2], Ra[i,2]), col=i)
  }
  mP = cbind(Ra[2:nsp,1], Ra[1:(nsp-1),2])
  mM = mY[,1]*(mK[,1]+1)/mY[,2]-mK[,2]
  
  points(Ra[1,1], Ra[1,2])   #L1 vertex
  points(Ra[2,1], Ra[2,2])   #L2 vertex
  points(Ra[2,1], Ra[1,2])   #Intersection vertex
  
  abline(a = xord, b = b1, col = 99)
  
  #Vectores de consumo
  for(i in 1:(nsp-1)){
    inter1 = mP[i,2]-mM[i]*mP[i,1]
    inter2 = mP[i,2]-mM[i+1]*mP[i,1]
    lines(c(0,1000000), c(inter1, inter1+mM[i]*1000000), col = i)
    lines(c(0,1000000), c(inter2, inter2+mM[i+1]*1000000), col = i+1)
  }
  #Interseccion de vectores de consumo y puntos de suministro
  xC1 = ((mM[1]*Ra[2,1])-Ra[1,2] - (b1*xab)+yab)/(mM[1] - b1)
  xC2 = ((mM[2]*Ra[2,1])-Ra[1,2] - (b1*xab)+yab)/(mM[2] - b1)
  yC1 = b1*(xC1 - xab) + yab
  yC2 = b1*(xC2 - xab) + yab
  points(xC1, yC1)
  points(xC2, yC2)
  
  #Region de coexistencia para el muestro (delimitada por ZNGI y vectores de consumo)
  lines(c(xC1, xC2), c(yC1, yC2))
  list("vcons" = mM, "pinter" = mP, "xord" = xord, "xC1" = xC1, "yC1" = yC1, "xC2" = xC2, "yC2" = yC2)
}


muestra = function(xord, intens=99){
  #Genera una muestra de puntos de la region de coexistencia 
  #Recibe:
  ### xord: coordenadas de las combinaciones de concentraciones de recursos dentro de la region de coexistencia
  ### intens: intensidad del muestreo (no. de puntos de abastecimiento sobre la linea de suministro)
  cbind(seq(0, xord, xord/intens), seq(xord, 0, -xord/intens))
}


minsoft = function(xx, yy, rho = 10){
  #Permite obtener el minimo para la resolucion de las ecuaciones diferenciales
  -log((exp(-rho*xx) + exp(-rho*yy)))/rho
}

convpar = function(par){
  #Nombra los parametros
  #Recibe:
  ### par: lista de parametros generados en la funcion hazparam()
  list("r1" = par$vr[1],
       "r2" = par$vr[2],
       "r3" = par$vr[3],
       "K11" = par$mK[1,1],
       "K21" = par$mK[2,1],
       "K31" = par$mK[3,1],
       "K12" = par$mK[1,2],
       "K22" = par$mK[2,2],
       "K32" = par$mK[3,2],
       "Y11" = par$mY[1,1],
       "Y21" = par$mY[2,1],
       "Y31" = par$mY[3,1],
       "Y12" = par$mY[1,2],
       "Y22" = par$mY[2,2],
       "Y32" = par$mY[3,2],
       "D" = par$vD[1],
       "Fr" = par$vF[1])
}


Tilman = function(t, state, parameters){
  #Ecuaciones del modelo de Tilman para tres especies y dos recursos
  #Recibe:
  ### t: tiempo para la resolucion de las ecuaciones
  ### state: estados iniciales de N1, N2, N3, R1 y R2
  ### parameters: valores generados para los parametros del modelo
  #Regresa:
  ### Lista con las dinamicas de las tres especies y los dos recursos
  with(as.list(c(state,parameters)),{
    dN1 = (r1*minsoft(R1/(K11+R1),R2/(K12+R2))-D)*N1                           #Dinamica de la especie 1
    dN2 = (r2*minsoft(R1/(K21+R1),R2/(K22+R2))-D)*N2                           #Dinamica de la especie 2
    dN3 = (r3*minsoft(R1/(K31+R1),R2/(K32+R2))-D)*N3                           #Dinamica de la especie 3
    dR1 = Fr*(A1-R1) - N1*r1*R1/((K11+R1)*Y11) - N2*r2*R1/((K21+R1)*Y21) - N3*r3*R1/((K31+R1)*Y31)     #Dinamica del recurso 1
    dR2 = Fr*(A2-R2) - N1*r1*R1/((K12+R2)*Y12) - N2*r2*R2/((K22+R2)*Y22) - N3*r3*R2/((K32+R2)*Y32)     #Dinamica del recurso 2
    list(c(dN1, dN2, dN3, dR1, dR2))
  })
}

hazdat = function(){
  #Une todas las funciones anteriores para obtener los valores resultantes de la resolucion de las ecuaciones en el equilibrio
  #Regresa:
  ### par: lista de valores de parametros
  ### npar: grafica de las ZNGI y vectores de consumo de las tres especies con los dos recursos, asi como la región de coexistencia
  ### sal: dinamicas de las tres especies y los dos recursos
  par(mfrow = c(3,3))
  sal = matrix(ncol=6, nrow=500)
  par = hazparam(3)
  npar = plotTilman(par$vr, par$vD, par$mK, par$mY)
  mues = muestra(npar$xord, 499)
  points(mues)
  for(i in 1:500){
    state = c(N1 = 0.1, N2 = 0.1, N3 = 0.1, R1 = mues[i,1], R2 = mues[i,2])
    parameters = c(convpar(par), A1=mues[i,1], A2=mues[i,2])
    times = seq(0, 500, 0.5)
    mod = ode(y = state, times = times, func = Tilman, parms = parameters)
    sal[i,] = mod[dim(mod)[1],]
  }
  
  plot(log(sal[,6]/sal[,5]), log(sal[,3]/sal[,2]), xlim=c(-10,10))
  plot(log(sal[,6]/sal[,5]), log(sal[,4]/sal[,2]), xlim=c(-10,10))
  plot(sal[,5], sal[,2])
  plot(sal[,5], sal[,3])
  plot(sal[,5], sal[,4])
  plot(sal[,6], sal[,2])
  plot(sal[,6], sal[,3])
  plot(sal[,6], sal[,4])
  list(par, npar, sal)
  
}


Ra_est = function(comp, umbral){
  #Hace la estimación de los valores de R*
  #Recibe:
  ### comp (list): lista con par, npar y sal
  ### umbral (int): valores umbrales para la estimacion
  #Regresa:
  ### Ras_mat (matrix: matriz con los R* estimados para cada par de especie-recurso con cada uno de los umbrales 
  Ras_mat = matrix(nrow = 3, ncol = 2)
  sal1 = comp[[3]]
  Ras_mat[1,1] = min(sal1[which(sal1[,2]>max(sal1[,2]*umbral)),5])
  Ras_mat[2,1] = min(sal1[which(sal1[,3]>max(sal1[,3]*umbral)),5])
  Ras_mat[3,1] = min(sal1[which(sal1[,4]>max(sal1[,4]*umbral)),5])
  Ras_mat[1,2] = min(sal1[which(sal1[,2]>max(sal1[,2]*umbral)),6])
  Ras_mat[2,2] = min(sal1[which(sal1[,3]>max(sal1[,3]*umbral)),6])
  Ras_mat[3,2] = min(sal1[which(sal1[,4]>max(sal1[,4]*umbral)),6])
  
  Ras_mat
}



explore = function(){
  #Ejecuta todas las funciones anteriores para simular el escenariio de una comunidad artificial por corrida
  ### Calcula error de estimacion respecto a los R* simulados y la correlacion entre estos valores con diferentes umbrales
  #Regresa:
  ### Lista con:
      # dens_max (vector): densidad maxima de cada especie
      # est_Ra (array): R* estimados utilizando diferentes valores umbrales 
      # obs_Ra (array): R* simulados para cada par de especie-recurso 
      # error_RAS (matrix): errores asociados a cada valor umbral
      # cor_coefs (vector): coeficientes de correlacion entre los R*
  comp = hazdat()
  dynamics = comp[[3]]
  dens_max = apply(dynamics, 2, max)[2:4]
  real_Ra = comp[[1]]$mRa
  
  umbrales = c(0.025, 0.01, 0.001, 0.00001)
  est_RAS = array(dim = c(3,2,4))             #nsp, nrec, no.umbrales + C2
  for (i in 1:4){
    est_RAS[,,i] = Ra_est(comp, umbrales[i])
  }
  
  
  error_RAS = array(dim = c(3,2,4))           #nsp, nrec, no.umbrales + C2
  for (i in 1:4){
    error_RAS[,,i] = (est_RAS[,,i]-real_Ra)/real_Ra
  }
  
  
  cor_coefs = c(cor(real_Ra)[1,2], cor(est_RAS[,,1])[1,2], cor(est_RAS[,,2])[1,2],
                cor(est_RAS[,,3])[1,2], cor(est_RAS[,,4])[1,2])
  
  list("dens_max" = dens_max, "real_Ra" = real_Ra, "est_RAS" = est_RAS, 
       "error_RAS" = error_RAS, "cor_coefs" = cor_coefs)
}


#Exploracion en 100 escenarios distintos
  
try_list = vector(mode = "list", length = 100)
for (i in 1:100){
  try_list[[i]] = explore()
}

saveRDS(try_list, "/Users/marco/Documents/Tesis/Databases/Try list.rds")
