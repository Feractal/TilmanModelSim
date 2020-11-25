library(deSolve)
library(easyNCDF)

#FUNCTIONS

hazparam = function(nsp){
  #Generates parameters for the resolution of DEs
  #Receive:
  # nsp: numbre of species
  #Return:
  # Parameter list
  D = runif(1)*2                             #Mortality rate
  r = D + runif(nsp)*2 + 1                   #Growth rate sp 1 (must be greater than D)
  Ra1 = sort(runif(nsp)*2)                   #Requeriment for resource 1
  Ra2 = sort(runif(nsp)*2,decreasing = T)    #Requeriment for resource 2
  K1 = Ra1*(r-D)/D                           #Half saturation constant for species i when limited by resource 1
  K2 = Ra2*(r-D)/D                           #Half saturation constant for species i when limited by resource 2
  m = sort(tan(runif(nsp)*(pi/2)))           
  Y2 = runif(nsp)*3                          #Individuals of sp i produced by unity of resource 2
  Y1 = ((K2+m)*Y2)/(K1+1)                    #Individuals of sp i produced by unity of resource 1
  list(vr = r, vD = rep(D,nsp), mK = cbind(K1,K2), mY = cbind(Y1,Y2), mRa = cbind(Ra1,Ra2))
}


plotTilman = function(vr,vD,mK,mY){
  #Generate graphs with ZNGI, consumption vectors and supply region
  #Receive:
  # vr: vector with growth rates
  # vD: vector with moratlity rates
  # mK: matrix with half saturation constants
  # mY: matrix with yield
  
  nsp = dim(mK)[1]
  mr = cbind(vr, vr)
  mD = cbind(vD, vD)
  Ra = mD*mK/(mr-mD)
  mr = mr[order(Ra[,1]),]
  mD = mD[order(Ra[,1]),]
  mY = mY[order(Ra[,1]),]
  mK = mK[order(Ra[,1]),]
  Ra = Ra[order(Ra[,1]),]
  
  #Distances between points and origin
  eles = sqrt(Ra[,1]^2 + Ra[,2]^2)
  inter = sqrt(Ra[2:nsp,1]^2 + Ra[1:(nsp-1),2]^2)
  distmax = c(eles, inter)[which.max(c(eles, inter))]
  
  #Supply gradient
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
  
  #Consumption vectors
  for(i in 1:(nsp-1)){
    inter1 = mP[i,2]-mM[i]*mP[i,1]
    inter2 = mP[i,2]-mM[i+1]*mP[i,1]
    lines(c(0,1000000), c(inter1, inter1+mM[i]*1000000), col = i)
    lines(c(0,1000000), c(inter2, inter2+mM[i+1]*1000000), col = i+1)
  }
  #Consumption vectors and supply points intersection
  xC1 = ((mM[1]*Ra[2,1])-Ra[1,2] - (b1*xab)+yab)/(mM[1] - b1)
  xC2 = ((mM[2]*Ra[2,1])-Ra[1,2] - (b1*xab)+yab)/(mM[2] - b1)
  yC1 = b1*(xC1 - xab) + yab
  yC2 = b1*(xC2 - xab) + yab
  points(xC1, yC1)
  points(xC2, yC2)
  
  #Supply region
  lines(c(xC1, xC2), c(yC1, yC2))
  list("vcons" = mM, "pinter" = mP, "xord" = xord, "xC1" = xC1, "yC1" = yC1, "xC2" = xC2, "yC2" = yC2)
}


muestra = function(xord, intens=99){
  #Generate sample from supply region
  #Receive:
  # xord: sample corrdinate
  # intens: sample intensity
  cbind(seq(0, xord, xord/intens), seq(xord, 0, -xord/intens))
}


minsoft = function(xx, yy, rho = 10){
  #soft minimum for DEs resolution
  -log((exp(-rho*xx) + exp(-rho*yy)))/rho
}

convpar = function(par){
  #Name the parameters
  #Receive:
  #par: list of parameters
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
       "D" = par$vD[1])
}


Tilman = function(t, state, parameters){
  #Set Tilman model equations for 3 species and 2 resources
  with(as.list(c(state,parameters)),{
    dN1 = (r1*minsoft(R1/(K11+R1),R2/(K12+R2))-D)*N1                           #Dynamics of sp 1
    dN2 = (r2*minsoft(R1/(K21+R1),R2/(K22+R2))-D)*N2                           #Dynamics of sp 2
    dN3 = (r3*minsoft(R1/(K31+R1),R2/(K32+R2))-D)*N3                           #Dynamics of sp 3
    dR1 = D*(A1-R1) - N1*r1*R1/((K11+R1)*Y11) - N2*r2*R1/((K21+R1)*Y21) - N3*r3*R1/((K31+R1)*Y31)     #Dynamics of resource 1
    dR2 = D*(A2-R2) - N1*r1*R1/((K12+R2)*Y12) - N2*r2*R2/((K22+R2)*Y22) - N3*r3*R2/((K32+R2)*Y32)     #Dynamics of resource 2
    list(c(dN1, dN2, dN3, dR1, dR2))
  })
}

hazdat = function(){
  #Join the functions to obtain parameters and dynamics for 3 species and 2 resources
  #Return:
  # par: parameters list
  # npar: graphs with ZNGI, consumption vectors and supply region
  # sal: dynamics of sp1, sp2, sp3, resource 1 and resource 2
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
  #Generates R* estimates under criterion 1 described in the methods section
  #Return:
  #Ras_mat: matrix with the R* estimated for each pair of species-resource
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
  #Explores the estimation errors with the R * values 
  #observed with criterion 1 and the correlations between these values 
  #with different thresholds.
  #Return:
      # List with:
      # dens_max: maximum densities of species
      # est_Ra: R* estimated using criterion 1
      # obs_Ra: "real" R* for each pair of species-resource
      # error_RAS: matrix with the error associated to each estimation of R* using different thresholds
      # cor_coefs: vector with correlation coefficients
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


#Explore in 100 scenarios
  
try_list = vector(mode = "list", length = 100)
for (i in 1:100){
  try_list[[i]] = explore()
}

saveRDS(try_list, "../Try list.rds")
