library(deSolve)
range = 5 # range of the lognormal distribution to use
MaxTotalDOCUptake = 2.42  # g/h per km stream, the maximum rate derived from  McLaughlin and Kaplan 2013
# this factor of 20 effectively extends the reach by to 20 km

BiofilmBacteriaFraction = .8; # for now set density of bacteria at 80% of their maximum possible
MaximumPossibleUptake = MaxTotalDOCUptake # same as above, at least for now

DOCin = function(t, subpool, poolDivisions) {
  segment = range / poolDivisions
  # this use of lognormal is not quite right
  # remember than proportion in a given pool have ln(k) is normally distributed
 
  if(t < 50) {
    mean = .5
  } else {
    mean = 1.7
  }
  mean = 1.7
  
  portion = plnorm(subpool * segment, mean = mean) - plnorm((subpool - 1) * segment, mean = mean)
  #print(segment)
  #print(portion)
  
  DOCinTotal = 211
  
  PoolDOCin = DOCinTotal * portion
  return(PoolDOCin)
}

k = function(subpool, poolDivisions){
  # this is also not quite right
  # because MaxTotalDOCUptake is not the max rate, but rather the max combined rate
  # of all those exponential decays
  lookup = seq(.1, .0034, length.out = poolDivisions + 1)
  k = lookup[subpool + 1]
  return(k)
}

pool =  function(t, y, params) {
  
  MaximumFeasibleUptakeRate = MaximumPossibleUptake  # per hour
  UptakeRateThisTimeStep = 0 # summation
  
  BacteriaDensity = y[PoolDivisions * 2 + 1]
  AlgaeDensity = y[PoolDivisions * 2 + 1]
  
  with(as.list(params), {
    differentials = c()
    TotalBacterialUptake = 0
    for(i in 1:PoolDivisions){
      index = (i-1)*2 + 1
      DOC = y[index]
      DOCoutflow = y[index+1]
      
      UptakeRate = k(i, PoolDivisions) * DOC ;
      if(UptakeRate + UptakeRateThisTimeStep > MaximumFeasibleUptakeRate){
        UptakeRate = MaximumFeasibleUptakeRate - UptakeRateThisTimeStep
      }
      UptakeRateThisTimeStep = UptakeRateThisTimeStep + UptakeRate
      
      #print(UptakeRateThisTimeStep)
      
      dDOC = DOCin(t, i, PoolDivisions) - UptakeRate * BacteriaDensity - DOCoutflow
      
      # Ad Hoc Algal contribution to DOC
      l = .05 #exhudation of labile DOC by algae, needs to be properly parameterized
      if( i < PoolDivisions / 10) { # lowest tenth are most labile
        dDOC = dDOC + (l * AlgaeDensity) / (PoolDivisions / 10)    # so we get a boost from algae here
      }
      
      dDOCoutflow = Dout * dDOC
      differentials[index] = dDOC
      differentials[index+1] = dDOCoutflow

      TotalBacterialUptake = TotalBacterialUptake + UptakeRate * BacteriaDensity
    }
    
    #
    # With the exception of K, these params are completely made up at this point
    # Needs proper parameterization
    #
    KBact = 10 # which is to say, we are allowing 10 kilometers of biofilm.
    BGE = .05 #  stating point: http://www.annualreviews.org/doi/abs/10.1146/annurev.ecolsys.29.1.503
    AlgaeInhibitsBact = .5 # not parameterized yet
    # bacteria need to grow based on DOC uptake
    dBacteriaDensity = BGE * TotalBacterialUptake * ( 1 - BacteriaDensity / KBact - AlgaeInhibitsBact * AlgaeDensity / KBact)
    
    KAlg = 10
    # these constants parameterized yet
    r = .1  # this is the photosynthesis production
    BactInhibitsAlgae = .5
    dAlgaeDensity = r * AlgaeDensity * ( 1 - AlgaeDensity / KAlg - BactInhibitsAlgae * BacteriaDensity / KAlg )
    
    differentials = c(differentials, dBacteriaDensity, dAlgaeDensity)
    
    
    #print(differentials)
    list(differentials)
  })
}

#
# DOCin - flow of DOC into the system, kg
# Dout - percentage flow out of the system, i.e. water flow/residence time

# 2 pools
# 2 different k
# assume the inflow to each pool is equal (this is not true)
# outflows are equal percentages (this will be true)
# also assume initial DOC is the same and 0 for both pools
PoolDivisions = 50
params = c(Dout = .2, PoolDivisions = PoolDivisions)

Initial.DOC = 0 # units ?
Initial.DOCoutflow = Initial.DOC * params['Dout']
Initial.BacteriaDensity = .1
Initial.AlgaeDensity = .2

one.week = 24*7
two.weeks = 24*14
max = 300
t = 1:max



state = rep(0, params['PoolDivisions'] * 2)
state = c(state, Initial.BacteriaDensity)
state = c(state, Initial.AlgaeDensity)
out = ode(y = state, times = t, func = pool, parms = params)

par(mfrow=c(2,3))
matplot(t, out[,2:3], type = "l", ylab = "DOC fraction 1", xlab = 'hours')
matplot(t, out[,params['PoolDivisions'] * 2+2], type = "l", ylab = "Bacteria Density", xlab = 'hours')
matplot(t, out[,params['PoolDivisions'] * 2+3], type = "l", ylab = "Algae Density", xlab = 'hours')
#matplot(t, sapply(1:PoolDivisions, function(subpool){out[,1+2*subpool]}))

#par(mfrow=c(1,2))
inflows = sapply(1:PoolDivisions, function(subpool){DOCin(1, subpool, PoolDivisions)})
outflows = sapply(1:PoolDivisions, function(subpool){out[,1+2*subpool][max-1]})
ylimit = max(inflows)
ks = -sapply(1:PoolDivisions, k, PoolDivisions);
plot(inflows~ks, type="l", ylim=c(0,ylimit), xlab='-k' )
plot(outflows~ks,type="l", ylim=c(0,ylimit), xlab='-k' )

#par(mfrow=c(1,1))
#matplot(t, out[,2:90], 2, type = "l")
#matplot(t, sapply(1:PoolDivisions, function(subpool){ sapply(t,function(t){ 
#  DOCin(t, subpool, PoolDivisions)})}))

#par(mfrow=c(2,2))
#t1 = sapply(1:PoolDivisions, function(subpool){out[1+subpool]})
