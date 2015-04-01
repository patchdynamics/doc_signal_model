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
  portion = plnorm(subpool * segment) - plnorm((subpool - 1) * segment)
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
  
  BiofilmBacteriaFraction = 1
  MaximumFeasibleUptakeRate = MaximumPossibleUptake  # per hour
  UptakeRateThisTimeStep = 0 # summation
  
  with(as.list(params), {
    differentials = c()
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
      
      dDOC = DOCin(t, i, PoolDivisions) - UptakeRate *  BiofilmBacteriaFraction - DOCoutflow
      dDOCoutflow = Dout * dDOC
      differentials[index] = dDOC
      differentials[index+1] = dDOCoutflow

    }
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
PoolDivisions = 100
params = c(Dout = .2, PoolDivisions = PoolDivisions)

Initial.DOC = 0 # units ?
Initial.DOCoutflow = Initial.DOC * params['Dout']
one.week = 24*7
two.weeks = 24*14
t = 1:two.weeks
max = one.week
t = 1:max


state = rep(0, params['PoolDivisions'] * 2)
out = ode(y = state, times = t, func = pool, parms = params)

par(mfrow=c(2,2))
matplot(t, out[,2:3], type = "l", ylab = "p", xlab = 'hours')
matplot(t, out[,4:5], type = "l", ylab = "p", xlab = 'hours')
#matplot(t, sapply(1:PoolDivisions, function(subpool){out[,1+2*subpool]}))

#par(mfrow=c(1,2))
inflows = sapply(1:PoolDivisions, function(subpool){DOCin(1, subpool, PoolDivisions)})
outflows = sapply(1:PoolDivisions, function(subpool){out[,1+2*subpool][max]})
ylimit = max(inflows)
ks = -sapply(1:PoolDivisions, k, PoolDivisions);
plot(inflows~ks, type="l", ylim=c(0,ylimit), xlab='-k' )
plot(outflows~ks,type="l", ylim=c(0,ylimit), xlab='-k' )

#par(mfrow=c(1,1))
#matplot(t, out[,2:90], 2, type = "l")


#par(mfrow=c(2,2))
#t1 = sapply(1:PoolDivisions, function(subpool){out[1+subpool]})
