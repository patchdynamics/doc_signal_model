library(deSolve)
range = 5 # range of the lognormal distribution to use
MaxTotalDOCUptake = 2.42 # g/h per km stream, the maximum rate derived from  McLaughlin and Kaplan 2013

DOCin = function(t, subpool, poolDivisions) {
  segment = range / poolDivisions
  # this use of lognormal is not quite right
  # remember than proportion in a given pool have ln(k) is normally distributed
  portion = plnorm(subpool * segment) - plnorm((subpool - 1) * segment)
  #print(segment)
  #print(portion)
  
  DOCinTotal = 317
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

  BiofilmDensity = 1;
    
  with(as.list(params), {
    differentials = c()
    for(i in 1:PoolDivisions){
      index = (i-1)*2 + 1
      DOC = y[index]
      DOCoutflow = y[index+1]
      dDOC = DOCin(t, i, PoolDivisions) - k(i, PoolDivisions) * DOC * BiofilmDensity - DOCoutflow
      dDOCoutflow = Dout * dDOC
      differentials[index] = dDOC
      differentials[index+1] = dDOCoutflow

    }
    print(differentials)
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
params = c(k=c(.5,.1), DOCin = 8, Dout = .2, PoolDivisions = PoolDivisions)

Initial.DOC = 0 # units ?
Initial.DOCoutflow = Initial.DOC * params['Dout']
one.week = 24*7
two.weeks = 24*14
t = 1:two.weeks
max = 40
t = 1:max


state = rep(0, params['PoolDivisions'] * 2)
out = ode(y = state, times = t, func = pool, parms = params)

par(mfrow=c(2,2))
matplot(t, out[,2:3], type = "l", ylab = "p", xlab = "time")
matplot(t, out[,4:5], type = "l", ylab = "p", xlab = "time")

#par(mfrow=c(1,2))
inflows = sapply(1:PoolDivisions, function(subpool){DOCin(1, subpool, PoolDivisions)})
outflows = sapply(1:PoolDivisions, function(subpool){out[,1+2*subpool][max]})
ylimit = max(inflows)
ks = -sapply(1:PoolDivisions, k, PoolDivisions);
plot(inflows~ks, type="l", ylim=c(0,ylimit) )
plot(outflows~ks,type="l", ylim=c(0,ylimit))

par(mfrow=c(1,1))
matplot(t, out[,2:90], 2, type = "l")


par(mfrom=c(2,2))
t1 = sapply(1:PoolDivisions, function(subpool){out[1+subpool]})
