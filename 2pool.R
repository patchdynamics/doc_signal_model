library(deSolve)

DOCin = function(t, subpool) {
  if(t < 80 || t > 120){
    return(8 / subpool)
  } else {
    return(20 / subpool)
  }
}

pool =  function(t, y, params) {
  DOC1 = y[1]
  DOCoutflow1 = y[2]
  DOC2 = y[3]
  DOCoutflow2 = y[4]
  
  with(as.list(params), {
    dDOC1 = DOCin(t, 1) / 2 - ( k1 * DOC1 ) - ( DOCoutflow1 ) 
    dDOC2 = DOCin(t, 2) / 2 - ( k2 * DOC2 ) - ( DOCoutflow2 ) 
    
    dDOCoutflow1 = Dout * dDOC1
    dDOCoutflow2 = Dout * dDOC2
    
    list(c(dDOC1, dDOCoutflow1, dDOC2, dDOCoutflow2))
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
params = c(k=c(.5,.1), DOCin = 8, Dout = .2)

Initial.DOC = 0 # units ?
Initial.DOCoutflow = Initial.DOC * params['Dout']
one.week = 24*7
two.weeks = 24*14
t = 1:two.weeks
t=1:200


state = c(Initial.DOC, Initial.DOCoutflow, Initial.DOC, Initial.DOCoutflow)
out = ode(y = state, times = t, func = pool, parms = params)

par(mfrow=c(2,2))
matplot(t, out[,2:3], type = "l", ylab = "p", xlab = "time")
matplot(t, out[,4:5], type = "l", ylab = "p", xlab = "time")

#par(mfrow=c(1,2))
inflows = c(DOCin(1, 1) / 2, DOCin(1,2) / 2)
outflows = c(out[,3][200], out[,5][200])
plot(inflows, ylim=c(0, 4), type="l")
plot(outflows, ylim=c(0, 4), type="l")
