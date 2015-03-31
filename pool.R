library(deSolve)

DOCin = function(t) {
  if(t < 80 || t > 120){
    return(8)
  } else {
    return(20)
  }
}

pool =  function(t, y, params) {
  DOC = y[1]
  with(as.list(params), {
    dDOC = DOCin(t) - ( k * DOC ) - ( Dout * DOC ) 
    d2DOCout = Dout * dDOC  # want the 2nd deriv (flow change) b/c DOCout accumlates
    list(c(dDOC, d2DOCout))
  })
}

#
# DOCin - flow of DOC into the system, kg
# Dout - percentage flow out of the system, i.e. water flow/residence time

params = c(k=.1, DOCin = 8, Dout = .6)
Initial.DOC = 100 # units ?
Initial.DOCout = 0
one.week = 24*7
two.weeks = 24*14
t = 1:two.weeks

state = c(Initial.DOC, Initial.DOCout)
out = ode(y = state, times = t, func = pool, parms = params)

matplot(t, out[,2:3], type = "l", ylab = "p", xlab = "time")