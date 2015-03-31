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
  DOCoutflow = y[2]
  with(as.list(params), {
    dDOC = DOCin(t) - ( k * DOC ) - ( DOCoutflow ) 
    dDOCoutflow = Dout * dDOC
    list(c(dDOC, dDOCoutflow))
  })
}

#
# DOCin - flow of DOC into the system, kg
# Dout - percentage flow out of the system, i.e. water flow/residence time

params = c(k=.1, DOCin = 8, Dout = .2)
Initial.DOC = 100 # units ?
Initial.DOCoutflow = Initial.DOC * params['Dout']
one.week = 24*7
two.weeks = 24*14
t = 1:two.weeks

state = c(Initial.DOC, Initial.DOCoutflow)
out = ode(y = state, times = t, func = pool, parms = params)

matplot(t, out[,2:3], type = "l", ylab = "p", xlab = "time")