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
    dDOC = -k * DOC + DOCin(t)
    return(list(dDOC))
  })
}

params = c(k=.1, DOCin = 8)
Initial.DOC = 100 # units ?
one.week = 24*7
two.weeks = 24*14
t = 1:two.weeks

out = ode(y = Initial.DOC, times = t, func = pool, parms = params)

matplot(t, out[,2], type = "l", ylab = "p", xlab = "time")