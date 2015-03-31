library(deSolve)

pool =  function(t, y, params) {
  DOC = y[1]
  with(as.list(params), {
    dDOC = -k * DOC + DOCin
    return(list(dDOC))
  })
}

params = c(k=.1, DOCin = 8)
Initial.DOC = 100 # units ?
t = 1:200

out = ode(y = Initial.DOC, times = t, func = pool, parms = params)

matplot(t, out[,2], type = "l", ylab = "p", xlab = "time")