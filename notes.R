sum(sapply(1:5, function(subpool){ DOCin(1, subpool, 5)  }))


out = sapply(1:5, function(subpool){ k(subpool, 5)  })

plot(plnorm(seq(0,1,.1), 1.7))

plot(dlnorm(seq(0,1.7*10,.1), 1.7) ~ seq(0,1.7*10,.1) )
plot(plnorm(seq(0,1.7*100,.1), 1.7) ~ seq(0,1.7*100,.1) )


plot(dlnorm(seq(0,1.7*100,.1), 1.7) ~ seq(0,1.7*100,.1) )


# use this to adjust mean of lognormal until it matches what they found in that paper
# understand where the K for semilabile ends on this scale

range = 100
out = sapply(1:range, function(subpool){ DOCin(1, subpool, range)  })
plot(out ~ lookup[0:100])
sum(out[97:100]) / sum(out)
sum(out[30:96]) / sum(out)

plot(dlnorm(seq(0,4.95,.05), mean= 1, sdlog=3) ~lookup )



sum(out[1:5]) / sum(out)
sum(out[6:42]) / sum(out)



sum(out[3:5]) / sum(out)
sum(out[6:50]) / sum(out)
# for instance, this would be coorrect of k = 5 is the upper limit on semilabile.


# so here is a distribution that makes sense at least
# 0.5, 1.9
# 1.9, 2.4
# .1 12
# 0.4 3
mean = 0.3
sdlog = 2.5
o = dlnorm(seq(0.0004,.8,.01), mean= mean, sdlog= sdlog)
plot(o~seq(0.0004,.8,.01), ylim=c(0,2) )
sum(o * .0079)

# and here is the shape like in the paper
plot(dlnorm(seq(0.0004,.8,.01), mean= mean, sdlog= sdlog)~log(seq(0.0004,.8,.01)) )


# now use that in the DOCin function
mini = 0.0004
maxi = .8
range = maxi - mini
DOCin2 = function(t, subpool, poolDivisions) {
  segment = range / poolDivisions
  # this use of lognormal is not quite right
  # remember than proportion in a given pool have ln(k) is normally distributed
  
  #if(t < 50) {
  #  mean = .5
  #} else {
  #  mean = 5
  #}
  #mean = 0.5 #+ .5 * sin(2 * pi * (t %% 120) / 120) # oscillate this
  #sdlog = 1.5
  
  
  # portion = plnorm(subpool * segment, mean = mean, sdlog = sdlog) - plnorm((subpool - 1) * segment, mean = mean, sdlog = sdlog)
  # trapezoid rule
  upper = (poolDivisions - subpool) * segment + mini
  lower = (poolDivisions - subpool -  1) * segment + mini
  portion = .5 * (dlnorm(upper, mean = mean, sdlog = sdlog) + dlnorm(lower, mean = mean, sdlog = sdlog) ) * segment
  
  #print(segment)
  #print(portion)
  
  DOCinTotal = 211  
  #DOCinTotal = 120  + 50 * sin(2 * pi * (t %% 120) / 120) # oscillate this
  
  PoolDOCin = DOCinTotal * portion
  return(PoolDOCin)
}


portion = function (subpool, poolDivisions){
  segment = range / poolDivisions
  print(segment)
  upper = subpool * segment + mini
  lower = (subpool - 1) * segment + mini
  portion = .5 * (dlnorm(upper, mean = mean, sdlog = sdlog) + dlnorm(lower, mean = mean, sdlog = sdlog) ) * segment
  return(portion)
}

poolDivisions = 100
portions = sapply(1:poolDivisions, function(subpool){ portion( subpool, poolDivisions)  })

out = sapply(1:poolDivisions, function(subpool){ DOCin2(1, subpool, poolDivisions)  })
lookup = seq(maxi, mini, length=100)
plot(out ~ lookup[0:100], ylim=c(0, 2))
sum(out[1:25]) / sum(out)
sum(out[25:58]) / sum(out)




