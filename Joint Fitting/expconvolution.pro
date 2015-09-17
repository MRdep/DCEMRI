function ExpConvolution, l, X
 
  n = n_elements(X)/2
  
  T = X[0:n-1] ; T is the first half go the X array, i.e. the time points

  A = X[n:*]   ; A is the second half of the X array i.e. the concentration values at each time point ; units of mM

 
  DT = T[1:n-1]-T[0:n-2] ; This is an array of dt between time points
  DA = A[1:n-1]-A[0:n-2] ; this is an array of da between conc values
  
  Z = l*DT ;Z is an array of K*dt
  
  E = exp(-Z); E is an array of  exponentials [exp-lDT0, exp-lDT1 etc
  E0 = 1-E ; E0 is an array [(1-exp-lDT0), (1-exp-lDT1)
  E1 = Z-E0
  
  Il = (A[0:n-2]*E0 + DA*E1/Z)/l
  
  Y = fltarr(n) ; Y is a double precision array with every element set to 0.0
  
  for i=0L,n-2 do Y[i+1] = E[i]*Y[i] + Il[i]
  
  return, Y
  
end
