function randomgauss,seed_temp,sigma,input

; Function to generate random deviate/s drawn from a gaussian distribution with parameter sigma. 
; sigma is the standard deviation of X and Y where R = sqrt(x^2 + Y^2), and R~R(M,s)
; 
; If we want to add noise to a previously noisefree signal such that the resulting signal is a sample from a rice distribution then:
; new = randomrice(old,s) 

s = size(input)

X = abs(sigma*randomn(seed_temp,s[1],s[2],s[3],s[4]) + input) ; adding gaussian noise to original data

return, X
end