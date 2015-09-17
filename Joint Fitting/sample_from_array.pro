function sample_from_array,array,no_samples,seed_no,no_replacement = no_replacement

; a function to take random samples from an array with or without replacement. To sample without replacement add the /no_replacement keyword. 
;print,array,no_samples
if n_elements(no_replacement) eq 0 then begin

sample_array = fltarr(no_samples)

s = n_elements(array)

index_array = floor(s*randomu(seed_no,no_samples))

sample_array = array[index_array]

endif

if n_elements(no_replacement) gt 0 then begin
 
sample_array = fltarr(no_samples)
index_array = fltarr(no_samples)
  
s = n_elements(array)

index_array[0] = floor(s*randomu(seed_no,1))


for i = 1, no_samples -1 do begin

index_array[i] = floor(s*randomu(seed_no,1))

index_1 = where(index_array[i] eq index_array[0:i-1],n)

while n gt 0 do begin

  index_array[i] = floor(s*randomu(seed_no,1))
 
  index_1 = where(index_array[i] eq index_array[0:i-1],n)

endwhile

endfor

sample_array = array[index_array]

endif

return, sample_array
 
end
