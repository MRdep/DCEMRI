PRO Tofts_ost_expC,timepoints, Result, data

  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON ost, T10_ost
 
  n_timepoints_ost = n_elements(timepoints) 
  
  Ct = FLTARR(n_timepoints_ost)
  data = FLTARR(n_timepoints_ost)
  
  Ktrans = result[0]
  ve = result[1]
  ost = result[2]
  
  timepoints_adjusted = timepoints-ost
     
  index = where(timepoints_adjusted lt 0) ; finds the timepoints that are negative
  
  timepoints_adjusted[index] = 0.0 ; sets the negative to zero

  adjusted_ca_for_ost = interpol(AIF[0:n_timepoints_ost-1],timepoints,timepoints_adjusted)

  X = [timepoints,adjusted_ca_for_ost]
  
  const = Ktrans/ve
     
  int1 = Expconvolution(const,X)

  Ct = Ktrans*int1 
 
  E1_0 = exp(-TR / T10_ost) ; no units
  R1_t = (1.0 / T10_ost) + (r * Ct) ; units of s-1
  T1_t = 1.0 / R1_t ;units of s
  E1_t = exp(-TR / T1_t) ; no units
  
  data = baseline * (1.0 - (cos(rad_dyn) * E1_0)) / (1.0 - E1_0) * (1.0 - E1_t) / (1.0 - (cos(rad_dyn) * E1_t))

END
