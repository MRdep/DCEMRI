pro Flash_func,FA2,p,data

 COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
 COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline

  T1 = p[1]
  E0 = exp(-TR/T1)
  S0 = p[0]

  data = S0*((sin(FA2)*(1-E0))/(1- cos(FA2)*E0)); FLASH EQUATION

end