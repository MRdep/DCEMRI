Pro TwoCXM_ifit_dc,timepoints, result, data

  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON ve_vp_constraints, constrain
  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
   
  t_factor = mean(timepoints[1:n_timepoints-1] - timepoints[0:n_timepoints-2])

  Fp = result[2]
  vp = result[3]
  ve = result[4]
  FE = result[5]
  ost = result[7]
  f = result[9]
   
  if constrain eq 'yes' then begin
  
  if ve + vp lt 1.0 then begin
  
  if ost gt 0 then  begin   
     
    AIF_time = timepoints - ost
    
    adj_AIF = interpol(AIF,timepoints,AIF_time) ; adjusting AIF for difference in ost between AIF and voxel.
  
  endif else begin ; if data is synthetic and has no offset time we don't need to adjust AIF
  
    adj_AIF = AIF
    
  endelse
   
  Tp_i = (Fp+FE)/vp
  Te_i = FE/ve
  Tb_i = Fp/vp
  
  K_p = 0.5*(Tp_i + Te_i + sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))
  K_m = 0.5*(Tp_i + Te_i - sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))
  E = (K_p - Tb_i)/(K_p - K_m)
    
  a1 = exp(-K_m*timepoints)
  a2 = exp(-K_p*timepoints)

  int1 = convol(adj_AIF,a1,center = 0,/edge_zero)*t_factor
  int2 = convol(adj_AIF,a2,center = 0,/edge_zero)*t_factor
 
  Ct = Fp*(E*int1 + (1.0-E)*int2)
  
  R1_t = Ct*r + 1.0/T10_v
  
  E1 = exp(-TR*R1_t)

  S0d_temp = baseline*((1- cos(f*rad_dyn)*exp(-TR/T10_v))/(sin(f*rad_dyn)*(1-exp(-TR/T10_v))))
  
  result[0] = S0v
  result[1] = T10_v
  result[6] = S0d_temp
  
  data = S0d_temp*sin(f*rad_dyn)*(1-E1)/(1-E1*cos(f*rad_dyn)) ; S0d is calculated in LS_startparams
  
  
  endif else begin
    
  data = fltarr(n_timepoints)
  
  endelse
  
  endif else begin ; if unconstrained fitting
    
    if ost gt 0 then  begin
    
      AIF_time = timepoints - ost
      
      adj_AIF = interpol(AIF,timepoints,AIF_time) ; adjusting AIF for difference in ost between AIF and voxel.
      
    endif else begin ; if data is synthetic and has no offset time we don't need to adjust AIF
    
      adj_AIF = AIF
      
    endelse
    
    Tp_i = (Fp+FE)/vp
    Te_i = FE/ve
    Tb_i = Fp/vp
    
    K_p = 0.5*(Tp_i + Te_i + sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))
    K_m = 0.5*(Tp_i + Te_i - sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))
    E = (K_p - Tb_i)/(K_p - K_m)
    
    a1 = exp(-K_m*timepoints)
    a2 = exp(-K_p*timepoints)
    
    int1 = convol(adj_AIF,a1,center = 0,/edge_zero)*t_factor
    int2 = convol(adj_AIF,a2,center = 0,/edge_zero)*t_factor
    
    Ct = Fp*(E*int1 + (1.0-E)*int2)
    
    R1_t = Ct*r + 1.0/T10_v
    
    E1 = exp(-TR*R1_t)
    
    S0d_temp = baseline*((1- cos(f*rad_dyn)*exp(-TR/T10_v))/(sin(f*rad_dyn)*(1-exp(-TR/T10_v))))
    
    result[0] = S0v
    result[1] = T10_v
    result[6] = S0d_temp

    data = S0d_temp*sin(f*rad_dyn)*(1-E1)/(1-E1*cos(f*rad_dyn))  

  endelse
  
 END
