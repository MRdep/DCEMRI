pro twoCXM_vfa_jfit_dc, x_meas, result, data3

  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON TwoCXM_params,E,K_m,K_p
  COMMON ve_vp_constraints, constrain
  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
    
  timepoints = x_meas[0:n_timepoints-1]
  
  rad_vfa = x_meas[n_timepoints:*]
  
  t_factor = mean(timepoints[1:n_timepoints-1] - timepoints[0:n_timepoints-2])

  S0v_temp = result[0]
  T10_v_temp = result[1]
  Fp = result[2]
  vp = result[3]
  ve = result[4]
  FE = result[5]
  S0d_temp = result[6] ; S0v is passed later on if required. Updates common block.
  ost = result[7]
  f = result[9]
  
  result[8] = (Fp*FE)/(Fp+FE); setting the value for Ktrans
    
  if constrain eq 'yes' then begin ; reject parameters where ve + vp > 1
  
    if ve + vp lt 1.0 then begin

      
        AIF_time = timepoints - ost
        
        adj_AIF = interpol(AIF,timepoints,AIF_time) ; adjusting AIF for difference in ost between AIF and voxel.
     
      
      Tp_i = (Fp + FE)/vp
      
      Te_i = FE/ve
      
      Tb_i = Fp/vp
      
      K_p = (Tp_i + Te_i + sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))/2
      K_m = (Tp_i + Te_i - sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))/2
      E = (K_p - Tb_i)/(K_p-K_m)
      
      ;X = [timepoints, adj_ca]
      
      ;calculating the time_series curve
      ;int1 = ExpConvolution(K_m,X)
      ;int2 = ExpConvolution(K_p,X)
      
      a1 = exp(-K_m*timepoints)
      a2 = exp(-K_p*timepoints)
      
      int1 = convol(adj_AIF,a1,center = 0,/edge_zero)*t_factor
      int2 = convol(adj_AIF,a2,center = 0,/edge_zero)*t_factor
     
      Ct = Fp*(E*int1 + (1.0-E)*int2)
      
      R1_t = Ct*r + 1.0/T10_v_temp
      E1 = exp(-TR*R1_t)
      
;      print, rad_dyn
;      print, rad_vfa     
            
      data1 = S0d_temp*sin(f*rad_dyn)*(1-E1)/(1-E1*cos(f*rad_dyn))

      data2 = S0v_temp*sin(f*rad_vfa)*(1-exp(-TR/T10_v_temp))/(1- cos(f*rad_vfa)*exp(-TR/T10_v_temp))   
   
      data3 = [data1,data2]
    
    endif else begin
    
      data3 = fltarr(n_timepoints + 3)
      
    endelse
    
  endif else begin
  

    AIF_time = timepoints - ost
      
    adj_AIF = interpol(AIF,timepoints,AIF_time) ; adjusting AIF for difference in ost between AIF and voxel.

    Tp_i = (Fp + FE)/vp
    
    Te_i = FE/ve
    
    Tb_i = Fp/vp
    
    K_p = (Tp_i + Te_i + sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))/2
    K_m = (Tp_i + Te_i - sqrt(( Tp_i + Te_i)*( Tp_i + Te_i) - 4*Te_i*Tb_i))/2
    E = (K_p - Tb_i)/(K_p-K_m)
    
    ;X = [timepoints, adj_ca]
    
    ;calculating the time_series curve
    ;int1 = ExpConvolution(K_m,X)
    ;int2 = ExpConvolution(K_p,X)
    
    a1 = exp(-K_m*timepoints)
    a2 = exp(-K_p*timepoints)
    
    int1 = convol(adj_AIF,a1,center = 0,/edge_zero)*t_factor
    int2 = convol(adj_AIF,a2,center = 0,/edge_zero)*t_factor
    
    Ct = Fp*(E*int1 + (1.0-E)*int2)
    
    R1_t = Ct*r + 1.0/T10_v_temp
    E1 = exp(-TR*R1_t)
        
    data1 = S0d_temp*sin(f*rad_dyn)*(1-E1)/(1-E1*cos(f*rad_dyn))
    
    data2 = S0v_temp*sin(f*rad_vfa)*(1-exp(-TR/T10_v_temp))/(1- cos(f*rad_vfa)*exp(-TR/T10_v_temp))      
   
    data3 = [data1,data2]
    
  endelse
  
end