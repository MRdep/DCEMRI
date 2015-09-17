;DCEMRIprocessor is a tool to fit tracer kinetic models to dynamic contrast enhanced MRI data. 

function DCEMRIprocessor_freeS0::Generate_maps, $
  fit_type,$ ; string. sequential:'ifit', joint:'jfit'
  data_type,$; string. 'phantom' or 'real'
  sa = sa,$; 1D array. arterial signal
  ca = ca,$; ID array. arterial contrast agent concentration
  update = update,$; logical. update reference dynamics and ref vfa (used for registration) after each Monte Carlo iteration. 
  n_vox_fit = n_vox_fit,$; integer. Number of voxels to sample during Monte Carlo experiements
  ostfORf = ostfORf,$; string. 'Free' or 'Fixed'. Should the ost be free or fixed during fitting. If free it is bound by +-5 seconds of the Tofts model estimate. 
  FAfORf = FAfORf,$; string. 'Free' or 'Fixed'. Should the flip angle error be estimated. Recommend to set as fixed. 
  aif_cap = aif_cap,$; integer. A user defined estimate of dynamic image in which enhancement occurs. At all timepoints before this ca is set to 0. 
  Hmt = Hmt,$; float. An estimate of the patient hematocrit. Used to convert sa to ca. 
  model = model,$; string. The tracer kinetic model to fit. The options are 'Tofts', 'ETofts' and 'XM' (two compartment exchange model)
  known_T1 = known_T1,$; float. A user defined estimate of precontrast T1. Useful if T1 mapping data does not exist or if thought to be unreliable. 
  known_S0v = known_S0v,$;float. Ditto for S0v.
  T1_dataType = T1_dataType,$; string. Type of T1 mapping data. At the moment only 'vfa'.
  baseline_cap = baseline_cap; integer. Similar to aif_cap. defines the last dynamic image up to which calculate the precontrast dynamic signal. 
  
  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON TwoCXM_params,E,K_m,K_p
  COMMON time, timepoints
  COMMON indices, index, xyz_indices
  COMMON signal_result, Res
  COMMON ost, T10_ost
  COMMON ve_vp_constraints, constrain
  COMMON iteration, j
  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
  COMMON sigma_info, fixed_sigma, equal_sigma, known_sigma ; fixed sigma will fixed sigmas to inputs from 'known sigma' or scan for sigma. equal sigma sets sigma_vfa = sigma_dyn.
  COMMON IRparams, little_tr, long_TR, rad_IR, ETL,n_centre
  
  ;defines how S0v and S0d should be fit
  equal_and_free = 1
  unequal_and_free = 0
  equal_by_baseline = 0
  
  COMPILE_OPT idl2
  
  device, decompose = 0 ; for linecolors to work, color decomposition needs to be set to zero.
  linecolors

  ;---------------------------------------- Data read-in from self.
  ; read in dynamic data
  ROI_timeseries = (*self.DCE_timeseries_p)*(*self.ROI_mask_p)/100.0 ; read in masked dce-time series
  TR = self.scan_params.TR
  alpha_dyn = self.scan_params.alpha_dyn
  rad_dyn = alpha_dyn*!dtor
  
  row= self.scan_params.row
  col= self.scan_params.column
  n_slices = self.scan_params.n_slices
  n_timepoints = self.scan_params.n_timepoints
  timepoints = *self.scan_params.timepoints_p
      
  ; contrast agent and arterial blood properties
  r = self.scan_params.relaxivity
  ;T10_a = 1480.0/1000; See Zhang Mag Res Med 70:1082-1086
  T10_a = 1200.0/1000 
  
  ;read in variable flip angle data
  if T1_datatype eq 'vfa' then begin
  vfa_images = *self.vfa_images_p; read in vfa images
  alpha_vfa = *self.scan_params.alpha_vfa_p
  rad_vfa = alpha_vfa*!dtor
  n_vfa = n_elements(rad_vfa)
  ref_vfa_images_LS = fltarr(col,row,n_slices,n_vfa)
  ref_vfa_images_ML = fltarr(col,row,n_slices,n_vfa)
  ref_vfa = fltarr(col,row,n_slices,n_vfa)
  endif
  
  noisefree_images = *self.noisefree_images_p
    
  ;----------------- Fitting parameters
  
  n_startpoints = 100 ; number of startpoints for fitting
  
  ;---------------------------------------- Defintion of reference storage arrays
  
  ref_timeseries_LS = fltarr(col,row,n_slices,n_timepoints); storage arrays for least squares fitting (assuming gaussian noise)
  
  ;---------------------------------------- ifit and jfit specific operations
  ; for sequential estimation we generate T1 map prior to tracer kinetic parameters.
  ; for joint estimation we also need to do this so that we can fit the Tofts model to estimate ost.
  ;--------------------------------------- Get AIF
  
  IF n_elements(sa) GT 0 THEN AIF = self->Convert_s2c(t10_a,sa, Hmt = Hmt, aif_cap = aif_cap)
  
  IF n_elements(ca) GT 0 THEN AIF = ca
  
  if aif_cap gt 0 then begin ; cap aif if needed. This puts zeroes in the AID up to the point where injection is thought to occur.
  AIF[0:aif_cap-1] = fltarr(aif_cap)
  endif
  
  ;-------------------------------------- Defining arrays to aid processing
  ;Here we convert our 4D arrays into 2D list type arrays by extracting the relevant information
  
  ROI_random_vox = fltarr(col,row,n_slices); an matrix to hold signal time series for randomly selected voxels.
  
  ;generating a mask which selects n_vox_fit random voxels from all voxels 
  if n_elements(n_vox_fit) gt 0 then begin ; add on if we don't want to fit to all roi voxels
  
    n_voxels = n_vox_fit
    
    ROI_timeseries = (*self.DCE_timeseries_p)*(*self.ROI_mask_p)/100.0
    
    index_all = where(ROI_timeseries[*,*,*] gt 0.)

    seed_for_sample = 1.0
        
    index = sample_from_array(index_all,n_vox_fit,seed_for_sample,/no_replacement) ; creating an index array to select n_vox_fit random voxels from the roi voxels.
 
    si = sort(index); indices are in random order. Need to be sorted
    index = index[si]
    ROI_random_vox[index] = 100.0

  endif else begin
    ;or if we want to fit to all voxels we skip the previous step
    
    index = where(ROI_timeseries[*,*,*] gt 0, n_voxels) ; Finding the linear indices and roi voxels
    
  endelse
  
  ; defining number of datapoints and dummy arrays for processing 
  
  if T1_datatype eq 'vfa' then begin
    if fit_type eq 'ifit' then  n_x = n_timepoints
    if fit_type eq 'jfit' then n_x = n_timepoints + n_vfa
    noisefree_array = fltarr(n_timepoints + n_vfa, n_voxels)
    vfa_array = fltarr(n_vfa,n_voxels)
  endif
  
  sv_array = fltarr(n_timepoints,n_voxels); an array of the measured sig-time curves
  sv_mod =  fltarr(n_x,n_voxels,n_startpoints)
  T1_array = fltarr(n_voxels) ; an array of T1 values. One element for each voxel in the roi.
  xyz_indices = findgen(3,n_voxels); a 3 by n_voxels array of x,y,z indices to hold the location of roi voxels.
  chisq_arr = fltarr(n_voxels,n_startpoints) ; an array to hold the chisq value for each voxel and startpoint.
  ost_array = fltarr(n_voxels) ; an array to hold ost if we want to set it beforehand. i.e. for bootstrapping.
  ;------------------------------------------------ Converting image arrays into 2D arrays where each row represent signal data for a voxel.
  
  z = 0 ; counter variable
  
  if n_elements(n_vox_fit) eq 0 then begin
    ROI_filler = ROI_timeseries
  endif else begin
    ROI_filler = roi_random_vox
  endelse
  
  FOR k = 0,n_slices-1 DO BEGIN
    FOR j = 0, row-1 DO BEGIN
      FOR i = 0, col-1 DO BEGIN
      
        IF ROI_filler[i,j,k] GT 0 THEN BEGIN
       ; print,i,j,k
          sv_array[*,z] = ROI_timeseries[i,j,k,*] ; filling the sig-time array with sig time curves from all ROI voxels
          
          if T1_datatype eq 'vfa' then vfa_array[*,z] = vfa_images[i,j,k,*]
          
          ;noisefree_array[*,z] = noisefree_images[i,j,k,*]
          
          xyz_indices[*,z] = [i,j,k]
          
          if n_elements(ost_map_old) gt 0 then ost_array[z] = ost_map_old[i,j,k] ; fills ost_array with passed values if needed
          z++ ; increment z
          
        ENDIF
        
      ENDFOR
    ENDFOR
  ENDFOR
  
  ;------------------------------------------ Defining storage arrays to hold maps

  Fp_map_LS = fltarr(col,row,n_slices)
  FE_map_LS = fltarr(col,row,n_slices)
  ve_map_LS = fltarr(col,row,n_slices)
  vp_map_LS = fltarr(col,row,n_slices)
  T1_map_LS = fltarr(col,row,n_slices)
  S0v_map_LS = fltarr(col,row,n_slices)
  S0d_map_LS = fltarr(col,row,n_slices)
  ost_map_LS = fltarr(col,row,n_slices)
  Ktrans_map_LS = fltarr(col, row, n_slices)
  FAerror_map_LS = fltarr(col,row,n_slices)
  chisq_map_LS = fltarr(col,row,n_slices)

  
  ;-------------------------------------------- Begin model fitting

  FOR i = 0, n_voxels-1 DO BEGIN ; fit the tracer kinetic (ifit) or tracer kinetic and SPGR models (jfit) to each voxel within the roi one at a time
    print, '------------------------'
    
    sv = sv_array[*,i] ; pick out signal-time curve for voxel i

    if T1_datatype eq 'vfa' then sv_vfa = vfa_array[*,i]; pick out signal-flip angle curves for voxel i

    noisefree_sig = noisefree_array[*,i]
    
    ;-------------------Defining the precontrast mean and dynamic noise
    
    cutoff = baseline_cap
    
    baseline = mean(sv[0:cutoff]); fixed for synthetic data. Used in sequential estimation to calculate S0d. 
    noise = stddev(sv[0:cutoff])
    SNR = baseline/noise
        
    ;---------------if we have no T1 estimate for a patient we use an estimate of T1 and S0v calculated from average of other patients
    
    if n_elements(known_T1) gt 0 then begin ; section of code for patients with no T1 estimate
    
    if T1_dataType eq 'vfa' then begin
    
    output = {ref_vfa:ref_vfa_images_LS,T1:known_T1,S0v:known_S0v}
    T10_ost = output.T1
    
    endif
    
    endif else begin ; for all other patients do this.
      
      if T1_dataType eq 'vfa' then output = self->SPGR_Fit(sv_vfa) ; fit SPGR model to vfa curves. Outputs T1 and S0v estiamtes. T1 estimate is used by both sequential and joint (ost calc)

      T10_ost = output.T1 ; set this T10_v estimate to the variable used to calculate ost.    
      
    endelse
        
    S0v = output.S0v ; for ifit this is not changed in the proceding code. For jfit it is changed.
    
    if fit_type eq 'ifit' then begin
    
      T1_map_LS[index[i]] = output.T1; for SE, assign Sov and T1 estimates to arrays and update common block
      S0v_map_LS[index[i]] = output.S0v
      T10_v = output.T1
      
      if T1_dataType eq 'vfa' then begin
      ref_vfa = output.ref_vfa
      if n_elements(known_T1) gt 0 then LS_mod_vfa = fltarr(n_vfa) else LS_mod_vfa = ref_vfa
      endif

    endif
        
    ;Defining the ost to be fed as an initial param
    if n_elements(ost_map_old) gt 0 then ost = ost_array[i] else ost = self->Get_ost_by_model(timepoints,sv,AIF) ;if ost_map_old supplied then ost should be taken from there and be fixed later.
    ;--------------------------------------------------- Fit model using least squares while varying vp
    
    if model eq 'Tofts' then n_startpoints = 1
    FOR j = 0, n_startpoints-1 DO BEGIN; ; fit 100 times while incrementing vp
    
      ; model dependent starting parameters
      start_params = self->Set_startparamsLS(fit_type,ost)
      
      if T1_datatype eq 'vfa' then sv_T1 = sv_vfa
        
      LS_output = self->Run_LSfit(fit_type,data_type,start_params,timepoints,sv,sv_T1,ostfORf = ostfORf,FAfORf = FAfORf, model = model, T1_datatype = T1_datatype);  output is in the form output = {modelled_sig:A, result:B, parinfo:C]
      
      sv_mod[*,i,j] = LS_output.modelled_signal ; this is just dyn for ifit and dyn + vfa for jfit

      LS_params = LS_output.parameters; LS_params is an array with all fitted parameters.
      
     
      chisq = LS_output.chisquared

      n = n_elements(LS_params)
      
      ;----------------------------------------------------- Choosing the best fit and filling LS maps
      
      chisq_arr[i,j] = chisq; for each voxel (i) there is a chisq value for each vp increment (j). We want to choose the fit with the lowest chisq.
      ;print, [LS_params,S0d,S0v,T10_v], chisq_arr[i,j]
      
      gt_zero_index = where(chisq_arr[i,*] gt 0) ; finding the indexes of the non-zero elements
      
      if gt_zero_index[0] gt -1 then chisq_no_zeroes = chisq_arr[i,gt_zero_index] else chisq_no_zeroes = -1
      
      tester1 = (chisq_arr[i,j] le min(chisq_no_zeroes)) && (chisq_arr[i,j] gt 0);
      
      IF tester1 EQ 1 THEN BEGIN;  if we are on the first startpoint and chisq gt 0
      
        S0v_map_LS[index[i]] = LS_params[0]
        T1_map_LS[index[i]] = LS_params[1]
        Fp_map_LS[index[i]] = LS_params[2]
        vp_map_LS[index[i]] = LS_params[3]
        ve_map_LS[index[i]] = LS_params[4]
        FE_map_LS[index[i]] = LS_params[5]
        S0d_map_LS[index[i]] = LS_params[6]
        ost_map_LS[index[i]] = LS_params[7]
        Ktrans_map_LS[index[i]] = LS_params[8]
        FAerror_map_LS[index[i]] = LS_params[9]
                
      ENDIF
      
      
    ENDFOR ; end of loop over n_startpoints
    
  
    gt_zero_index1 = where(chisq_arr[i,*] gt 0) ; extracting the modelled curve from the best fit.
    min_chisq = min(chisq_arr[i,gt_zero_index1],index_min)
    LS_mod_sig = sv_mod[*,i,index_min]
    
    ;calculate the normalised chisq for all timepoints.
    chisq_map_LS[index[i]] = min_chisq/n_timepoints

    ref_timeseries_LS[xyz_indices[0,i],xyz_indices[1,i],xyz_indices[2,i],*] = LS_mod_sig[0:n_timepoints-1] ; filling the reference time series array for registration.
    
    ;------------------------------------------------- Plotting LS fits

    IF fit_type eq 'ifit' THEN BEGIN

      LS_mod_dyn = LS_mod_sig[0:n_timepoints-1]
      
      if T1_datatype eq 'vfa' then begin 
      ;ref_vfa_images_LS[xyz_indices[0,i],xyz_indices[1,i],xyz_indices[2,i],*] = LS_mod_vfa
      plot, [rad_vfa, rad_vfa[2] + timepoints],[sv_vfa,sv],psym = 4, xtitle = 'Time (minutes)', ytitle = 'Signal intensity'
      oplot, [rad_vfa, rad_vfa[2]+ timepoints], [LS_mod_vfa,LS_mod_dyn] ,color = 7
      wait, 0.00001
      endif

      xyouts, 400,220, 'S0v = ' + strtrim(string(S0v_map_LS[index[i]]),1) ,/Device,CHARSIZE = 2
      xyouts, 400,200, 'S0d = ' + strtrim(string(S0d_map_LS[index[i]]),1) ,/Device,CHARSIZE = 2
      xyouts, 400,180, 'T1 = ' + strtrim(string(T1_map_LS[index[i]]),1) + 's',/Device,CHARSIZE = 2
      if model eq 'XM' then xyouts, 400,160, 'Fp = ' + strtrim(string(Fp_map_LS[index[i]]),1) + ' ml/min/ml' ,/Device,CHARSIZE = 2
      if model eq 'ETofts' || model eq 'Tofts' then xyouts, 400,160, 'Ktrans = ' + strtrim(string(Ktrans_map_LS[index[i]]),1) + 'ml-1' ,/Device,CHARSIZE = 2
      xyouts, 400,140, 'vp = ' + strtrim(string(vp_map_LS[index[i]]),1) + ' ml/ml',/Device,CHARSIZE = 2
      xyouts, 400,120, 've = ' + strtrim(string(ve_map_LS[index[i]]),1) + ' ml/ml',/Device,CHARSIZE = 2
      if model eq 'XM' then xyouts, 400,100,'FE = ' + strtrim(string(FE_map_LS[index[i]]),1) + ' ml/min/ml',/Device,CHARSIZE = 2
      xyouts, 400,80,'ost = ' + strtrim(string(ost_map_LS[index[i]]),1) + ' min',/Device,CHARSIZE = 2
      
      if model eq 'XM' then print, 'LS params','   S0v   ',S0v_map_LS[index[i]],'   T1   ',T1_map_LS[index[i]],'   Fp   ',Fp_map_LS[index[i]], '   vp   ',vp_map_LS[index[i]],'   ve   ',ve_map_LS[index[i]],'   FE   ',FE_map_LS[index[i]],'   S0d   ',S0d_map_LS[index[i]],'ost',ost_map_LS[index[i]]
      if model eq 'UM' then print, 'LS params','   S0v   ',S0v_map_LS[index[i]],'   T1   ',T1_map_LS[index[i]],'   Fp   ',Fp_map_LS[index[i]], '   vp   ',vp_map_LS[index[i]],'   FE   ',FE_map_LS[index[i]],'   S0d   ',S0d_map_LS[index[i]],'ost',ost_map_LS[index[i]]
      if model eq 'Tofts' || model eq 'ETofts' then  print, 'LS params','   S0v   ',S0v_map_LS[index[i]],'   T1   ',T1_map_LS[index[i]],'   vp   ',vp_map_LS[index[i]],'   ve   ',ve_map_LS[index[i]],'   Ktrans   ',Ktrans_map_LS[index[i]],'   S0d   ',S0d_map_LS[index[i]],'ost',ost_map_LS[index[i]]
    ENDIF 
    
    IF fit_type eq 'jfit' THEN BEGIN
         
      LS_mod_dyn = LS_mod_sig[0:n_timepoints - 1]
           
      if T1_datatype eq 'vfa' then begin
      LS_mod_vfa = LS_mod_sig[n_timepoints:*]
      ref_vfa_images_LS[xyz_indices[0,i],xyz_indices[1,i],xyz_indices[2,i],*] = LS_mod_vfa
      plot, [rad_vfa, timepoints],[sv_vfa,sv],psym = 4, xtitle = 'Time (minutes)', ytitle = 'Signal intensity'
      oplot, [rad_vfa, timepoints], [LS_mod_sig[n_timepoints:*],LS_mod_sig[0:n_timepoints-1]],color = 7
      wait, 0.000001
      endif
      
      if T1_datatype eq 'SR' then begin
        LS_mod_SR = LS_mod_sig[n_timepoints:*]
        ref_SR_images_LS[xyz_indices[0,i],xyz_indices[1,i],xyz_indices[2,i],*] = LS_mod_SR
        plot, [Inv_times/10., Inv_times[4]/10. + timepoints[0:80]],[sv_SR,sv[0:80]],psym = 4, xtitle = 'Time (minutes)', ytitle = 'Signal intensity'
        oplot, [Inv_times/10., Inv_times[4]/10. + timepoints[0:80]], [LS_mod_SR,LS_mod_dyn[0:80]] ,color = 7  
        wait, 0.000001    
      endif

      xyouts, 400,220, 'S0v = ' + strtrim(string(S0v_map_LS[index[i]]),1) ,/Device,CHARSIZE = 2
      xyouts, 400,200, 'S0d = ' + strtrim(string(S0d_map_LS[index[i]]),1) ,/Device,CHARSIZE = 2
      xyouts, 400,180, 'T1 = ' + strtrim(string(T1_map_LS[index[i]]),1) + 's',/Device,CHARSIZE = 2
       if model eq 'XM' then xyouts, 400,160, 'Fp = ' + strtrim(string(Fp_map_LS[index[i]]),1) + ' ml/min/ml' ,/Device,CHARSIZE = 2
       if model eq 'Tofts' || model eq 'ETofts'  then xyouts, 400,160, 'Ktrans = ' + strtrim(string(Ktrans_map_LS[index[i]]),1) + ' ml-1' ,/Device,CHARSIZE = 2
      xyouts, 400,140, 'vp = ' + strtrim(string(vp_map_LS[index[i]]),1) + ' ml/ml',/Device,CHARSIZE = 2
      xyouts, 400,120, 've = ' + strtrim(string(ve_map_LS[index[i]]),1) + ' ml/ml',/Device,CHARSIZE = 2
       if model eq 'XM' then xyouts, 400,100,'FE = ' + strtrim(string(FE_map_LS[index[i]]),1) + ' ml/min/ml',/Device,CHARSIZE = 2
      xyouts, 400,80,'ost = ' + strtrim(string(ost_map_LS[index[i]]),1) + ' min',/Device,CHARSIZE = 2
      
      if model eq 'XM' then print, 'LS params','   S0v   ',S0v_map_LS[index[i]],'   T1   ',T1_map_LS[index[i]],'   Fp   ',Fp_map_LS[index[i]], '   vp   ',vp_map_LS[index[i]],'   ve   ',ve_map_LS[index[i]],'   FE   ',FE_map_LS[index[i]],'   S0d   ',S0d_map_LS[index[i]],'ost',ost_map_LS[index[i]], 'FAerror ',FAerror_map_LS[index[i]]
      if model eq 'Tofts' || model eq 'ETofts'  then  print, 'LS params','   S0v   ',S0v_map_LS[index[i]],'   T1   ',T1_map_LS[index[i]],'   Ktrans   ',Ktrans_map_LS[index[i]], '   vp   ',vp_map_LS[index[i]],'   ve   ',ve_map_LS[index[i]],'   S0d   ',S0d_map_LS[index[i]],'ost',ost_map_LS[index[i]], 'FAerror ',FAerror_map_LS[index[i]]
      
    ENDIF
    
endfor
    
end
   
function DCEMRIprocessor_freeS0::Set_startparamsLS,fit_type,ost

  COMMON iteration, j
  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
  ;start params for 2CXM
  
  Fp_int = 0.5 ; default 0.5
  vp_int = 0.0000001 + 0.01*j ; if this is set to zero then funny things happen
  ve_int = 0.2 ;default 0.2
  FE_int = 0.2 ; be careful if using 2CUM as we get discontinuities in IRF when Fp ~ FE, default 0.2
  S0d_int = 5000.0; default ; 5000
  if fit_type eq 'ifit' then S0d_int = baseline*((1- cos(rad_dyn)*exp(-TR/T10_v))/(sin(rad_dyn)*(1-exp(-TR/T10_v))))
  S0v_int = 5000.0
  T10_v_int = 1.0
  ost_int = ost ; this will be the estimate from fitting to Tofts model or passed using ost_map_old
  Ktrans = 0.2
  f = 1.

  start_params = [S0v_int,T10_v_int,Fp_int,vp_int,ve_int,FE_int,S0d_int,ost_int, Ktrans,f]
  
  return, start_params
  
end
;----------------------------------------------------------------------------------------------------------------------------------------------
function DCEMRIprocessor_freeS0::Set_LScontraints, model, fit_type, data_type, start_params, ostfORf = ostfORf, FAfORf = FAfORf

  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
  
  
  result = start_params; S0v,T10,Fp,vp,ve,FE,S0d,ost, Ktrans, f

  parinfo = replicate({limited:[0,0],limits:[0.D,0],fixed:0, tied: ''},n_elements(result)) ; parinfo is an array of structures
  
  parinfo[0].limited[0] = 1    ;set lower limit for S0v
  parinfo[0].limits[0] = 0.D
  parinfo[0].limited[1] = 1    ;set upper limit for S0v
  parinfo[0].limits[1] = 80000.D ; usually 20000 for simulated and cervix data

  if equal_by_baseline eq 1 then parinfo[0].fixed = 1
  
  parinfo[1].limited[0] = 1    ;set lower limit for T10
  parinfo[1].limits[0] = 0.D
  parinfo[1].limited[1] = 1    ;set lower limit for T10
  parinfo[1].limits[1] = 10.D
    
  parinfo[2].limited[0] = 1   ;set lower limit for Fp
  parinfo[2].limits[0] = 0.D
  if model eq 'ETofts' || model eq 'Tofts' then parinfo[2].fixed = 1.
  
  parinfo[3].limited[0] = 1    ;set lower limit for vp
  parinfo[3].limits[0] = 0.D
  parinfo[3].limited[1] = 1    ;set upper limit for vp
  parinfo[3].limits[1] = 1.D
  parinfo[3].fixed = 1
  
  parinfo[4].limited[0] = 1    ;set lower limit for ve
  parinfo[4].limits[0] = 0.D
  parinfo[4].limited[1] = 1    ;set upper limit for ve
  parinfo[4].limits[1] = 1.D
  if model eq 'UM' then parinfo[4].fixed = 1.
  
  parinfo[5].limited[0] = 1    ;set lower limit for FE
  parinfo[5].limits[0] = 0.D
  
  if model eq 'ETofts' || model eq 'Tofts' then parinfo[5].fixed = 1.
  
  parinfo[6].limited[0] = 1    ;set lower limit for S0d
  parinfo[6].limits[0] = 0.D
  parinfo[6].limited[1] = 1    ;set lower limit for S0d
  parinfo[6].limits[1] = 80000.D
  
   parinfo[8].limited[0] = 1    ;set lower limit for Ktrans
  parinfo[8].limits[0] = 0.D
  parinfo[8].limited[1] = 1    ;set lower limit for Ktrans
  parinfo[8].limits[1] = 10.D
  
  if equal_by_baseline eq 1 then parinfo[6].fixed = 1
  if equal_and_free eq 1 then parinfo[0].tied = 'p[6]'
  
  ;fix parameter to reduce searching unnessasarily.
  if fit_type eq 'ifit' then parinfo[0].fixed = 1 ; this param is still passed in ifit.
  if fit_type eq 'ifit' then parinfo[1].fixed = 1 ; this param is still passed in ifit.
  if fit_type eq 'ifit' then parinfo[6].fixed = 1 ; this param is still passed in ifit.
  
  ;constraints for ost calculation.
  
  if ostfORf eq 'free' then begin
  
    temp_ost = start_params[7]
    
    parinfo[7].limited[0] = 1 ;set lower limit for ost
    
    ;temp_ost - 5.0/60 is the estimated ost - 5/60th of a minute. Is this a wide enough range?
    if temp_ost - 5.0/60 lt 0  then parinfo[7].limits[0] = 0 else parinfo[7].limits[0] = temp_ost - 5.0/60
    
    parinfo[7].limited[1] = 1;set upper limit for ost
    parinfo[7].limits[1] = temp_ost + 5.0/60
    
  endif else begin ;this will fix ost in fitting. The intial value will be set to tht calculated from fitting to Tofts or that passed by ost_old.
  
    parinfo[7].fixed = 1
    
  endelse
   
  if model eq 'XM' || model eq 'UM' then parinfo[8].fixed = 1
  
  ;flip angle correction
  if FAfORf eq 'free' then begin
  parinfo[9].limited[0] = 1    ;set lower limit for f
  parinfo[9].limits[0] = 0.5D
  parinfo[9].limited[1] = 1    ;set lower limit for f
  parinfo[9].limits[1] = 1.5D
  if fit_type eq 'ifit' then parinfo[9].fixed = 1
  endif else begin
  parinfo[9].fixed = 1
  endelse
    
  return, {parinfo:parinfo, result:result}
  
end



;----------Run_LSfit
;PURPOSE:
;Method to fit the 2CXM to dynamic timeseries data using sequential/joint estimation
;INPUTS:
;fit_type -> string, sequential estimation ('ifit') or joint estimation ('ifit')
;data_type -> string, tells the program the type of data it is analysing ('phantom' or 'real'). Has implications regarding offset time estimation.
;start_params -> array, of intial parameter estimates for fitting. start_params = [Fp,vp,ve,FE,T10_v,S0,ost]
;timepoints -> array, timepoint vector for the dynamic time series.
;sv -> array, dce signal vector for the voxel.
;rad_vfa -> array flip angles used for T1 mapping.
;sv_vfa -> array, vfa signal vector for the voxel.
;ostfORf -> string array, tells us allow ost to be a free parameter or fixedduring fitting.
;model = model -> string, '2CXM' or '2CUM'. Can easily add more.

function DCEMRIprocessor_freeS0::Run_LSfit,fit_type,data_type,start_params,timepoints,sv,sv_T1,ostfORf = ostfORf,FAfORf = FAfORf,model = model,T1_datatype = T1_datatype

  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON sigma_info, fixed_sigma, equal_sigma, known_sigma ; fixed sigma will fixed sigmas to inputs from 'known sigma' or scan for sigma. equal sigma sets sigma_vfa = sigma_dyn.
  
  LS_constraints = self->Set_LScontraints(model, fit_type, data_type, start_params, ostfORf = ostfORf, FAfORf = FAfORf)
  
  param_estimates = LS_constraints.result
    
  parinfo = LS_constraints.parinfo
    
  if fit_type eq 'ifit' then begin
   
    x_meas = timepoints
    y_meas = sv
    weights = findgen(n_elements(x_meas))
    weights[*] = 1.0
     
    if model eq 'UM' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='twoCUM_ifit_dc',parinfo=parinfo,/noderivative,/quiet)  
    if model eq 'XM' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='twoCXM_ifit_dc',parinfo=parinfo,/noderivative,/quiet)
    if model eq 'Tofts' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='Tofts_ifit_dc',parinfo=parinfo,/noderivative,/quiet)
    if model eq 'ETofts' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='ETofts_ifit_dc',parinfo=parinfo,/noderivative,/quiet)
  
    result = param_estimates ; may need to adjust for Tofts and ET
    
  endif
  
  if fit_type eq 'jfit' then begin
  
    if T1_datatype eq 'vfa' then begin
      y_meas= [sv,sv_T1]
      x_meas = [timepoints, rad_vfa]
    endif
    
    if T1_datatype eq 'SR' then begin
      y_meas= [sv,sv_T1]
      x_meas = [timepoints, Inv_times]
    endif
    
    if T1_datatype eq 'IR' then begin
      y_meas= [sv,sv_T1]
      x_meas = [timepoints, Inv_times]
    endif
        
    weights = findgen(n_elements(x_meas))
    
    if known_sigma.sigmaDYN eq known_sigma.sigmaVFA then begin ; logic
    
      weights[*] = 1.0
     ;weights[0:10] = 5.0
      
    endif else begin
    
      weights[0:n_timepoints-1] = 1.0/(known_sigma.sigmaDYN^2)
      weights[n_timepoints:*] = 1.0/(known_sigma.sigmaVFA^2) 

    endelse
    
    if model eq 'XM' then begin
    if T1_datatype eq 'vfa' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='twoCXM_vfa_jfit_dc',parinfo=parinfo, /noderivative, /quiet)
    if T1_datatype eq 'SR' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='twoCXM_SR_jfit_dc',parinfo=parinfo, /noderivative, /quiet)
    ; function to do IR jointly is not written yet. 
    ;if T1_datatype eq 'IR' then curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='twoCXM_IR_jfit_dc',parinfo=parinfo, /noderivative, /quiet)
    endif

    if model eq 'Tofts' then begin
    curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='Tofts_vfa_jfit_dc',parinfo=parinfo, /noderivative, /quiet)
    endif  
    
    if model eq 'ETofts' then begin
    curve = mpcurvefit(x_meas,y_meas,weights,param_estimates,chisq=chisq, function_name='ETofts_vfa_jfit_dc',parinfo=parinfo, /noderivative, /quiet)
    endif
    
    result = param_estimates
    
    ; we minimise the joint chisq, but select best vp using only the dynamic chisq?
    ;chisq = total((curve[0:n_timepoints-1] - y_meas[0:n_timepoints-1])^2) 
    
  endif
    
  output = {modelled_signal:curve,parameters:result,chisquared:chisq, parameter_info:parinfo}
  
  return, output
  
end

;-------------------------------------------------------------------------------------------------------------------
;Method to perform registration of the dynamic timeseries with a reference volume using translational registration.
;-------------------------------------------------------------------------------------------------------------------

function DCEMRIprocessor_freeS0::RIGID_reg_dynamics, objective_function, optimiser,reg_params_k = reg_params_k, ref_timeseries = ref_timeseries

  compile_opt idl2
  
  COMMON registration_params,ref_timepoint,dce_timepoint,reg_timepoint,roi,pre_vector, obj_fcn
  COMMON stepper, step_size,boundary
  
  Obj_fcn = objective_function
  
  dce_timeseries = *self.dce_timeseries_p ; Get dce_timeseries from storage
  ;----------------------------- import images and masks from storage
  
  if n_elements(ref_timeseries) eq 0 then begin ; if ref timeseries not supplied
    ref_timeseries = *self.ref_timeseries_p
  endif else begin
    ref_timeseries = ref_timeseries
  endelse
  
  s = size(ref_timeseries)
  
  roi_fourD = *self.roi_mask_p; get roi mask from storage
  
  roi = roi_fourD[*,*,*,0]/100.0
  
  n_timepoints = s[4]
  
  reg_timeseries = fltarr(s[1],s[2],s[3],s[4]) ; generate an array for the registered timeseries
  
  ;-----------------Generate a "pre vector" for each voxel in roi -------------------------
  
  roi_index = where(roi_fourD[*,*,*] gt 0, n_voxels) ; find the roi voxels
  
  pre_vector = fltarr(n_voxels,4) ; [xi,yi,zi,1]
  
  n_startpoints  = 1
  
  xf_i = fltarr(3,n_startpoints)
  xf_f = fltarr(3,n_startpoints)
  of_min_i = fltarr(n_startpoints)
  of_min_f = fltarr(n_startpoints)
  min_xf = fltarr(4,s[4])
  reg_params = fltarr(3,s[4])
  
  n = 0
  for i = 0, s[1] -1 do begin
    for j = 0, s[2] -1 do begin
      for k = 0, s[3] -1 do begin
      
        if ref_timeseries[i,j,k] gt 0 then begin
        
          pre_vector[n,0] = i
          pre_vector[n,1] = j
          pre_vector[n,2] = k
          pre_vector[n,3] = 1.0
          
          n++
          
        endif
        
      endfor
    endfor
  endfor
  
  ;exclude voxels that fail to fit properly (i.e. duds in ref timeseries)
  index_nan = where(finite(ref_timeseries,/NAN))

  ref_timeseries[index_nan] = 0.0
  
  
  for i = 0, n_timepoints - 1 do begin
  
    ref_timepoint = ref_timeseries[*,*,*,i] ;selects the ith timepoint of the reference vol (just roi voxels)
    dce_timepoint = DCE_timeseries[*,*,*,i] ; selects the ith timepoint of the shifted image for reg (full vol)
    for j = 0, n_startpoints -1 do begin
    
      reg_timepoint = fltarr(s[1],s[2],s[3])
      
      CASE optimiser OF
      
        'sa_then_powell': BEGIN
        
          ; start by finding an estimate of global minimum using simulated annealing
          print, 'timepoint ', i
          nstep = 1000
          step_size = 0.01 ; was 0.1 changed to 0.01 on 24/03/2014
          boundary = 5.0
          start_state = [0.0,0.0,0.0]
          best_params = sa(nstep = nstep, start_state = start_state)
          
          xf_i[*,j] = best_params.best_state_tab ; best reg params after sa
          of_min_i[j] = best_params.best_fitness_tab
          
          
          ; refine estimate with POWELL
          print, 'sa estimate', xf_i[*,j]
          
          dx = xf_i[0,j]
          dy = xf_i[1,j]
          dz = xf_i[2,j]
          
          X = [dx,dy,dz]
          
          ftol = 1.0e-15
          
          xi = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] ; intial direction vectors for search
          
          powell,X,xi,ftol,fmin,'register_fcn'
          
          xf_f[*,j] = X ; best reg params after sa and powell
          
          print, 'powell refinement', xf_f[*,j]
          
          of_min_f[j] = fmin
          
        END
        
        'sa': BEGIN
        
          ; start by finding an estimate of global minimum using simulated annealing
          print, 'timepoint ', i
          nstep = 1000
          step_size = 0.01
          boundary = 5.0
          start_state = [0.0,0.0,0.0]
          best_params = sa(nstep = nstep, start_state = start_state)
          
          xf_f[*,j] = best_params.best_state_tab ; best reg params after sa
          of_min_f[j] = best_params.best_fitness_tab
          print, 'sa estimates',xf_f
          
          
        END
        
        'powell': BEGIN
        
          ;------------just use powell
          
          dx = 0.0
          dy = 0.0
          dz = 0.0
          
          X = [dx,dy,dz]
          
          ftol = 1.0e-15
          
          xi = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] ; intial direction vectors for search
          
          powell,X,xi,ftol,fmin,'register_fcn'
          
          xf_f[*,j] = X
          of_min_f[j] = fmin
          print, 'powell estimate', xf_f
        END
        
        
      ENDCASE
      
    endfor ; end of loop over n_startpoints
    
    ; ------------------- find the transformation parameters with the minimum SAD (from all startpoints)
    
    minimum_of = min(of_min_f,min_index)
    
    min_xf[*,i] = [minimum_of,xf_f[0,min_index],xf_f[1,min_index],xf_f[2,min_index]]
    
    print, 'objective function minimum', minimum_of
    
    ;--------------------Calculating the transformation vectors -------------------------------
    reg_params[*,i] = min_xf[1:3,i]
    
    trans = [[1,0,0, reg_params[0,i]],[0,1,0,reg_params[1,i]],[0,0,1,reg_params[2,i]],[0,0,0,1]]
    
    pre_vector_allvox = fltarr(s[1]*s[2]*s[3],4) ; performing the transformation on the volume.
    
    c = 0
    
    for l = 0, s[1] -1 do begin
      for m = 0, s[2] -1 do begin
        for n = 0, s[3] -1 do begin
        
          pre_vector_allvox[c,0] = l
          pre_vector_allvox[c,1] = m
          pre_vector_allvox[c,2] = n
          pre_vector_allvox[c,3] = 1.0
          
          c++
          
        endfor
      endfor
    endfor
    
    post_vector = trans##pre_vector_allvox
    
    ;--------------------------- Preparing arrays for interpolate function
    
    interp_x = 2*pre_vector_allvox[*,0] - post_vector[*,0]
    interp_y = 2*pre_vector_allvox[*,1] - post_vector[*,1]
    interp_z = 2*pre_vector_allvox[*,2] - post_vector[*,2]
    
    interp = interpolate(dce_timepoint[*,*,*],interp_x,interp_y,interp_z,missing = 0.0)
    
    num = n_elements(interp)
    
    for k = 0, num-1 do begin
    
      reg_timepoint[pre_vector_allvox[k,0],pre_vector_allvox[k,1],pre_vector_allvox[k,2]] = interp[k]
      
    endfor
    
    reg_timeseries[*,*,*,i] = reg_timepoint
    
  endfor; end of loop over n_timepoints
  
 ; display_series, image_array = [dce_timeseries, ref_timeseries, reg_timeseries], matrix = [s[1],s[2]], ns = s[3], nt = s[4], scale = 1500
  
  *self.dce_timeseries_p = reg_timeseries ; set the new dce_timeseries to the registered timeseries
  ;
  reg_output = {reg_timeseries:reg_timeseries,ref_timeseries:ref_timeseries,reg_params:reg_params}
  ;
  return, reg_output
  ;
end

;------------------------------------------------------------------------------------------------
;Method to performs translational registration of fa images to the average dce timeseries volume.
;------------------------------------------------------------------------------------------------

function DCEMRIprocessor_freeS0::RIGID_reg_fa2dyn, objective_function, optimiser

  compile_opt idl2
  
  COMMON registration_params,ref_timepoint,dce_timepoint,reg_timepoint,roi,pre_vector, obj_fcn
  COMMON stepper, step_size,boundary
  
  
  Obj_fcn = objective_function
  
  ;----------------------------- import images and masks from storage
  
  dce_timeseries = *self.dce_timeseries_p ; Get dce_timeseries from storage
  vfa_images = *self.O_vfa_images_p ; Get fa images from storage
  
  roi_fourD = *self.roi_mask_p; get roi mask from storage
  roi = roi_fourD[*,*,*,0]; get roi mask from storage
  roi_timeseries = dce_timeseries*roi_fourD/100.0 ; 100 accounts for the fact that roi mask is between 0 and 100 not 0 and 1
  
  row = self.scan_params.row
  col = self.scan_params.column
  n_slices = self.scan_params.n_slices
  n_timepoints = self.scan_params.n_timepoints
  
  ref_vol = fltarr(col,row, n_slices)
  
  ;generating the ref vol which is the mean to all dynamic timepoints
  ;  for i = 0, n_timepoints-1 do begin
  ;
  ;    ref_vol[*,*,*] = roi_timeseries[*,*,*,i] + ref_vol[*,*,*]
  ;
  ;  endfor
  ;
  ;  ref_vol = ref_vol/n_timepoints
  
  ref_vol = (roi_timeseries[*,*,*,0]+roi_timeseries[*,*,*,1] + roi_timeseries[*,*,*,2] + roi_timeseries[*,*,*,3] + roi_timeseries[*,*,*,4]+ roi_timeseries[*,*,*,5]+ roi_timeseries[*,*,*,6]+ roi_timeseries[*,*,*,7])/8
  
  s = size(ref_vol)
  
  ;-----------------------------  make new image arrays
  
  n_vfa = 3
  
  reg_params = fltarr(3,n_vfa)
  
  reg_fa = fltarr(s[1],s[2],s[3],n_vfa) ; generate an array for the registered timeseries
  
  ; -----------------Generate a "pre vector" for each voxel in roi -------------------------
  
  roi_index = where(ref_vol gt 0, n_voxels)
  
  pre_vector = fltarr(n_voxels,4) ; [xi,yi,zi,1]
  
  n_startpoints  = 1
  
  xf_i = fltarr(3,n_startpoints)
  xf_f = fltarr(3,n_startpoints)
  of_min_i = fltarr(n_startpoints)
  of_min_f = fltarr(n_startpoints)
  min_xf = fltarr(4,n_vfa)
  
  
  n = 0
  for i = 0,s[1] -1 do begin
    for j = 0, s[2] -1 do begin
      for k = 0, s[3] -1 do begin
      
        ; if ref_timeseries[i,j,k,0] gt 0 then begin
        if ref_vol[i,j,k] gt 0 then begin
        
          pre_vector[n,0] = i
          pre_vector[n,1] = j
          pre_vector[n,2] = k
          pre_vector[n,3] = 1.0
          
          n++
          
        endif
        
      endfor
    endfor
  endfor
  
  ;  ;------------------------------ perform the optimisation
  
  print, 'number',n_elements(reg_params_k)
  
  
  for i = 0, n_vfa - 1 do begin
  
    dce_timepoint = vfa_images[*,*,*,i]; selecting an fa image to register
    
    for j = 0, n_startpoints -1 do begin
    
      reg_fa1 = fltarr(s[1],s[2],s[3]) ; an array to store a single registered fa image
      
      
      
      ;-------------------------renaming of volumes to work with common block
      
      ref_timepoint = ref_vol
      reg_timepoint = reg_fa1
      
      ;------------------------------------------------------------------------
      
      ; start by finding an estimate of global minimum using SIMULATED ANNEALING
      
      Case optimiser of
      
        'sa_then_powell': begin
        
          nstep = 1000
          start_state = [0.0,0.0,0.0]
          step_size = 0.1
          boundary = 10.0
          best_params = sa(nstep = nstep, start_state = start_state)
          
          
          xf_i[*,j] = best_params.best_state_tab
          of_min_i[j] = best_params.best_fitness_tab
          
          ;------------ refine estimate with POWELL and sinc integration
          print, xf_i[*,j]
          
          dx = xf_i[0,j]
          dy = xf_i[1,j]
          dz = xf_i[2,j]
          
          X = [dx,dy,dz]
          
          ftol = 1.0e-15
          
          xi = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] ; intial direction vectors for search
          
          powell,X,xi,ftol,fmin,'register_fcn'
          
          xf_f[*,j] = X
          of_min_f = fmin
          
          
        End
        
        
        'sa': begin
        
          nstep = 1000
          start_state = [0.0,0.0,0.0]
          step_size = 0.01
          boundary = 10.0
          best_params = sa(nstep = nstep, start_state = start_state)
          
          
          xf_f[*,j] = best_params.best_state_tab
          of_min_f[j] = best_params.best_fitness_tab
          
          
          
          
        End
        
        'powell': begin
        
          ;------------refine estimate with POWELL
          dx = 0.0
          dy = 0.0
          dz = 0.0
          
          X = [dx,dy,dz]
          
          ftol = 1.0e-15
          
          xi = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] ; intial direction vectors for search
          
          powell,X,xi,ftol,fmin,'register_fcn'
          
          xf_f[*,j] = X
          of_min_f = fmin
          
        End
        
        
      Endcase
      
    endfor ; end of loop over n_startpoints
    ;display_series, image_array = [filtered_ref_timepoint,filtered_roi_timepoint], matrix = [size_image[1], size_image[2]], ns = size_image[3],nt = 1
    ; calculating the transformed timeseries
    
    ; ------------------- find the transformation parameters with the minimum SAD (from all startpoints)
    
    ; print output
    minimum_of = min(of_min_f,min_index)
    
    min_xf[*,i] = [minimum_of,xf_f[0,min_index],xf_f[1,min_index],xf_f[2,min_index]]
    
    print, 'registration params', min_xf[1:3,i]
    print, 'objective function minimum', minimum_of
    ;--------------------Calculating the post transformation vectors -------------------------------
    reg_params[*,i] =  min_xf[1:3,i]
    
    trans = [[1,0,0,reg_params[0,i]],[0,1,0,reg_params[1,i]],[0,0,1,reg_params[2,i]],[0,0,0,1]]
    
    post_vector = trans##pre_vector
    
    ;---------------------------Preparing arrays for interpolate function
    
    interp_x = 2*pre_vector[*,0] - post_vector[*,0]
    interp_y = 2*pre_vector[*,1] - post_vector[*,1]
    interp_z = 2*pre_vector[*,2] - post_vector[*,2]
    
    interp = interpolate(dce_timepoint[*,*,*],interp_x,interp_y,interp_z,missing = -1)
    
    num = n_elements(interp)
    
    for k = 0, num-1 do begin
    
      if interp[k] gt 0 then begin
      
        reg_timepoint[pre_vector[k,0],pre_vector[k,1],pre_vector[k,2]] = interp[k]
        
      endif
      
    endfor
    
    reg_fa[*,*,*,i] = reg_timepoint
    
    
  endfor ; end of loop over n_timepoints
  
  ;set the newly created registered timeseries as fa_images
  *self.vfa_images_p = reg_fa
  
  reg_output = {reg_fa:reg_fa,reg_params:reg_params}
  
  return, reg_output
  ;display_series, image_array = [reg_fa3, vfa_images], matrix = [size_image[1], size_image[2]], ns = size_image[3],nt = 3, scale = 1500
  
end


;-----------------------------------------------------------
;Method to convert MR signal to contrast agent concentration
;-----------------------------------------------------------

function DCEMRIprocessor_freeS0::Convert_s2c,t10,sy,Hmt = Hmt, aif_cap = aif_cap

  ; Method to convert a signal->concentration
  
  ; Checked for accuracy on 12/06/2013. passed.
  
  compile_opt idl2
  
  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  
  E = exp(-TR/t10)

  Q = (sin(rad_dyn)*(1.0-E))/(1.0- cos(rad_dyn)*E) ; ; rad_dyn is the flip angle of the DCE-MRI time series not of the flip angle scans used for T1 mapping
  
  ;baseline = mean(sy[0:9]) ; take a mean of the precontrast datapoints. Maybe incorporate the BAT into this calc?
  
  S0_AIF = mean(sy[0:aif_cap-1])/Q; really dodgy. Need to vary from patient to patient 

  ;A1 = alog((sy*cos(rad_dyn)/S0_AIF*sin(rad_dyn)) + 1.0)
  
  A2 = alog((sy*cos(rad_dyn) - S0_AIF*sin(rad_dyn))/((sy - S0_AIF*sin(rad_dyn))))

  ;R1_1 = A1/TR
  ;cy_1 = (R1_1-1.0/t10)/r ; a array of concentration values
  
  ;R1_2 = A2/TR
  ;cy_2 = (R1_2-1.0/t10)/r ; a array of concentration values
  
;  cgplot,cy_1, color = 'red'
;  cgoplot,cy_2

  R1 = A2/TR
  cy = (R1-1.0/t10)/r ; a array of concentration values
  a = plot(cy)
  if n_elements(Hmt) gt 0 then begin
  
    cp = cy/(1-Hmt)
    
  endif
  a = plot(cp, /overplot)
  print, 'Finished conversion of AIF signal to AIF concentration'
  return, cp
  
end


;--------------------------------------------------------------------------------------
;Method to estimates the bolus arrival time by fitting the Tofts model to the concentration time curve
;and varying the arrival time as a free parameter
;------------------------------------------------------------------------------------------------------

function DCEMRIprocessor_freeS0::GET_ost_by_model,timepoints,sv,AIF

  ; computes the offset time between the AIF and tissue conc BATs by fitting the Tofts model with the offset time as a free param.
  ;COMMON signal_result, Res
  compile_opt idl2
  
  n_timepoints_ost = n_elements(timepoints)/3.
  sv_ost = sv[0:n_timepoints_ost-1]
  timepoints_ost = timepoints[0:n_timepoints_ost-1]
  
  
  Ktrans = 0.5 ; volume transfer constant
  ve_k = 0.25 ; EES volume
  ost = 5.0/60 ; offset time between the signal time curve and the AIF. In minutes.
 
  Result1 = [Ktrans,ve_k,ost]

  parinfo = replicate({limited:[0,0],limits:[0.D,0]},n_elements(Result1)) ; parinfo is an array of structures
  
  parinfo[0].limited[0] = 1   ;set lower limit for Ktrans
  parinfo[0].limits[0] = 0.D
  parinfo[1].limited[0] = 1    ;set lower limit for ve
  parinfo[1].limits[0] = 0.D
  parinfo[1].limited[1] = 1    ;set upper limit for ve
  parinfo[1].limits[1] = 1.D
  parinfo[2].limited[0] = 1    ;set lower limit for ost
  parinfo[2].limits[0] = 0.D
  parinfo[2].limited[1] = 1    ;set upper limit for ost
  parinfo[2].limits[1] = 25.0/60
  
  Result = Result1
  
  weights = fltarr(n_timepoints_ost)
  weights[*] = 1
  
  Res = mpcurvefit(timepoints_ost,sv_ost,weights,Result,chisq=chisq,status=status,function_name='Tofts_ost_expC',parinfo=parinfo,/noderivative,/quiet); LS fitting using ExpConvolution
 
  ost = Result[2]

  return, ost
  
end

;-----------------------------------------------------
;Method to compute a T1 map given a set of VFA images
;----------------------------------------------------

function DCEMRIprocessor_freeS0::SPGR_Fit,sv_vfa

  compile_opt idl2
  
  COMMON fit_params, TR,rad_dyn,rad_vfa,Inv_times,r,T10_v,S0d,S0v,baseline,n_timepoints,AIF,adj_AIF
  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
  
  weights = [1.0,1.0,1.0] ; setting the uncertainty on each measurement of signal to be the same
  Result1 = [5000.0,1.5] ; intialising the values [S0,T1]
  n_vfa = n_elements(rad_vfa)
  
  parinfo = replicate({limited:[0,0],limits:[0.D,0],value:0.D, fixed:0}, n_elements(Result1)) ; defining the variable parinfo
    
  parinfo[0].limited[0] = 1    ;set lower limit for S0_vfa
  parinfo[0].limits[0] = 0.D
  parinfo[0].limited[1] = 1    ;set upper limit for S0_vfa
  parinfo[0].limits[1] = 80000.D
  
  if equal_by_baseline eq 1 then begin
  parinfo[0].fixed = 1   
  endif
  
  parinfo[1].limited[0] = 1    ;set lower limit for T1
  parinfo[1].limits[0] = 0.D
  parinfo[1].limited[1] = 1    ;set upper limit for T1
  parinfo[1].limits[1] = 10.D ; usually 4.0 for simulated data
  
  S = sv_vfa
  
  P = Result1

  ref_vfa = MPCURVEFIT(rad_vfa,S,weights,P,Chisq = chisq, FUNCTION_NAME= 'FLASH_func',PARINFO = parinfo, /noderivative, /quiet) ; fit the function FLASH_func to the S vs alpha time datapoint
  
  S0v_temp = P[0]
  
  T1_temp = P[1]
  
  output = {T1:T1_temp, S0v:S0v_temp, ref_vfa:ref_vfa}
  
  return,output
  
end

pro DCEMRIprocessor_freeS0::Get_data,      $
  TR=TR,                                   $
  long_TR=long_TR,                                   $
  little_TR=little_TR,                                   $
  alpha_dyn = alpha_dyn,                   $
  alpha_vfa = alpha_vfa,                   $
  alpha_IR = alpha_IR,                   $
  alpha_SR = alpha_SR,                   $
  Inv_times = Inv_times, $
  row = row,                               $
  column = column,                         $
  n_timepoints = n_timepoints,             $
  timepoints = timepoints,                 $
  n_slices = n_slices,                     $
  patient_filename = patient_filename,     $
  DCE_timeseries = DCE_timeseries,         $
  O_DCE_timeseries = O_DCE_timeseries,     $
  reg_timeseries = reg_timeseries,         $
  vfa_images = vfa_images,                   $
  IR_images = IR_images,                   $
  SR_images = SR_images,                   $
  noisefree_images = noisefree_images,      $
  O_vfa_images = O_vfa_images,               $
  O_SR_images = O_SR_images,                $
  reg_vfa_images = reg_vfa_images,           $
  ROI_mask = ROI_mask,                     $
  T1_map = T1_map,                         $
  temp_res = temp_res,                     $
  relaxivity = relaxivity,                 $
  GT_maps = GT_maps,                       $
  temporal_resolution = temporal_resolution
  
  compile_opt idl2
  
  if arg_present(TR)                        then          self.scan_params.TR = TR
  if arg_present(long_TR)                   then          self.scan_params.long_TR = long_TR
  if arg_present(little_tr)                 then          self.scan_params.little_tr = little_tr
  if arg_present(Inv_times)                 then          self.scan_params.Inv_times_p = ptr_new(Inv_times) ; flip angle values for flipa gnle images.
  if arg_present(alpha_dyn)                 then          self.scan_params.alpha_dyn = alpha_dyn ; flip angle values for dynamics
  if arg_present(alpha_vfa)                 then          self.scan_params.alpha_vfa_p = ptr_new(alpha_vfa) ; flip angle values for flipa gnle images.
  if arg_present(alpha_IR)                  then          self.scan_params.alpha_IR = alpha_IR ; flip angle values for IR t1 measurements
  if arg_present(alpha_SR)                  then          self.scan_params.alpha_SR = alpha_SR ; flip angle values  for SR t1 measurements
  if arg_present(row)                       then          self.scan_params.row = row
  if arg_present(column)                    then          self.scan_params.column = column
  if arg_present(n_slices)                  then          self.scan_params.n_slices = n_slices
  if arg_present(n_timepoints)              then          self.scan_params.n_timepoints = n_timepoints
  if arg_present(relaxivity)                then          self.scan_params.relaxivity = relaxivity
  if arg_present(patient_filename)          then          self.patient_filename = patient_filename
  if arg_present(temporal_resolution)       then          self.scan_params.temporal_resolution = temporal_resolution
  if arg_present(timepoints)                then          self.scan_params.timepoints_p = ptr_new(timepoints)
  if arg_present(DCE_timeseries)            then          self.DCE_timeseries_p = ptr_new(DCE_timeseries) ; this assigns the pointer DCE_timeseries_p to point at DCEtimeseries
  if arg_present(O_DCE_timeseries)          then          self.O_DCE_timeseries_p = ptr_new(O_DCE_timeseries)
  if arg_present(reg_timeseries)            then          self.reg_timeseries_p = ptr_new(reg_timeseries)
  if arg_present(vfa_images)                then          self.vfa_images_p = ptr_new(vfa_images)
  if arg_present(IR_images)                 then          self.IR_images_p = ptr_new(IR_images)
  if arg_present(SR_images)                 then          self.SR_images_p = ptr_new(SR_images)
  if arg_present(O_vfa_images)              then          self.O_vfa_images_p = ptr_new(O_vfa_images)
  if arg_present(O_SR_images)               then          self.O_SR_images_p = ptr_new(O_SR_images)
  if arg_present(reg_vfa_images)            then          self.reg_vfa_images_p = ptr_new(reg_vfa_images)
  if arg_present(reg_SR_images)             then          self.reg_SR_images_p = ptr_new(reg_SR_images)
  if arg_present(ROI_mask)                  then          self.ROI_mask_p = ptr_new(ROI_mask)
  if arg_present(T1_map)                    then          self.T1_map_p = ptr_new(T1_map)
  if arg_present(noisefree_images)          then          self.noisefree_images_p = ptr_new(noisefree_images)
  
  
  ref_vfa_images_p = ptr_new(/allocate)
  ref_timeseries_p = ptr_new(/allocate)
  
  print, 'Finished data readin'
end

;------------------------------------------------------------------------------------------------------------------------
pro DCEMRIprocessor_freeS0::cleanup
  compile_opt idl2
  
  ;freeing memory used for pointers
  ; memory used for objetct is freed manually after the user has finished using the class
  
  if ptr_valid(self.DCE_timeseries_p) then  ptr_free, self.DCE_timeseries_p
  if ptr_valid(self.O_DCE_timeseries_p) then  ptr_free, self.O_DCE_timeseries_p
  if ptr_valid(self.O_vfa_images_p) then  ptr_free, self.O_vfa_images_p
  if ptr_valid(self.O_SR_images_p) then  ptr_free, self.O_SR_images_p
  if ptr_valid(self.ref_timeseries_p) then  ptr_free, self.ref_timeseries_p
  if ptr_valid(self.reg_timeseries_p) then  ptr_free, self.reg_timeseries_p
  if ptr_valid(self.ref_vfa_images_p) then  ptr_free, self.ref_vfa_images_p
   if ptr_valid(self.ref_IR_images_p) then  ptr_free, self.ref_IR_images_p
  if ptr_valid(self.ref_SR_images_p) then  ptr_free, self.ref_SR_images_p
  if ptr_valid(self.reg_vfa_images_p) then  ptr_free, self.reg_vfa_images_p
  if ptr_valid(self.reg_SR_images_p) then  ptr_free, self.reg_SR_images_p
  if ptr_valid(self.roi_mask_p) then  ptr_free, self.roi_mask_p
  ;if ptr_valid(self.roi_random_vox_mask_p) then  ptr_free, self.roi_random_vox_mask_p
  if ptr_valid(self.T1_map_p) then  ptr_free, self.T1_map_p
  if ptr_valid(self.vfa_images_p) then  ptr_free, self.vfa_images_p
  if ptr_valid(self.IR_images_p) then  ptr_free, self.IR_images_p
  if ptr_valid(self.SR_images_p) then  ptr_free, self.SR_images_p
  if ptr_valid(self.scan_params.timepoints_p) then  ptr_free, self.scan_params.timepoints_p
  if ptr_valid(self.noisefree_images_p) then  ptr_free, self.noisefree_images_p
  
  print,'cleanup complete'
  
end

function DCEMRIprocessor_freeS0::Init

  compile_opt idl2
  
  self.ref_vfa_images_p = ptr_new(/allocate)
  self.ref_IR_images_p = ptr_new(/allocate)
  self.ref_SR_images_p = ptr_new(/allocate)
  self.ref_timeseries_p = ptr_new(/allocate)
  self.reg_timeseries_p = ptr_new(/allocate)
  self.reg_vfa_images_p = ptr_new(/allocate)
  self.reg_SR_images_p = ptr_new(/allocate)  
  ;self.ROI_random_vox_mask_p = ptr_new(/allocate)
  self.Est_maps_p = ptr_new(/allocate)
  self.scan_params.alpha_vfa_p = ptr_new(/allocate)
  self.noisefree_images_p = ptr_new(/allocate)
  print,'init complete'
  return,1
  
end

;-------------------------------------------------------------------------------------------------------------------------

pro DCEMRIprocessor_freeS0_define

  compile_opt idl2
  
  
  ;  DCEMRIprocessor scan params structure
  ;
  void  = {Scan_parameters, $ ; scan parameters structure key {
    TR : 0.0,             $ ;   int ;
    long_TR : 0.0,             $ ;   int ;  
    little_tr : 0.0,             $ ;   int ;    
    Inv_times_p :ptr_new(), $ 
    TE  : 0.0,      $ ;   float
    alpha_dyn    : 0,      $ ;   float
    alpha_vfa_p:ptr_new(), $; float array
    alpha_IR    : 0,      $
    alpha_SR    : 0,      $
    row    : 0,           $ ;   int
    column : 0 ,       $ ; int
    relaxivity: 0.0,                $; float
    n_timepoints: 0,     $
    timepoints_p: ptr_new(), $
    temporal_resolution: 0.0, $
    n_slices: 0}; may need to change for different datasets
    
    
  ;DCEMRIprocessor class structure :
  
  class = {DCEMRIprocessor_freeS0,            $
    fit_type    : '',                                 $ ;  A string containing the specified fitting technique
    fit_method  : '',                                 $ ;  A string containing the specified fit method
    fit_domain  : '',                                 $ ;  A string containing the specified fit_domain
    scan_params: {scan_parameters},                           $ ; A stucture of scan parameters - filled by calling the IDL variables containing the image data
    DCE_timeseries_p: ptr_new() ,           $ ; saves memory for a pointer to reference DCE-MRI timeseiress
    O_DCE_timeseries_p: ptr_new(),        $ ; an array to hold the original DCE-MRI timeseries before any registration is performed
    O_SR_timeseries_p: ptr_new(),        $ ; an array to hold the original DCE-MRI timeseries before any registration is performed
    ref_timeseries_p: ptr_new(),             $ ; A blank array to store the modelled (reference) timeseries (for registration + bootstrapping purposes)
    reg_timeseries_p: ptr_new(),        $
    vfa_images_p : ptr_new() ,                 $ ; a ptr to the flip angle array [FA1, FA2, FA3] if T1 map has
    IR_images_p : ptr_new() ,                 $ ; a ptr to the IR array [IR1, IR2, IR3,IR4,IR5]
    SR_images_p : ptr_new() ,                 $ ; a ptr to the SR array [IR1, IR2, IR3,IR4,IR5]
    O_vfa_images_p:ptr_new(),                    $
    O_SR_images_p:ptr_new(),                    $
    ref_vfa_images_p: ptr_new(),               $; A blank array to store the modelled (reference) fa images (for bootstrapping purposes only)
    ref_IR_images_p: ptr_new(),               $; A blank array to store the modelled (reference) fa images (for bootstrapping purposes only)
    ref_SR_images_p: ptr_new(),               $; A blank array to store the modelled (reference) fa images (for bootstrapping purposes only)
    reg_vfa_images_p:ptr_new(),        $
    reg_SR_images_p:ptr_new(),        $
    ROI_mask_p: ptr_new() ,                   $ ; a ptr to the ROI mask
    noisefree_images_p: ptr_new() ,                   $ ; noisefree
    ;ROI_random_vox_mask_p: ptr_new() ,                   $
    T1_map_p: ptr_new(),                  $ ; a ptr to the T1 map if single estimation is used
    Est_maps_p: ptr_new(),                    $
    patient_filename : ''}
  ;
  
end

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

