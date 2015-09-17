function run_2CXM_phantom,$
  snr_dyn = snr_dyn,$
  snr_vfa = snr_vfa,$
  T1_datatype = T1_datatype,$
  model = model, $
  MultiThreadIndex= MultiThreadIndex,$
  n_repeats = n_repeats,$
  n_vox_fit = n_vox_fit

  !EXCEPT=0
  
  ;PURPOSE: wrapper for dcemriProcessor enabling Monte Carlo simulations.
  
  COMMON ve_vp_constraints, constrain
  COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
  COMMON sigma_info, fixed_sigma, equal_sigma, known_sigma
  
  if MultiThreadIndex eq 1 then begin
    fit_type = 'ifit'
    aif_type = '1'
  endif
  
  if MultiThreadIndex eq 2 then begin
    fit_type = 'jfit'
    aif_type = '1'
  endif
  
  if MultiThreadIndex eq 3 then begin
    fit_type = 'ifit'
    aif_type = '2'
  endif
  
  if MultiThreadIndex eq 4 then begin
    fit_type = 'jfit'
    aif_type = '2'
  endif
  
  if MultiThreadIndex eq 5 then begin
    fit_type = 'ifit'
    aif_type = '3'
  endif
  
  if MultiThreadIndex eq 6 then begin
    fit_type = 'jfit'
    aif_type = '3'
  endif
  
  constrain = 'yes'
  ostfORf = 'free';ostfORf is a setting that allows ost to vary in the actual model fit as a free parameter
  FAfORf = 'fixed';if free then fits for flip angle error
  mean_baseline = 163.0
  sigma_struct = {fixed_sigma:1,equal_sigma:0,known_sigma:{sigmaDYN:mean_baseline/5.,sigmaVFA:mean_baseline/11.2}}
  
  fixed_sigma = sigma_struct.fixed_sigma
  equal_sigma = sigma_struct.equal_sigma
  known_sigma = sigma_struct.known_sigma
  
  ; define imaging parameters
  TR =  3.2/1000
  TE = 1.2/1000
  alpha_dyn = 25.
  alpha_fa = [5.,10.,35.]
  temporal_resolution = 2.5
  r = 4.5
  
  ;------------------------------ Read in phantom data
  ;dynamics
  dyn_array = read_analyze('Dynamics/', n_timepoints = 96); reads in the synthetic dynamic time series
  s = size(dyn_array)
  column = s[1]
  row = s[2]
  n_slices = s[3]
  n_timepoints = s[4]
  
  ;vfa images
  nfa = 3; number of unique flip angles
  
  Fa_array = fltarr(s[1],s[2],s[3],nfa); define a 4D array to hold flip angle images
  
  filename_hdr5_0 =  'VFA\fa_5.hdr'
  filename_hdr10_0 =  'VFA\fa_10.hdr'
  filename_hdr35_0 =  'VFA\fa_35.hdr'
  
  image_obj = Obj_New('IDLffanalyze'); define an analyze object to read in analyze format files
  
  image_obj->READFILE,filename_hdr5_0
  Fa_array[*,*,*,0] = image_obj->GETPIXELDATA()
  
  image_obj->READFILE,filename_hdr10_0
  Fa_array[*,*,*,1] = image_obj->GETPIXELDATA()
  
  image_obj->READFILE,filename_hdr35_0
  Fa_array[*,*,*,2] = image_obj->GETPIXELDATA()
  
  ;reads in ROI mask
  image_obj->READFILE, 'ROImask.hdr'
  roi_mask = image_obj->GETPIXELDATA()
  index = where(roi_mask gt 0)
  roi_fourD = fltarr(column,row,n_slices,n_timepoints)
  roi_mask[index] = 100.0
  
  for i = 0, n_timepoints -1 do begin
  
    roi_fourD[*,*,*,i] = roi_mask
    
  endfor
  
  tumour_vox = dyn_array[*,*,*,0:6]*roi_fourD[*,*,*,0:6]/100.0
  
  ;----------------------------------AIFs
  ;accurate AIF
  if aif_type eq '1' then begin
    ;AIF from randomly chosen patient from bladder data
    restore, 'PhantomAIF.dat'
    ca =  fit_ext
    aif_cap = 4 ; this sets the first four timepoints to zero.
    baseline_cap = 5.0; check this
    timepoints = findgen(n_timepoints)*2.5/60; timepoints in minutes
  endif
  
  ;overestimated AIF (scale x1.5)
  if aif_type eq '2' then begin
    ;AIF from randomly chosen patient from bladder data
    restore, 'PhantomAIF.dat'
    ca =  fit_ext*1.5
    aif_cap = 4 ; this sets the first four timepoints to zero.
    baseline_cap = 5.0; check this
    timepoints = findgen(n_timepoints)*2.5/60; timepoints in minutes
  endif
  
  ;underestimated AIF (scale x0.5)
  if aif_type eq '3' then begin
    ;AIF from randomly chosen patient from bladder data
    restore, 'PhantomAIF.dat'
    ca =  fit_ext*0.5
    aif_cap = 4 ; this sets the first four timepoints to zero.
    baseline_cap = 5.0; check this
    timepoints = findgen(n_timepoints)*2.5/60 ; timepoints in minutes
  endif
  
  ;-------------------------------- Add noise to images
  
  sigmaDYN = mean_baseline/snr_dyn
  dyn_array_noisy = randomgauss(1,sigmaDYN,dyn_array); adds gaussian noise
  sigmaVFA = mean_baseline/snr_vfa
  Fa_array_noisy = randomgauss(2,sigmaVFA,Fa_array); adds gaussian noise
  
  actual_sigma = {sigmaDYN:sigmaDYN, sigmaVFA:sigmaVFA}
  noisefree_images = fltarr(s[1],s[2],s[3],n_timepoints + nfa)
  noisefree_images[*,*,*,0:n_timepoints-1] = dyn_array
  noisefree_images[*,*,*,n_timepoints:*] = fa_array
  
  ;----------------------------- fit T1 mapping and tracer kinetic model to data using 'in house' DCEMRIprocessor code
  
  DCEMRIprocessor_freeS0_define; defines all variables needed for the DCEMRIprocessor object
  processor = obj_new('DCEMRIprocessor_freeS0')
  
  processor->GET_DATA,                    $ ; input data to dcemriProcessor object
    TR = TR,                              $
    alpha_dyn = alpha_dyn,                $
    alpha_vfa = alpha_fa,                 $
    n_timepoints = n_timepoints,          $
    n_slices = n_slices,                  $
    row = row,                            $
    column = column,                      $
    relaxivity = r,                       $
    timepoints = timepoints,              $
    DCE_timeseries = dyn_array_noisy,     $
    O_DCE_timeseries = dyn_array_noisy,   $
    vfa_images = fa_array_noisy,          $
    O_vfa_images = fa_array_noisy,        $
    ROI_mask = ROI_fourD,                 $
    noisefree_images = noisefree_images
    
  Fp_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1); define arrays to hold parameter estimates
  FE_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  vp_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  ve_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  T1_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  S0d_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  S0v_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  ost_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  FA_array_repeats = fltarr(s[1],s[2],s[3],n_repeats+1)
  
    ; make a folder in directory to save parameter estimates 
    CD, CURRENT=c
    string =  c + '\' + fit_type + '_snrDYN' + strtrim(string(snr_dyn,FORMAT = '( 2(I5))'),1) + '_snrVFA' + strtrim(string(snr_vfa,FORMAT = '( 2(I5))'),1) + '_' + aif_type
    File_MKdir, string

    ; main function to fit models. Parameter estimates and fits are outputed in a structure
    ; fits to n_vox_fit voxels once
    output_struct = processor->Generate_maps(fit_type,'phantom',ca = ca,/update,T1_datatype = 'vfa',model = model, ostfORf = ostfORf,FAfORf = FAfORf, aif_cap = aif_cap,baseline_cap = baseline_cap, n_vox_fit = n_vox_fit) ; fits model for the first time. Outputs a structure containing parameter maps for each parameter.
    
    ref_timeseries = output_struct.ref_timeseries
    ref_vfa_images = output_struct.ref_vfa_images
   
    ; repeat fits to n_vox_fit voxels (Monte Carlo experiement)
    for i = 0, n_repeats do begin; For each repeat maps are recorded into arrays below.
    
      Fp_array_repeats[*,*,*,i] = output_struct.Fp_map
      FE_array_repeats[*,*,*,i] = output_struct.FE_map
      vp_array_repeats[*,*,*,i] = output_struct.vp_map
      ve_array_repeats[*,*,*,i] = output_struct.ve_map
      T1_array_repeats[*,*,*,i] = output_struct.T1_map
      S0d_array_repeats[*,*,*,i] = output_struct.S0d_map
      S0v_array_repeats[*,*,*,i] = output_struct.S0v_map
      ost_array_repeats[*,*,*,i] = output_struct.ost_map
      FA_array_repeats[*,*,*,i] =  output_struct.flip_angle_map
      
      if i le n_repeats - 1 then begin

        seed_d = i+2
        seed_vfa = i+3
        ; adds different sample of noise to noise-free images
        dyn_array_noisy = randomgauss(seed_d,sigmaDYN,dyn_array); adds a new sample of noise to the noisefree data
        Fa_array_noisy = randomgauss(seed_vfa,sigmaVFA,Fa_array)
        ;updates images in DCEMRIprocessor object
        processor->GET_DATA, $
          DCE_timeseries = dyn_array_noisy,$
          vfa_images = fa_array_noisy ; updates the dce and vfa images in dcemriprocessor.
        ;refits to new images  
        output_struct = processor->Generate_maps(fit_type,'phantom',T1_datatype = 'vfa',ca = ca,n_vox_fit = n_vox_fit, model = model, ostfORf = ostfORf,aif_cap = aif_cap, baseline_cap = baseline_cap); repeat fitting to new noisy data
        
      endif
      
    endfor
    
    ;saves parameter maps
    save,Fp_array_repeats,filename = string + '\Fp.dat'
    save,FE_array_repeats,filename = string + '\FE.dat'
    save,vp_array_repeats,filename = string + '\vp.dat'
    save,ve_array_repeats,filename = string + '\ve.dat'
    save,T1_array_repeats,filename = string + '\T1.dat'
    save,S0d_array_repeats,filename = string + '\S0d.dat'
    save,S0v_array_repeats,filename = string + '\S0v.dat'
    save,ost_array_repeats,filename = string + '\ost.dat'
      
  return, 1
  obj_destroy,processor
  
END


