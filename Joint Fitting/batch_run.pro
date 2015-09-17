pro batch_run

;This file runs the script 'run_2cxm_phantom' in batch mode. It uses multithreading to access multiple cores to increase speed.

!EXCEPT=0
 
 mean_baseline = 163.0
 
 COMMON S0_info, equal_and_free, unequal_and_free, equal_by_baseline
 COMMON sigma_info, fixed_sigma, equal_sigma, known_sigma
 
 ;parameters relating to how M0v and M0d are varied during fitting
 equal_and_free = 1
 unequal_and_free = 0
 equal_by_baseline = 0
 
 ;multithread indices. 6 jobs, sequential and joint estimation for 3 AIF conditions  
 index = [1,2,3,4,5,6]

 T1_dataType = 'vfa'; could adapt to fit to IR or SR T1 mapping data
 model = 'XM'; could change to Tofts or E.Tofts 
 n_repeats = 100; number of Monte Carlo repeats
 n_vox_fit = 504; number of voxels to sample
 snr_dyn = 5.0; required SNR of dynamic images 
 snr_vfa = 11.2; required SNR of VFA images -> snr_dyn*sqrt(5) to account for 5 signal averages

split_for,0,n_elements(index)-1,nsplit = 8,command = [$
'tag = ""',$
'a = run_2CXM_phantom(snr_dyn = snr_dyn, snr_vfa = snr_vfa, T1_dataType = T1_dataType,model = model,multithreadindex = index[i],n_repeats = n_repeats,n_vox_fit = n_vox_fit)'],$
varnames=['snr_dyn','snr_vfa','T1_datatype','model','index','n_repeats','n_vox_fit']

end