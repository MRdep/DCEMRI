for i=5,5 do begin &$
tag = "" &$
a = run_2CXM_phantom(snr_dyn = snr_dyn, snr_vfa = snr_vfa, T1_dataType = T1_dataType,model = model,multithreadindex = index[i],tag = tag,n_repeats = n_repeats,n_vox_fit = n_vox_fit) &$
endfor
