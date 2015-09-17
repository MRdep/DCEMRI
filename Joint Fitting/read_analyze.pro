function read_analyze, directory, n_timepoints = n_timepoints

;Reads in a phantom dynamic timeseries. 

  image_obj = OBJ_NEW('IDLffanalyze') ; read in info and pixels for the whole tumour mask
  
  image_obj->READFILE, directory + '\dyn_0.hdr'
  img = image_obj->GETPIXELDATA()
  
  size_img = size(img)
  img_4D = fltarr(size_img[1],size_img[2],size_img[3],n_timepoints)
  dir = directory
  
  for i = 0, n_timepoints -1 do begin

    filename_hdr =  dir + '\dyn_' + strtrim(strcompress(string(i)),1) + '.hdr'
    image_obj->READFILE,filename_hdr
    img_4D[*,*,*,i] = image_obj->GETPIXELDATA()

  endfor
 
  obj_destroy, image_obj
  return, img_4D
  
end
