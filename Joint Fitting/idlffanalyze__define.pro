;
;  File:  IDLffAnalyze__define.pro
;
;  Defines the IDLffAnalyze object which can be used to read/write
;  ANALYZE 7.5 format files.
;
;  RTK, 23-Nov-2004
;  Last update:  09-Dec-2004
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  Interprets Analyze 7.5 files as defined in "ANALYZE 7.5 File Format" PDF
;  with (some) extensions for SPM99 and SPM2.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;
;  Private:
;

;-----------------------------------------------------------
;  IDLffAnalyze::isNum
;+
;  Assert that a value is a number of a specific type and
;  within a specific range.
;
;  @returns True (1) if all assertions pass, (0) otherwise.
;
;  @param n {in}{type=any}
;    The number (or array) which will be tested.  If an
;    array, each element of the array must pass the test.
;
;  @keyword POSITIVE {in}{type=boolean}{optional}
;    Number must be positive (and not zero unless /ZERO also set)
;
;  @keyword NEGATIVE {in}{type=boolean}{optional}
;    Number must be negative (and not zero unless /ZERO also set)
;
;  @keyword ZERO {in}{type=boolean}{optional}
;    Number must be exactly zero or adds zero to the allowed range
;    if /POSITIVE or /NEGATIVE set.
;
;  @keyword LIMITS {in}{type=2-element vector}{optional}
;    If present, represents the lower and upper limits allowed for n.
;    If set, other keywords except TYPE are ignored.
;
;  @keyword INTEGER {in}{type=boolean}{optional}
;    If set input must be an integer (IDL types 1,2,3,12,13,14, or 15).
;
;  @keyowrd FLOATING {in}{type=boolean}{optional}
;    If set input must be floating point real (IDL types 4 and 5).
;
;  @keyword TYPE {in}{type=integer}{optional}
;    If present, also asserts that n is of the given IDL type.  Overrides
;    /INTEGER or /FLOATING.
;
;@private
;-
function IDLffAnalyze::isNum, n, POSITIVE=pos, NEGATIVE=neg, ZERO=zero,   $
  LIMITS=lim, TYPE=type, INTEGER=integer,  $
  FLOATING=floating
  compile_opt idl2, HIDDEN
  
  ;
  ;  Assert that n is a number and of the proper type
  ;
  t = size(n, /TYPE)
  
  if where(t eq [0,7,8,10,11]) ne -1 then  $
    return, 0
    
  if n_elements(type) ne 0 then begin
    if t ne type[0] then  return, 0
  endif else begin
    z = keyword_set(integer) + 2*keyword_set(floating)
    case z of
      0: ; no other type keywords, ignore
      1: if (where(t eq [1,2,3,12,13,14,15]) eq -1) then  return, 0  ; INTEGER
      2: if (t ne 4) and (t ne 5) then  return, 0                    ; FLOATING
      3: message, 'Incompatible keywords: /INTEGER and /FLOATING.'
    endcase
  endelse
  
  ;
  ;  Use limits or other keywords
  ;
  if ~keyword_set(lim) then begin
    ;
    ;  Assert that n is positive, negative, or zero, etc.
    ;
    z = keyword_set(pos) + 2*keyword_set(neg) + 4*keyword_set(zero)
    
    case z of
      0: ; no keywords set, ignored
      1: if (min(n) le 0.0) then  return, 0         ; positive only
      2: if (max(n) ge 0.0) then  return, 0         ; negative only
      3: if (where(n eq 0.0) ne -1) then return, 0  ; pos or neg, anything but zero
      4: if (where(n ne 0.0) ne -1) then return, 0  ; zero only
      5: if (min(n) lt 0.0) then  return, 0         ; pos or zero
      6: if (max(n) gt 0.0) then  return, 0         ; neg or zero
      7: ; any value, ignore these keywords
    endcase
    
  endif else begin
    ;
    ;  Assert that all elements are in [lim[0], lim[1]]
    ;
    if (n_elements(lim) ne 2) then  $
      message, 'LIMITS must be a two-element numeric vector.'
    if (where(size(lim,/TYPE) eq [0,7,8,10,11]) ne -1) then  $
      message, 'LIMITS must be a numeric vector.'
      
    lo = lim[0]
    hi = lim[1]
    if (lo ge hi) then  $
      message, 'LIMITS must be of the form [lo,hi] with lo > hi.'
      
    if (min(n) lt lo) or (max(n) gt hi) then  $
      return, 0
  endelse
  
  ;  All assertions pass
  return, 1
end


;-----------------------------------------------------------
;  IDLffAnalyze::Swap
;+
;  @param n {in}{type=numeric}  The value to be byte swapped.
;
;@private
;-
pro IDLffAnalyze::Swap, n
  compile_opt idl2, HIDDEN
  
  case size(n,/TYPE) of
    2: byteorder, n, /SSWAP      ; short
    3: byteorder, n, /LSWAP      ; long
    4: byteorder, n, /LSWAP      ; float
    5: byteorder, n, /L64SWAP    ; double
    6: byteorder, n, /L64SWAP    ; complex
    ELSE: ; nothing, leave as is
  endcase
end


;-----------------------------------------------------------
;  IDLffAnalyze::SwapHeader
;+
;  Byte swap a HEADER structure.
;
;  @param h {in}{type=HEADER structure}
;    The structure to byte swap.
;
;@private
;-
pro IDLffAnalyze::SwapHeader, h
  compile_opt idl2, HIDDEN
  
  ;  Must be a header structure
  if (size(h,/TYPE) ne 8) then  $
    message, 'Argument must be a HEADER structure.'
  if (size(h,/SNAME) ne 'HEADER') then  $
    message, 'Argument must be a HEADER structure.'
    
  ;
  ;  Only swap numeric values, type determines how it is swapped.
  ;
  q = h.sizeof_hdr    &  self->Swap,q  &  h.sizeof_hdr = q     ;   int sizeof_hdr;
  q = h.extents       &  self->Swap,q  &  h.extents = q        ;   int extents;
  q = h.session_error &  self->Swap,q  &  h.session_error = q  ;   short int session_error;
  q = h.dim           &  self->Swap,q  &  h.dim = q            ;   short int dim[8];
  q = h.unused1       &  self->Swap,q  &  h.unused1 = q        ;   short int unused1;
  q = h.datatype      &  self->Swap,q  &  h.datatype = q       ;   short int datatype;
  q = h.bitpix        &  self->Swap,q  &  h.bitpix = q         ;   short int bitpix;
  q = h.dim_un0       &  self->Swap,q  &  h.dim_un0 = q        ;   short int dim_un0;
  q = h.pixdim        &  self->Swap,q  &  h.pixdim = q         ;   float pixdim[8];
  q = h.vox_offset    &  self->Swap,q  &  h.vox_offset = q     ;   float vox_offset;
  q = h.roi_scale     &  self->Swap,q  &  h.roi_scale = q      ;   float roi_scale;
  q = h.funused1      &  self->Swap,q  &  h.funused1 = q       ;   float funused1;  // scalefactor in SPM
  q = h.funused2      &  self->Swap,q  &  h.funused2 = q       ;   float funused2;  // dcoff in SPM2
  q = h.cal_max       &  self->Swap,q  &  h.cal_max = q        ;   float cal_max;
  q = h.cal_min       &  self->Swap,q  &  h.cal_min = q        ;   float cal_min;
  q = h.compressed    &  self->Swap,q  &  h.compressed = q     ;   int compressed;
  q = h.verified      &  self->Swap,q  &  h.verified = q       ;   int verified;
  q = h.glmax         &  self->Swap,q  &  h.glmax = q          ;   int glmax;
  q = h.glmin         &  self->Swap,q  &  h.glmin = q          ;   int glmin;
  q = h.orient        &  self->Swap,q  &  h.orient = q         ;   char orient;
  q = h.views         &  self->Swap,q  &  h.views = q          ;   int views;
  q = h.vols_added    &  self->Swap,q  &  h.vols_added = q     ;   int vols_added;
  q = h.start_field   &  self->Swap,q  &  h.start_field = q    ;   int start_field;
  q = h.field_skip    &  self->Swap,q  &  h.field_skip = q     ;   int field_skip;
  q = h.omax          &  self->Swap,q  &  h.omax = q           ;   int omax;
  q = h.omin          &  self->Swap,q  &  h.omin = q           ;   int omin;
  q = h.smax          &  self->Swap,q  &  h.smax = q           ;   int smax;
  q = h.smin          &  self->Swap,q  &  h.smin = q           ;   int smin;
  
  ;
  ;  Originator is a special case.  Swap as if it were being used by
  ;  SPM since it appears that little else actually uses this field.
  ;
  ;  char originator[10]; // image_origin in SPM
  ;
  w = h.originator
  t = w[0] &  w[0] = w[1]  &  w[1] = t
  t = w[2] &  w[2] = w[3]  &  w[3] = t
  t = w[4] &  w[4] = w[5]  &  w[5] = t
  h.originator = w
end


;-----------------------------------------------------------
;  IDLffAnalyze::IDLtoAnalyzeType
;+
;  @returns The Analyze type for a given IDL type, -1 if
;    Analyze does not support the IDL type.
;
;  @param itype {in}{type=integer}
;    The IDL type code.
;
;@private
;-
function IDLffAnalyze::IDLtoAnalyzeType, itype
  compile_opt idl2, HIDDEN
  
  idx = where(itype eq indgen(16))
  return, ([0,2,4,8,16,64,32,-1,-1,-1,-1,-1,-1,-1,-1,-1])[idx[0]]
end


;-----------------------------------------------------------
;  IDLffAnalyze::AnalyzetoIDLType
;+
;  @returns The IDL type for a given Analyze type.
;
;  @param atype {in}{type=integer}
;    The Analyze type code.
;
;@private
;-
function IDLffAnalyze::AnalyzetoIDLType, atype
  compile_opt idl2, HIDDEN
  
  idx = (where(atype eq [0,1,2,4,8,16,64,32]))[0]
  return, ([0,0,1,2,3,4,5,6])[idx]
end


;-----------------------------------------------------------
;  IDLffAnalyze::IDLtoAnalyzeBitPix
;+
;  @returns The Analyze bits per pixel for a given IDL type, -1 if
;    Analyze does not support the IDL type.
;
;  @param itype {in}{type=integer}
;    The IDL type code.
;
;@private
;-
function IDLffAnalyze::IDLtoAnalyzeBitPix, itype
  compile_opt idl2, HIDDEN
  
  idx = where(itype eq indgen(16))
  return, ([0,8,16,32,32,64,64,-1,-1,-1,-1,-1,-1,-1,-1,-1])[idx[0]]
end


;-----------------------------------------------------------
;  IDLffAnalyze::SystemByteOrder
;+
;  @returns True (1) if the system is big-endian, false (0)
;    if the system is little-endian.
;
;@private
;-
function IDLffAnalyze::SystemByteOrder
  compile_opt HIDDEN
  
  b = (a=123)
  byteorder, a, /SWAP_IF_LITTLE_ENDIAN
  return, (b ne a) ? 0 : 1
end


;-----------------------------------------------------------
;  IDLffAnalyze::Destroy
;+
;  Free memory used by this object.
;
;@private
;-
pro IDLffAnalyze::Destroy
  compile_opt idl2, HIDDEN
  
  ;  remove any image data
  if ptr_valid(self.imgdata) then  ptr_free, self.imgdata
end


;
;  These functions (superficially) validate property assignments:
;

;-----------------------------------------------------------
;  IDLffAnalyze::UpdateAll
;+
;  This is a dangerous function.  *No* checking whatsoever
;  is performed for consistent values.
;
;  @param all {in}{type=header struct} Structure containing
;             the entire Analyze header.
;
;@private
;-
pro IDLffAnalyze::UpdateAll, all
  compile_opt idl2, HIDDEN
  
  ; Must be a header structure
  if (size(all,/SNAME) ne 'HEADER') then  $
    message, 'ALL must be a header structure.'
    
  self.hdr = all
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateDim
;+
;  @param dim {in}{type=vector} Dimension vector.
;
;@private
;-
pro IDLffAnalyze::UpdateDim, dim
  compile_opt idl2, HIDDEN
  
  ; Must be integer
  if (~self->isNum(dim, /INTEGER)) then  $
    message, 'DIM must be integer.'
  if (n_elements(dim) ne 8) then  $
    message, 'DIM must be an 8 element vector.'
    
  self.hdr.dim = dim
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateGlmin
;+
;  @param glmin {in}{type=int32} Minimum pixel value for entire file.
;
;@private
;-
pro IDLffAnalyze::UpdateGlmin, glmin
  compile_opt idl2, HIDDEN
  
  ; Must be an integer
  if (~self->isNum(glmin, /INTEGER)) then  $
    message, 'GLMIN must be an integer.'
  if (n_elements(glmin) ne 1) then  $
    message, 'GLMIN must be a scalar.'
    
  self.hdr.glmin = glmin  ; IDL converts to LONG if necessary
  ; w/o warning if precision is lost.
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateGlmax
;+
;  @param glmax {in}{type=int32} Maximum pixel value for entire file.
;
;@private
;-
pro IDLffAnalyze::UpdateGlmax, glmax
  compile_opt idl2, HIDDEN
  
  ; Must be an integer
  if (~self->isNum(glmax, /INTEGER)) then  $
    message, 'GLMAX must be an integer.'
  if (n_elements(glmax) ne 1) then  $
    message, 'GLMAX must be a scalar.'
    
  self.hdr.glmax = glmax
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdatePixdim
;+
;  @param pixdim {in}{type=vector} Pixel dimensions vector.
;
;@private
;-
pro IDLffAnalyze::UpdatePixdim, pixdim
  compile_opt idl2, HIDDEN
  
  ;  Must be floating point
  if (~self->isNum(pixdim, /FLOATING)) then  $
    message, 'PIXDIM must be a 32-bit float.'
  if (n_elements(pixdim) ne 8) then  $
    message, 'PIXDIM must be an 8 element vector.'
    
  ; Assign w/o regard for meaningfulness
  self.hdr.pixdim = pixdim
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateVoxUnits
;+
;  @param vox_units {in}{type=character*4} Spatial voxel units of measure.
;
;@private
;-
pro IDLffAnalyze::UpdateVoxUnits, vox_units
  compile_opt idl2, HIDDEN
  
  ; Must be a single string
  if (size(vox_units,/TYPE) ne 7) then  $
    message, 'VOX_UNITS must be a string.'
  if (n_elements(vox_units) ne 1) then  $
    message, 'VOX_UNITS must be a single string.'
    
  ; Chop if too long
  if (strlen(vox_units) le 4) then begin
    self.hdr.vox_units = byte('    ')
    self.hdr.vox_units = byte(vox_units)
  endif else begin
    self.hdr.vox_units = byte(strmid(vox_units,0,4))
  endelse
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateCalUnits
;+
;  Mindlessly chops input that is too long and adds trailing
;  spaces to input that is too short.
;
;  @param cal_units {in}{type=character*8} Calibration units.
;
;@private
;-
pro IDLffAnalyze::UpdateCalUnits, cal_units
  compile_opt idl2, HIDDEN
  
  ; Must be a single string
  if (size(cal_units,/TYPE) ne 7) then  $
    message, 'CAL_UNITS must be a string.'
  if (n_elements(cal_units) ne 1) then  $
    message, 'CAL_UNITS must be a single string.'
    
  ; Chop if too long
  if (strlen(cal_units) le 8) then begin
    self.hdr.cal_units = byte('        ')
    self.hdr.cal_units = byte(cal_units)
  endif else begin
    self.hdr.cal_units = byte(strmid(cal_units,0,8))
  endelse
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateDatatype
;+
;  N.B. This function expects the Analyze datatype codes, not
;       IDL datatype codes.  See self->DT_XXX() functions.
;
;  @param datatype {in}{type=int16} Data type code for image pixels.
;
;@private
;-
pro IDLffAnalyze::UpdateDatatype, datatype
  compile_opt idl2, HIDDEN
  
  ; Must be numeric and scalar
  if (~self->isNum(datatype, /POSITIVE)) then  $
    message, 'DATATYPE must be numeric.'
  if (n_elements(datatype) ne 1) then  $
    message, 'DATATYPE must be a scalar.'
    
  ; Must also be a valid code
  if (where(datatype eq [0,1,2,4,8,16,32,64,128,255]) eq -1) then $
    message, 'DATATYPE must be a valid Analyze 7.5 datatype code.'
    
  self.hdr.datatype = datatype
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateBitpix
;+
;  @param bitpix {in}{type=int16} Bits per pixel.
;
;@private
;-
pro IDLffAnalyze::UpdateBitpix, bitpix
  compile_opt idl2, HIDDEN
  
  ; Must be numeric and scalar
  if (~self->isNum(bitpix, /POSITIVE)) then  $
    message, 'BITPIX must be numeric.'
  if (n_elements(bitpix) ne 1) then  $
    message, 'BITPIX must be a scalar.'
    
  ; Must also be 1, 8, 16, 32, or 64
  if (where(bitpix eq [1,8,16,32,64]) eq -1) then $
    message, 'BITPIX must be 1, 8, 16, 32, or 64.'
    
  self.hdr.bitpix = bitpix
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateVoxOffset
;+
;  @param vox_offset {in}{type=float32}  Offset, in bytes,
;                    to the start of voxel data in the .img file.
;
;@private
;-
pro IDLffAnalyze::UpdateVoxOffset, vox_offset
  compile_opt idl2, HIDDEN
  
  ;  Must be a positive or zero integer scalar
  if (~self->isNum(vox_offset, /POS, /ZERO)) then  $
    message, 'VOX_OFFSET must be a number >= 0.'
  if (n_elements(vox_offset) ne 1) then  $
    message, 'VOX_OFFSET must be a scalar.'
    
  self.hdr.vox_offset = vox_offset
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateCalMin
;+
;  @param cal_min {in}{type=float32} Minimum calibration value.
;
;@private
;-
pro IDLffAnalyze::UpdateCalMin, cal_min
  compile_opt idl2, HIDDEN
  
  ;  Must be a scalar number
  if (~self->isNum(cal_min)) then  $
    message, 'CAL_MIN must be a number.'
  if (n_elements(cal_min) ne 1) then  $
    message, 'CAL_MIN must be a scalar.'
    
  self.hdr.cal_min = cal_min
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateCalMax
;+
;  @param cal_max {in}{type=float32} Maximum calibration value.
;
;@private
;-
pro IDLffAnalyze::UpdateCalMax, cal_max
  compile_opt idl2, HIDDEN
  
  ;  Must be a scalar number
  if (~self->isNum(cal_max)) then  $
    message, 'CAL_MAX must be a number.'
  if (n_elements(cal_max) ne 1) then  $
    message, 'CAL_MAX must be a scalar.'
    
  self.hdr.cal_max = cal_max
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateOrient
;+
;
;  N.B. This is the primitive approach to updating the
;       orientation parameter.  It simply checks for values
;       in range and mindlessly sets the value given.
;
;  Use IDLffAnalyze::SetOrientation instead unless one is very
;  certain of what he is doing.
;
;  @param orient {in}{type=byte} Slice orientation indicator.
;
;@private
;-
pro IDLffAnalyze::UpdateOrient, orient
  compile_opt idl2, HIDDEN
  
  ; Must be an integer scalar
  if (~self->isNum(orient, /INTEGER)) then  $
    message, 'ORIENT must be an integer.'
  if (n_elements(orient) ne 1) then  $
    message, 'ORIENT must be a scalar.'
    
  self.hdr.orient = orient
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateScaleFactor
;+
;  @param scal {in}{type=float32} Pixel scale factor (SPM2, SPM99)
;
;@private
;-
pro IDLffAnalyze::UpdateScaleFactor, scal
  compile_opt idl2, HIDDEN
  
  ; Must be a numeric scalar
  if (~self->isNum(scal)) then  $
    message, 'SCALEFACTOR must be numeric.'
  if (n_elements(scal) ne 1) then  $
    message, 'SCALEFACTOR must be a scalar.'
    
  self.hdr.funused1 = scal  ; IDL will not change the type of funused1
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateDCoff
;+
;  @param dcoff {in}{type=float32} Image intensity zero intercept. (SPM2)
;
;@private
;-
pro IDLffAnalyze::UpdateDCoff, dcoff
  compile_opt idl2, HIDDEN
  
  ; Must be a numeric scalar
  if (~self->isNum(dcoff)) then  $
    message, 'DC offset must be numeric.'
  if (n_elements(dcoff) ne 1) then  $
    message, 'DC offset must be a scalar'
    
  self.hdr.funused2 = dcoff
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateImageOrigin
;+
;  @param im {in}{type=int16 vector} Image origin (SPM2, SPM99)
;
;@private
;-
pro IDLffAnalyze::UpdateImageOrigin, im
  compile_opt idl2, HIDDEN
  
  ; Must be numeric and of three elements
  if (~self->isNum(im, TYPE=2)) then  $
    message, 'IMAGE_ORIGIN must be a short integer vector.'
  if (n_elements(im) ne 3) then  $
    message, 'IMAGE_ORIGIN must be a three element vector.'
    
  ;
  ;  Since image_origin is overlayed on originator (10 element, char)
  ;  the values have to be assigned manually according with the
  ;  proper endian-ness.
  ;
  w = self.hdr.originator
  if (self.endian eq 1) then begin
    w[0] = byte(im[0]/256)  &  w[1] = byte(im[0] mod 256)  ; big-endian
    w[2] = byte(im[1]/256)  &  w[3] = byte(im[1] mod 256)
    w[4] = byte(im[2]/256)  &  w[5] = byte(im[2] mod 256)
  endif else begin
    w[1] = byte(im[0]/256)  &  w[0] = byte(im[0] mod 256)  ; little-endian
    w[3] = byte(im[1]/256)  &  w[2] = byte(im[1] mod 256)
    w[5] = byte(im[2]/256)  &  w[4] = byte(im[2] mod 256)
  endelse
  self.hdr.originator = w
end


;-----------------------------------------------------------
;  IDLffAnalyze::UpdateEndian
;+
;  @param endian {in}{type=integer}  Endian, 1 (big), 0 (little).
;
;@private
;-
pro IDLffAnalyze::UpdateEndian, endian
  compile_opt idl2, HIDDEN
  
  ;  Must be 0 or 1
  if (~self->isNum(endian, /INTEGER, LIMITS=[0,1])) then  $
    message, 'ENDIAN must be 0 (little) or 1 (big).'
    
  self.endian = endian
end


;----------------------------------------------------------
;  IDLffAnalyze::ReadTRSImageFile
;+
;  @returns  The properly read and byte swapped image data
;            for the given pathname.
;
;  @param f {in}{type=string}
;    Pathname for the image file (.img).
;
;  @param h {in}{type=instance of HEADER structure}
;    The header values for this file.
;
;  @param endian {in}{type=boolean}
;    File byte order, 1=big, 0=little.
;
;  @keyword FLIP {in}{type=boolean}{optional}
;    If set, flip the proper axis of the image data.
;
;@private
;-
function IDLffAnalyze::ReadTRSImageFile, f, h, endian, FLIP=flip
  compile_opt idl2, HIDDEN
  
  ;  Create the image data array of the proper type
  ans = make_array(h.dim[1:4], TYPE=self->AnalyzetoIDLType(h.datatype))
  
  ;  Ready the file
  openr, u, f, /GET_LUN
  pix = make_array(h.dim[1], TYPE=self->AnalyzetoIDLType(h.datatype))
  
  ;
  ;  Read the file
  ;
  if (~keyword_set(flip)) then begin
    for v = 0, h.dim[4]-1 do begin          ; normal
      for z = 0, h.dim[3]-1 do begin
        for y = 0, h.dim[2]-1 do begin
          readu, u, pix
          ans[*,y,z,v] = pix
        endfor
      endfor
    endfor
  endif else begin
    for v = 0, h.dim[4]-1 do begin          ; flipped
      for z = 0, h.dim[3]-1 do begin
        for y = h.dim[2]-1,0,-1 do begin
          readu, u, pix
          ans[*,y,z,v] = pix
        endfor
      endfor
    endfor
  endelse
  
  free_lun, u
  
  ;  change endian if necessary
  if (endian ne self->SystemByteOrder()) then  $
    self->Swap, ans
    
  return, ans
end


;----------------------------------------------------------
;  IDLffAnalyze::ReadCORImageFile
;+
;  @returns  The properly read and byte swapped image data
;            for the given pathname.
;
;  @param f {in}{type=string}
;    Pathname for the image file (.img).
;
;  @param h {in}{type=instance of HEADER structure}
;    The header values for this file.
;
;  @param endian {in}{type=boolean}
;    File byte order, 1=big, 0=little.
;
;  @keyword FLIP {in}{type=boolean}{optional}
;    If set, flip the proper axis of the image data.
;
;@private
;-
function IDLffAnalyze::ReadCORImageFile, f, h, endian, FLIP=flip
  compile_opt idl2, HIDDEN
  
  ;  Create the image data array of the proper type
  ans = make_array(h.dim[1], h.dim[3], h.dim[2], h.dim[4], TYPE=self->AnalyzetoIDLType(h.datatype))
  
  ;  Ready the file
  openr, u, f, /GET_LUN
  pix = make_array(h.dim[1], TYPE=self->AnalyzetoIDLType(h.datatype))
  
  ;
  ;  Read the file
  ;
  if (~keyword_set(flip)) then begin
    for v = 0, h.dim[4]-1 do begin          ; normal
      for y = 0, h.dim[3]-1 do begin
        for z = 0, h.dim[2]-1 do begin
          readu, u, pix
          ans[*,y,z,v] = pix
        endfor
      endfor
    endfor
  endif else begin
    for v = 0, h.dim[4]-1 do begin          ; flipped
      for y = 0, h.dim[3]-1 do begin
        for z = h.dim[2]-1,0,-1 do begin
          readu, u, pix
          ans[*,y,z,v] = pix
        endfor
      endfor
    endfor
  endelse
  
  free_lun, u
  
  ;  change endian if necessary
  if (endian ne self->SystemByteOrder()) then  $
    self->Swap, ans
    
  return, ans
end


;----------------------------------------------------------
;  IDLffAnalyze::ReadSAGImageFile
;+
;  @returns  The properly read and byte swapped image data
;            for the given pathname.
;
;  @param f {in}{type=string}
;    Pathname for the image file (.img).
;
;  @param h {in}{type=instance of HEADER structure}
;    The header values for this file.
;
;  @param endian {in}{type=boolean}
;    File byte order, 1=big, 0=little.
;
;  @keyword FLIP {in}{type=boolean}{optional}
;    If set, flip the proper axis of the image data.
;
;@private
;-
function IDLffAnalyze::ReadSAGImageFile, f, h, endian, FLIP=flip
  compile_opt idl2, HIDDEN
  
  ;  Create the image data array of the proper type
  ans = make_array(h.dim[3], h.dim[1], h.dim[2], h.dim[4], TYPE=self->AnalyzetoIDLType(h.datatype))
  
  ;  Ready the file
  openr, u, f, /GET_LUN
  pix = make_array(h.dim[1], TYPE=self->AnalyzetoIDLType(h.datatype))
  
  ;
  ;  Read the file
  ;
  if (~keyword_set(flip)) then begin
    for v = 0, h.dim[4]-1 do begin        ; normal
      for x = 0, h.dim[3]-1 do begin
        for z = 0, h.dim[2]-1 do begin
          readu, u, pix
          ans[x,*,z,v] = pix
        endfor
      endfor
    endfor
  endif else begin
    for v = 0, h.dim[4]-1 do begin        ; flipped
      for x = 0, h.dim[3]-1 do begin
        for z = h.dim[2]-1,0,-1 do begin
          readu, u, pix
          ans[x,*,z,v] = pix
        endfor
      endfor
    endfor
  endelse
  
  free_lun, u
  
  ;  change endian if necessary
  if (endian ne self->SystemByteOrder()) then  $
    self->Swap, ans
    
  return, ans
end


;----------------------------------------------------------
;  IDLffAnalyze::ReadImageFile
;+
;  @returns The image data for the given file.
;
;  @param f {in}{type=string}
;    The pathname of the image file.
;
;  @param h {in}{type=HEADER structure}
;    The header file structure associated with this image data.
;
;  @param orient {in}{type=integer}
;    Orientation to use when reading the file.
;
;  @keyword ENDIAN {in}{type=integer}
;    Byte order of the image data.
;
;@private
;-
function IDLffAnalyze::ReadImageFile, f, h, orient, ENDIAN=endian
  compile_opt idl2, HIDDEN
  
  ;
  ;  Add a .img extension if not already present.
  ;
  if (strmatch(f,'*.IMG',/FOLD_CASE) eq 1) then begin
    fname = f  ; already ends in .img
  endif else begin
    if (strmatch(f,'*.HDR',/FOLD_CASE) eq 1) then begin
      fname = strmid(f,0,strlen(f)-4) + '.img'
    endif else begin
      fname = f + '.img'
    endelse
  endelse
  
  ;  Verify that there is a data file to read.
  if (~file_test(fname,/READ)) then  $
    message, 'Unable to read file ['+fname+']'
    
  ;  Check for invalid dimension values
  if (h.dim[0] ne 4) then  h.dim[0] = 4  ; always 4D
  if (h.dim[3] eq 0) then  h.dim[3] = 1  ; with one slice
  if (h.dim[4] eq 0) then  h.dim[4] = 1  ; and one volume
  if (h.dim[1] eq 0) then  message, 'Invalid header file used, cols == 0.'
  if (h.dim[2] eq 0) then  message, 'Invalid header file used, rows == 0.'
  
  ;  Load according to given orientation
  q = h.dim[1:3]
  case orient of
    0:  ans = self->ReadTRSImageFile(fname, h, endian)
    1:  BEGIN
      ans = self->ReadCORImageFile(fname, h, endian)
      h.dim[2:3] = [q[2],q[1]]
    END
    2:  BEGIN
      ans = self->ReadSAGImageFile(fname, h, endian)
      h.dim[1:3] = [q[2],q[0],q[1]]
    END
    ;
    ;  Rare formats:
    ;
    3:  ans = self->ReadTRSImageFile(fname, h, endian, /FLIP)
    4:  ans = self->ReadCORImageFile(fname, h, endian, /FLIP)
    5:  ans = self->ReadSAGImageFile(fname, h, endian, /FLIP)
    ELSE:  message, 'ORIENT has an illegal value.'
  endcase
  
  ; Return image data
  return, ans
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteSagittalFlipFile
;+
;  Write the data to a single file in sagittal flipped order.
;
;  @param f {in}{type=string}  Output base name.
;
;@private
;-
pro IDLffAnalyze::WriteSagittalFlipFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Need this before bytes get swapped
  ;
  xmax = self.hdr.dim[1]
  nvol = self.hdr.dim[4]
  
  ;
  ;  Swap bytes if necessary
  ;
  swapped = 0
  if (self.endian ne self->SystemByteOrder()) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata   ; and data itself
    swapped = 1
  endif
  
  ;
  ;  XYZ for dim and pixdim reflect how the data is
  ;  stored, therefore, these need to be adjusted as well:
  ;    XYZ -> YZX
  ;
  orig = self.hdr.dim[1:3]
  porig= self.hdr.pixdim[1:3]
  self.hdr.dim[1:3]    = self.hdr.dim[[2,3,1]]
  self.hdr.pixdim[1:3] = self.hdr.pixdim[[2,3,1]]
  
  ;  Header
  openw, u, f+'.hdr', /GET_LUN
  writeu, u, self.hdr
  free_lun, u
  
  ; Image data
  openw, u, f+'.img', /GET_LUN
  
  for v=0, nvol-1 do begin
    for x = 0, xmax-1 do begin
      w = reverse(reform((*self.imgdata)[x,*,*,v]),2)
      writeu, u, w
    endfor
  endfor
  
  free_lun, u
  
  ;
  ;  Re-arrange the dim and pixdim vectors
  ;
  self.hdr.dim[1:3]    = orig
  self.hdr.pixdim[1:3] = porig
  
  ;
  ;  Swap back, if necessary
  ;
  if (swapped eq 1) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
  endif
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteCoronalFlipFile
;+
;  Write the data to a single file in coronal flipped order.
;
;  @param f {in}{type=string}  Output base name.
;
;@private
;-
pro IDLffAnalyze::WriteCoronalFlipFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Need this before bytes get swapped
  ;
  ymax = self.hdr.dim[2]
  nvol = self.hdr.dim[4]
  
  ;
  ;  Swap bytes if necessary
  ;
  swapped = 0
  if (self.endian ne self->SystemByteOrder()) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
    swapped = 1
  endif
  
  ;
  ;  XYZ for dim and pixdim reflect how the data is
  ;  stored, therefore, these need to be adjusted as well:
  ;    XYZ -> XZY
  ;
  self.hdr.dim[2:3]    = self.hdr.dim[[3,2]]
  self.hdr.pixdim[2:3] = self.hdr.pixdim[[3,2]]
  
  ;  Header
  openw, u, f+'.hdr', /GET_LUN
  writeu, u, self.hdr
  free_lun, u
  
  ; Image data
  openw, u, f+'.img', /GET_LUN
  
  for v=0, nvol-1 do begin
    for y = 0, ymax-1 do begin
      w = reverse(reform((*self.imgdata)[*,y,*,v]),2)
      writeu, u, w
    endfor
  endfor
  
  free_lun, u
  
  ;
  ;  Re-arrange the dim and pixdim vectors
  ;
  self.hdr.dim[2:3]    = self.hdr.dim[[3,2]]
  self.hdr.pixdim[2:3] = self.hdr.pixdim[[3,2]]
  
  ;
  ;  Swap back, if necessary
  ;
  if (swapped eq 1) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
  endif
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteTransverseFlipFile
;+
;  Write the data to a single file in transverse flipped order.
;
;  @param f {in}{type=string}  Output base name.
;
;@private
;-
pro IDLffAnalyze::WriteTransverseFlipFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Get limit prior to any byte swap
  ;
  zmax = self.hdr.dim[3]
  nvol = self.hdr.dim[4]
  
  ;
  ;  Swap bytes if necessary
  ;
  swapped = 0
  if (self.endian ne self->SystemByteOrder()) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
    swapped = 1
  endif
  
  ;  Header
  openw, u, f+'.hdr', /GET_LUN
  writeu, u, self.hdr
  free_lun, u
  
  ; Image data
  openw, u, f+'.img', /GET_LUN
  
  ; Output images with Y (rows) flipped
  for v=0, nvol-1 do begin
    for z = 0, zmax-1 do begin
      w = reverse(reform((*self.imgdata)[*,*,z,v]),2)
      writeu, u, w
    endfor
  endfor
  
  free_lun, u
  
  ;
  ;  Swap back, if necessary
  ;
  if (swapped eq 1) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
  endif
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteSagittalFile
;+
;  Write the data to a single file in sagittal order.
;
;  @param f {in}{type=string}  Output base name.
;
;@private
;-
pro IDLffAnalyze::WriteSagittalFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Need this before bytes get swapped
  ;
  xmax = self.hdr.dim[1]
  nvol = self.hdr.dim[4]
  
  ;
  ;  Swap bytes if necessary
  ;
  swapped = 0
  if (self.endian ne self->SystemByteOrder()) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
    swapped = 1
  endif
  
  ;
  ;  XYZ for dim and pixdim reflect how the data is
  ;  stored, therefore, these need to be adjusted as well:
  ;    XYZ -> YZX
  ;
  orig = self.hdr.dim[1:3]
  porig= self.hdr.pixdim[1:3]
  self.hdr.dim[1:3]    = self.hdr.dim[[2,3,1]]
  self.hdr.pixdim[1:3] = self.hdr.pixdim[[2,3,1]]
  
  ;  Header
  openw, u, f+'.hdr', /GET_LUN
  writeu, u, self.hdr
  free_lun, u
  
  ; Image data
  openw, u, f+'.img', /GET_LUN
  
  for v=0, nvol-1 do begin
    for x = 0, xmax-1 do begin
      w = reform((*self.imgdata)[x,*,*,v])
      writeu, u, w
    endfor
  endfor
  
  free_lun, u
  
  ;
  ;  Re-arrange the dim and pixdim vectors
  ;
  self.hdr.dim[1:3]    = orig
  self.hdr.pixdim[1:3] = porig
  
  ;
  ;  Swap back, if necessary
  ;
  if (swapped eq 1) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
  endif
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteCoronalFile
;+
;  Write the data to a single file in coronal order.
;
;  @param f {in}{type=string}  Output base name.
;
;@private
;-
pro IDLffAnalyze::WriteCoronalFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Need this before bytes get swapped
  ;
  ymax = self.hdr.dim[2]
  nvol = self.hdr.dim[4]
  
  ;
  ;  Swap bytes if necessary
  ;
  swapped = 0
  if (self.endian ne self->SystemByteOrder()) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
    swapped = 1
  endif
  
  ;
  ;  XYZ for dim and pixdim reflect how the data is
  ;  stored, therefore, these need to be adjusted as well:
  ;    XYZ -> XZY
  ;
  self.hdr.dim[2:3]    = self.hdr.dim[[3,2]]
  self.hdr.pixdim[2:3] = self.hdr.pixdim[[3,2]]
  
  ;  Header
  openw, u, f+'.hdr', /GET_LUN
  writeu, u, self.hdr
  free_lun, u
  
  ; Image data
  openw, u, f+'.img', /GET_LUN
  
  for v = 0, nvol-1 do begin
    for y = 0, ymax-1 do begin
      w = reform((*self.imgdata)[*,y,*,v])
      writeu, u, w
    endfor
  endfor
  
  free_lun, u
  
  ;
  ;  Re-arrange the dim and pixdim vectors
  ;
  self.hdr.dim[2:3]    = self.hdr.dim[[3,2]]
  self.hdr.pixdim[2:3] = self.hdr.pixdim[[3,2]]
  
  ;
  ;  Swap back, if necessary
  ;
  if (swapped eq 1) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
  endif
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteTransverseFile
;+
;  Write the data to a single file in transverse order.
;  This is the fastest since the data is already stored
;  this way.
;
;  @param f {in}{type=string}  Output base name.
;
;@private
;-
pro IDLffAnalyze::WriteTransverseFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Swap bytes if necessary
  ;
  swapped = 0
  if (self.endian ne self->SystemByteOrder()) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
    swapped = 1
  endif
  
  ;  Header
  openw, u, f+'.hdr', /GET_LUN
  writeu, u, self.hdr
  free_lun, u
  
  ; Image data
  openw, u, f+'.img', /GET_LUN
  writeu, u, *self.imgdata
  free_lun, u
  
  ;
  ;  Swap back, if necessary
  ;
  if (swapped eq 1) then begin
    h = self.hdr
    self->SwapHeader, h
    self.hdr = h
    self->Swap, *self.imgdata
  endif
end


;----------------------------------------------------------
;  IDLffAnalyze::WriteSingleFile
;+
;  Output the current data as a single Analyze file using
;  orient to properly orient the data.
;
;  @param f {in}{type=string}  Output pathname.
;
;@private
;-
pro IDLffAnalyze::WriteSingleFile, f
  compile_opt idl2, HIDDEN
  
  ;
  ;  Respond to the orient field
  ;
  case self.hdr.orient of
    0:  self->WriteTransverseFile, f
    1:  self->WriteCoronalFile, f
    2:  self->WriteSagittalFile, f
    ;
    ;  These are rare, probably shouldn't be supported:
    ;
    3:  self->WriteTransverseFlipFile, f
    4:  self->WriteCoronalFlipFile, f
    5:  self->WriteSagittalFlipFile, f
    ELSE:  message, 'ORIENT has an invalid value.'
  endcase
end


;-----------------------------------------------------------
;  IDLffAnalyze::WriteMultipleFiles
;+
;  Output the current data as a series of individual image
;  files.  This is always done in transverse order (ie, as
;  a set of transverse images) regardless of the value of
;  orient in the header.
;
;  @param f {in}{type=string}  Output base file name.
;
;@private
;-
pro IDLffAnalyze::WriteMultipleFiles, f
  compile_opt idl2, HIDDEN
  
  ;  Set up for byte swapping
  swap = 0
  if (self.endian ne self->SystemByteOrder()) then  $
    swap = 1
    
  ;
  ;  Loop through all axial slices outputing
  ;  each to a new pair of files (.hdr, .img)
  ;  with a header that reflects the updated
  ;  dimensions and proper byte swapping applied.
  ;
  nvol = self.hdr.dim[4]
  for vol=0, nvol-1 do begin
    for slice=0, self.hdr.dim[3]-1 do begin
      ;
      ;  Get the image data and construct an appropriate header
      ;
      img = reform((*self.imgdata)[*,*,slice,vol])
      h = {header}
      h.sizeof_hdr = long(348)
      h.extents    = long(16384)
      h.regular    = byte('r')
      h.dim = self.hdr.dim
      h.dim[3:4] = [1,1]  ; one slice, one volume
      h.pixdim[3:4] = [0.,0.]
      h.orient = 0
      h.glmin = min(img)
      h.glmax = max(img)
      h.bitpix= self.hdr.bitpix
      h.datatype = self.hdr.datatype$
      bname = (nvol eq 1) ? f+strtrim(string(slice, format='(I04)'),2)  $
        : f+strtrim(string(vol),2)+'_'+strtrim(string(slice, format='(I04)'),2)
      hname = bname+'.hdr'
      iname = bname+'.img'
      
      ;
      ;  Swap if necessary
      ;
      if (swap eq 1) then begin
        self->SwapHeader, h
        self->Swap, img
      endif
      
      ;
      ;  Write the header and image data
      ;
      openw, u, hname, /GET_LUN
      writeu, u, h
      free_lun, u
      
      openw, u, iname, /GET_LUN
      writeu, u, img
      free_lun, u
    endfor
  endfor
end



;
;  Public:
;

;-----------------------------------------------------------
;  IDLffAnalyze::GetProperty
;+
;
;-
pro IDLffAnalyze::GetProperty,   $
  ALL=all, DIM=dim, GLMIN=glmin, GLMAX=glmax, PIXDIM=pixdim,    $
  VOX_UNITS=vox_units, CAL_UNITS=cal_units, DATATYPE=datatype,  $
  BITPIX=bitpix, VOX_OFFSET=vox_offset, CAL_MIN=cal_min,        $
  CAL_MAX=cal_max, ORIENT=orient, SCALEFACTOR=scalefactor,      $
  IMAGE_ORIGIN=image_origin, ENDIAN=endian, DCOFF=dcoff
  compile_opt idl2
  
  if arg_present(endian)       then      endian= self.endian
  if arg_present(all)          then         all= self.hdr
  if arg_present(dim)          then         dim= self.hdr.dim
  if arg_present(glmin)        then       glmin= self.hdr.glmin
  if arg_present(glmax)        then       glmax= self.hdr.glmax
  if arg_present(pixdim)       then      pixdim= self.hdr.pixdim
  if arg_present(vox_units)    then   vox_units= self.hdr.vox_units
  if arg_present(cal_units)    then   cal_units= self.hdr.cal_units
  if arg_present(datatype)     then    datatype= self.hdr.datatype
  if arg_present(bitpix)       then      bitpix= self.hdr.bitpix
  if arg_present(vox_offset)   then  vox_offset= self.hdr.vox_offset
  if arg_present(cal_min)      then     cal_min= self.hdr.cal_min
  if arg_present(cal_max)      then     cal_max= self.hdr.cal_max
  if arg_present(orient)       then      orient= self.hdr.orient
  if arg_present(scalefactor)  then scalefactor= self.hdr.funused1
  if arg_present(dcoff)        then       dcoff= self.hdr.funused2
  if arg_present(image_origin) then begin
    w = self.hdr.originator
    if (endian eq 1) then begin
      image_origin = [fix(256*w[0]+w[1]),     $  ; big-endian
        fix(256*w[2]+w[3]),     $
        fix(256*w[4]+w[5])]
    endif else begin
      image_origin = [fix(256*w[1]+w[0]),     $  ; little-endian
        fix(256*w[3]+w[2]),     $
        fix(256*w[5]+w[4])]
    endelse
  endif
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetNCols
;+
;  @returns The number of pixel columns.
;-
function IDLffAnalyze::GetNCols
  compile_opt idl2
  return, self.hdr.dim[1]
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetNRows
;+
;  @returns The number of pixel rows.
;-
function IDLffAnalyze::GetNRows
  compile_opt idl2
  return, self.hdr.dim[2]
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetNSlices
;+
;  @returns The number of slices.
;-
function IDLffAnalyze::GetNSlices
  compile_opt idl2
  return, self.hdr.dim[3]
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetNVolumes
;+
;  @returns The number of volumes or time-points.
;-
function IDLffAnalyze::GetNVolumes
  compile_opt idl2
  return, self.hdr.dim[4]
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetOrientation
;+
;  Return the orientation field value as a string and/or
;  a number.  The orientation string maps XYZ axes, which are
;  always in radiological, left-handed order (ala Analyze) to
;  the final direction of the axis.  For example, normal axial
;  (TRS) images, in radiological order, follow this pattern:
;
;  XYZ -> R-L, P-A, I-S, "LAS" which means that X runs from
;  (patient) right to left, posterior to anterior (back to front)
;  and inferior to superior (feet to head) if the patient is in
;  a supine position (on his back, face up).
;
;  Orient is used to properly place the data read from the file into
;  memory or when the file is written to disk.
;
;  A "-" in front of the returned order string means "flipped".
;  These are rare in the wild and not recommended for use.
;
;  @keyword TEXT {in}{type=boolean}{optional}
;    If set, return the orientation as a string, not a number.
;-
function IDLffAnalyze::GetOrientation, TEXT=text
  compile_opt idl2
  
  if ~keyword_set(text) then begin
    return, self.hdr.orient
  endif else begin
    return, (["TRS","COR","SAG","-TRS","-COR","-SAG"])[self.hdr.orient]
  endelse
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetVoxelDims
;+
;  Get the pixdim dimensions.  These are the actual physical
;  dimensions in mm or ms (assumed units).  Provides a
;  shortcut to GetProperty, PIXDIM=p and removes the need to
;  remember which elements go with which label.
;
;  @keyword VOXEL_WIDTH {out}{type=double}{optional}
;    If present, the voxel width (along readout, for example).
;
;  @keyword VOXEL_HEIGHT {out}{type=double}{optional}
;    The voxel height (eg, along phase axis) in mm (assumed).
;
;  @keyword SLICE_THICKNESS {out}{type=double}{optional}
;    Slice thickness (as opposed to slice spacing) in mm (assumed).
;
;  @keyword TIME_POINTS {out}{type=double}{optional}
;    Time points at which volumes were acquired.  Really number
;    of volumes (corresponding to dim[4]) which would make this
;    perhaps the interval between acquisitions, in ms (assumed).
;-
pro IDLffAnalyze::GetVoxelDims, VOXEL_WIDTH=width, VOXEL_HEIGHT=height,  $
  SLICE_THICKNESS=thick, TIME_POINTS=timep
  compile_opt idl2
  
  if arg_present(width)  then  width = self.hdr.pixdim[1]
  if arg_present(height) then  height= self.hdr.pixdim[2]
  if arg_present(thick)  then  thick = self.hdr.pixdim[3]
  if arg_present(timep)  then  timep = self.hdr.pixdim[4]
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetProperty
;+
;  Allows direct setting of header properties.  The most important
;  of these properties are set by SetPixelData, SetVoxelDimensions
;  or SetBitPix.  With this function it is possible to create
;  an invalid header file that does not reflect the actual data in
;  the .img file.  Caveat emptor!
;-
pro IDLffAnalyze::SetProperty,  $
  ALL=all, DIM=dim, GLMIN=glmin, GLMAX=glmax, PIXDIM=pixdim,    $
  VOX_UNITS=vox_units, CAL_UNITS=cal_units, DATATYPE=datatype,  $
  BITPIX=bitpix, VOX_OFFSET=vox_offset, CAL_MIN=cal_min,        $
  CAL_MAX=cal_max, ORIENT=orient, SCALEFACTOR=scalefactor,      $
  IMAGE_ORIGIN=image_origin, DCOFF =dcoff, ENDIAN=endian
  compile_opt idl2
  
  ;
  ;  Update elements using Analyze names
  ;
  if (n_elements(all) ne 0)          then  self->UpdateAll, all
  if (n_elements(dim) ne 0)          then  self->UpdateDim, dim
  if (n_elements(glmin) ne 0)        then  self->UpdateGlmin, glmin
  if (n_elements(glmax) ne 0)        then  self->UpdateGlmax, glmax
  if (n_elements(pixdim) ne 0)       then  self->UpdatePixdim, pixdim
  if (n_elements(vox_units) ne 0)    then  self->UpdateVoxUnits, vox_units
  if (n_elements(cal_units) ne 0)    then  self->UpdateCalUnits, cal_units
  if (n_elements(datatype) ne 0)     then  self->UpdateDatatype, datatype
  if (n_elements(bitpix) ne 0)       then  self->UpdateBitpix, bitpix
  if (n_elements(vox_offset) ne 0)   then  self->UpdateVoxOffset, vox_offset
  if (n_elements(cal_min) ne 0)      then  self->UpdateCalMin, cal_min
  if (n_elements(cal_max) ne 0)      then  self->UpdateCalMax, cal_max
  if (n_elements(orient) ne 0)       then  self->UpdateOrient, orient
  if (n_elements(scalefactor) ne 0)  then  self->UpdateScaleFactor, scalefactor
  if (n_elements(dcoff) ne 0)        then  self->UpdateDcoff, dcoff
  if (n_elements(image_origin) ne 0) then  self->UpdateImageOrigin, image_origin
  if (n_elements(endian) ne 0)       then  self->UpdateEndian, endian
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetVoxelDims
;+
;  Set the pixdim dimensions.  These are the actual physical
;  dimensions in mm or ms (assumed units).
;
;  @keyword VOXEL_WIDTH {in}{type=double}{optional}
;    The voxel width (eg, along read axis) in mm (assumed).
;
;  @keyword VOXEL_HEIGHT {in}{type=double}{optional}
;    The voxel height (eg, along phase axis) in mm (assumed).
;
;  @keyword SLICE_THICKNESS {in}{type=double}{optional}
;    Slice thickness (as opposed to slice spacing) in mm (assumed).
;
;  @keyword TIME_POINTS {in}{type=double}{optional}
;    Time points at which volumes were acquired.  Really number
;    of volumes (corresponding to dim[4]) which would make this
;    perhaps the interval between acquisitions, in ms (assumed).
;-
pro IDLffAnalyze::SetVoxelDims, VOXEL_WIDTH=width, VOXEL_HEIGHT=height,  $
  SLICE_THICKNESS=thick, TIME_POINTS=timep
  compile_opt idl2
  
  p = self.hdr.pixdim
  
  if (n_elements(width) ne 0) then begin
    if (~self->isNum(width, /POS)) then  $
      message, 'VOXEL_WIDTH must be positive.'
    p[1] = width
  endif
  
  if (n_elements(height) ne 0) then begin
    if (~self->isNum(height, /POS)) then  $
      message, 'VOXEL_HEIGHT must be positive.'
    p[2] = height
  endif
  
  if (n_elements(thick) ne 0) then begin
    if (~self->isNum(thick, /POS)) then  $
      message, 'SLICE_THICK must be positive.'
    p[3] = thick
  endif
  
  if (n_elements(timep) ne 0) then begin
    if (~self->isNum(timep, /POS)) then  $
      message, 'TIME_POINTS must be positive.'
    p[4] = timep
  endif
  
  self.hdr.pixdim = p
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetBitPix
;+
;  Set the bits per pixel field.
;
;  @param bitpix {in}{type=integer}
;    Number of bits per pixel of data.
;-
pro IDLffAnalyze::SetBitPix, bitpix
  compile_opt idl2
  
  ;  Assert that bitpix is a positive number
  if ~self->isNum(bitpix, /POS) then  $
    message, 'Bits per pixel must be a positive integer.'
    
  ;  Must be one of the following
  if where(bitpix eq [1,8,16,32,64]) eq -1 then $
    message, 'Bits per pixel must be 1, 8, 16, 32, or 64.'
    
  ;  Assign a valid number
  self.hdr.bitpix = bitpix
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetOrientation
;+
;  Set the output orientation.  Internally, data is stored and
;  accessed in transverse, LAS order.
;
;  Value  String
;  0      TRS
;  1      COR
;  2      SAG
;  3      -TRS
;  4      -COR
;  5      -SAG
;
;  The orientation string maps XYZ axes, which are
;  always in radiological, left-handed order (ala Analyze) to
;  the final direction of the axis.  For example, normal axial
;  (TRS) images, in radiological order, follow this pattern:
;
;  XYZ -> R-L, P-A, I-S, "LAS" which means that X runs from
;  (patient) right to left, posterior to anterior (back to front)
;  and inferior to superior (feet to head) if the patient is in
;  a supine position (on his back, face up).
;
;  Orient is used to properly place the data read from the file into
;  memory or when the file is written to disk.
;
;  A "-" in front of the order string means "flipped".
;  These are rare in the wild and not recommended for use.
;-
pro IDLffAnalyze::SetOrientation, ord
  compile_opt idl2
  
  ;  Must be a scalar
  if (n_elements(ord) ne 1) then  $
    message, 'ORIENT must be a scalar or single string.'
    
  ;  Can be a string
  if (size(ord,/TYPE) eq 7) then begin
    ans = where(strupcase(ord) eq ['TRS','COR','SAG','-TRS','-COR','-SAG'])
    if (ans  eq -1) then  message, 'Illegal string given for ORIENT.'
    self.hdr.orient = ans
  endif else begin
    if (~self->isNum(ord, LIMITS=[0,5])) then   $
      message, 'ORIENT must be between 0 and 5, inclusive.'
    self.hdr.orient = ord
  endelse
end


;
;  Get/Set pixel data methods:
;

;-----------------------------------------------------------
;  IDLffAnalyze::QueryPixelData
;+
;  @returns True (1) if there is valid image data, false (0)
;       otherwise.
;-
function IDLffAnalyze::QueryPixelData
  return, ptr_valid(self.imgdata)
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetPixelData
;+
;  @returns The pixel data, or query image data status.
;
;  @keyword POINTER {in}{type=boolean}{optional}
;    If set, return a pointer to the image data.  This
;    is against "good" object-oriented design but these
;    files might be very large and direct access through a
;    pointer will not result in copying megabytes of data.
;-
function IDLffAnalyze::GetPixelData, POINTER=pointer
  compile_opt idl2
  
  if ~ptr_valid(self.imgdata) then  $
    message, 'No image data available.'
    
  if ~keyword_set(pointer) then begin
    return, *self.imgdata ; this might require a lot of memory
  endif else begin
    return, self.imgdata  ; memory friendly, but dangerous
  endelse
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetVolume
;+
;  @returns A copy of the desired volume.
;
;  @param v {in}{type=integer}  Volume (time point) number.
;-
function IDLffAnalyze::GetVolume, v
  compile_opt idl2
  
  ;  There must be some data
  if ~ptr_valid(self.imgdata) then  $
    message, 'No image data available.'
    
  ;  Validate the volume number
  if ~self->isNum(v, /POS, /ZERO, /INTEGER) then  $
    message, 'Volume (time point) number must be an integer >=0.'
  if (v ge self.hdr.dim[4]) then  $
    message, 'Volume (time point) index out of bounds.'
    
  return, reform((*self.imgdata)[*,*,*,v])
end


;-----------------------------------------------------------
;  IDLffAnalyze::GetImage
;+
;  @returns  The selected image, from the given volume and
;            orientation.
;
;  @param i {in}{type=integer}
;      The image number, zero-based.
;
;  @param v0 {in}{type=integer}{optional}
;      The source volume (time point), zero-based.  Defaults
;      to 0 if not given.
;
;  @keyword AXIAL {in}{type=boolean}{optional}
;      Select an axial image, default.
;
;  @keyword TRS {in}{type=boolean}{optional}
;      Synonym for AXIAL.
;
;  @keyword TRANSVERSE {in}{type=boolean}{optional}
;      Synonym for AXIAL.
;
;  @keyword CORONAL {in}{type=boolean}{optional}
;      Select a coronal image.
;
;  @keyword SAGITTAL {in}{type=boolean}{optional}
;      Select a sagittal image.
;-
function IDLffAnalyze::GetImage, i, v0,  AXIAL=axial, TRS=trs,  $
  TRANSVERSE=transverse, CORONAL=coronal, SAGITTAL=sagittal
  compile_opt idl2
  
  ;  There must be some image data
  if ~ptr_valid(self.imgdata) then  $
    message, 'No image data available.'
    
  ;  Keywords are mutually exclusive
  z = keyword_set(axial) + 2*keyword_set(trs) + 4*keyword_set(transverse) +  $
    8*keyword_set(coronal) + 16*keyword_set(sagittal)
    
  if (where(z eq [0,1,2,4,8,16]) eq -1) then  $
    message, 'Incompatible keywords specified.'
    
  ;  Volume (time point) number must exist
  v = (n_params() eq 2) ? v0 : 0
  if (v ge self.hdr.dim[4]) then  $
    message, 'Volume (time point) out of range.'
    
  ;  Get the desired image
  msg = 'Image number out of range.'
  
  case z of
    16: BEGIN
      if (i ge self.hdr.dim[1]) then  message, msg  ; sagittal
      img = reform((*self.imgdata)[i,*,*,v])
    END
    8: BEGIN
      if (i ge self.hdr.dim[2]) then  message, msg  ; coronal
      img = reform((*self.imgdata)[*,i,*,v])
    END
    ELSE: BEGIN
      if (i ge self.hdr.dim[3]) then  message, msg  ; transverse
      img = (*self.imgdata)[*,*,i,v]
    END
  endcase
  
  return, img
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetPixelData
;+
;  Set the pixel data and associated header elements using
;  the input array.  Will set glmin, glmax, datatype, bitpix,
;  and dim but does not set orient.
;
;  @param d {in}{type=numeric array, 2, 3, or 4D}
;    The data to store.
;
;  @keyword NO_COPY {in}{type=boolean}{optional}
;    If set the original array input will be destroyed
;    once the data has been assigned to the heap.  Ie,
;    the NO_COPY keyword to ptr_new().
;-
pro IDLffAnalyze::SetPixelData, d, NO_COPY=no_copy0
  compile_opt idl2
  
  no_copy = keyword_set(no_copy0) ? no_copy0 : 0
  
  ;
  ;  Must be a numeric array that is supported by the
  ;  Analyze 7.5 "standard".
  ;
  if (where(size(d,/TYPE) eq [1,2,3,4,5,6]) eq -1) then  $
    message, 'Input array not of a type supported by Analyze 7.5.'
    
  ;  Must be 2, 3, or 4 dimensional
  dims = size(d,/DIM)
  ndims = n_elements(dims)
  if (ndims lt 2) or (ndims gt 4) then  $
    message, 'Input array must be of 2, 3 or 4 dimensions.'
    
  ;  Remove any existing data
  if ptr_valid(self.imgdata) then  $
    ptr_free, self.imgdata
    
  ;
  ;  Reset the header structure to default values, including orient
  ;  (all data set this way is "transverse")
  ;
  self.hdr = self.default
  
  ;  Set the dimensions
  self.hdr.dim[0] = 4                           ; always say it is 4-dimensional
  self.hdr.dim[1] = dims[0]                     ; X direction (columns)
  self.hdr.dim[2] = dims[1]                     ; Y direction (rows)
  self.hdr.dim[3] = (ndims ge 3) ? dims[2] : 1  ; Z direction (slices), 2D is one slice
  self.hdr.dim[4] = (ndims eq 4) ? dims[3] : 1  ; volumes (time points), 3D is one time point
  
  ;  Set the datatype
  self.hdr.datatype = self->IDLtoAnalyzeType(size(d,/TYPE))
  
  ;  Set bits per pixel
  self.hdr.bitpix = self->IDLtoAnalyzeBitPix(size(d,/TYPE))
  
  ;  Set the glmin and glmax
  self.hdr.glmin = min(d)
  self.hdr.glmax = max(d)
  
  ;  Set endian to the order this system uses
  self.endian = self->SystemByteOrder()
  
  ;  Set the data itself
  self.imgdata = ptr_new(d, NO_COPY=no_copy)
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetVolume
;+
;  Update an existing volume (time point).  Datatype and
;  dimensions must match the existing volume.
;
;  @param v {in}{type=integer}
;    The volume (time point) number.
;
;  @param vol {in}{type=3D array}
;    The new volume data.
;-
pro IDLffAnalyze::SetVolume, v, vol
  compile_opt idl2
  
  ;  Is there any data?
  if ~ptr_valid(self.imgdata) then  $
    message, 'No image data available.'
    
  ;  Validate v, volume number
  if ~self->isNum(v, /POS, /ZERO, /INTEGER) then  $
    message, 'Volume (time point) number must be an integer >=0.'
  if (v ge self.hdr.dim[4]) then  $
    message, 'Volume (time point) index out of bounds.'
    
  ;  Validate vol, the new volume data
  if n_elements(size(vol,/DIM)) ne 3 then  $
    message, 'Input volume must be 3 dimensions.'
    
  ;  What about the datatype?
  if (self->IDLtoAnalyzeType(size(vol,/TYPE)) ne self.hdr.datatype) then  $
    message, 'Data type of input must be the same as existing image data.'
    
  ;  Now for the dimensions
  if ~array_equal(size(vol,/DIM),self.hdr.dim[1:3]) then  $
    message, 'Input volume must match the dimensions of the existing image data.'
    
  ;  Assign the new image data
  (*self.imgdata)[*,*,*,v] = vol
end


;-----------------------------------------------------------
;  IDLffAnalyze::SetImage
;+
;  Update an existing image in a particular volume with a
;  particular orientation.  Dimensions and datatype must
;  match (exactly, no behind-the-scenes promotion of types).
;
;  @param i {in}{type=integer}
;    Image number to update.
;
;  @param v {in}{type=integer}
;    Volume which holds the image.
;
;  @param img {in}{type=matrix}
;    The new image data.
;
;  @keyword AXIAL {in}{type=boolean}{optional}
;      Select an axial image, default.
;
;  @keyword TRS {in}{type=boolean}{optional}
;      Synonym for AXIAL.
;
;  @keyword TRANSVERSE {in}{type=boolean}{optional}
;      Synonym for AXIAL.
;
;  @keyword CORONAL {in}{type=boolean}{optional}
;      Select a coronal image.
;
;  @keyword SAGITTAL {in}{type=boolean}{optional}
;      Select a sagittal image.
;-
pro IDLffAnalyze::SetImage, i, v, img, AXIAL=axial,  TRS=trs,  $
  TRANSVERSE=transverse, CORONAL=coronal, SAGITTAL=sagittal
  compile_opt idl2
  
  ;  There must be some image data
  if ~ptr_valid(self.imgdata) then  $
    message, 'No image data available.'
    
  ;  Keywords are mutually exclusive
  z = keyword_set(axial) + 2*keyword_set(trs) + 4*keyword_set(transverse) +  $
    8*keyword_set(coronal) + 16*keyword_set(sagittal)
    
  if (where(z eq [0,1,2,4,8,16]) eq -1) then  $
    message, 'Incompatible keywords specified.'
    
  ;  Volume (time point) number must exist
  if ~self->isNum(v, /POS, /ZERO, /INTEGER) then  $
    message, 'Volume (time point) number must be an integer >=0.'
  if (v ge self.hdr.dim[4]) then  $
    message, 'Volume (time point) out of range.'
    
  ;  Check the datatype
  if (self->IDLtoAnalyzeType(size(img,/TYPE)) ne self.hdr.datatype) then  $
    message, 'Data type of input must be the same as existing image data.'
    
  ;  Assign the desired image
  msg = 'Image number out of range.'
  msg2 = 'Input image dimensions do not match existing data.'
  d = size(img,/DIM)
  
  case z of
    16: BEGIN
      if (i ge self.hdr.dim[1]) then  message, msg  ; sagittal
      if (d[0] ne self.hdr.dim[2]) and   $
        (d[1] ne self.hdr.dim[3]) then  message, msg2
      (*self.imgdata)[i,*,*,v] = img
    END
    8: BEGIN
      if (i ge self.hdr.dim[2]) then  message, msg  ; coronal
      if (d[0] ne self.hdr.dim[1]) and   $
        (d[1] ne self.hdr.dim[3]) then  message, msg2
      (*self.imgdata)[*,i,*,v] = img
    END
    ELSE: BEGIN
      if (i ge self.hdr.dim[3]) then  message, msg  ; transverse
      if (d[0] ne self.hdr.dim[1]) and   $
        (d[1] ne self.hdr.dim[2]) then  message, msg2
      (*self.imgdata)[*,*,i,v] = img
    END
  endcase
end


;
;  File I/O methods:
;

;----------------------------------------------------------
;  IDLffAnalyze::ReadHeaderFile
;+
;  @returns  An instance of {header} for the given .hdr file.
;
;  @param f {in}{type=string}
;    Pathname to the header file.
;
;  @keyword ENDIAN {out}{type=integer}
;    Set to the byte order of the header (1=big, 0=little).
;-
function IDLffAnalyze::ReadHeaderFile, f, ENDIAN=endian
  compile_opt idl2
  
  ;  Filename must really exist with read permission
    if (~file_test(f,/READ)) then  $
    message, 'Unable to read file ['+f+']'
    
  ;  File must be 384 bytes long, no extended headers
  info = file_info(f)
  if (info.size ne 348) then  $
    message, 'Non-standard header files not supported, length > 384 bytes.'
    
  ;  Read the structure.
  ans = {header}
  openr, u, f, /GET_LUN
  readu, u, ans
  free_lun, u
  
  ;
  ;  Determine file byte order and set local ENDIAN property.
  ;  Used by ReadImageFile() to swap bytes properly.
  ;
  if (ans.sizeof_hdr eq 348) then begin
    ; system and header using same byte order
    endian = self->SystemByteOrder()
  endif else begin
    ; system and header using opposite byte order
    endian = ~self->SystemByteOrder()  ; store true order
    self->SwapHeader, ans  ; swap to get to local system order
  endelse
  
  ;  Return the structure
  return, ans
end


;-----------------------------------------------------------
;  IDLffAnalyze::ReadFile
;+
;  Reads an Analyze 7.5 file.
;
;  @param f {in}{type=string vector}{optional}
;    Name or names of the files to load.  If not given
;    an open dialog box is used.
;
;  @keyword ORIENTATION {in}{type=string or integer}{optional}
;    Orientation to use when reading the file, overrides any
;    specified in the header.  If not given use the value in
;    the header.  Integer 0..5 or 'TRS','COR','SAG' with '-'
;    in front for flipped (values 3..5).
;
;  @keyword _EXTRA {in}{type=any}{optional}
;    Extra keywords for dialog_pickfile.
;-
pro IDLffAnalyze::ReadFile, f, ORIENTATION=orient0, _EXTRA=extra
  compile_opt idl2
  
  ;
  ;  Set up the orientation, 0..5 or -1 for header
  ;
  orient = -1
  if (n_elements(orient0) ne 0) then begin
    if (n_elements(orient0) ne 1) then  $
      message, 'Orientation must be a string or integer scalar.'
    if (size(orient0,/TYPE) eq 7) then begin
      idx = where(strupcase(orient0) eq ['TRS','COR','SAG','-TRS','-COR','-SAG'])
      if (idx eq -1) then  $
        message, 'Orientation must be a valid string.'
      orient = idx[0]
    endif else begin
      if (~self->isNum(orient0, /INTEGER, LIMITS=[0,5])) then  $
        message, 'Orientation must be an integer, [0,5] inclusive.'
      orient = fix(orient0)
    endelse
  endif
  
  ; If f given?
  if (n_params() eq 0) then begin
    f = dialog_pickfile(/MUST_EXIST, FILTER='*.hdr', _EXTRA=extra)
    if f[0] eq '' then return  ; user cancel
    f = f[0]
  endif else begin
    ;  Validate f
    if (size(f,/TYPE) ne 7) then  $
      message, 'Argument must be a string.'
    if (n_elements(f) ne 1) then  $
      message, 'Argument must be a single string.'
  endelse
  
  ;  Read the header file
  h = self->ReadHeaderFile(f, ENDIAN=endian)
  
  ;  Set orientation properly
  if (orient eq -1) then orient = h.orient
  
  ;  Now read the image file
  img = self->ReadImageFile(f, h, orient, ENDIAN=endian)
  
  ;  All read successfully, update necessary fields
  if ptr_valid(self.imgdata) then  ptr_free, self.imgdata
  self.imgdata = ptr_new(img,/NO_COPY)
  self.hdr = h
  self.hdr.orient = orient  ; orientation actually used
  self.endian = endian
end


;-----------------------------------------------------------
;  IDLffAnalyze::WriteFile
;+
;  Output the image data to a file or set of files.
;
;  @param f {in}{type=string}{optional}
;    The output file name.  Base name only, no extensions.
;
;  @keyword MULTIPLE_FILES {in}{type=boolean}{optional}
;    If set will create multiple .hdr/.img pairs for each
;    image in the data set using axial ordering (ie, axial
;    images, ignores the orient field).
;
;  @keyword _EXTRA {in}{type=any}{optional}
;    Extra keywords for dialog_pickfile.
;-
pro IDLffAnalyze::WriteFile, f, MULTIPLE_FILES=mult, _EXTRA=extra
  compile_opt idl2
  
  ; Must have something to write
  if ~ptr_valid(self.imgdata) then  $
    message, 'No image data available.'
    
  ;
  ;  Get a filename
  ;
  if (n_params() eq 1) then begin
    if (size(f,/TYPE) ne 7) then  $  ; f given, validate it
      message, 'Output filename must be a string.'
    if (n_elements(f) ne 1) then  $
      message, 'Output filename must be a single string.'
  endif else begin
    ; f not present, use a dialog box
    if n_elements(dp) ne 0 then begin
      f = dialog_pickfile(_EXTRA=extra)
    endif else begin
      f = dialog_pickfile(_EXTRA=extra)
    endelse
    if f[0] eq '' then return
    f = f[0]  ; dialog_pickfile returns a vector
  endelse
  
  ;
  ;  Decide whether or not to output a single file or
  ;  multiple files.
  ;
  if ~keyword_set(mult) then begin
    self->WriteSingleFile, f
  endif else begin
    self->WriteMultipleFiles, f
  endelse
end



;
;  Lifecycle methods:
;

;-----------------------------------------------------------
;  IDLffAnalyze::Init
;+
;  Create an instance of this class.
;
;  @returns  True (1) if the class successfully created.
;
;  @param f {in}{type=string}{optional}
;    If present, must be a string giving the name of a
;    header file.  This file (and its associated .img
;    file which must be in the same directory) are loaded.
;
;  @keyword ORIENTATION {in}{type=integer or string}{optional}
;    Orientation to use when reading the data, if f is given.
;
;  @keyword DATA {in}{type=array}{optional}
;    If present, the data from this array is used to
;    set up the object.  This overrides any given argument.
;-
function IDLffAnalyze::Init, f, DATA=data, ORIENTATION=orient
  compile_opt idl2
  
  ;
  ;  Catch any errors so memory isn't wasted
  ;
  CATCH, err
  if (err ne 0) then begin
    catch, /CANCEL
    self->Destroy  ; release any allocated memory
    return, 0      ; object creation failed
  endif
  
  ;
  ;  Set certain constant values, save as default
  ;
  self.hdr.sizeof_hdr = long(348)
  self.hdr.extents    = long(16384)
  self.hdr.regular    = byte('r')
  
  self.default = self.hdr
  
  ;  Use the system byte order by default
  self.endian = self->SystemByteOrder()
  
  if keyword_set(data) then begin
    self->SetPixelData, data
  endif else begin
    if (n_elements(f) eq 0) then begin
      ; no arguments at all...
    endif else begin
      self->ReadFile, f, ORIENTATION=orient
    end
  endelse
  
  ;  have a valid object
  return, 1
end


;-----------------------------------------------------------
;  IDLffAnalyze::Cleanup
;+
;  Destory the object.
;-
pro IDLffAnalyze::Cleanup
  compile_opt idl2
  
  self->Destroy  ; separate destroy method so it can be
  ; called from Init if necessary
end


;-----------------------------------------------------------
;  IDLffAnalyze__define
;+
;  Define the class and other structures.
;-
pro idlffanalyze__define
  compile_opt idl2  ;  N.B. matters here!
  
  ;
  ;  These structures define the contents of the .hdr file:
  ;
  
  ;
  ;  Analyze 7.5 header structure - incorporates the
  ;    header_key, image_dimension, and data_history structures.
  ;
  void  = {header,                         $ ; struct header_key {
    sizeof_hdr : 348,             $ ;   int sizeof_hdr;
    data_type  : bytarr(10),      $ ;   char data_type[10];
    db_name    : bytarr(18),      $ ;   char db_name[18];
    extents    : 16384,           $ ;   int extents;
    session_error : 0s,           $ ;   short int session_error;
    regular    : byte('r'),       $ ;   char regular;
    hkey_un0   : 0b,              $ ;   char hkey_un0; }
    $ ;
    $ ;  Analyze 7.5 image_dimension structure:
    $ ;                                        ; struct image_dimension {
    dim : fix(intarr(8),TYPE=2),  $ ;   short int dim[8];
    vox_units : bytarr(4),        $ ;   char vox_units[4];
    cal_units : bytarr(8),        $ ;   char cal_units[4];
    unused1   : 0s,               $ ;   short int unused1;
    datatype  : 0s,               $ ;   short int datatype;
    bitpix    : 0s,               $ ;   short int bitpix;
    dim_un0   : 0s,               $ ;   short int dim_un0;
    pixdim    : fltarr(8),        $ ;   float pixdim[8];
    vox_offset: FLOAT(0.0),       $ ;   float vox_offset;
    roi_scale : FLOAT(0.0),       $ ;   float roi_scale;
    funused1  : FLOAT(0.0),       $ ;   float funused1;  // scalefactor in SPM
    funused2  : FLOAT(0.0),       $ ;   float funused2;  // dcoff in SPM2
    cal_max   : FLOAT(0.0),       $ ;   float cal_max;
    cal_min   : FLOAT(0.0),       $ ;   float cal_min;
    compressed: 0,                $ ;   int compressed;
    verified  : 0,                $ ;   int verified;
    glmax     : 0,                $ ;   int glmax;
    glmin     : 0,                $ ;   int glmin; }
    $ ;
    $ ;  Analyze 7.5 data_history structure:
    $ ;                                        ; struct data_history {
    descrip    : bytarr(80),      $ ;   char descrip[80];
    aux_file   : bytarr(24),      $ ;   char aux_file[24];
    orient     : 0b,              $ ;   char orient;
    originator : bytarr(10),      $ ;   char originator[10]; // image_origin in SPM
    generated  : bytarr(10),      $ ;   char generated[10];
    scannum    : bytarr(10),      $ ;   char scannum[10];
    patient_id : bytarr(10),      $ ;   char patient_id[10];
    exp_date   : bytarr(10),      $ ;   char exp_date[10];
    exp_time   : bytarr(10),      $ ;   char exp_time[10];
    hist_un0   : bytarr(3),       $ ;   char hist_un0[3];
    views      : 0,               $ ;   int views;
    vols_added : 0,               $ ;   int vols_added;
    start_field: 0,               $ ;   int start_field;
    field_skip : 0,               $ ;   int field_skip;
    omax       : 0,               $ ;   int omax;
    omin       : 0,               $ ;   int omin;
    smax       : 0,               $ ;   int smax;
    smin       : 0                $ ;   int smin;
  }                                  ; }
  
  ;
  ;  Class instance vars:
  ;
  class = { IDLffAnalyze,            $
    endian : 0,              $ ;  1=big-endian, 0=little-endian (read/write)
    hdr    : {header},       $ ;  header data
    default: {header},       $ ;  default header values
    imgdata: ptr_new()       $ ;  actual image data
  }
end


;
;  End idlffanalyze__define.pro
;
