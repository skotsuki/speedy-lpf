dset ./Mask_GLB1.00.bin
title land/sea mask (global)
UNDEF -9999
 XDEF          360 LINEAR   -179.5000       1.000000
 YDEF          150 LINEAR   -59.50000       1.000000    
zdef 1 linear 1 1 
tdef 1 LINEAR 00:00Z01Jan1971 1yr
vars 1
mask 0 99 (0:sea  1:land  2:inland water)
endvars
