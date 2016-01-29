- Explanation of possible Output files

a) If --savecoeff, then only the samples of the coefficients will be saved and no peaks. This is very fast, but the output files can be huge.

b) If --savecoeff --savemeancoeff, then only the mean ODF coefficients (across the samples) will be saved for each voxel and again no peaks. This is very fast and also requires less memory/storage.

c) If --savemeancoeff, then the mean ODF coefficients (across the samples) will be saved for each voxel, along with the ODF peaks.

d) If none of the above is set, only the ODF peaks are saved.






- Multi-shell data assumptions

The current implementation of qboot can estimate multi-shell ODFs, but assume the following for the data:

a) Three-shell data are assumed. The bvalues should form an arithmetic progression, e.g. 1000, 2000, 3000 or 2000, 4000, 6000.

b) There is no assumption on the number of directions in each shell. Each shell can have its own directions and number of directions. The minimum number of directions in a shell is dictated by the number of coefficients estimated, i.e. if lmax=4 => Num_of_coeff=15 => All shells should have at least 15 directions. This is because interpolation of the data occurs, when: i) different number of directions exist in different shells or ii) same number of directions, but corresponding directions between shells (e.g. dir1_shell1 vs dir1_shell2) differ more than 1 degree. In this case each shell is fitted with SH of the same order, which is by default 10, unless the number of directions in the shell with the less points is not enough (then order is decreased).

c) It is assumed that all data from each shell are grouped together and shells are one after the other in data, bvels and bvals. So, if 3 directions are available for 3 shells these should appear as:
    dir1_shell1, dir2_shell1, dir3_shell1,     dir1_shell2,dir2_shell2, dir3_shell2,      dir1_shell3,dir2_shell3, dir3_shell3  

 

