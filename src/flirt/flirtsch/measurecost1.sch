# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
clear U
setrow UA  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
measurecost 7 UA 0 0 0 0 0 0 abs
printparams U

