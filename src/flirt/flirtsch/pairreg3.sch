# 2mm scale
setscale 2
setoption smoothing 2
setoption paramsubset 6  1 0 0 0 0 0 0 0 0 0 0 0  0 1 0 0 0 0 0 0 0 0 0 0  0 0 1 0 0 0 0 0 0 0 0 0  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0  0 0 0 0 0 1 0 0 0 0 0 0 
clear U
clear UC
clear UD
clear UE
clear UF
# the line below is the dummy identity - to be useful must have an init
setrow UC 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
copy U UD
setoption boundguess 1
optimise 12 UD:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2
copy U UE
optimise 12 UE:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
sort U
copy U UF
# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
clear U
# also try the identity transform as a starting point at this resolution
clear UG
copy UF:1 UG
setrow UG  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UG:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
sort U

