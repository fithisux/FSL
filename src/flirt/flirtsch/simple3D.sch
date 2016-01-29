# 8mm scale
setscale 8
setoption smoothing 8
clear UU
clear U
setrow UU 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UU:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
copy U UA
# 4mm scale
setscale 4
setoption smoothing 4
clear U
clear UB
# optimise best 1 candidate (pre and post 8mm optimisations)
clear U
optimise 7 UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
sort U
copy U UB
# 2mm scale
setscale 2
setoption smoothing 2
clear U
clear UC
clear UD
clear UE
clear UF
# remeasure costs at this scale
measurecost 7 UB 0 0 0 0 0 0 rel
sort U
copy U UC
clear U
optimise 7  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
copy U UD
setoption boundguess 1
if MAXDOF > 7
 clear U
if MAXDOF > 7
 optimise 9  UD:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
copy U UE
if MAXDOF > 9
 clear U
if MAXDOF > 9
 optimise 12 UE:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2
sort U
copy U UF
# 1mm scale
setscale 1
setoption smoothing 1
setoption boundguess 1
clear U
# also try the identity transform as a starting point at this resolution
setrow UF  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1
sort U

