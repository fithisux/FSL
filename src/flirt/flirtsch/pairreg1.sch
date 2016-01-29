# 8mm scale
setscale 8
setoption smoothing 8
clear S
clear P
search
# 4mm scale
setscale 4
setoption smoothing 4
clear U
clear UA 
clear UB
clear US
clear UP
# remeasure costs at this scale
measurecost 7 S 0 0 0 0 0 0 rel
copy U US
clear U
measurecost 7 P 0 0 0 0 0 0 rel
copy U UP
dualsort US UP
# optimise best 3 candidates (pre and post 8mm optimisations)
clear U
optimise 7 US:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UP:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
# also try the identity transform as a starting point at this resolution
clear UQ
setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1
optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
sort U
copy U UA
# select best 4 optimised solutions and try perturbations of these
clear U
copy UA:1-4 U
optimise 7 UA:1-4  1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4 -1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   1.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0  -1.0   0.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   0.0   1.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   0.0  -1.0   0.0   0.0   0.0   0.0  rel 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0   0.1  abs 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0  -0.1  abs 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0   0.2  abs 4
optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0  -0.2  abs 4
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
