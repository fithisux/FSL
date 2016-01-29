fslsurface.o: fslsurface.cc fslsurface.h fslsurface_structs.h
fslsurfacefns.o: fslsurfacefns.cc fslsurfacefns.h fslsurface_structs.h \
 fslsurface.h fslsurface_structs.h
fslsurfacegl.o: fslsurfacegl.cc fslsurfacegl.h fslsurface.h \
 fslsurface_structs.h
fslsurfaceio.o: fslsurfaceio.cc fslsurface.h fslsurface_structs.h \
 fslsurface_structs.h fslsurfaceio.h
run_test.o: run_test.cc fslsurface.h fslsurface_structs.h
