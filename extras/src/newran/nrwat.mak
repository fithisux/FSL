everything:   tryrand.exe


.cpp.obj:
	      wpp386/bt=nt  $*.cpp

OBJ_T = tryrand.obj tryrand1.obj tryrand2.obj tryrand3.obj &
  tryrand4.obj hist.obj newran.obj extreal.obj myexcept.obj


tryrand.exe:    $(OBJ_T)
	        echo file tryrand.obj,tryrand1.obj,tryrand2.obj   > link.lnk
	        echo file tryrand3.obj,tryrand4.obj,hist.obj     >> link.lnk
	        echo file newran.obj,extreal.obj,myexcept.obj    >> link.lnk
	        echo name tryrand.exe                            >> link.lnk
		echo SYSTEM nt                                   >> link.lnk
	        wlink @link.lnk


newranxx = include.h newran.h boolean.h myexcept.h

tryrand.obj:    $(newranxx) tryrand.cpp 

tryrand1.obj:   $(newranxx) tryrand1.cpp 

tryrand2.obj:   $(newranxx) tryrand2.cpp 

tryrand3.obj:   $(newranxx) tryrand3.cpp 

tryrand4.obj:   $(newranxx) tryrand4.cpp 

hist.obj:       $(newranxx) hist.cpp 

newran.obj:     $(newranxx) newran.cpp 

extreal.obj:    include.h boolean.h extreal.h myexcept.h extreal.cpp 

myexcept.obj:   include.h boolean.h myexcept.h myexcept.cpp 
