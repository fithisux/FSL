# $Id: rules.mk,v 1.30 2006/03/27 08:48:16 steve Exp $
all: 

help:
	@echo " make                Rebuild project targets";
	@echo " make debug          Rebuild targets with debugging options";
	@echo " make commit         Send changes to the CVS repository";
	@echo " make update         Download changes that others may have commited";
	@echo "                     to your working copy, ie. freshen your copy";
	@echo " make tagstable      Tag the project as stable in the CVS repository";
	@echo " make stable-branch  Create and tag a branch as stable for this project in the CVS repository";
	@echo "";
	@echo " make install        Install into your local FSLDEVDIR"; 
	@echo " make install-fmrib  Install includes, lib, scripts and binaries to FSLDIR";

debug:
	${MAKE} all "DBGFLAGS=${MACHDBGFLAGS}" "OPTFLAGS="

nodebug:
	${MAKE} all "DBGFLAGS="

static:
	${MAKE} all "USEDCSTATICFLAGS=${CSTATICFLAGS}" "USEDCXXSTATICFLAGS=${CXXSTATICFLAGS}"

initprj:
	@${MAKE} help

distclean: clean
	${RM} -f /tmp/fslgrot ${XFILES} ${TESTXFILES}

clean:
	${RM} -f /tmp/fslgrot *.o *.a *.exe core depend.mk

install:
	@${MAKE} "DESTDIR=${FSLDEVDIR}" master-install-script

install-fmrib: 
	@numunst=`cvs status | grep "Sticky Tag" | grep -v -i stable | wc -l` ; export numunst ; \
	if [ $$numunst -ge 1 ] ; then \
	  echo " " ; \
	  echo "ERROR :: Project does not have stable tag" ; \
	  echo "Can only use install-fmrib after: cvs co -r stable projname" ; \
	else  \
	nummod=`cvs status | grep Status | grep -i modified | wc -l` ; export nummod ; \
	if [ $$nummod -ge 1 ] ; then \
	  echo "Cannot install as some files are locally modified!" ; \
	  echo "All files must be tagged stable and checked in (cvs ci) first" ; \
	 else \
	  ${MAKE} "DESTDIR=${FSLDIR}" "FSLDEVDIR=${FSLDIR}" master-install-script ; \
	fi ; fi

install-stable:
	@${MAKE} help

master-install-script:
	@if [ "X${PROJNAME}X" = XX ] ; then \
		echo " " ; \
		echo "No PROJNAME defined in the Makefile" ; \
		echo "    ... aborting install" ; \
		echo " " ; \
		exit 4 ; \
	fi;
	@${MAKE} all
	@${MAKE} install-non-src-files

install-non-src-files:
	@${MAKE} exeinstall
	@if [ "${FSLINFMRIB}" = "YES" ] ; then \
		${MAKE} XFILES="${FXFILES}" SCRIPTS="${FSCRIPTS}" RUNTCLS="" RUNAVWS="" exeinstall ; \
	fi;
	@${MAKE} install-non-src-exe-files

install-non-src-exe-files:
	@${MAKE} hdrinstall
	@${MAKE} libinstall
	@${MAKE} tclinstall
	@${MAKE} docinstall
	@${MAKE} refdocinstall
	@${MAKE} postinstallscript

postinstallscript:
	@echo " "


maninstall:
	@echo " "
	@echo "Installing man pages"
	@if [ ! -d ${dest_MANDIR} ] ; then \
		${MKDIR} ${dest_MANDIR} ; \
		${CHMOD} g+w ${dest_MANDIR} ; \
		${MKDIR} ${dest_MANDIR}/man1 ; \
		${CHMOD} g+w ${dest_MANDIR}/man1 ; \
		${MKDIR} ${dest_MANDIR}/man2 ; \
		${CHMOD} g+w ${dest_MANDIR}/man2 ; \
		${MKDIR} ${dest_MANDIR}/man3 ; \
		${CHMOD} g+w ${dest_MANDIR}/man3 ; \
		${MKDIR} ${dest_MANDIR}/man4 ; \
		${CHMOD} g+w ${dest_MANDIR}/man4 ; \
		${MKDIR} ${dest_MANDIR}/man5 ; \
		${CHMOD} g+w ${dest_MANDIR}/man5 ; \
		${MKDIR} ${dest_MANDIR}/man6 ; \
		${CHMOD} g+w ${dest_MANDIR}/man6 ; \
		${MKDIR} ${dest_MANDIR}/man7 ; \
		${CHMOD} g+w ${dest_MANDIR}/man7 ; \
		${MKDIR} ${dest_MANDIR}/man8 ; \
		${CHMOD} g+w ${dest_MANDIR}/man8 ; \
		${MKDIR} ${dest_MANDIR}/mann ; \
		${CHMOD} g+w ${dest_MANDIR}/mann ; \
		${MKDIR} ${dest_MANDIR}/manl ; \
		${CHMOD} g+w ${dest_MANDIR}/manl ; \
		fi
	@for mansec in n l 1 2 3 4 5 6 7 8 ; do \
	   for manfile in doc/*.$$mansec refdoc/*.$$mansec verylongdummyname ; do \
		if [ -f $$manfile ] ; then \
			${INSTALL} -m 0664 $$manfile ${dest_MANDIR}/man$$mansec ; \
			echo ${INSTALL} -m 0664 $$manfile ${dest_MANDIR}/man$$mansec ; \
		fi \
	   done \
	done
	@if ${CHMOD} -R ug+w ${dest_MANDIR} ; then \
		echo "${CHMOD} -R ug+w ${dest_MANDIR}" ; \
	fi

docinstall:
	@echo " "
	@echo "Installing (userguide) documents"
	@if [ ! -d ${dest_DOCDIR} ] ; then \
		${MKDIR} ${dest_DOCDIR} ; \
		${CHMOD} g+w ${dest_DOCDIR} ; \
		fi
	@${RM} -r -f /tmp/fslgrot ${dest_DOCDIR}/${PROJNAME}
	@${MKDIR} ${dest_DOCDIR}/${PROJNAME}
	@${CHMOD} g+w ${dest_DOCDIR}/${PROJNAME}
	@for docfile in doc/* verylongdummyname ; do \
		if [ -f $$docfile -o -d $$docfile ] ; then \
			${CP} -R -f $$docfile ${dest_DOCDIR}/${PROJNAME}/ ; \
			echo ${CP} -R -f $$docfile ${dest_DOCDIR}/${PROJNAME}/ ; \
		fi \
	done
	@if ${CHMOD} -R ug+w ${dest_DOCDIR}/${PROJNAME}/ ; then \
		echo "${CHMOD} -R ug+w ${dest_DOCDIR}/${PROJNAME}/" ; \
	fi

refdocinstall:
	@echo " "
	@echo "Installing reference documents"
	@if [ ! -d ${dest_REFDOCDIR} ] ; then \
		${MKDIR} ${dest_REFDOCDIR} ; \
		${CHMOD} g+w ${dest_REFDOCDIR} ; \
		fi
	@${RM} -r -f /tmp/fslgrot ${dest_REFDOCDIR}/${PROJNAME}
	@${MKDIR} ${dest_REFDOCDIR}/${PROJNAME}
	@${CHMOD} g+w ${dest_REFDOCDIR}/${PROJNAME}
	@for refdocfile in refdoc/* verylongdummyname ; do \
		if [ -f $$refdocfile -o -d $$refdocfile ] ; then \
			${CP} -R -f $$refdocfile ${dest_REFDOCDIR}/${PROJNAME}/ ; \
			echo ${CP} -R -f $$refdocfile ${dest_REFDOCDIR}/${PROJNAME}/ ; \
		fi \
	done
	@if ${CHMOD} -R ug+w ${dest_REFDOCDIR}/${PROJNAME}/ ; then \
		echo "${CHMOD} -R ug+w ${dest_REFDOCDIR}/${PROJNAME}/" ; \
	fi

xml2doc:
	@echo " "
	@echo "Making documents from xml source"
	@for xmlfile in doc/*.xml refdoc/*.xml verylongdummyname ; do \
	   if [ -f $$xmlfile ] ; then \
		${FSLML} -h $$xmlfile ; \
		${FSLML} -t $$xmlfile ; \
		${FSLML} -m $$xmlfile ; \
		${FSLML} -l $$xmlfile ; \
		echo "${FSLML} $$xmlfile" ; \
	   fi \
	done

tclinstall:
	@echo " "
	@echo "Installing tcl scripts"
	@if [ ! -d ${dest_TCLDIR} ] ; then \
		${MKDIR} ${dest_TCLDIR} ; \
		${CHMOD} g+w ${dest_TCLDIR} ; \
		fi
	@install_yn=no; \
	for tclfile in ${TCLFILES} verylongdummyname ; do \
		if [ -f $$tclfile ] ; then \
			install_yn=yes; \
			${INSTALL} -m 0775 $$tclfile ${dest_TCLDIR}/ ; \
			echo ${INSTALL} -m 0775 $$tclfile ${dest_TCLDIR}/ ; \
		fi \
	done; \
	if [ $$install_yn = "yes" ] ; then \
		( cd ${dest_TCLDIR} ; echo 'auto_mkindex . *.tcl' | ${TCLSH} ) \
	fi

exeinstall:
	@echo " "
	@echo "Installing stable binaries"
	@if [ ! -d ${dest_BINDIR} ] ; then \
		${MKDIR} ${dest_BINDIR} ; \
		${CHMOD} g+w ${dest_BINDIR} ; \
		fi
	@for exefile in ${XFILES} ${SCRIPTS} verylongdummyname ; do \
		if [ -f $$exefile ] ; then \
			${INSTALL} -m 0775 $$exefile ${dest_BINDIR}/ ; \
			echo ${INSTALL} -m 0775 $$exefile ${dest_BINDIR}/ ; \
		fi \
	done
	@for lntarget in ${RUNTCLS} verylongdummyname ; do \
		if [ $$lntarget != verylongdummyname ] ; then \
			if [ `uname` = Darwin -o X`uname | grep CYGWIN`X != XX ] ; then \
				lntarget=$${lntarget}_gui ; \
			fi ; \
			cd ${dest_BINDIR} ; ${RM} -f /tmp/fslgrot $$lntarget ; \
			ln -s Runtcl $$lntarget ; \
			echo ln -s Runtcl $$lntarget ; \
		fi \
	done
	@for lntarget in ${RUNAVWS} verylongdummyname ; do \
		if [ $$lntarget != verylongdummyname ] ; then \
			cd ${dest_BINDIR} ; ${RM} -f /tmp/fslgrot $$lntarget ; \
			ln -s runavw $$lntarget ; \
			echo ln -s runavw $$lntarget ; \
		fi \
	done

libinstall:
	@echo " "
	@echo "Installing stable library archives"
	@if [ ! -d ${dest_LIBDIR} ] ; then \
		${MKDIR} ${dest_LIBDIR} ; \
		${CHMOD} g+w ${dest_LIBDIR} ; \
		fi
	@for libfile in ${AFILES} verylongdummyname ; do \
		if [ -f $$libfile ] ; then \
			${INSTALL} -m 0664 $$libfile ${dest_LIBDIR} ; \
			echo  ${INSTALL} -m 0664 $$libfile ${dest_LIBDIR} ; \
			if [ `uname` = Darwin -o X`uname | grep CYGWIN`X != XX ] ; then \
				${RANLIB} ${dest_LIBDIR}/$$libfile ; \
			fi ; \
		fi \
	done

hdrinstall:
	@echo " "
	@echo "Installing stable header files"
	@if [ ! -d ${dest_INCDIR} ] ; then \
		${MKDIR} ${dest_INCDIR} ; \
		${CHMOD} g+w ${dest_INCDIR} ; \
		fi
	@if [ ! -d ${dest_INCDIR}/${PROJNAME} ] ; then \
		${MKDIR} ${dest_INCDIR}/${PROJNAME} ; \
		${CHMOD} g+w ${dest_INCDIR}/${PROJNAME} ; \
		fi
	@for hdrfile in ${HFILES} verylongdummyname ; do \
		if [ -f $$hdrfile ] ; then \
			${INSTALL} -m 0664 $$hdrfile ${dest_INCDIR}/${PROJNAME}/ ; \
			echo ${INSTALL} -m 0664 $$hdrfile ${dest_INCDIR}/${PROJNAME}/ ; \
		fi \
	done

showlibs:
	@echo ${LIBS}

tagstable: commit
	@if [ "X${PROJNAME}X" = XX ] ; then \
		echo " " ; \
		echo "No PROJNAME defined in the Makefile" ; \
		echo "    ... aborting install" ; \
		echo " " ; \
		exit 4 ; \
	fi;
	cvs tag -F stable 
	cvs tag -F stable-`date +%d%b%Y`


stable-branch: 
	dt=`date +%d%b%Y` ; export dt ; \
	bnchnm=stable-branch-$$dt; export bnchnm ; \
	if [ ! -f branchname.log ] ; then \
	  echo "Initial creation of branchname log file" > branchname.log ; \
	  cvs add branchname.log ; \
	  cvs ci -m 'Added branchname.log for stable checkin purposes' ; \
	fi  ; \
	cvs tag -b $$bnchnm  ; \
	cvs update -r $$bnchnm  ; \
	echo $$bnchnm >> branchname.log  ; \
	cvs ci -m 'Updated branchname.log'  ; \
	cvs admin -N stable-branch:$$bnchnm  ; \
	cvs tag -F stable ; \
	cvs tag -F stable-$$dt  ; \
	cvs tag -F last-merge  ; \
	cvs update -A 
        

insertcopyright:
	${FSLCONFDIR}/common/insertcopyright * */*	

update:
	cvs -q -r update


commit:
	cvs -q commit

depend: 
	${RM} -f /tmp/fslgrot depend.mk
	${MAKE} depend.mk

depend.mk:
	@echo Building dependency file depend.mk
	@for srcfile in *.c *.cc *.cpp *.inc *.hpp verylongdummyname ; do \
		if [ -f $$srcfile ] ; then \
			${CC} ${DEPENDFLAGS} ${AccumulatedIncFlags} $$srcfile >> depend.mk ; \
		else \
			touch depend.mk; \
		fi \
	done

include depend.mk
