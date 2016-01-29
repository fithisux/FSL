#!/bin/sh
# FSL Host type identifier
#  - Used by several scripts to identify the host type for build environment

# Written by Duncan Mortimer
#  FMRIB Analysis Group, University of Oxford

# SHCOPYRIGHT


#### Identify machine ####
machtype=""

system_type=`uname -s`
case ${system_type} in
    SunOS)
        AWK=nawk
        ;;
    *)
        AWK=awk
        ;;
esac

gcc_version=`gcc -dumpversion 2> /dev/null`
if [ $? -eq 0 ]; then
    # GCC is installed
    gcc_version=`echo $gcc_version | $AWK -F . '{ printf "%s.%s", $1,$2 }'`
    gcc_host=`gcc -dumpmachine`

    gcc_cpu_type=`echo $gcc_host | $AWK -F - '{ printf "%s", $1; }'`
    gcc_os_vendor=`echo $gcc_host | $AWK -F - '{ printf "%s", $2; }'`
    gcc_os_name=`echo $gcc_host | $AWK -F - '{ printf "%s", $3; }'`

    case ${system_type} in
        Darwin)
            # We are going to build a universal (ppc32, ppc64, x86_32, x86_64)
            #    combined binary, so we only need 1 host type
            os_release=`uname -r | $AWK -F . '{ printf "%s", $1 }'`
	    if [ $os_release -ge 12 ]; then
		# Using LLVM/clang from now on
		llvm_vstr=`cc -v 2>&1`
		llvm_v=`echo ${llvm_vstr} | awk '{ print $4 }'`
		os_v=`echo ${llvm_vstr} | awk '{ print $11 }' | awk -F. '{ print $1 }' | sed 's/x86_64-//'`
		machtype=${os_v}-llvm${llvm_v}
            elif [ $os_release -ge 8 ]; then
                # This is 10.4 so Universal builds possible
                gcc_host="${gcc_os_vendor}-${gcc_os_name}"
                # Note, for BSDish platforms, gcc_os_name includes a version
		machtype=${gcc_host}-gcc${gcc_version}
            else
	        # Prior to Tiger and gcc4, things were a mess...
		gcc_host="${gcc_cpu_type}-${gcc_os_vendor}"
		machtype=${gcc_host}-gcc${gcc_version}
                # Note, gcc_os_vendor in this case is actually gcc_os_name above
	    fi
	    ;;
        Linux)
	    if [ "${gcc_cpu_type}" = "x86_64" -o "${gcc_cpu_type}" = "x86-64" ]; then
		gcc_host="${gcc_os_name}_64"
	    elif [ `echo ${gcc_cpu_type} | grep i.86` ]; then
		gcc_host="${gcc_os_name}_32"
	    else
		gcc_host="${gcc_cpu_type}-${gcc_os_vendor}-${gcc_os_name}"
            fi
	    machtype=${gcc_host}-gcc${gcc_version}
	    ;;
	*)
	    gcc_host="${gcc_cpu_type}-${gcc_os_vendor}-${gcc_os_name}"
	    machtype=${gcc_host}-gcc${gcc_version}
	    ;;
    esac
    
fi
echo ${machtype}
exit 0

