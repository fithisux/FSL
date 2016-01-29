# Auto-configure options to fix TCL directory access to 64bit NFS volumes from 32bit hosts

# Written by Duncan Mortimer

cflags="${cflags} -DHAVE_STRUCT_DIRENT64=1"
