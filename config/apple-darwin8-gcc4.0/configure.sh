# Auto-configure options for Mac OS X Universal build

# Written by Duncan Mortimer

macosx_universal_opts="-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch ppc64 -arch i386 -arch x86_64"
cflags="${cflags} ${macosx_universal_opts}"
cxxflags="${cxxflags} ${macosx_universal_opts}"
ldflags="${ldflags} -Wl,-search_paths_first ${macosx_universal_opts}"
configure_opts="${configure_opts} --disable-dependency-tracking"
