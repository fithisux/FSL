# Auto-configure options for Mac OS X Universal build

# Written by Duncan Mortimer

macosx_universal_opts="-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch ppc -arch i386"
cflags="${cflags} ${macosx_universal_opts} -mmacosx-version-min=10.4"
cxxflags="${cxxflags} ${macosx_universal_opts} -mmacosx-version-min=10.4"
ldflags="${ldflags} -Wl,-search_paths_first ${macosx_universal_opts}"
configure_opts="${configure_opts} --disable-dependency-tracking"

