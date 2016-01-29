# Auto-configure options for Mac OS X Universal build

# Written by Duncan Mortimer

macosx_opts="-arch x86_64 -mmacosx-version-min=10.9"
cflags="${cflags} ${macosx_opts}"
cxxflags="${cxxflags} ${macosx_opts}"
ldflags="${ldflags} -Wl,-search_paths_first ${macosx_opts}"
configure_opts="${configure_opts} --disable-dependency-tracking"
