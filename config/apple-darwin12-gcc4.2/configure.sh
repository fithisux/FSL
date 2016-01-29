# Auto-configure options for Mac OS X Universal build

# Written by Duncan Mortimer

macosx_universal_opts="-arch x86_64 -mmacosx-version-min=10.8"
cflags="${cflags} ${macosx_universal_opts}"
cxxflags="${cxxflags} ${macosx_universal_opts}"
ldflags="${ldflags} -Wl,-search_paths_first ${macosx_universal_opts}"
configure_opts="${configure_opts} --disable-dependency-tracking"
