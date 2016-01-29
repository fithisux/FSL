# Auto-configure options for Mac OS X Universal build

# Written by Duncan Mortimer

macosx_universal_opts="-arch x86_64"
cflags="${cflags} ${macosx_universal_opts} -mmacosx-version-min=10.8"
cxxflags="${cxxflags} ${macosx_universal_opts} -mmacosx-version-min=10.8"
ldflags="${ldflags} -Wl,-search_paths_first ${macosx_universal_opts}"
configure_opts="${configure_opts} --disable-dependency-tracking --x-includes=/usr/X11/include --x-libraries=/usr/X11/lib"

