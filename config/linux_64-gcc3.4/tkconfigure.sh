# Auto-configure options for Centos4_64 Build

# Written by Matthew Webster
ldflags="${ldflags} -L/usr/X11R6/lib64/"
configure_opts="${configure_opts} --x-libraries=/usr/X11R6/lib64"

