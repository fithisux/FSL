include ${FSLCONFDIR}/default.mk

PROJNAME = giftiio

LIBTOOL=libtool

ifeq ($(FSLMACHTYPE),apple-darwin8-gcc4.0)
	ARCHFLAGS = -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -I/usr/X11R6/include/
	ARCHLDFLAGS = -Wl,-search_paths_first -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -L/Developer/SDKs/MacOSX10.4u.sdk/usr/X11R6/lib/
endif

CFLAGS += -std=c99
ARCHFLAGS += -DHAVE_ZLIB

gifti_tool: gifti_tool.o gifti_io.o gifti_xml.o
	$(CC) -o $@ gifti_tool.o gifti_io.o gifti_xml.o \
	      $(LDFLAGS) -lexpat -lniftiio -lznz -lz -lm

gifti_test: gifti_test.o gifti_io.o gifti_xml.o
	$(CC) -o $@ gifti_test.o gifti_io.o gifti_xml.o \
	      $(LDFLAGS) -lexpat -lniftiio -lznz -lz -lm

libgiftiio.a: gifti_test.o gifti_io.o gifti_xml.o
	${AR} -r libgiftiio.a gifti_test.o gifti_io.o gifti_xml.o

all: gifti_tool gifti_test libgiftiio.a


