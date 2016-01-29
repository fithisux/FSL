PROJECT(FSL)

FIND_PATH( FSL_DIR etc/fslconf/fsl.sh
	   /usr/local/fsl )

INCLUDE_DIRECTORIES( ${FSL_DIR}/extras/include/zlib )
ADD_DEFINITIONS ( -DHAVE_ZLIB )

SUBDIRS(znzlib niftiio fslio newmat libprob cprob miscmaths utils)
