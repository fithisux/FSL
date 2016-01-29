FIND_PATH( FSL_DIR etc/fslconf/fsl.sh
	   /usr/local/fsl )

INCLUDE_DIRECTORIES( ${FSL_DIR}/extras/include/zlib )

SUBDIRS(storage fslview)
