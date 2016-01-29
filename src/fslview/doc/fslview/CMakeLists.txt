PROJECT(FSLVIEW)

FILE(GLOB docs *.html *.png *.css fsl-bg.jpg fslview.adp)
INSTALL(FILES ${docs} DESTINATION doc/fslview)
