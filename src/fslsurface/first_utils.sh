#!/bin/sh
#   Copyright (C) 2012 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 5.0 (c) 2012, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/9564.
export LC_ALL=C
useBvars=0;
doVertexAnalysis=0;
infile=""
reffile=""
outname=""
imagepath=""
designame=""
DOF=0;
DEBUG=0;

#if 0 
reconMode=0;


while [ $# -gt 0 ] ; do 
	if [ $1 = "--vertexAnalysis" ] ; then 
		doVertexAnalysis=1
		shift
	elif [ $1 = "--usebvars" ] ; then 
		useBvars=1
		shift
	elif [ $1 = "--useReconNative" ] ; then 
		reconMode=1
		shift
	elif [ $1 = "--useReconMNI" ] ; then 
		reconMode=2
		shift
	elif  [ $1 = "-i" ] ; then 
		infile=$2
		shift 2
    elif  [ $1 = "-r" ] ; then 
        reffile=$2
        shift 2
	elif  [ $1 = "-d" ] ; then 
		designame=$2
		shift 2
	elif  [ $1 = "-o" ] ; then 
		outname=$2
		shift 2
	elif  [ $1 = "-a" ] ; then 
		imagepath="-a $2"
		shift 2
	elif  [ $1 = "--useRigidAlign" ] ; then 
		DOF=`echo "${DOF}+6" | bc`
		shift 1
	elif  [ $1 = "--useScale" ] ; then 
		DOF=`echo "${DOF}+1" | bc`
		shift 1
	elif  [ $1 = "--debug" ] ; then 
        echo "debug"
		DEBUG=1
		shift 1
	else
		echo "Unrecognized option $1"
		exit 1;
	fi

done
#Done running analysis
	if [ $doVertexAnalysis -eq 1 ] ;then 
#if [ $useBvars -eq 1 ] ; then
		
			if [ -f  ${outname}_surfaces_aligned.txt ] ; then 
				rm ${outname}_surfaces_aligned.txt
			fi
		   
			echo "Reconstruct Data from bvars... $imagepath"
            if [ $useBvars = 1 ] ; then 
                        if [ $reconMode -eq 0 ] ; then 
                            echo "Please choose the space to reconstruct surfaces into"
                            exit 1
                        elif [ $reconMode -eq 1 ] ; then 
                            echo "do reconMode"
            #./fslsurface_utils --mni_template=${FSLDIR}/data/standard/MNI152_T1_1mm --doReconSurfacesFromBvars  $imagepath -s $infile -o $outname 
                             ${FSLDIR}/bin/fslsurfacemaths -reconAllFromBvarsAndSave $infile $outname 
                                reffile="${outname}_modelmean.gii"
                         elif [ $reconMode -eq 2 ] ; then 
                            echo "reecon mni"
            #./fslsurface_utils --mni_template=${FSLDIR}/data/standard/MNI152_T1_1mm --doReconSurfacesFromBvars --reconMNISpace $imagepath -s $infile -o $outname 

                          fi 
            else
                echo "copy file"
                cp $infile ${outname}_list.txt


            fi
			for surf in `cat ${outname}_list.txt` ; do 
				echo "Registering surface ${surf}..." 
#				echo ${fsurf}_2_mni_dof${DOF}.gii
				fsurf=`basename $surf .gii`
                fsurf=`basename $surf .vtk`
				dsurf=`dirname $surf`
			#	echo ${fsurf}_to_mni_dof${DOF}.gii

				if [ ${DOF} -ge 6 ] ; then 
					#./build/Debug/fslsurface_utils --doLeastSquaresReg -s $surf -r ${outname}_mean_mni.gii -o ${fsurf}_2_mni_dof${DOF}.gii -d ${DOF}
                    fsurf=`basename $surf .gii`
                    ${FSLDIR}/bin/fslsurfacemaths $surf -reg_lq ${reffile} ${DOF} ${dsurf}/${fsurf}_to_mni_dof${DOF}.mat  ${dsurf}/${fsurf}_to_mni_dof${DOF}.gii					

                    echo "${SURFDIR}/fslsurfacemaths $surf -reg_lq ${outname}_modelmean.gii ${DOF} ${fsurf}_to_mni_dof${DOF}.mat  ${fsurf}_to_mni_dof${DOF}.gii"
#./fslsurface_utils --doLeastSquaresReg -s $surf -r ${outname}_mean_mni.gii -o ${dsurf}/${fsurf}_to_mni_dof${DOF}.gii -d ${DOF}
			#		echo "DOF ${DOF}"
					echo ${dsurf}/${fsurf}_to_mni_dof${DOF}.gii >> ${outname}_surfaces_aligned.txt
				else	
			#	echo "dof0 $surf ${dsurf}/${fsurf}_to_mni_dof0.gii"
					mv $surf  ${dsurf}/${fsurf}_to_mni_dof0.gii
					echo  ${dsurf}/${fsurf}_to_mni_dof${DOF}.gii >> ${outname}_surfaces_aligned.txt
				fi
			
			done
echo "run vetex MVGLM "


			 ${FSLDIR}/bin/fslsurfacemaths -vertexMVGLM ${outname}_surfaces_aligned.txt ${designame} none ${outname}_F.gii
			echo "done fslSurfaceutils vertex analysis"
#			./build/Debug/fslsurface_utils  --doVertexGLMfit -i ${outname}_surfaces_aligned.txt
#	./fslsurface_utils  --doVertexGLMfit --designame=${designame} -s ${outname}_surfaces_aligned.txt -o ${outname}_F
		
            echo DEBUG $DEBUG
			if [ $DEBUG = 0 ] ;then
                echo clean
                if [ $useBvars = 1 ] ; then
                                rm `cat ${outname}_surfaces_aligned.txt` `cat ${outname}_list.txt`
                fi
            fi
		
#fi
	fi
	
	
	
