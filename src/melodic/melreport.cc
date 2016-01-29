/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
    independent components
   
    melreport.cc - report generation

    Christian F. Beckmann, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */

#include "newimage/newimageall.h"
#include "utils/log.h"
#include "melreport.h"
#include "melhlprfns.h"
#include "miscmaths/miscprob.h"

namespace Melodic{
  
	void MelodicReport::IC_rep(MelGMix &mmodel, int cnum, int dim, Matrix ICstats){

		if( bool(opts.genreport.value()) ){
    	addlink(mmodel.get_prefix()+".html",num2str(cnum));
			IChtml.setDir(report.getDir(),mmodel.get_prefix()+".html");

      {//start IC page
				IChtml << "<HTML><HEAD><link REL=stylesheet TYPE=text/css href=file:" +
					(string) getenv("FSLDIR") +"/doc/fsl.css>" << endl
					<< "<style type=\"text/css\">OBJECT { width: 100% }</style>"
	       	<< "<TITLE>FSL</TITLE></HEAD>" << endl
	  			<< "<IFRAME  height=" << int(melodat.get_numfiles()/30 + 1)*50 
					<<"px width=100% src=nav.html frameborder=0></IFRAME><BR>"<< endl
	       	<< "<Center>" << endl;
			
				if(cnum>1)
	  			IChtml << "<a href=\"" << string("IC_")+num2str(cnum-1)
					<<".html\">&#60;</a>&nbsp;-&nbsp;";
				else
	  			IChtml << "<a href=\"" << string("IC_")+num2str(melodat.get_mix().Ncols())
					<<".html\">&#60;</a>&nbsp;-&nbsp;";
	
				if(cnum<dim)
	  			IChtml << "<a href=\"" << string("IC_")+num2str(cnum+1)
					<<".html\">&#62;</a>";
				else 
	  			IChtml << "<a href=\"" << string("IC_")+num2str(1)
					<<".html\">&#62;</a>";
				IChtml << "<p><H3>MELODIC Component " << num2str(cnum)
				<< "<br></b></H3><hr><p>" << endl;
			}
      {//output IC stats
    		if(ICstats.Storage()>0&&ICstats.Nrows()>=cnum){
	  			IChtml << fixed << setprecision(2) << std::abs(ICstats(cnum,1)) << " % of explained variance";
	  			if(ICstats.Ncols()>1)
	    			IChtml << "; &nbsp; &nbsp; " << std::abs(ICstats(cnum,2)) << " % of total variance";
					if(ICstats.Ncols()>2&&opts.addsigchng.value()){
	    			IChtml << "<p>" <<endl;
	    			IChtml << " &nbsp; &nbsp; " << ICstats(cnum,3) << " % signal change (pos peak voxel); &nbsp; &nbsp;" << ICstats(cnum,4) << "% signal change (peak neg. voxel)" << endl ;
	  			}
	  			IChtml << "<hr><p>" << endl;
				}
      }
      
      volume4D<float> tempVol;       
      if(mmodel.get_threshmaps().Storage()>0&&
	    (mmodel.get_threshmaps().Ncols() == mmodel.get_data().Ncols()))
	    {//Output thresholded IC map	
				tempVol.setmatrix(mmodel.get_threshmaps().Row(1),melodat.get_mask());
				volume<float> map1;
				volume<float> map2;
				map1 = threshold(tempVol[0],float(0.0), tempVol[0].max());
				map2 = threshold(tempVol[0],tempVol[0].min(), float(0.0));
	
				volume<float> newvol;

				miscpic newpic;
	
				float map1min = std::max((map1 + binarise(tempVol[0],tempVol[0].min(), 
						    	float(0.0)) * map1.max()).min(),float(0.001));
			  float map1max = std::max(map1.max(),float(0.001));
				float map2min = std::min(map2.min(),float(-0.001));
				float map2max = std::min((map2 + binarise(tempVol[0],float(0.0), 
						  		tempVol[0].max()) * map2.min()).max(),float(-0.001));
	
    		newpic.overlay(newvol, melodat.get_bg(), map1, map2, 
		       		melodat.get_bg().percentile(0.02),
		       		melodat.get_bg().percentile(0.98),
		       		map1min, map1max, map2max, map2min, 
		       		0, 0);
				char instr[10000];
	
				sprintf(instr," ");
				strcat(instr,axials_instr.c_str());
				strcat(instr,string(report.appendDir(mmodel.get_prefix()+
					     	"_thresh.png")).c_str());
				newpic.set_title(string("Component No. "+num2str(cnum)+
								" - thresholded IC map  ") + mmodel.get_infstr(1));
				newpic.set_cbar(string("ysb"));
	
				if((std::abs(map1.max()-map1.min())>0.01) || (std::abs(map2.max()-map2.min())>0.01))
	  			newpic.slicer(newvol, instr); 
				else
	  			newpic.slicer(newvol, instr);
	  		IChtml << "<a href=\""+mmodel.get_prefix()+"_MM.html\">";
	    	IChtml << "<img BORDER=0 SRC=\""+mmodel.get_prefix()
	    		+"_thresh.png\" ALT=\"MMfit\"></A><p>" << endl;
      }		
			
      {//plot time course
    	IChtml << "<H3> Temporal mode </H3><p>" << endl <<endl;     	
		miscplot newplot;
		
			Matrix tmptc = melodat.get_Tmodes(cnum-1).Column(1).t();

			newplot.col_replace(0,0xFF0000);

			newplot.add_label(string("IC ")+num2str(cnum)+" time course");
			//add GLM OLS fit
			if(melodat.Tdes.Storage()>0 &&
		  melodat.glmT.get_beta().Nrows() == melodat.Tdes.Ncols()){
				tmptc &= melodat.glmT.get_beta().Column(cnum).t() * melodat.Tdes.t();
				newplot.add_label("full model fit");
			}

			//add deviation around time course
			if(melodat.get_Tmodes(cnum-1).Ncols()>1 && opts.varplots.value()){
				Matrix tmp = stdev(melodat.get_Tmodes(cnum-1).Columns(2,melodat.get_Tmodes(cnum-1).Ncols()).t(),1);
				tmptc &= melodat.get_Tmodes(cnum-1).Column(1).t()+tmp;
				tmptc &= melodat.get_Tmodes(cnum-1).Column(1).t()-tmp;
				newplot.add_label("std error across subjects");
				newplot.col_replace(tmptc.Nrows()-1,0x808080);
				newplot.col_replace(tmptc.Nrows()-2,0x808080);				
			}
		
	    if(opts.tr.value()>0.0)
				newplot.add_xlabel(string("Time (seconds); TR = ")+
				float2str(opts.tr.value(),0,2,0)+" s");
			else 
		  	newplot.add_xlabel(string("Time (TRs)"));
			
			newplot.add_ylabel("Normalised Response");
			newplot.set_yrange(tmptc.Row(1).Minimum()-0.05*(tmptc.Row(1).Maximum() - 
				tmptc.Row(1).Minimum()),tmptc.Row(1).Maximum()+
				0.05*(tmptc.Row(1).Maximum()-tmptc.Row(1).Minimum()));
			newplot.grid_swapdefault();
			
			newplot.set_xysize(750,150);
			
	        newplot.timeseries(tmptc,
			  report.appendDir(string("t")+num2str(cnum)+".png"),
			  string("Timecourse No. ")+num2str(cnum), 
			  opts.tr.value(),150,12,1,false);

			if(melodat.get_Tmodes(cnum-1).Ncols()>1)
				tmptc &= melodat.get_Tmodes(cnum-1).Columns(2,melodat.get_Tmodes(cnum-1).Ncols()).t();
	     write_ascii_matrix(report.appendDir(string("t")
				+num2str(cnum)+".txt"),tmptc.t());
	     IChtml << "<A HREF=\"" << string("t")
	  	  +num2str(cnum)+".txt" << "\"> ";
			 IChtml << "<img BORDER=0 SRC=\"" 
	  	  +string("t")+num2str(cnum)+".png\"></A><p>" << endl;	
	
			if(melodat.get_numfiles()>1 && melodat.explained_var.Storage()>0 
				&& melodat.explained_var.Ncols()>=cnum && opts.varvals.value())
				IChtml << "Rank-1 approximation of individual time courses explains " 
				<< std::abs(melodat.explained_var(cnum)) << "% of variance.<p>" << endl;
			}//time series plot

	 		if(!opts.pspec.value())
	   	{//plot frequency  
    		miscplot newplot;
	    	RowVector empty(1);
	 			empty = 0.0;
				int fact = int(std::pow(10.0,int(std::log10(float(melodat.get_Tmodes(0).Nrows())))));
				if(opts.logPower.value())
					newplot.add_ylabel(string("log-Power"));
				else
					newplot.add_ylabel(string("Power"));
		
				Matrix fmixtc = calc_FFT(melodat.get_Tmodes(cnum-1).Column(1), opts.logPower.value());
			  
				newplot.set_Ylabel_fmt("%.0f");
				newplot.set_yrange(0.0,1.02*fmixtc.Maximum());
				newplot.grid_swapdefault();
				newplot.set_xysize(750,150);
				
				if(opts.tr.value()>0.0){
					newplot.add_xlabel(string("Frequency (in Hz / ")+num2str(fact)+ " )");
	  			newplot.timeseries(empty | fmixtc.t(),
			  	report.appendDir(string("f")+
						num2str(cnum)+".png"),
			     	string("Powerspectrum of timecourse"),
			     	fact/(opts.tr.value()*melodat.get_Tmodes(0).Nrows()), 150,0,2);
				}else{	
					newplot.add_xlabel(string("Frequency (in cycles); ")
					+"frequency(Hz)=cycles/("
			    +num2str(melodat.get_Tmodes(0).Nrows())
			    +"* TR); period(s)=("
			    +num2str(melodat.get_Tmodes(0).Nrows())
			    +"* TR)/cycles"
					);
	  			newplot.timeseries(fmixtc.t(),
			    report.appendDir(string("f")+num2str(cnum)+".png"),
			     	string("Powerspectrum of timecourse"));
				}
				write_ascii_matrix(report.appendDir(string("f")
			 		+num2str(cnum)+".txt"), fmixtc);
				IChtml << "<A HREF=\"" << string("f")
	  			+num2str(cnum)+".txt" << "\"> ";
				IChtml << "<img BORDER=0 SRC=\"" 
	  			+string("f")+num2str(cnum)+".png\"></A><p>" << endl;
      }//frequency plot
   		{//add T-mode GLM F-stats for full model fit & contrasts
				if(melodat.Tdes.Storage() > 0 &&
					melodat.glmT.get_beta().Nrows() == melodat.Tdes.Ncols()){
							IChtml << " <TABLE border=1 bgcolor=ffffff cellpadding=5>" <<
								"<CAPTION><EM> <b>GLM (OLS) on time series </b></EM></CAPTION>" << endl
								<< "<TR valign=middle><TH ><EM>GLM &beta;'s</EM> <TH> <EM> F-test on <br> full model fit </em>";
							if(melodat.Tcon.Storage() > 0)
								IChtml << "<TH ><EM>Contrasts</EM>"<<endl;
							IChtml << "<TR><TD><TABLE border=0><TR><TD align=right>" << endl; 
							for(int ctr=1;ctr <= melodat.Tdes.Ncols();ctr++)
								IChtml << " PE(" <<num2str(ctr)+"): <br>" << endl;
							IChtml << "<TD align=right>" << endl;
							for(int ctr=1;ctr <= melodat.Tdes.Ncols();ctr++)
								IChtml << melodat.glmT.get_beta().Column(cnum).Row(ctr) << "<br>" <<endl;
							IChtml << "</TABLE>" <<
								" <TD align=center> F = "<< melodat.glmT.get_f_fmf().Column(cnum) << 
								" <BR> dof1 = " << melodat.Tdes.Ncols() << "; dof2 = " 
								<< melodat.glmT.get_dof() << "<BR>" <<endl;
							if(melodat.glmT.get_pf_fmf().Column(cnum).AsScalar() < 0.05)
								IChtml << fixed << setprecision(5) <<"<b> p < " << melodat.glmT.get_pf_fmf().Column(cnum) <<
								"<BR> (uncorrected for #comp.)<b></TD>" << endl;
							else
								IChtml << fixed << setprecision(5) << " p < " << 
								melodat.glmT.get_pf_fmf().Column(cnum) << 
								"<BR> (uncorrected for #comp.)</TD>" << endl;
						if(melodat.Tcon.Storage() > 0	&&
						melodat.Tdes.Ncols() == melodat.Tcon.Ncols()){
							IChtml << fixed << setprecision(2) << "<TD><TABLE border=0><TR><TD align=right>" <<endl;
							for(int ctr=1; ctr <= melodat.Tcon.Nrows() ; ctr++)
								IChtml << "COPE(" << num2str(ctr) << "): <br>" << endl;
							IChtml << "<td align=right>" << endl;
							for(int ctr=1; ctr <= melodat.Tcon.Nrows() ; ctr++)
								IChtml <<" z = <BR>" <<endl;
							IChtml << "<td align=right>" << endl;						
							for(int ctr=1; ctr <= melodat.Tcon.Nrows() ; ctr++)
								IChtml << melodat.glmT.get_z().Column(cnum).Row(ctr) <<";<BR>" <<endl;
							IChtml << "<td align=right>" << endl;
							for(int ctr=1; ctr <= melodat.Tcon.Nrows() ; ctr++)
								if(melodat.glmT.get_p().Column(cnum).Row(ctr).AsScalar() < 0.05)
									IChtml << fixed << setprecision(5) << "<b> p < " << melodat.glmT.get_p().Column(cnum).Row(ctr) << 
									"</b><BR>" << endl;
								else
									IChtml << fixed << setprecision(5) <<" p < " << melodat.glmT.get_p().Column(cnum).Row(ctr) << 
									"<BR>" << endl;
							IChtml << "</TABLE></td></tr>" << endl;
					  }
					}
				IChtml << "</TABLE><p>" << endl;
			}
  
      if(cnum <= (int)melodat.get_Smodes().size())
	    {//plot subject mode 
			  
	  		Matrix smode;
	  		smode = melodat.get_Smodes(cnum-1);
	
	  		if(smode.Nrows() > 1){
	  	  	IChtml << "<hr><H3> Sessions/Subjects mode </H3><p>" << endl <<endl;
		    	miscplot newplot;

					//add GLM OLS fit
					//newplot.setscatter(smode,2);

					if(melodat.Sdes.Storage() > 0&&
					melodat.glmS.get_beta().Nrows() == melodat.Sdes.Ncols()){
						smode |= melodat.Sdes * melodat.glmS.get_beta().Column(cnum);
						newplot.add_label(string("IC ")+num2str(cnum)+" subject/session-mode");
						newplot.add_label("full model fit");
					}
					newplot.grid_swapdefault();
					newplot.set_Ylabel_fmt("%.2f");
					newplot.add_xlabel(" Subject Number");
					

					newplot.set_xysize(750,150);
	      	newplot.timeseries(smode.t(), 
			    	report.appendDir(string("s")+num2str(cnum)+".png"),
			      string("Subject/Session mode No. ") + num2str(cnum));
					newplot.clear_xlabel();
					newplot.clear_labels();
					newplot.set_xysize(120,200);
					newplot.set_minmaxscale(1.1);

					newplot.boxplot((Matrix)smode.Column(1),
			    	report.appendDir(string("b")+num2str(cnum)+".png"),
			        string("Subject/Session mode"));
	      	write_ascii_matrix(report.appendDir(string("s")
		 				+num2str(cnum)+".txt"),  smode);
	      	IChtml << "<A HREF=\"" << string("s")
	        	+num2str(cnum)+".txt" << "\"> ";
	      	IChtml << "<img BORDER=0 SRC=\"" 
	      		+string("s")+num2str(cnum)+".png\"></A>" << endl;
					IChtml << "<A HREF=\"" << string("s")
	        	+num2str(cnum)+".txt" << "\"> ";
	      	IChtml << "<img BORDER=0 SRC=\"" 
	      		+string("b")+num2str(cnum)+".png\"></A><p>" << endl;
	    	}
   			{//add S-mode GLM F-stats for full model fit & contrasts
			  	if(melodat.Sdes.Storage() > 0 &&
					melodat.glmS.get_beta().Nrows() == melodat.Sdes.Ncols()){
							IChtml << " <TABLE border=1 bgcolor=ffffff cellpadding=5>" <<
								"<CAPTION><EM> <b>GLM (OLS) on subject/session-mode </b></EM></CAPTION>" << endl
								<< "<TR valign=middle><TH ><EM>GLM &beta;'s</EM> <TH> <EM> F-test on <br> full model fit </em>";
							if(melodat.Scon.Storage() > 0)
								IChtml << "<TH ><EM>Contrasts</EM>"<<endl;
							IChtml << "<TR><TD><TABLE border=0><TR><TD align=right>" << endl; 
							for(int ctr=1;ctr <= melodat.Sdes.Ncols();ctr++)
								IChtml << " PE(" <<num2str(ctr)+"): <br>" << endl;
							IChtml << "<TD align=right>" << endl;
							for(int ctr=1;ctr <= melodat.Sdes.Ncols();ctr++)
								IChtml << melodat.glmS.get_beta().Column(cnum).Row(ctr) << "<br>" <<endl;
							IChtml << "</TABLE>" <<
								" <TD align=center> F = "<< melodat.glmS.get_f_fmf().Column(cnum) << 
								" <BR> dof1 = " << melodat.Sdes.Ncols() << "; dof2 = " 
								<< melodat.glmS.get_dof() << "<BR>" <<endl;
							if(melodat.glmS.get_pf_fmf().Column(cnum).AsScalar() < 0.05)
								IChtml << fixed << setprecision(5) <<"<b> p < " << melodat.glmS.get_pf_fmf().Column(cnum) <<
								"<BR> (uncorrected for #comp.)<b></TD>" << endl;
							else
								IChtml << fixed << setprecision(5) << " p < " << 
								melodat.glmS.get_pf_fmf().Column(cnum) << 
								"<BR> (uncorrected for #comp.)</TD>" << endl;
				  if(melodat.Scon.Storage() > 0 	&& melodat.Sdes.Storage() > 0 &&
							melodat.Sdes.Ncols() == melodat.Scon.Ncols()){
							IChtml << fixed << setprecision(2) << "<TD><TABLE border=0><TR><TD align=right>" <<endl;
							for(int ctr=1; ctr <= melodat.Scon.Nrows() ; ctr++)
								IChtml << "COPE(" << num2str(ctr) << "): <br>" << endl;
							IChtml << "<td align=right>" << endl;
							for(int ctr=1; ctr <= melodat.Scon.Nrows() ; ctr++)
								IChtml <<" z = <BR>" <<endl;
							IChtml << "<td align=right>" << endl;						
							for(int ctr=1; ctr <= melodat.Scon.Nrows() ; ctr++)
								IChtml << melodat.glmS.get_z().Column(cnum).Row(ctr) <<";<BR>" <<endl;
							IChtml << "<td align=right>" << endl;
							for(int ctr=1; ctr <= melodat.Scon.Nrows() ; ctr++)
								if(melodat.glmS.get_p().Column(cnum).Row(ctr).AsScalar() < 0.05)
									IChtml << fixed << setprecision(5) << "<b> p < " << melodat.glmS.get_p().Column(cnum).Row(ctr) << 
									"</b><BR>" << endl;
								else
									IChtml << fixed << setprecision(5) <<" p < " << melodat.glmS.get_p().Column(cnum).Row(ctr) << 
									"<BR>" << endl;
							IChtml << "</TABLE></td></tr>" << endl;
					}
				}
				IChtml << "</TABLE><p>" << endl;
				}
	    }//subject mode plot
   
      if(mmodel.get_threshmaps().Storage()>0&&
	 			(mmodel.get_threshmaps().Ncols() == mmodel.get_data().Ncols())&&
	 			(mmodel.get_threshmaps().Nrows()>1))
	    {//Output other thresholded IC map
	  
	  	for(int tctr=2; tctr<=mmodel.get_threshmaps().Nrows(); tctr++){
	    	tempVol.setmatrix(mmodel.get_threshmaps().Row(tctr),melodat.get_mask());
	    	volume<float> map1;
	    	volume<float> map2;
	    	map1 = threshold(tempVol[0],float(0.0), tempVol[0].max());
	    	map2 = threshold(tempVol[0],tempVol[0].min(), float(0.0));

	    	volume<float> newvol; 
	    	miscpic newpic;
	    
	    	float map1min = (map1 + binarise(tempVol[0],tempVol[0].min(), 
					     	float(0.0)) * map1.max()).min();
	    	float map1max = map1.max();
	    	float map2min = map2.min();
	    	float map2max = (map2 + binarise(tempVol[0],float(0.0), 
					     	tempVol[0].max()) * map2.min()).max();
	    
	    	//cerr << endl << map1min << " " << map1max << endl
	    	//  << map2min << " " << map2max << endl;
	    
	    	//	    if(map1.max()-map1.min()>0.01)
	    	newpic.overlay(newvol, melodat.get_bg(), map1, map2, 
			   	melodat.get_bg().percentile(0.02),
			   	melodat.get_bg().percentile(0.98),
			   	map1min, map1max, map2max, map2min, 
			   	0, 0);
	    
	    	char instr[10000];
	    
	    	sprintf(instr," ");
	    	strcat(instr,axials_instr.c_str());
	    	strcat(instr,string(report.appendDir(mmodel.get_prefix()+"_thresh"+
						 	num2str(tctr)+".png")).c_str());
	    	newpic.set_title(string("Component No. "+num2str(cnum)+
				    	" - thresholded IC map ("+
				    	num2str(tctr)+")  ")+
			     		mmodel.get_infstr(tctr));
	    	newpic.set_cbar(string("ysb"));
	    	//cerr << instr << endl;
	    	newpic.slicer(newvol, instr); 
	    
	    	IC_rep_det(mmodel, cnum, dim);
	    	IChtml << "<a href=\""+mmodel.get_prefix()+"_MM.html\">";
	    	IChtml << "<img BORDER=0 SRC=\""+mmodel.get_prefix()
	      		+"_thresh"+num2str(tctr)+".png\" ALT=\"MMfit\"></A><p>" << endl;
	  	}
      }
      
      { //finish IC page
	    	IChtml<< "<HR><FONT SIZE=1>This page produced automatically by "
	      	<< "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC</A> Version "  
	      	<< version << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	      	<< "FMRIB Software Library</A>.</FONT></Center>" << endl
	      		<< "</BODY></HTML>" << endl;
      } //finish IC page
      IC_rep_det(mmodel, cnum, dim);
		}	
	}

  void MelodicReport::IC_rep_det(MelGMix &mmodel, int cnum, int dim){
    if( bool(opts.genreport.value()) ){

      {//start IC2 page
				IChtml2.setDir(report.getDir(),mmodel.get_prefix()+"_MM.html");
				IChtml2 << "<HTML><HEAD><link REL=stylesheet TYPE=text/css href=file:" +
					(string) getenv("FSLDIR") +"/doc/fsl.css>" << endl
					<< "<style type=\"text/css\">OBJECT { width: 100% }</style>"
	       	<< "<TITLE>FSL</TITLE></HEAD>" << endl
					<< "<IFRAME  height="<< int(melodat.get_numfiles()/30 + 1)*50 
					<<"px width=100% src=nav.html frameborder=0></IFRAME><p>"<< endl	
					<< "<CENTER>";
				if(cnum>1)
	  			IChtml2 << "<b><a href=\"" << string("IC_")+num2str(cnum-1)
					<<"_MM.html\">&#60;</a>&nbsp;-&nbsp;";
				else
					IChtml2 << "<b><a href=\"" << string("IC_")+num2str(melodat.get_mix().Ncols())
					<<"_MM.html\">&#60;</a>&nbsp;-&nbsp;";
			//	IChtml << "<a href=\"00index.html\">&nbsp;index&nbsp;</a>" ;
	
				if(cnum<dim)
	  			IChtml2 << "<a href=\"" << string("IC_")+num2str(cnum+1)
					<<"_MM.html\">&#62;</a>";
				else 
					IChtml2 << "<a href=\"" << string("IC_")+num2str(1)
					<<"_MM.html\">&#62;</a>";
				IChtml2 << "<p><H3>Component " << num2str(cnum)
				<< " Mixture Model fit <br></b></H3><hr><p>" << endl;
			}

      volume4D<float> tempVol; 

      if(melodat.get_IC().Storage()>0)
			{//Output raw IC map

	  		//	tempVol.setmatrix(melodat.get_IC().Row(cnum),
	  		//	  melodat.get_mask());
	  		tempVol.setmatrix(mmodel.get_data(),
			  	melodat.get_mask());
	  		volume<float> map1;
	  		volume<float> map2;
	  		map1 = threshold(tempVol[0],float(0.0), 
			   	tempVol[0].max());
	  		map2 = threshold(tempVol[0],tempVol[0].min(), 
			   	float(-0.0));

	  		volume<float> newvol; 
	  		miscpic newpic;

	  		//	float map1min = (map1 + binarise(tempVol[0],tempVol[0].min(), 
	  		//		float(0.0)) * map1.max()).robustmin();
	  		float map1max = map1.percentile(0.99);
	  		float map2min = map2.percentile(0.01);
	  		//float map2max = (map2 + binarise(tempVol[0],float(0.0), 
	  		//		tempVol[0].max()) * map2.min()).robustmax();
	
	  		newpic.overlay(newvol, melodat.get_bg(), map1, map2, 
			 		float(0.0),
			 		float(0.0),
			 		float(0.01), map1max, float(-0.01), map2min, 
			 		0, 0);

	  		char instr[10000];
       
	  		sprintf(instr," ");
	  		strcat(instr,axials_instr.c_str());
	  		strcat(instr,string(report.appendDir(mmodel.get_prefix()+
		 			".png")).c_str());
	  		newpic.set_title(string("Component No. "+num2str(cnum)+
				  " - raw Z transformed IC map (1 - 99 percentile)"));
	  		newpic.set_cbar(string("ysb"));
		
	  		newpic.slicer(newvol, instr);
			}
      IChtml2 << "<a href=\""+mmodel.get_prefix()+".html\">";
      IChtml2 << "<img BORDER=0 SRC=\""+ mmodel.get_prefix()+
				".png\"><A><p>" << endl;

      if(mmodel.get_probmap().Storage()>0&&
	 			(mmodel.get_probmap().Ncols() == mmodel.get_data().Ncols())&&
	 			(mmodel.get_probmap().Nrows() == mmodel.get_data().Nrows()))
			{//Output probmap  
	  		tempVol.setmatrix(mmodel.get_probmap(),melodat.get_mask());
	  
	  		volume<float> map;
	  		map = tempVol[0];
      
	  		volume<float> newvol; 
	  		miscpic newpic;

	  		newpic.overlay(newvol, melodat.get_bg(), map, map, 
			 		melodat.get_bg().percentile(0.02),
			 		melodat.get_bg().percentile(0.98),
			 		float(0.1), float(1.0), float(0.0), float(0.0),
			 		0, 0);

	  		char instr[10000];

	  		sprintf(instr," ");
	  		strcat(instr,"-l render1 ");
	  		strcat(instr,axials_instr.c_str());
	  		strcat(instr,string(report.appendDir(mmodel.get_prefix()+
					"_prob.png")).c_str());      
	  		newpic.set_title(string("Component No. "+num2str(cnum)+
				  " - Mixture Model probability map"));
      
	  		newpic.set_cbar(string("y"));
	  		newpic.slicer(newvol, instr); 

	  		IChtml2 << "<a href=\""+mmodel.get_prefix()+".html\">";
	  		IChtml2 << "<img BORDER=0 SRC=\""+ mmodel.get_prefix()+
	    		"_prob.png\">" << endl;
	
	  		IChtml2 << "</A><p>" << endl;
			}

			RowVector dat = mmodel.get_data().Row(1);
			if(dat.Maximum()>dat.Minimum())
      {//Output GGM/GMM fit
				miscplot newplot;

				if(mmodel.get_type()=="GGM"){
	  		newplot.add_label("IC map histogram");
	  		newplot.add_label("full GGM fit");
	  		newplot.add_label("background Gaussian");
	  		newplot.add_label("Gamma distributions");
	  		newplot.gmmfit(mmodel.get_data().Row(1),
			 		mmodel.get_means(),
			 		mmodel.get_vars(),
			 		mmodel.get_pi(),
			 		report.appendDir(mmodel.get_prefix()+"_MMfit.png"),
			 		string(mmodel.get_prefix() +
					" Gaussian/Gamma Mixture Model("+num2str(mmodel.mixtures())+") fit"),
			 		true, float(0.0), float(0.0)); 
   
			}
				else{
	  		newplot.add_label("IC map histogram");
	  		newplot.add_label("full GMM fit");
	  		newplot.add_label("individual Gaussians");
	  		newplot.gmmfit(mmodel.get_data().Row(1),
			 		mmodel.get_means(),
			 		mmodel.get_vars(),
			 		mmodel.get_pi(),
			 		report.appendDir(mmodel.get_prefix()+"_MMfit.png"),
			 		string(mmodel.get_prefix() +
					" Gaussian Mixture Model("+num2str(mmodel.mixtures())+") fit"), 
			 		false, float(0.0), float(2.0));
   
			}

			//	IChtml2 << "<A HREF=\"" +mmodel.get_prefix()+"_MMfit_detail.png\"> ";
				IChtml2 << "<img BORDER=0 SRC=\""+ mmodel.get_prefix()+
	  			"_MMfit.png\"><p>" << endl;
      } //GGM/GMM plot
      
      {//MM parameters
				IChtml2 << "<br> &nbsp;" << mmodel.get_prefix() 
					<< " Mixture Model fit <br>" << endl
					<< "<br> &nbsp; Means :  " << mmodel.get_means() << endl
					<< "<br> &nbsp;  Vars  :  " << mmodel.get_vars()  << endl
					<< "<br> &nbsp;  Prop. :  " << mmodel.get_pi()    << endl;
	    }

      { //finish IC2 page
				IChtml2<< "<HR><FONT SIZE=1>This page produced automatically by "
	       	<< "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC</A> Version " 
	       	<< version << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	       	<< "FMRIB Software Library</A>.</FONT></CENTER>" << endl
	       		<< "</BODY></HTML>" << endl;
      } //finish IC2 page
    }
  }

  void MelodicReport::IC_simplerep(string prefix, int cnum, int dim){
    if( bool(opts.genreport.value()) ){
      addlink(prefix+".html",num2str(cnum));
      IChtml.setDir(report.getDir(),prefix+".html");
      
      {//start IC page
	
				IChtml << "<HTML> " << endl
	       	<< "<TITLE>MELODIC Component " << num2str(cnum)
	       	<< "</TITLE>" << endl
	       	<< "<BODY BACKGROUND=\"file:" << getenv("FSLDIR") 
	    		<< "/doc/images/fsl-bg.jpg\">" << endl 
	     		<< "<hr><CENTER><H1>MELODIC Component " << num2str(cnum)
	       	<< "</H1>"<< endl;
	
				if(cnum>1)
	  			IChtml << "<a href=\"" << string("IC_")+num2str(cnum-1)
		 			<<".html\">previous</a>&nbsp;-&nbsp;";
	
				IChtml << "<a href=\"00index.html\">&nbsp;index&nbsp;</a>" ;
	
				if(cnum<dim)
	  			IChtml << "&nbsp;-&nbsp;<a href=\"" << string("IC_")+num2str(cnum+1)
		 				<<".html\">next</a><p>";
	
					IChtml << "<hr><p>" << endl;
      }

      volume4D<float> tempVol; 

      if(melodat.get_IC().Storage()>0)
			{//Output raw IC map

	  		tempVol.setmatrix(melodat.get_IC().Row(cnum),
			    melodat.get_mask());
	  		volume<float> map1;
	  		volume<float> map2;
	  		map1 = threshold(tempVol[0],float(0.0), 
			   	tempVol[0].max());
	  		map2 = threshold(tempVol[0],tempVol[0].min(), 
			   	float(-0.0));

	  		volume<float> newvol; 
	  		miscpic newpic;

	  		//	float map1min = (map1 + binarise(tempVol[0],tempVol[0].min(), 
	  		//		float(0.0)) * map1.max()).robustmin();
	  		float map1max = map1.percentile(0.99);
	  		float map2min = map2.percentile(0.01);
	  		//float map2max = (map2 + binarise(tempVol[0],float(0.0), 
	  		//		tempVol[0].max()) * map2.min()).robustmax();
	
	  		newpic.overlay(newvol, melodat.get_bg(), map1, map2, 
			 		float(0.0),
			 		float(0.0),
			 		float(0.01), map1max, float(-0.01), map2min, 
			 		0, 0);

	  		char instr[10000];
	
	  		sprintf(instr," ");
	  		strcat(instr,axials_instr.c_str());
	  		strcat(instr,string(report.appendDir(prefix+
					".png")).c_str());
	  		newpic.set_title(string("Component No. "+num2str(cnum)+
				  " - raw Z transformed IC map (1 - 99 percentile)"));
	  		newpic.set_cbar(string("ysb"));
	
	  		newpic.slicer(newvol, instr);
			}
     
      IChtml << "<img BORDER=0 SRC=\""+ prefix+
				".png\"><p>" << endl;
	
      {//plot time course
				miscplot newplot;
	
				if(opts.tr.value()>0.0)
	  			newplot.timeseries(melodat.get_Tmodes(cnum-1).t(),
			     	report.appendDir(string("t")+
				     	num2str(cnum)+".png"),
			     		string("Timecourse (in seconds); TR = ")+
			     		float2str(opts.tr.value(),0,2,0)+" s", 
			     		opts.tr.value(),150,4,1);
				else
	  			newplot.timeseries(melodat.get_Tmodes(cnum-1).t(),
			     	report.appendDir(string("t")+
					   	num2str(cnum)+".png"),
			     		string("Timecourse (in TRs)"));
				write_ascii_matrix(report.appendDir(string("t")
			 		+num2str(cnum)+".txt"),
			   	melodat.get_Tmodes(cnum-1));
				IChtml << "<A HREF=\"" << string("t")
	  			+num2str(cnum)+".txt" << "\"> ";
				IChtml << "<img BORDER=0 SRC=\"" 
	  			+string("t")+num2str(cnum)+".png\"></A><p>" << endl;
      }//time series plot
      
      {//plot frequency 
				miscplot newplot;
				int fact = int(std::pow(10.0,
					int(std::log10(float(melodat.get_Tmodes(0).Nrows())))));
	
				if(opts.tr.value()>0.0)
	  			newplot.timeseries(melodat.get_fmix().Column(cnum).t(),
			     	report.appendDir(string("f")+
					  num2str(cnum)+".png"),
			     	string("FFT of timecourse (in Hz / ") +
			     	num2str(fact)+")",
			     	fact/(opts.tr.value()*melodat.get_Tmodes(0).Nrows()),
			     	150,0,2);
				else
	  			newplot.timeseries(melodat.get_fmix().Column(cnum).t(),
			     	report.appendDir(string("f")+
			   		num2str(cnum)+".png"),
			     	string(string("FFT of timecourse (in cycles); ")
				    +"frequency(Hz)=cycles/("
				    +num2str(melodat.get_Tmodes(0).Nrows())
				    +"* TR); period(s)=("
				    +num2str(melodat.get_Tmodes(0).Nrows())
				    +"* TR)/cycles"));
				write_ascii_matrix(report.appendDir(string("f")
			 		+num2str(cnum)+".txt"),
			   	melodat.get_Tmodes(cnum-1));
				IChtml << "<A HREF=\"" << string("f")
	  			+num2str(cnum)+".txt" << "\"> ";
				IChtml << "<img BORDER=0 SRC=\"" 
	  			+string("f")+num2str(cnum)+".png\"></A><p>" << endl;
      }//frequency plot
 
      { //finish IC page
				IChtml<< "<HR><FONT SIZE=1>This page produced automatically by "
	      	<< "<A HREF=\"http://www.fmrib.ox.ac.uk/fsl/melodic/index.html\"> MELODIC</A> Version " 
	      	<< version << " - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL - "
	      	<< "FMRIB Software Library</A>.</FONT>" << endl
	      		<< "</BODY></HTML>" << endl;
      } //finish IC page 
    } 
  }

  void MelodicReport::PPCA_rep(){
    
    {//plot time course
      report << "<p><hr><b>PCA estimates </b> <p>" << endl;
 
      Matrix what;
      miscplot newplot;

      what  = melodat.get_EV();
      what &= melodat.get_EVP();

      newplot.add_label("ordered Eigenvalues");
      newplot.add_label("% of expl. variance");

      if(melodat.get_PPCA().Storage()>0){
				what = what.Columns(1,melodat.get_PPCA().Nrows());
				if(opts.allPPCA.value()&&melodat.get_PPCA().Ncols()==7){
					what &= melodat.get_PPCA().Columns(3,7).t();
					newplot.add_label("Laplace");
					newplot.add_label("BIC");
					newplot.add_label("MDL");
					newplot.add_label("RRN");
					newplot.add_label("AIC");
				}else{
					what &= melodat.get_PPCA().Column(1).t();
					newplot.add_label("dim. estimate");
				}
      }
			
			newplot.set_Ylabel_fmt("%.2f");
			newplot.add_xlabel("Number of included components");
			newplot.set_yrange(0.0,1.02);
      newplot.grid_swapdefault();
      newplot.timeseries(what,
			 	report.appendDir("EVplot.png"),
			 	string("Eigenspectrum Analysis"), 
			 	0,450,4,0);

      report << "<img BORDER=0 SRC=\"EVplot.png\"><p>" << endl;
    }//time series plot
  }
	
	void MelodicReport::Smode_rep(){
	  if(melodat.get_Smodes().size()>0){
		report << "<p><hr><b>TICA Subject/Session modes </b> <br>" << endl;
		miscplot newplot;
		report << "Boxplots show the relative response amplitudes across the "
			<< "session/subject domain (" << melodat.get_numfiles() 
			<< " input files). Components have been sorted in decreasing order of "
			<< " the median response per component. <br><br>";
			
		outMsize("Smode.at(0)", melodat.get_Smodes().at(0));
		Matrix allmodes = melodat.get_Smodes().at(0);
		for(int ctr = 1; ctr < (int)melodat.get_Smodes().size();++ctr)
			allmodes |= melodat.get_Smodes().at(ctr);
	
		outMsize("allmodes", allmodes);
		newplot.add_xlabel("Component No.");
		newplot.add_ylabel("");
		if(allmodes.Ncols()<100)
			newplot.set_xysize(100+30*allmodes.Ncols(),300);
		else
			newplot.set_xysize(1200,300);
		newplot.boxplot(allmodes,report.appendDir(string("bp_Smodes.png")),
  			string("Subject/Session modes"));
  		report << "<img BORDER=0 SRC=\"bp_Smodes.png\"><p>" << endl;
      }
	}
}
