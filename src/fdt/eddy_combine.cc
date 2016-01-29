//     eddy_combine.cc : In case JAC resampling is used in the eddy output, then no pairs of phase encoding directions are combined. 
//                       This is performed by this tool that combines (by taking the average of) pairs of phase encoding directions.
//                       Input text files provide the information needed to find out which volumes to combine.
//                       Pos_SeriesVol and Neg_SeriesVol are text files produced by DiffPreprocPipeline.sh    
//     Stamatios Sotiropoulos, FMRIB Analysis Group, 2012

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "newimage/newimageall.h"

using namespace NEWIMAGE;

void print_usage(const string& progname) {
  cout << endl;
  cout << "Usage: eddy_combine <PosData> <Posbvals> <Posbvecs> <Pos_SeriesVol> <NegData> <Negbvals> <Negbvecs> <Neg_SeriesVol> <outputpath> <onlymatched_flag>" << endl;
  cout << "If onlymatched_flag=1, then only the volumes that have matched Pos/Neg pairs will be included in the output" << endl;
  cout << endl;
}


//Copy a single volume from one 4D image to another 
void CopyVolume(const volume4D<float>& In,volume4D<float>& Out, const int InVolume, const int OutVolume){
  for(int k=0;k<In.zsize();k++)
    for(int j=0;j<In.ysize();j++)	    
      for(int i=0;i<In.xsize();i++)
	Out.value(i,j,k,OutVolume)=In.value(i,j,k,InVolume);
} 

//Copy the mean of two volumes volume from one 4D image to another
void CopyMeanVolume(const volume4D<float>& In1,const volume4D<float>& In2, volume4D<float>& Out, const int InVolume1, const int InVolume2,const int OutVolume){
  for(int k=0;k<In1.zsize();k++)
    for(int j=0;j<In1.ysize();j++)	    
      for(int i=0;i<In1.xsize();i++)
	Out.value(i,j,k,OutVolume)=0.5*(In1.value(i,j,k,InVolume1)+In2.value(i,j,k,InVolume2));
} 


int main(int argc,char *argv[])
{
  try {
    string progname=argv[0];
    if (argc !=11){ 
	print_usage(progname);
	return 1; 
    }
    //Read Input 
    volume4D<float> PosData, NegData;
    Matrix PosVols, NegVols, Posbvals, Posbvecs, Negbvals, Negbvecs;

    read_volume4D(PosData,string(argv[1]));
    Posbvals=read_ascii_matrix(string(argv[2]));
    Posbvecs=read_ascii_matrix(string(argv[3]));
    PosVols=read_ascii_matrix(string(argv[4]));
    read_volume4D(NegData,string(argv[5]));
    Negbvals=read_ascii_matrix(string(argv[6]));
    Negbvecs=read_ascii_matrix(string(argv[7]));
    NegVols=read_ascii_matrix(string(argv[8]));
    string oname=string(argv[9]);
    int flag=atoi(argv[10]);

    //Initialise Output    
    int tsize=0;  int xsize(PosData.xsize());  int ysize(PosData.ysize());  int zsize(PosData.zsize());
    if (flag==0)
      tsize=Posbvals.Ncols()+Negbvals.Ncols()-PosVols.Column(1).Sum(); //number of total input volumes minus the matched volumes
    else
      tsize=PosVols.Column(1).Sum(); //number of matched volumes
    volume4D<float> Output_Data(xsize,ysize,zsize,tsize);
    Matrix Output_bvals(1,tsize); Matrix Output_bvecs(3,tsize);
    ColumnVector CombinedVols(tsize); CombinedVols=0.0;
    Output_Data.copyproperties(PosData); //copy header info

    int OutVolIndex=0; int PosVolIndex=0; int NegVolIndex=0;
    for (int l=1; l<=PosVols.Nrows(); l++){ //For each series pair
      for (int offset=0; offset<PosVols(l,1); offset++){ //Handle corresponding volumes
	CopyMeanVolume(PosData,NegData,Output_Data,PosVolIndex,NegVolIndex,OutVolIndex);
	Output_bvals(1,OutVolIndex+1)=Posbvals(1,PosVolIndex+1);
	Output_bvecs.Column(OutVolIndex+1)=Posbvecs.Column(PosVolIndex+1);
	CombinedVols(OutVolIndex+1)=1;  //Indicate in this array that the respective volume is a combination of two input volumes
	PosVolIndex++; 	NegVolIndex++; OutVolIndex++;
      }
      for (int offset=0; offset<(PosVols(l,2)-PosVols(l,1)); offset++){ //Handle non-corresponding Pos volumes
	if (flag==0){  //Include in the Output if requested
	  CopyVolume(PosData,Output_Data,PosVolIndex,OutVolIndex);
	  Output_bvals(1,OutVolIndex+1)=Posbvals(1,PosVolIndex+1);
	  Output_bvecs.Column(OutVolIndex+1)=Posbvecs.Column(PosVolIndex+1);
	  OutVolIndex++;
	}
	PosVolIndex++; 	
      }
      for (int offset=0; offset<(NegVols(l,2)-NegVols(l,1)); offset++){ //Handle non-corresponding Neg volumes
	if (flag==0){ //Include in the Output if requested
	  CopyVolume(NegData,Output_Data,NegVolIndex,OutVolIndex);
	  Output_bvals(1,OutVolIndex+1)=Negbvals(1,NegVolIndex+1);
	  Output_bvecs.Column(OutVolIndex+1)=Negbvecs.Column(NegVolIndex+1);
	  OutVolIndex++;
	}
	NegVolIndex++; 	
      }
    }

  string onameData=oname+"/data";
  save_volume4D(Output_Data,onameData);
  onameData=oname+"/bvals";
  write_ascii_matrix(Output_bvals,onameData);
  onameData=oname+"/bvecs";
  write_ascii_matrix(Output_bvecs,onameData);
  if (flag==0){   //If the output file has a mixture of combined (flag=1) and uncombined (flag=0) volumes, store this info in a text file
    onameData=oname+"/dataCombinedFlag.txt";
    write_ascii_matrix(CombinedVols,onameData);
  }
  
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  }  catch(Exception &e) {
    exit(EXIT_FAILURE);
  } 
}


