/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "tsplotcode.h"
#include "filemanager.h"
#include "newmatio.h"
#include "modelfit.h"
#include "tracker.h"
#include "storage/timeseries.hpp"
#include <utility>
using namespace std;
#include <math.h>

struct ModelFit::Implementation {

  Implementation(const QString& fd):
    m_featDir(fd),  m_nevs(0), m_nftests(0),
    m_higherLevel(false)
  {
    m_curCope=1;
    m_curPE=0;

    COPE_PE=true;
    // this variable indicates Cope or PE(i) curve. True indicates Cope curve
    readModel(m_featDir + "/design.mat");
    for(unsigned short i=0; i<m_nevs; i++)
      m_contrastList.push_back(std::make_pair(QString("PE"+QString::number(i+1)), false));
    if(FileManager::checkFileExists(string(m_featDir.ascii()) + "/design.con"))
      readConFile(m_featDir + "/design.con");
    if(FileManager::checkFileExists(string(m_featDir.ascii()) + "/design.fts"))
      readFtsFile(m_featDir + "/design.fts");
    loadFiltFuncData(m_featDir + "/filtered_func_data");
    loadPeData(m_featDir + "/stats/pe");

    m_higherLevel = FileManager::checkFileExists(string(m_featDir.ascii()) + "/design.lev");
  }

  // multiply the whole matrix xI with the scalar value betaI; refer to newmat
  Matrix xI_BetaI(Matrix xI, float betaI)
  {  
    return xI * betaI;
  }

  void  readModel(const QString & filePath)
  {
    FILE *fp;

    if((fp=fopen(filePath.latin1(),"r")))
      {
	m_nevs = findKey(fp,"/NumWaves").toInt();
	m_npts = findKey(fp,"/NumPoints").toInt();

	m_designMatrix.ReSize(m_npts, m_nevs);

	findKey(fp,"/Matrix");

	double dataVal;

	for(unsigned int j=1;j<=m_npts;j++)
	  {
	    for(unsigned int i=1;i<=m_nevs;i++)
	      {
		fscanf(fp,"%lf", &dataVal);
		m_designMatrix(j,i) = dataVal;
	      }
	  }
	fclose(fp);
      }
    else
      {
	throw Exception("Failed to open design.mat file!");
      }
  }

  void readFtsFile(const QString & filepath)
  {
    ifstream is(filepath);
    if(!is)
      throw std::ios::failure("Couldn't open design f-test file!");
    
    char buffer[1024];

    is.getline(buffer, 1024);
    std::string s(buffer);

    // s should contain "%! VEST-Waveform file

    int n(0);

    while(s != string("/NumContrasts"))
      {
	is >> s;
	is >> n;
	if(s == string("/NumContrasts"))
	  m_nftests = n;
      }
  }

  void readConFile(const QString & filePath)
  {
    TRACKER("ModelFit::Implementation::readConFile(const QString & filePath)");

    FILE *fp;

    if((fp=fopen(filePath.latin1(),"r")))
      {
	unsigned int nevs = findKey(fp,"/NumWaves").toInt();
	m_ncontrasts      = findKey(fp,"/NumContrasts").toInt();

	if(nevs != m_nevs)
	  throw Exception("Inconsitent number of EVs in design.con file!");
	m_copeVectors.ReSize(m_ncontrasts, m_nevs);

	double dataVal;

	for(unsigned int j=1;j<=m_ncontrasts;j++)
	  {
	    MESSAGE(QString("Looking for /ContrastName%1").arg(j));
	    QString name = findKey(fp,QString("/ContrastName%1").arg(j));
	    MESSAGE(QString("/ContrastName%1 = %2").arg(j).arg(name));
	    if(!name.isEmpty()) {
	      MESSAGE("Appending named contrast");
	      m_contrastList.push_back(std::make_pair(QString("COPE%1 (%2)").arg(j).arg(name), true));
	    } else {
	      MESSAGE("Appending unnamed contrast");
	      m_contrastList.push_back(std::make_pair(QString("COPE%1").arg(j), true));
	    }
	  }
	
	findKey(fp,"/Matrix");

	for(unsigned int jj=1;jj<=m_ncontrasts;jj++)
	    for(unsigned int i=1;i<=m_nevs;i++)
	      {
		fscanf(fp,"%lf",&dataVal);
		m_copeVectors(jj,i) = dataVal;
	      }
	fclose(fp);
      }
    else
      {
        throw Exception("Failed to open contrast file!");
      }
  }

  QString findKey(FILE *fd, const QString & keyName )
  {
    char charArray[100];
    QString valueStr;
    bool keyFound(false);
  
    fseek(fd,0,SEEK_SET);
  
    while(!keyFound && fgets(charArray, 100, fd))
      {
	QString lineStr(charArray);
	if(lineStr.find(keyName) != -1)
	  {
	    valueStr = lineStr.remove(0,keyName.length());
	    keyFound = true;
	  }
      }	  
    //file pointer left at pos just after key
    return valueStr.simplifyWhiteSpace();
  }

  void loadPeData(const QString & filePath)
  {  
    TRACKER("ModelFit::Implementation::loadPeData");
  
    for(unsigned int n = 1; n <= m_nevs; n++)
      {
	QString peFile = filePath + QString::number(n);
	Image::Handle peData = Image::load((const char *)peFile);
	m_peData.push_back(peData->getVolume(0));
      }
  }

  void loadFiltFuncData(const QString & filePath)
  {  
    TRACKER("ModelFit::Implementation::loadFiltFuncData");

    m_filtFuncData = Image::load((const char *)filePath);
  }

  ColumnVector getPeColumnVector(short x, short y, short z)
  {  
    TRACKER("ModelFit::Implementation::getPeColumnVector");

    unsigned int nevs = m_peData.size();
    ColumnVector peCv(nevs);

    for(unsigned int n = 0; n < nevs; ++n)
      {
	peCv(n+1) = m_peData[n]->value(x,y,z);
      }

    return peCv;
  }

  QString m_featDir;

  Matrix m_designMatrix;
  Matrix m_copeVectors;
  std::vector<Volume::Handle> m_peData;
  Image::Handle m_filtFuncData;

  TimeSeries::Handle m_timeSeriesModel;
  
  //QStringList m_contrList;
  std::vector< std::pair<QString, bool> > m_contrastList;
  
  unsigned int m_npts, m_nevs, m_ncontrasts, m_nftests, m_curCope, m_curPE;
  
  std::list< ModelFitObserver *> m_modelFitObservers;
  
  bool COPE_PE;

  bool m_higherLevel;
};

ModelFit::Handle ModelFit::create(const QString& fd)
{
  Handle dst;
  try {
    ModelFit::Handle newModel(new ModelFit(fd));
    if(newModel)
      dst = newModel;
  } catch(...) {
    throw;
  }
  return dst;
}

ModelFit::ModelFit(const QString& fd):   
  m_impl(new ModelFit::Implementation(fd))
{
}

ModelFit::~ModelFit(){}

QString& ModelFit::featDir() const { return m_impl->m_featDir; }

TimeSeries::Handle ModelFit::getDataTimeSeries(short x, short y, short z)
{  
  TRACKER("ModelFit::getDataTimeSeries");

  return m_impl->m_filtFuncData->getTimeSeries(x,y,z);
}

Image::Handle ModelFit::getFilteredFuncImage()
{
  TRACKER("ModelFit::getFilteredFuncImage()");

  return m_impl->m_filtFuncData;
}

TimeSeries::Handle ModelFit::fullModel(short x, short y, short z, float mean)
{
  Matrix colVector(m_impl->m_npts, 1), result(m_impl->m_npts, 1);
  colVector=0.0; result=0.0;
  
  for(unsigned int i=1; i<=m_impl->m_nevs; i++)  // full model mean + x1beta1 + x2beta2 + x3beta3
    {
      colVector = m_impl->m_designMatrix.Column(i);
    
      result += m_impl->xI_BetaI(colVector, m_impl->m_peData[i-1]->value(x, y, z)); 
    }
  
  if(!m_impl->m_higherLevel)
    result += mean;
  
  TimeSeries::Handle pePlotCurve = TimeSeriesD::create(x, y, z, m_impl->m_npts);
  
  for(unsigned int n = 1; n <= m_impl->m_npts; n++)
    pePlotCurve->setValue(n-1, (float)result(n, 1));
    
  return pePlotCurve;
}

TimeSeries::Handle ModelFit::CopeCurve(short x, short y, short z, float mean)
{
  TimeSeries::Handle copeCurve = TimeSeriesD::create(x, y, z, m_impl->m_npts);
  
  Matrix colVector(m_impl->m_npts, 1), result(m_impl->m_npts, 1);
  result = 0.0;
  
  for(unsigned int i=1; i<=m_impl->m_nevs; i++)
    {
      colVector = m_impl->m_designMatrix.Column(i);
    
      result += (m_impl->xI_BetaI(colVector, m_impl->m_peData[i-1]->value(x, y, z)) * 
		 m_impl->m_copeVectors(curFit()-m_impl->m_nevs+1, i));
    }
  
  if(!m_impl->m_higherLevel)
    result += mean;
  
  for(unsigned int n = 1; n <= m_impl->m_npts; n++)
    copeCurve->setValue(n-1, (float)result(n, 1));
    
  return copeCurve;
}

//!
//! @brief Evaluates the curve PE(i) where i is equal to the value of variable m_curPE
//!        m_curPE is in sync with the PE selected surrently
//!        if PE1 is selected the curve plotted will be PE(1) and similarly for other PEs
//! @return TimeSeries::Handle which can be appended to @ref CurveDataList which stores
//!         a list of curves that are to be plotted
//!
TimeSeries::Handle ModelFit::peCurve(short x, short y, short z, float mean)
{
  TimeSeries::Handle peCurve = TimeSeriesD::create(x, y, z, m_impl->m_npts);
  
  Matrix result(m_impl->m_npts, 1);
  
  result = (m_impl->xI_BetaI(m_impl->m_designMatrix.Column(m_impl->m_curPE+1),
			     m_impl->m_peData[m_impl->m_curPE]->value(x, y, z))) + mean;
  
  for(unsigned int n = 1; n <= m_impl->m_npts; n++)
    peCurve->setValue(n-1, (float)result(n, 1));
    
  return peCurve;
  
}

TimeSeries::Handle ModelFit::perCentChange(short x, short y, short z, float mean)
{
  TimeSeries::Handle perCentChange = TimeSeriesD::create(x, y, z, m_impl->m_npts);
  
  Matrix result(m_impl->m_npts, 1);
  
  result = (m_impl->m_designMatrix.Column(m_impl->m_curPE+1) * 100) / mean;
  
  for(unsigned int n = 1; n <= m_impl->m_npts; n++)
    perCentChange->setValue(n-1, (float)result(n, 1));
    
  return perCentChange;
}

//!
//! @brief Checks whether the directory passed to it is a valid FEAT directory or not.
//! @param path The path to be tested directory
//! @return true if input is a .feat directory otherwise returns False
//!
bool ModelFit::isFeatDir(const QString& path)
{
  bool valid(false);
  // get the directory immediately preceding filename and check if it is .feat
  // This guards against opening an image from a subdirectory within .feat directory
  // and detecting that directory as a current feat directory
  if((path.section("/", -2, -2).findRev(".feat"))>-1) {
    valid = true;
  }
   
  return valid;
}

QString ModelFit::getConName(unsigned int i) const
{
  if( (i < 0) || (i >= numFits()) )
    throw Exception("Invalid contrast index");

  return m_impl->m_contrastList[i].first;
}

unsigned int ModelFit::curFit(void) const
{
  return m_impl->m_curPE;
}

void ModelFit::curFit(unsigned int i)
{
  if( (i < 0) || (i >= numFits()) )
    throw Exception("Invalid contrast index");

  m_impl->m_curPE=i;
  
  notify();
}

void ModelFit::attach(ModelFitObserver *o)
{
  m_impl->m_modelFitObservers.remove(o);
  m_impl->m_modelFitObservers.push_back(o);
}

void ModelFit::detach(ModelFitObserver *o)
{
  m_impl->m_modelFitObservers.remove(o);
}

struct Update
{
public:
  Update(ModelFit *m): m_modelFit(m) {} 
    
  void operator ()(ModelFitObserver *m_o)
  {
    m_o->update(m_modelFit);
  }
private:
  ModelFit *m_modelFit;
};

void ModelFit::notify()
{
  for_each(m_impl->m_modelFitObservers.begin(), m_impl->m_modelFitObservers.end(), Update(this));
}

unsigned int ModelFit::numEVs(void) const
{
  return m_impl->m_nevs;
}

unsigned int ModelFit::numContrasts(void) const
{
  return m_impl->m_ncontrasts;
}

unsigned int ModelFit::numFtests(void) const
{
  return m_impl->m_nftests;
}

unsigned int ModelFit::numFits(void) const
{
  return m_impl->m_ncontrasts + m_impl->m_nevs;
}

bool ModelFit::copePe(void) const
{
  return m_impl->m_contrastList[curFit()].second;
}
