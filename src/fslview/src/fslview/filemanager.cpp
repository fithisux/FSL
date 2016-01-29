/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "tracker.h"

#include "lookuptable.h"
#include "filemanager.h"

#include <cstdio>
#include <ctype.h>

#include <iostream>
#include <fstream>
#include <string>
#include <limits>

#include <sstream>

#include <qobject.h>
#include <qstringlist.h>
#include <qfile.h>
#include <qdir.h>
#include <qdom.h>

using namespace std;

FileManager::FileManager()
{ 
}

FileManager::~FileManager()
{
}

struct Output
{
  Output(ostream& os) : m_os(os) {}
  void operator()(const BaseCluster::Handle& h)
  {
    m_os << h << std::endl;
  }
  ostream& m_os;
};

//! @brief Read clusters from cluster_*.txt files in a FEAT directory
//! @param filename The filename in which the clusters should be found
//! @param clusters A handle for the @ref ClusterDataList which will hold all the data
void FileManager::readClusters(const std::string& filename, ClusterList& clusters)
{
  STATIC_TRACKER("FileManager::readClusters(const std::string&, ClusterList&)");
  clusters.clear();

  MESSAGE("Reading " + filename);
  std::ifstream cf(filename.c_str());

  if(!cf)
    throw std::ios::failure("FileManager::readClusters failed to open file: " + filename);

  std::list<string> headers;
				// Accumulate a list of header names
  string header;
  getline(cf, header);
  istringstream ss(header);
  string name;
  while(getline(ss, name, '\t'))
    {
      //      cout << name << endl;
      headers.push_back(name);
    }

  				// ...read the remaining cluster data.
  BaseCluster::Handle c = Cluster::create(headers);
  cf >> c;
  while(!cf.eof())
    {
      //     std::cout << c << std::endl;
      clusters.push_back(c);

      c = Cluster::create(headers);
      cf >> c;
   }
}

//! @brief Read clusters from cluster_*.txt files in a FEAT directory
//! @param filename The filename in which the clusters should be found
//! @param clusters A handle for the @ref ClusterDataList which will hold all the data
void FileManager::readTalairachClusters(const std::string& filename, ClusterList& clusters)
{
  STATIC_TRACKER("FileManager::readTalairachClusters(const std::string&, ClusterList&)");
  clusters.clear();

  MESSAGE("Reading " + filename);
  std::ifstream cf(filename.c_str());

  if(!cf)
    throw std::ios::failure("FileManager::readClusters failed to open file: " + filename);

  std::list<string> headers;
				// Accumulate a list of header names
  string header;
  getline(cf, header);
  istringstream ss(header);
  string name;
  while(getline(ss, name, '\t'))
    {
      //      cout << name << endl;
      headers.push_back(name);
    }

				// ...read the remaining cluster data.
  BaseCluster::Handle c = TalairachCluster::create(headers);
  cf >> c;
  while(!cf.eof())
    {
      //      std::cout << c << std::endl;
      clusters.push_back(c);

      c = TalairachCluster::create(headers);
      cf >> c;
    }
}

//! @brief Read a LUT style LookUpTable from file
//!
//! @param filename Name of the file in which the lookuptable can
//!        be found.
//! @param lut Pointer to a LookUpTable to be populated from the file
// void  FileManager::readLutFile(const std::string& filename, LookUpTable* lut)
// {
//   vector<double> red, green, blue;
//   lut_actions a(red, green, blue);
//   lut_grammar lg(a);

//   parse(lg, filename.c_str());

//   vector<double>::const_iterator r = red.begin();
//   vector<double>::const_iterator g = green.begin();
//   vector<double>::const_iterator b = blue.begin();
	
//   lut->allocateMemory(red.size());

//   for(;r != red.end();)
//     {
//       cout << *r << "," << *g << "," << *b << endl;
//       ++r; ++g; ++b;
//     }
// }

void  FileManager::readLutFile(const std::string& filename, LookUpTable* lut)
{
  float red, blue, green;
  char chLine[64] = "";
  int linenumber = 0,count = 0;
  
  FILE* fp = fopen(filename.c_str(),"r");
  // FILE* out = fopen("render3t.ldt","w");
  if(!fp)
    throw std::ios::failure("Couldn't open lut file!");

  fseek(fp,0,SEEK_SET);

  // Should start %!VEST-LUT

  fgets(chLine,64,fp);
  if(!strcmp(chLine, "%!VEST-LUT"))
    {
      fclose(fp);
      throw std::ios::failure("File is not a valid lut file! Should start: %!VEST-LUT");
    }

  while(fgets(chLine,64,fp)){if(strncmp(chLine,"<-color{",8) == 0) ++count;} 

  lut->allocateMemory(count);

  fseek(fp,0,SEEK_SET);

  while(fgets(chLine,64,fp))
  {
    if(strncmp(chLine,"<-color{",8) == 0)
      {
      sscanf(chLine + 8,"%f , %f , %f",&red,&green,&blue);
      lut->pushValue((unsigned char)(red*255.0),(unsigned char)(green*255.0), (unsigned char)(blue*255.0),
                     linenumber);
      //    fprintf(out,"\"%f,%f,%f\",\n",red,green,blue);
      ++linenumber;
      }    
  }
  fclose(fp);
  // fclose(out);
}  

//! @brief Read an RGB style LookUpTable from file
//!
//! @param filename Name of the file in which the lookuptable can
//!        be found.
//! @param lut Pointer to a LookUpTable to be populated from the file
void  FileManager::readRgbFile(const std::string& filename, LookUpTable* lut)
{
  int  index, red, blue, green,
       idxVal(0),idxValMax(0);
  char name[64]   = "";
  char chLine[128] = "";
  char idxStr[10] = "";
  char firstChar;

  FILE* fp = fopen(filename.c_str(),"r");
  if(!fp)
    throw std::ios::failure("Couldn't open rgb file!");

  fseek(fp,0,SEEK_SET);

  fgets(chLine,64,fp);
  if(!strcmp(chLine, "%!VEST-LUT"))
    {
      fclose(fp);
      throw std::ios::failure("File is not a valid rgb file!");
    }

  while(fgets(chLine,128,fp))
    {
      sscanf(chLine," %c",&firstChar);

      if(isdigit(firstChar))
        {      
          sscanf(chLine," %s",idxStr);
          idxVal = atoi(idxStr);
          if (idxVal > idxValMax)idxValMax = idxVal;
        }

    }

  lut->allocateMemory(idxValMax + 1);

  fseek(fp,0,SEEK_SET);

  while(fgets(chLine,128,fp))
  {
    sscanf(chLine," %c",&firstChar);
    if(isdigit(firstChar))
      {
        sscanf(chLine,"%i %s %i %i %i",&index,name,&red,&green,&blue);
        lut->pushValue((unsigned char)red,(unsigned char)green,(unsigned char)blue, std::string(name), index);
      }
  }

  fclose(fp);
}  

LutElement readColourNode( QDomNode &node )
{
  QStringList colours( QStringList::split(",", node.firstChild().nodeValue()) );
  
  LutElement elem( node.toElement().attribute( "index", "-1" ).toUInt(),
		   node.toElement().attribute( "label", "Unknown" ) );
  switch( colours.size() )
    {
    case 3:
      elem.setColours( colours[0].toUInt(),
		       colours[1].toUInt(),
		       colours[2].toUInt() );
      break;
    case 4:
      elem.setColours( colours[0].toUInt(),
		       colours[1].toUInt(),
		       colours[2].toUInt(),
		       colours[3].toUInt() );
      break;
    default:
      break;
    }

  return elem;
}

//! @brief Read an LML style LookUpTable from file
//!
//! @param filename Name of the file in which the lookuptable can
//!        be found.
//! @param lut Pointer to a LookUpTable to be populated from the file
void FileManager::readLMLFile(const std::string& filename, LookUpTable* lut)
{
    // read the XML file and create DOM tree
    QFile opmlFile( filename );
    if ( !opmlFile.open( IO_ReadOnly ) )
      throw std::ios::failure( QObject::tr("FileManager", "Cannot open file %1").arg(filename) );

    QString errorMsg, errorLine;
    QDomDocument domTree;
    if ( !domTree.setContent( &opmlFile, &errorMsg, &errorLine ) ) {
      opmlFile.close();
      throw std::ios::failure( QObject::tr("FileManager", "Parsing error for file %1\n%2\nLine:%3").arg(filename).arg(errorMsg).arg(errorLine) );
    }
    opmlFile.close();

    // get the header information from the DOM
    QDomElement root = domTree.documentElement();
    QDomNode node;
    node = root.firstChild();
    while ( !node.isNull() ) {
      if ( node.isElement() && node.nodeName() == "header" ) {
	QDomNode currentNode = node.toElement().firstChild();
	while ( !currentNode.isNull() )
	  if ( currentNode.isElement() ) {
            // case for the different header entries
            if ( currentNode.nodeName() == "name" ) {
	      QDomText textChild = currentNode.firstChild().toText();
	      if ( !textChild.isNull() ) {
				// Process header->name here
		lut->setLutName(textChild.data());
	      }
	    } else if ( currentNode.nodeName() == "type" ) {
	      QDomText textChild = currentNode.firstChild().toText();
	      if ( !textChild.isNull() ) {
				// Process header->name here
		if(textChild.data() == "Indexed") 
		  lut->isIndexLut(true);
	      }
	    } else if ( currentNode.nodeName() == "colour_type" ) {
	    }
	    currentNode = currentNode.nextSibling();
	  }
	break;
      }
      node = node.nextSibling();
    }
    node = root.firstChild();
    while ( !node.isNull() ) {
      if ( node.isElement() && node.nodeName() == "volumes" ) {
	QDomNode currentNode = node.toElement().firstChild();
	while ( !currentNode.isNull() )
	  if ( currentNode.isElement() ) {
            if ( currentNode.nodeName() == "label" ) {
	      int n    = currentNode.toElement().attribute( "index" ).toInt();
	      string l = currentNode.firstChild().toText().data();
	      cout << n << " " << l << endl;
	      lut->setVolumeName(n, l);
            }
	    currentNode = currentNode.nextSibling();
	  }
	break;
      } else if ( node.isElement() && node.nodeName() == "data" ) {
	lut->allocateMemory(0);
	QDomNode currentNode = node.toElement().firstChild();
	while ( !currentNode.isNull() )
	  if ( currentNode.isElement() ) {
            if ( currentNode.nodeName() == "colour" ) {
	      LutElement e(readColourNode(currentNode));
	      lut->pushValue(e);
            }
	    currentNode = currentNode.nextSibling();
	  }
	break;
      }
      node = node.nextSibling();
    }
}

Image::Handle readImage(const string& atlasdir, QDomNode& imagesNode, const string& tag)
{
  QDomNode node(imagesNode.toElement().firstChild());
  Image::Handle im;
  string filename;

  while ( !node.isNull() ) {
    if ( node.isElement() && node.nodeName() == tag ) {
      QDomText textChild = node.firstChild().toText();
      if ( !textChild.isNull() )
	filename = string(textChild.data().ascii());
      im = Image::load(atlasdir + "/" + filename, false);
    }
    node = node.nextSibling();
  }
  return im;
}

Atlas::Handle FileManager::readXMLAtlas(const string& atlasdir, 
					const string& filename)
{
  string fullpath(atlasdir + "/" + filename);
  QFile opmlFile( fullpath );
  if ( !opmlFile.open( IO_ReadOnly ) )
    throw 
      std::ios::failure( QObject::tr("Couldn't open file %1").arg(fullpath) );
  QString errorMsg, errorLine;
  QDomDocument domTree;
  if ( !domTree.setContent( &opmlFile, &errorMsg, &errorLine ) ) {
    opmlFile.close();
    throw 
      std::ios::failure( QObject::tr("Parsing error for file %1\n%2\nLine:%3")
			 .arg(filename).arg(errorMsg).arg(errorLine) );
  }
  opmlFile.close();


  // get the header information from the DOM
  QDomElement root = domTree.documentElement();
  QDomNode node;
  
  string atlasname("Unset");
  Atlas::Type type(Atlas::Unknown);

  Atlas::ImageStore images, summaries;

  node = root.firstChild();
  while ( !node.isNull() ) {
    if ( node.isElement() && node.nodeName() == "header" ) {
      QDomNode currentNode = node.toElement().firstChild();
      while ( !currentNode.isNull() )
	if ( currentNode.isElement() ) {
	  // case for the different header entries
	  if ( currentNode.nodeName() == "name" ) {
	    QDomText textChild = currentNode.firstChild().toText();
	    if ( !textChild.isNull() ) {
	      atlasname = string(textChild.data().ascii());
	    }
	  } else if ( currentNode.nodeName() == "type" ) {
	    QDomText textChild = currentNode.firstChild().toText();
	    if ( !textChild.isNull() ) {
	      if(textChild.data() == "Label") 
		type = Atlas::Label;
	      else
		type = Atlas::Probabilistic;
	    }
	  } else if ( currentNode.nodeName() == "images" ) {
	    images.push_back(readImage(atlasdir, currentNode, "imagefile"));
	    summaries.push_back(readImage(atlasdir, currentNode, "summaryimagefile"));
	  }
	  currentNode = currentNode.nextSibling();
	}
      break;
    }
    node = node.nextSibling();
  }

  Atlas::Handle atlas;
      
  switch(type)
    {
    case Atlas::Probabilistic:
      atlas = ProbabilisticAtlas::create(images, summaries, atlasname); 
      break;
    case Atlas::Label:
      atlas = LabelAtlas::create(images, summaries, atlasname);
      break;
    default:
      break;
    }

  node = root.firstChild();
  while ( !node.isNull() ) {
    if ( node.isElement() && node.nodeName() == "data" ) {
      QDomNode currentNode = node.toElement().firstChild();
      while ( !currentNode.isNull() )
	if ( currentNode.isElement() ) {
	  if ( currentNode.nodeName() == "label" ) {
	    int n    = currentNode.toElement().attribute( "index" ).toInt();
	    string l = currentNode.firstChild().toText().data();
	    atlas->addLabel(n, l);
	    int x = currentNode.toElement().attribute( "x", "0" ).toInt();
	    int y = currentNode.toElement().attribute( "y", "0" ).toInt();
	    int z = currentNode.toElement().attribute( "z", "0" ).toInt();
	    int v = currentNode.toElement().attribute( "v", "0" ).toInt();
	    string ref = currentNode.toElement().attribute( "ref", "" );
	    atlas->addCentre(n, x, y, z, v);
	    //atlas->addReference(n, ref);
	  }
	  currentNode = currentNode.nextSibling();
	}
      break;
    }
    node = node.nextSibling();
  }

  return atlas;
}

bool FileManager::checkFileExists(const std::string& path)
{
  bool exists(false);
  std::ifstream is(path.c_str());
  if(is)
    exists = true;
  return exists;
}

vector<string> FileManager::getFilenames(const std::string& path, const std::string& glob)
{
  QDir d(path, glob);

  vector<string> result;

  for(unsigned int i = 0; i < d.count(); ++i)
    result.push_back(d[i]);

  return  result;
}
