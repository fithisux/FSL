/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#include "clusterbrowser.h"
#include "clusterdata.h"
#include "filemanager.h"

#include <iostream>
#include <sstream>

#include <qcheckbox.h>
#include <qlistview.h>
#include <qcombobox.h>

class ClusterListItem: public QListViewItem
{
public:
  ClusterListItem(QListView* v, Cluster::Handle c, TalairachCluster::Handle t, bool s):
    QListViewItem(v),
    m_showTalairach(s), m_cluster(c), m_talCluster(t)
  {
    refresh();
  }

  void refresh()
  {
    if( m_cluster->initialised() ) {
      setText( 0, m_cluster->inqIndex());
      setText( 1, m_cluster->inqSize());
      setText( 2, m_cluster->inqP());
      setText( 3, m_cluster->inqMinusLog10P()); 
      setText( 4, m_cluster->inqMaxZ());
      setText(11, m_cluster->inqMaxCOPE());
      setText(15, m_cluster->inqMeanCOPE());
   } else {
      setText( 0, m_talCluster->inqIndex());
      setText( 1, m_talCluster->inqSize());
      setText( 2, m_talCluster->inqP());
      setText( 3, m_talCluster->inqMinusLog10P());
      setText( 4, m_talCluster->inqMaxZ());
      setText(11, m_talCluster->inqMaxCOPE());
      setText(15, m_talCluster->inqMeanCOPE());
    }

    BaseCluster::Handle temp;
    if(m_showTalairach)
      temp = m_talCluster;
    else temp = m_cluster;

    setText( 5, temp->inqMaxZx());
    setText( 6, temp->inqMaxZy());
    setText( 7, temp->inqMaxZz());
    setText( 8, temp->inqMaxCOGx());
    setText( 9, temp->inqMaxCOGy());
    setText(10, temp->inqMaxCOGz());

    setText(12, temp->inqMaxCOPEx());
    setText(13, temp->inqMaxCOPEy());
    setText(14, temp->inqMaxCOPEz());
  }

  int compare(QListViewItem *i, int col, bool ascending) const
  {
    return key(col, ascending).toFloat() - i->key(col, ascending).toFloat();
  }

  Cluster::Handle getCluster()  const { return m_cluster; }
  TalairachCluster::Handle getTCluster() const { return m_talCluster; }

private:

  bool m_showTalairach;
  Cluster::Handle m_cluster; 
  TalairachCluster::Handle m_talCluster;
};

ClusterBrowser::ClusterBrowser(QWidget* parent, Image::Handle i,
			       Cursor::Handle c, ModelFit::Handle m):
  ClusterBrowserBase(parent, 0, WDestructiveClose), m_currentSelection(0),
  m_imageInfo(i->getInfo()), m_cursor(c), m_model(m)
{
  try {
    m_statComboBox->clear();

    for(unsigned int i = 1; i <= m_model->numContrasts(); ++i) {
      ostringstream basename;
      basename << m_model->featDir() << '/' << "cluster_zstat" << i;

      ClusterList clusters, tclusters;
      if( FileManager::checkFileExists(basename.str() + ".txt") )
	FileManager::readClusters(basename.str() + ".txt", clusters);
      if( FileManager::checkFileExists(basename.str() + "_std.txt") )
	FileManager::readTalairachClusters(basename.str() + "_std.txt", tclusters);
      else if( FileManager::checkFileExists(basename.str() + "_tal.txt") )
	FileManager::readTalairachClusters(basename.str() + "_tal.txt", tclusters);

      ClusterListPair cp(std::make_pair(clusters, tclusters));

      ostringstream name;
      name << "zstat" << i;
      m_clusterTables.push_back(std::make_pair(name.str(), cp));
      m_statComboBox->insertItem(name.str());
    }

    for(unsigned int i = 1; i <= m_model->numFtests(); ++i) {
      std::ostringstream basename;
      basename << m_model->featDir() << '/' << "cluster_zfstat" << i;

      ClusterList clusters, tclusters;
      if(FileManager::checkFileExists(basename.str() + ".txt") )
	FileManager::readClusters(basename.str() + ".txt", clusters);
      if( FileManager::checkFileExists(basename.str() + "_std.txt") )
	FileManager::readTalairachClusters(basename.str() + "_std.txt", tclusters);
      else if( FileManager::checkFileExists(basename.str() + "_tal.txt") )
	FileManager::readTalairachClusters(basename.str() + "_tal.txt", tclusters);

      ClusterListPair cp(std::make_pair(clusters, tclusters));
      
      ostringstream name;
      name << "zstatf" << i;
      m_clusterTables.push_back(std::make_pair(name.str(), cp));
      m_statComboBox->insertItem(name.str());
    }
    
    selectStatistic(0);

  } catch (const std::ios::failure& e) {
    throw ClusterBrowser::Exception(std::string("ClusterBrowser::ClusterBrowser Couldn't initialise\n") + e.what());
  } catch (...) {
    throw;
  }
}

void ClusterBrowser::showTalairach(bool s)
{
  m_showTalairach = s;
  selectStatistic(m_currentSelection);
}

void ClusterBrowser::selectStatistic(int n)
{
  m_currentSelection = n;
  m_clusterListView->clear();
  for(int i = 0; i < m_clusterListView->columns(); ++i) {
    m_clusterListView->setColumnWidth(i, 0);
    m_clusterListView->setColumnWidthMode(i, QListView::Maximum);
  }
#if (QT_VERSION < 0x030200)
  m_clusterListView->setSorting(1, false);
#else
  m_clusterListView->setSortColumn(1);
  m_clusterListView->setSortOrder(Qt::Descending);
#endif

  ClusterTable& t(m_clusterTables.at(n));
  ClusterListPair& p(t.second);
  ClusterList& clusters(p.first);
  ClusterList& tclusters(p.second);

  ClusterList::iterator ti = tclusters.begin();
  ClusterList::iterator ci = clusters.begin();
  
  while( (ci != clusters.end()) ||
	 (ti != tclusters.end()) ) {

    bool allowPlain(false);
    Cluster::Handle c = Cluster::create();     
    if( ci != clusters.end() ) {
      c = boost::dynamic_pointer_cast<Cluster>(*ci);
      ++ci;
      allowPlain = true;
    }

    bool allowTal(false); 
    TalairachCluster::Handle t = TalairachCluster::create();
    if( ti != tclusters.end() ) {
      t = boost::dynamic_pointer_cast<TalairachCluster>(*ti);
      ++ti;
      allowTal = true;
    }

    bool displayTal( (m_showTalairach && allowTal) || !allowPlain );
    new ClusterListItem(m_clusterListView, c, t, displayTal);

    if( !(allowPlain && allowTal) ) {
      m_talairachCheckBox->setDisabled(true);
      if( allowTal )
	m_talairachCheckBox->setChecked(true);
    } else
      m_talairachCheckBox->setEnabled(true);
  }
}

void ClusterBrowser::clusterListSelectionChanged(QListViewItem *item)
{
  if(ClusterListItem * cl = dynamic_cast<ClusterListItem*>(item)) {

    if( !m_talairachCheckBox->isChecked() ) {
      Cluster::Handle cluster  = cl->getCluster();
      cluster->setCursorToMaxZ(m_cursor);
    } else {
      TalairachCluster::Handle tcluster = cl->getTCluster();
      tcluster->setCursorToMaxZ(m_imageInfo, m_cursor);
    }
  }
}

void ClusterBrowser::closeEvent(QCloseEvent* e)
{
  emit windowClose(e);
}
