/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "overlayinfodialog.h"
#include "imagedisplaysetting.h"
#include <qapplication.h>
#include <qlabel.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qvalidator.h>

#include "tracker.h"
#include "assistantclient.h"

//#include "fileopen.xpm"

class AddLutToCombo {
public:
  AddLutToCombo(QComboBox* c,LookUpTable::Handle l):m_cb(c),m_curLut(l),m_index(0)
  {
    m_cb->clear();
  }
  void operator()(LookUpTable::Handle l)
  { 
    if(l->isVisible()) {
      m_cb->insertItem(l->inqLutName().c_str());
      if(l == m_curLut)
	m_cb->setCurrentItem(m_index);
      ++m_index;
    }
  }
 
private:
  QComboBox* m_cb;
  LookUpTable::Handle m_curLut;
  int m_index;
};


class AddModImageToCombo {
public:
  AddModImageToCombo(QComboBox* c,Image::Handle i):
    m_cb(c),m_curModImg(i),m_index(0)
  {
    m_cb->clear();
    m_cb->insertItem("None");++m_index;
    
  }
  void operator()(Image::Handle i)
  { 
    m_cb->insertItem(i->getInfo()->inqImageName().c_str());  
    if(m_curModImg.get())
      {
        if(i == m_curModImg){m_cb->setCurrentItem(m_index);}
      }
    else
      {
        m_cb->setCurrentItem(0);
      }
    ++m_index;
  }
 
private:
  QComboBox* m_cb;
  Image::Handle m_curModImg;
  int m_index;
};

OverlayInfoDialog::OverlayInfoDialog(QWidget* w,OverlayList::Handle l,ImageGroup::Handle i):
  OverlayInfoDialogBase(w,0,FALSE), m_overlayList(l),m_imageGroup(i),m_blockOverlayListUpdate(false)
{  
  m_overlayList->attach(this); 

  synchronizeDialog();
  #if (QT_VERSION >= 300)  
  QRegExp regexp("[a-zA-Z0-9_+.-]*");  
  QRegExpValidator *m_nameValidator = new QRegExpValidator(regexp,this);
  m_overlayEdit->setValidator(m_nameValidator);
  #endif
}

OverlayInfoDialog::~OverlayInfoDialog()
{  
  m_overlayList->detach(this);
}

void OverlayInfoDialog::update(const OverlayList* l, OverlayListMsg msg)
{
  if(!m_blockOverlayListUpdate)
    synchronizeDialog();
}

void OverlayInfoDialog::help()
{
  AssistantClient::getInstance()->showPage("imageinfo.html");
}

void OverlayInfoDialog::intentChanged(int i)
{
  MetaImage::Handle mi(m_overlayList->getActiveMetaImage());

  if(mi) {
    mi->getInfo()->setPurpose(ImageIntent::Code(i));
  }
}

void OverlayInfoDialog::synchronizeDialog()
{   
  MetaImage::Handle metaImage = m_overlayList->getActiveMetaImage();

  if(!metaImage) {
    clearDialog();
  } else {
    ImageDisplaySetting::Handle disp = metaImage->getDs();
    ImageInfo::Handle           info = metaImage->getImage()->getInfo();

    QString voxStr(QString("%1 x %2 x %3").arg(info->inqX())
		   .arg(info->inqY())
		   .arg(info->inqZ()));  
    QString dimStr(QString("%1 x %2 x %3 mm").arg(info->inqXDim())
		   .arg(info->inqYDim())
		   .arg(info->inqZDim()));
    QString volStr(QString("%1").arg(info->inqNumVolumes()));

    
    QString bitsStr(QString("%1 (%2 bpp)").arg(info->inqDtAsString())
		   .arg(info->inqBitsPerVoxel()));
   
    m_overlayEdit->blockSignals(true);
    m_overlayEdit->setText(info->inqImageName().c_str());   
    m_overlayEdit->blockSignals(false);
    m_fileNameLabel->setText(info->inqFileName().c_str());
    m_voxLabel->setText(voxStr);
    m_dimLabel->setText(dimStr);
    m_volLabel->setText(volStr);
    m_bppLabel->setText(bitsStr);
    
    m_dtiCombo->blockSignals(true);
    m_lutCombo->blockSignals(true);
    m_negLutCombo->blockSignals(true);
    m_negativeLuts->blockSignals(true);
    m_modCombo->blockSignals(true);

    std::for_each(m_imageGroup->beginLutList(),
                  m_imageGroup->endLutList(),
                  AddLutToCombo(m_lutCombo,disp->inqLookUpTable()));
    m_lutCombo->setCurrentText(disp->inqLookUpTable()->inqLutName());

    std::for_each(m_imageGroup->beginLutList(),
		  m_imageGroup->endLutList(),
		  AddLutToCombo(m_negLutCombo,disp->inqSecondaryLookUpTable()));
    if(disp->inqSecondaryLookUpTable())
      m_negLutCombo->setCurrentText(disp->inqSecondaryLookUpTable()->inqLutName());

    std::for_each(m_imageGroup->begin(),
                  m_imageGroup->end(),
                  AddModImageToCombo(m_modCombo,disp->inqModImage()));

    if(metaImage) 
      {            
        int mode = disp->inqDtiDisplay();
        switch (mode)
          {
          case DtiDisplay(None):
            m_dtiCombo->setCurrentItem(0);
            m_lutCombo->setEnabled(true);
	    m_negativeLuts->setChecked(disp->inqUseSecondaryLookUpTable());
	    m_negLutCombo->setEnabled(disp->inqUseSecondaryLookUpTable());
            m_modCombo->setEnabled(false);
            break;
          case DtiDisplay(Lines):
	    m_dtiCombo->setCurrentItem(1);           
            m_lutCombo->setEnabled(false);
	    m_negLutCombo->setEnabled(false);
	    m_negativeLuts->setChecked(false);
            m_modCombo->setEnabled(false);
            break;
          case DtiDisplay(RGB):
            m_dtiCombo->setCurrentItem(2);            
            m_lutCombo->setEnabled(false);
	    m_negLutCombo->setEnabled(false);
	    m_negativeLuts->setChecked(false);
            m_modCombo->setEnabled(true);
	    break;
	  case DtiDisplay(LinesRGB):
	    m_dtiCombo->setCurrentItem(3);           
            m_lutCombo->setEnabled(false);
	    m_negLutCombo->setEnabled(false);
	    m_negativeLuts->setChecked(false);
            m_modCombo->setEnabled(false);
            break;
          }
	m_intentCombo->setCurrentItem(info->inqPurpose());
      }
    
    m_dtiCombo->blockSignals(false);   
    m_lutCombo->blockSignals(false);   
    m_negLutCombo->blockSignals(false);   
    m_negativeLuts->blockSignals(false);
    m_modCombo->blockSignals(false);   
  }
}

void OverlayInfoDialog::clearDialog()
{
  m_overlayEdit->setText("");
  m_lutCombo->clear();
  m_voxLabel->clear();
  m_dimLabel->clear();
  m_volLabel->clear();
  m_bppLabel->clear();
  m_modCombo->clear();
}

void OverlayInfoDialog::modComboChanged(int n)
{  
  Image::Handle img;
  
  if(n){img = m_imageGroup->getImage(n - 1);}
  m_overlayList->setModImage(img);
}

void OverlayInfoDialog::lutComboChanged(int n)
{
  LookUpTable::Handle lut = m_imageGroup->getLut(n);
  m_overlayList->setLookUpTable(lut);
}

void OverlayInfoDialog::slutBoxChecked(bool state)
{
  LookUpTable::Handle lut;

  m_overlayList->getActiveMetaImage()->getDs()->setUseSecondaryLookUpTable(state);
  m_negLutCombo->setEnabled(state);
  m_overlayList->notify(OverlayListMsg(LookUpTable));
}

void OverlayInfoDialog::slutComboChanged(int n)
{
  LookUpTable::Handle lut = m_imageGroup->getLut(n);
  m_overlayList->setSecondaryLookUpTable(lut);
}

void OverlayInfoDialog::dtiComboChanged(int n)
{
  MetaImage::Handle metaImage = m_overlayList->getActiveMetaImage();

  if(metaImage.get()) 
    {
      if (metaImage->getInfo()->isDtiCompatible())
        {
	  metaImage->getInfo()->setDtiImage(true);
          switch(n)
            {
            case 0:
              metaImage->getDs()->setDtiDisplay(DtiDisplay(None));
              m_lutCombo->setEnabled(true);
              m_modCombo->setEnabled(false);
              m_overlayList->notify(OverlayListMsg(DtiMode));
              break;
            case 1:
              metaImage->getDs()->setDtiDisplay(DtiDisplay(Lines)); 
              m_lutCombo->setEnabled(false);              
              m_modCombo->setEnabled(true);
              m_overlayList->notify(OverlayListMsg(DtiMode));
              break;
            case 2:   
              metaImage->getDs()->setDtiDisplay(DtiDisplay(RGB));
              m_lutCombo->setEnabled(false);              
              m_modCombo->setEnabled(true);
	      metaImage->getDs()->inqBriCon()->setRange(0.0, 1.0);
              m_overlayList->notify(OverlayListMsg(LookUpTable));
              break;
	    case 3:
              metaImage->getDs()->setDtiDisplay(DtiDisplay(LinesRGB)); 
              m_lutCombo->setEnabled(false);              
              m_modCombo->setEnabled(true);
              m_overlayList->notify(OverlayListMsg(DtiMode));
	      break;
            }
        }
      else
        {
          emit message("Warning: Image is not a valid DTI image.",2000);
          m_dtiCombo->blockSignals(true);
          m_dtiCombo->setCurrentItem(0);   
          m_dtiCombo->blockSignals(false);
        }
    }
}

void OverlayInfoDialog::lutButtonPressed()
{
  QString fn = QFileDialog::getOpenFileName( QString::null, "LUTs (*.lut *.rgb *.lml)", this );
  if ( !fn.isEmpty() ) {
    QApplication::setOverrideCursor(Qt::waitCursor);
    emit message( QString("Loading lookup table.... %1").arg(fn), 2000 );
    LookUpTable::Handle lookUpTable = LookUpTable::load((const char *)fn);
    m_imageGroup->addLookUpTable( lookUpTable );

    QApplication::restoreOverrideCursor();

  }  else {
    emit message( "Loading aborted", 2000 );
  }
}

void OverlayInfoDialog::overlayTextChanged( const QString & newName ) 
{ 
  MetaImage::Handle metaImage = m_overlayList->getActiveMetaImage();

  if(metaImage && !newName.isEmpty()) 
    {
//       if (metaImage->inqReadOnly()) {
// 	emit message("Warning: Image must be unlocked before it can be renamed.",2000);    
// 	m_overlayEdit->blockSignals(true);   
// 	m_overlayEdit->setText(metaImage->getImage()->getInfo()->inqImageName().c_str());   
// 	m_overlayEdit->blockSignals(false);	
//       } else { 
// 	metaImage->setImageName(newName.latin1());
//       }
      metaImage->setImageName(newName.latin1());
      m_blockOverlayListUpdate = true;
      m_imageGroup->notify(ImageGroup::NameChange);
      m_blockOverlayListUpdate = false;
    }

}

