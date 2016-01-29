/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32)
#pragma warning (disable:4786)
#endif

#include "cursorwidget.h"
#include "imagedisplaysetting.h"
#include "overlaylist.h"

#include <qspinbox.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qlabel.h>
#include <qvalidator.h>
#include <qtooltip.h>

#include <iostream>

using namespace std;

//#define DEBUGGING
#include "tracker.h"

class VoxBox: public QSpinBox
{
 public:
  VoxBox(QWidget *parent): QSpinBox(parent)
  {
    setFont(QFont("Arial", 10));
    setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)1, (QSizePolicy::SizeType)3,
				0, 0, sizePolicy().hasHeightForWidth() ) );   
  }
};

class ValBox: public QLineEdit
{
 public:
  ValBox(QWidget *parent): QLineEdit(parent)
  {
    setFont(QFont("Arial", 10));
    setAlignment(AlignLeft);
    setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)1, (QSizePolicy::SizeType)3,
				0, 0, sizePolicy().hasHeightForWidth() ) );   
    setMinimumWidth( QFontMetrics(font()).width(QString("-0.00000e-00")) );
    setMaximumWidth( QFontMetrics(font()).width(QString("-0.00000e-00")) );
    setReadOnly(true);
  }
};

class MMBox: public QLineEdit
{
 public:
  MMBox(QWidget *parent): QLineEdit(parent)
  {
    setFont(QFont("Arial", 10));
    setAlignment(AlignLeft);
    setSizePolicy( QSizePolicy( (QSizePolicy::SizeType)1, (QSizePolicy::SizeType)3,
				0, 0, sizePolicy().hasHeightForWidth() ) );   
    setMinimumWidth( QFontMetrics(font()).width(QString("-0000.00")) );
    setMaximumWidth( QFontMetrics(font()).width(QString("-0000.00")) );  
  }
};

CursorWidget::CursorWidget(QWidget *parent, const Cursor::Handle& c, OverlayList::Handle ol):
  CursorWidgetBase(parent), m_cursor(c), m_overlayList(ol)
{
  setFont(QFont("Arial", 8));
  m_cursor->attach(this);
  m_overlayList->attach(this);
                
  connect(m_xVoxBox, SIGNAL(valueChanged(int)), SLOT(voxBoxChanged(int)));
  connect(m_yVoxBox, SIGNAL(valueChanged(int)), SLOT(voxBoxChanged(int)));
  connect(m_zVoxBox, SIGNAL(valueChanged(int)), SLOT(voxBoxChanged(int)));

  connect(m_xMmBox,SIGNAL(lostFocus()), SLOT(mmBoxChanged()));          
  connect(m_yMmBox,SIGNAL(lostFocus()), SLOT(mmBoxChanged()));
  connect(m_zMmBox,SIGNAL(lostFocus()), SLOT(mmBoxChanged()));
  connect(m_xMmBox,SIGNAL(returnPressed()), SLOT(mmBoxChanged()));          
  connect(m_yMmBox,SIGNAL(returnPressed()), SLOT(mmBoxChanged()));
  connect(m_zMmBox,SIGNAL(returnPressed()), SLOT(mmBoxChanged()));
  
  connect(m_volumeBox, SIGNAL(valueChanged(int)), SLOT(internalVolumeValueChanged(int)));

  setInputValidators();

  m_valBoxState=true;

  show();

  update(m_cursor);
}

void CursorWidget::internalVolumeValueChanged(int v)
{
  m_overlayList->getActiveMetaImage()->getDs()->setCurrentVolume(v);
  //  m_cursor->setVolume(v);
  emit volumeValueChanged(v);
}

CursorWidget::~CursorWidget()
{
  m_cursor->detach(this);
  m_overlayList->detach(this);
}

void CursorWidget::setValBoxState(bool state)
{
  m_valBoxState=state;
  m_valBox->setEnabled(m_valBoxState);
}

void CursorWidget::update(const Cursor::Handle& c)
{
  TRACKER("CursorWidget::update(const Cursor::Handle& c)");

  MESSAGE( QString("c = %1, %2, %3, %4").arg(c->inqX()).arg(c->inqY()).arg(c->inqZ()).arg(c->inqV()) );

  blockBoxSignals(true);
 
  Image::Handle image = m_overlayList->getActiveMetaImage()->getImage();

  string coordsysstring("Error!");
  switch(image->getInfo()->inqCoordSystem()) 
    {
    case ImageCoordSystem::Unknown: coordsysstring = "Unknown"; break;
    case ImageCoordSystem::ScannerAnatomical: coordsysstring = "Scanner Anatomical"; break;
    case ImageCoordSystem::AlignedAnatomical: coordsysstring = "Aligned Anatomical"; break;
    case ImageCoordSystem::Talairach: coordsysstring = "Talairach"; break;
    case ImageCoordSystem::MNI_152: coordsysstring = "MNI_152"; break;
    default:
      break;
    }
  m_xformDescription->setText(tr("Coordinate space: %1").arg(coordsysstring));

  if(c) {
    short radiogX(c->inqX());
    if(image->getInfo()->inqNoDimensions())
      {
	m_xMmBox->setText("");
	m_yMmBox->setText("");
	m_zMmBox->setText("");
      }
    else
      {
	float x(0), y(0), z(0);
	
	if(!image->getInfo()->isStoredRadiological())
	  radiogX = image->getInfo()->inqX()-1-radiogX;
	
	image->getInfo()->voxToMMCoord(radiogX, c->inqY(), c->inqZ(),
				       x, y, z);
	m_xMmBox->setText(tr("%1").arg(x, 3, 'f', 2));  
	m_yMmBox->setText(tr("%1").arg(y, 3, 'f', 2));  
	m_zMmBox->setText(tr("%1").arg(z, 3, 'f', 2));
      }
    
    m_xVoxBox->setValue(radiogX);
    m_yVoxBox->setValue(c->inqY());
    m_zVoxBox->setValue(c->inqZ());
  }
  setVolumeValue(m_overlayList->getActiveMetaImage()->getDs()->inqCurrentVolume());

  updateValBox();

  blockBoxSignals(false);
}

void CursorWidget::update(const OverlayList* ol, OverlayListMsg msg)
{
  TRACKER("CursorWidget::update(const OverlayList* ol, OverlayListMsg msg)");
  blockSignals(true);
  blockBoxSignals(true);

  if(msg == OverlayListMsg(Select))  { update(m_cursor); }
  if(msg == OverlayListMsg(DtiMode)) updateValBox();
  if(msg == OverlayListMsg(DtiMode) || msg == OverlayListMsg(Select) || msg == OverlayListMsg(Visibility))
    {
      bool state(ol->getActiveMetaImage()->inqVisibility());
      setValBoxState(state);
      MetaImage::Handle mi = m_overlayList->getActiveMetaImage();
      m_volumeBox->setEnabled(mi->getImage()->getInfo()->inqNumVolumes() > 1);
      m_volumeBox->setRange(0, mi->getImage()->getInfo()->inqNumVolumes() - 1);
      m_volumeBox->setValue(mi->getDs()->inqCurrentVolume());

      updateValBox();
    }
  blockBoxSignals(false);
  blockSignals(false);
}

void CursorWidget::voxBoxChanged(int v)
{
  TRACKER("CursorWidget::voxBoxChanged");
  blockBoxSignals(true);
  
  Image::Handle image = m_overlayList->getActiveMetaImage()->getImage();

  if(image->getInfo()->inqNoDimensions())
    {
      m_xMmBox->setText("");
      m_yMmBox->setText("");
      m_zMmBox->setText("");
    }
  else
    {
      float x(0), y(0), z(0);
      image->getInfo()->voxToMMCoord(m_xVoxBox->value(), m_yVoxBox->value(), m_zVoxBox->value(), x, y, z);
      m_xMmBox->setText(tr("%1").arg(x, 3, 'f', 2));
      m_yMmBox->setText(tr("%1").arg(y, 3, 'f', 2));
      m_zMmBox->setText(tr("%1").arg(z, 3, 'f', 2));
    }

  m_cursor->detach(this);
  short radiogX(m_xVoxBox->value());
  if(!image->getInfo()->isStoredRadiological())
    radiogX = image->getInfo()->inqX()-1-radiogX;
  m_cursor->setCursor(radiogX,m_yVoxBox->value(),m_zVoxBox->value());
  m_cursor->attach(this);

  updateValBox();

  blockBoxSignals(false);
}

void CursorWidget::mmBoxChanged()
{
  TRACKER("CursorWidget::mmBoxChanged");
       
  Image::Handle image = m_overlayList->getActiveMetaImage()->getImage();
  
  if(image->getInfo()->inqNoDimensions())
    {
      m_xMmBox->setText("");
      m_yMmBox->setText("");
      m_zMmBox->setText("");
    }
  else
    { 
      blockBoxSignals(true);

      int pos;  
  
      QString xText = m_xMmBox->text();
      QString yText = m_yMmBox->text();
      QString zText = m_zMmBox->text();

      if(QValidator::Intermediate == m_xBoxValidator->validate(xText,pos))
        m_xMmBox->setText(fixMmBoxVal(m_xBoxValidator,xText));
      if(QValidator::Intermediate == m_yBoxValidator->validate(yText,pos))
        m_yMmBox->setText(fixMmBoxVal(m_yBoxValidator,yText));
      if(QValidator::Intermediate == m_zBoxValidator->validate(zText,pos))
        m_zMmBox->setText(fixMmBoxVal(m_zBoxValidator,zText));

      float xMm = m_xMmBox->text().toFloat();
      float yMm = m_yMmBox->text().toFloat();
      float zMm = m_zMmBox->text().toFloat(); 

      short xCur(0), yCur(0), zCur(0);
      image->getInfo()->mmToVoxCoord(xMm ,yMm, zMm , xCur, yCur, zCur);
  
      m_xVoxBox->setValue(xCur);
      m_yVoxBox->setValue(yCur);
      m_zVoxBox->setValue(zCur);

      short radiogX(xCur);
      if(!image->getInfo()->isStoredRadiological())
	radiogX = image->getInfo()->inqX()-1-radiogX;
      m_cursor->detach(this);
      m_cursor->setCursor(radiogX,yCur,zCur);
      m_cursor->attach(this);

      updateValBox();

      blockBoxSignals(false);
    }
}

void CursorWidget::blockBoxSignals(bool state)
{
  m_xVoxBox->blockSignals(state);
  m_yVoxBox->blockSignals(state);
  m_zVoxBox->blockSignals(state);  
  m_xMmBox->blockSignals(state);
  m_yMmBox->blockSignals(state);
  m_zMmBox->blockSignals(state);
  m_volumeBox->blockSignals(state);
}

void CursorWidget::updateValBox()
{
  Image::Handle image;
  ImageDisplaySetting::Handle ds;

  if(MetaImage::Handle mi = m_overlayList->getActiveMetaImage()) {
    image = mi->getImage();
    ds    = mi->getDs();
  } else {
    image = m_overlayList->getMainMetaImage()->getImage();    
    ds    = m_overlayList->getMainMetaImage()->getDs(); 
  }
 
  //  int vol = std::min(m_cursor->inqV(), short (image->getInfo()->inqNumVolumes()-1));
  int vol = ds->inqCurrentVolume();
  float i = image->getVolume(vol) -> value(m_cursor->inqX(),m_cursor->inqY(),m_cursor->inqZ());
  
  if(ds->inqDtiDisplay() == DtiDisplay(None) && m_valBoxState) {      
    m_valBox->setEnabled(true);
    QString valStr;
    valStr.setNum(i,'g');
    
    if (valStr.length() > 9){valStr.setNum(i,'g',2);}
    
    QToolTip::remove(m_valBox);
    QToolTip::add(m_valBox, tr("Voxel value: %1").arg(i,8));
    
    m_valBox->setText(valStr); 
    m_valBox->repaint();
  } else {
    m_valBox->setEnabled(false);
    m_valBox->clear();
  }
}

void CursorWidget::setInputValidators()
{
  TRACKER("CursorWidget::setInputValidators");
  Image::Handle image = m_overlayList->getActiveMetaImage()->getImage();
  
  m_xVoxBox->setRange(0,image->getInfo()->inqX() - 1);
  m_yVoxBox->setRange(0,image->getInfo()->inqY() - 1);
  m_zVoxBox->setRange(0,image->getInfo()->inqZ() - 1);

  m_volumeBox->setEnabled(image->getInfo()->inqNumVolumes() > 1);
  m_volumeBox->setRange(0, image->getInfo()->inqNumVolumes() - 1);

  if(!image->getInfo()->inqNoDimensions())
    {
 
      m_xBoxValidator = new QDoubleValidator(this);
      m_yBoxValidator = new QDoubleValidator(this);
      m_zBoxValidator = new QDoubleValidator(this);
    
      m_xMmBox->setValidator(m_xBoxValidator);
      m_yMmBox->setValidator(m_yBoxValidator);
      m_zMmBox->setValidator(m_zBoxValidator); 

      float xMax, yMax, zMax, xMin, yMin, zMin;
      image->getInfo()->voxToMMCoord((image->getInfo()->inqX() - 1),
				     (image->getInfo()->inqY() - 1),
				     (image->getInfo()->inqZ() - 1) ,
				     xMax, yMax, zMax);
      image->getInfo()->voxToMMCoord(0, 0, 0, xMin, yMin, zMin);

      //Check that xMin < xMax otherwise Qt Validators are troublesome, need inverting.

      if(xMin < xMax){m_xBoxValidator->setRange(xMin,xMax,2);}
      else           {m_xBoxValidator->setRange(xMax,xMin,2);}
  
      if(yMin < yMax){m_yBoxValidator->setRange(yMin,yMax,2);}
      else           {m_yBoxValidator->setRange(yMax,yMin,2);}
  
      if(zMin < zMax){m_zBoxValidator->setRange(zMin,zMax,2);}
      else           {m_zBoxValidator->setRange(zMax,zMin,2);}
    }
}

QString CursorWidget::fixMmBoxVal(QDoubleValidator* v,QString & s)
{
  TRACKER("CursorWidget::fixMmBoxVal()");
  float val = s.toFloat();
  
  if(val > v->top())         {val = v->top();}
  else if(val < v->bottom()) {val = v->bottom();}

  QString retStr(tr("%1").arg(val,3,'f',2));
  return retStr;
}

void CursorWidget::setVolumeValue(int v)
{
  TRACKER("CursorWidget::setVolumeValue(int v)");
  MESSAGE( QString("v = %1").arg(v) );

  m_volumeBox->blockSignals(true);
  m_volumeBox->setValue(v);
  m_volumeBox->blockSignals(false);
}

