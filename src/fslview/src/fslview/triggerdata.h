/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2005 University of Oxford  */

/*  CCOPYRIGHT */

#if !defined(TRIGGERDATA_H)
#define TRIGGERDATA_H


// @brief Stores a single sequence of FEAT trigger events
//
// FEAT trigger events are a list of trigger epochs and a single common duration value.
// Use this class to store trigger sequences for peri-stimulus plots etc.
class TriggerData
{
public:
  typedef boost::shared_ptr<TriggerData> Handle;

  static Handle create();

private:
  TriggerData();

  struct Implementation;
  auto_ptr<Implementation> m_impl;
};
#endif
