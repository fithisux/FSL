//
// C++ Interface: command
//
// Description: Command Design Strategy for implementing Menu commands;
// V Rama Aravind, 9/11/04
//
// Author: Rama Aravind Vorray <rama@fmrib.ox.ac.uk>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __COMMAND_H
#define __COMMAND_H

//! @brief Abstract Command class
//! @see Design Patterns
//! @author V Rama Aravind
class Command
{
public:
  virtual ~Command(){}

  //! @brief Implement this method to provide behaviour
  virtual void execute(void) = 0;

protected:
  Command(){}
};

#endif
