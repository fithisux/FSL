#if !defined(PROPERTIES_H)
#define PROPERTIES_H

#include <boost/shared_ptr.hpp>

class Properties {
public:
  typedef boost::shared_ptr< Properties > Handle;

  ~Properties();
  
  bool inqAskCreate4dMask() const;
  bool inqCreate4dMask() const;
  void setAskCreate4dMask(bool);
  void setCreate4dMask(bool);

  static Handle create();
private:
  Properties();

  struct Implementation;  
  const std::auto_ptr< Implementation > m_impl;
};

#endif
