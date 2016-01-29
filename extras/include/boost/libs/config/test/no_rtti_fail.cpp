//  This file was automatically generated on Mon Aug 25 10:51:51 2008
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-4.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.//
//  Revision $Id: no_rtti_fail.cpp,v 1.1.1.1 2015/02/27 16:50:26 mwebster Exp $
//


// Test file for macro BOOST_NO_RTTI
// This file should not compile, if it does then
// BOOST_NO_RTTI should not be defined.
// See file boost_no_rtti.ipp for details

// Must not have BOOST_ASSERT_CONFIG set; it defeats
// the objective of this file:
#ifdef BOOST_ASSERT_CONFIG
#  undef BOOST_ASSERT_CONFIG
#endif

#include <boost/config.hpp>
#include "test.hpp"

#ifdef BOOST_NO_RTTI
#include "boost_no_rtti.ipp"
#else
#error "this file should not compile"
#endif

int main( int, char *[] )
{
   return boost_no_rtti::test();
}

