
// Copyright Aleksey Gurtovoy 2000-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: largest_int.cpp,v 1.1.1.2 2015/02/27 16:50:37 mwebster Exp $
// $Date: 2015/02/27 16:50:37 $
// $Revision: 1.1.1.2 $

#include <boost/mpl/aux_/largest_int.hpp>
#include <boost/mpl/aux_/test.hpp>


MPL_TEST_CASE()
{
    MPL_ASSERT(( is_same< mpl::aux::largest_int<bool,bool>::type, bool > ));
    MPL_ASSERT(( is_same< mpl::aux::largest_int<bool,char>::type, char > ));
    MPL_ASSERT(( is_same< mpl::aux::largest_int<char,bool>::type, char > ));
    MPL_ASSERT(( is_same< mpl::aux::largest_int<int,unsigned>::type, unsigned > ));
    MPL_ASSERT(( is_same< mpl::aux::largest_int<unsigned,long>::type, long > ));
}
