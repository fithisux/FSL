
// Copyright Aleksey Gurtovoy 2000-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: logical.cpp,v 1.1.1.2 2015/02/27 16:50:37 mwebster Exp $
// $Date: 2015/02/27 16:50:37 $
// $Revision: 1.1.1.2 $

#include <boost/mpl/logical.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/aux_/test.hpp>

struct unknown;

using mpl::true_;
using mpl::false_;

MPL_TEST_CASE()
{
    MPL_ASSERT(( mpl::and_< true_,true_ > ));
    MPL_ASSERT_NOT(( mpl::and_< false_,true_ > ));
    MPL_ASSERT_NOT(( mpl::and_< true_,false_ > ));
    MPL_ASSERT_NOT(( mpl::and_< false_,false_ > ));
    MPL_ASSERT_NOT(( mpl::and_< false_,unknown > ));
    MPL_ASSERT_NOT(( mpl::and_< false_,unknown,unknown > ));

    MPL_ASSERT(( mpl::or_< true_,true_ > ));
    MPL_ASSERT(( mpl::or_< false_,true_ > ));
    MPL_ASSERT(( mpl::or_< true_,false_ > ));
    MPL_ASSERT_NOT(( mpl::or_< false_,false_ > ));
    MPL_ASSERT(( mpl::or_< true_,unknown > ));
    MPL_ASSERT(( mpl::or_< true_,unknown,unknown > ));

    MPL_ASSERT_NOT(( mpl::not_< true_ > ));
    MPL_ASSERT(( mpl::not_< false_ > ));
}
