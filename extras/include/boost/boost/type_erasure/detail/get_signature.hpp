// Boost.TypeErasure library
//
// Copyright 2012 Steven Watanabe
//
// Distributed under the Boost Software License Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// $Id: get_signature.hpp,v 1.1.1.1 2015/02/27 15:20:28 mwebster Exp $

#ifndef BOOST_TYPE_ERASURE_DETAIL_GET_SIGNATURE_HPP_INCLUDED
#define BOOST_TYPE_ERASURE_DETAIL_GET_SIGNATURE_HPP_INCLUDED

#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/remove_pointer.hpp>

namespace boost {
namespace type_erasure {
namespace detail {

template<class Concept>
struct get_signature {
    BOOST_TYPEOF_NESTED_TYPEDEF_TPL(nested, &Concept::apply)
    typedef typename boost::remove_pointer<
        typename nested::type
    >::type type;
};

}
}
}

#endif
