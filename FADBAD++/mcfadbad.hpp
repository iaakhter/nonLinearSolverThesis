// Copyright (C) 2009-2016 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.


#ifndef MC__MCFADBAD_HPP
#define MC__MCFADBAD_HPP
#include "interval.hpp"

#include "fadbad.h"


namespace fadbad
{

template <> struct Op<mc::Interval>
{
   typedef double Base;
   typedef mc::Interval I;
   static Base myInteger( const int i ) { return Base(i); }
   static Base myZero() { return myInteger(0); }
   static Base myOne() { return myInteger(1);}
   static Base myTwo() { return myInteger(2); }
   static double myPI() { return mc::PI; }
   static I myPos( const I& x ) { return  x; }
   static I myNeg( const I& x ) { return -x; }
   template <typename U> static I& myCadd( I& x, const U& y ) { return x+=y; }
   template <typename U> static I& myCsub( I& x, const U& y ) { return x-=y; }
   template <typename U> static I& myCmul( I& x, const U& y ) { return x*=y; }
   template <typename U> static I& myCdiv( I& x, const U& y ) { return x/=y; }
   static I myInv( const I& x ) { return mc::inv( x ); }
   static I mySqr( const I& x ) { return mc::pow( x, 2 ); }
   template <typename X, typename Y> static I myPow( const X& x, const Y& y ) { return mc::pow( x, y ); }
   static I myCheb( const I& x, const unsigned n ) { return mc::cheb( x, n ); }
   static I mySqrt( const I& x ) { return mc::sqrt( x ); }
   static I myLog( const I& x ) { return mc::log( x ); }
   static I myExp( const I& x ) { return mc::exp( x ); }
   static I mySin( const I& x ) { return mc::sin( x ); }
   static I myCos( const I& x ) { return mc::cos( x ); }
   static I myTan( const I& x ) { return mc::tan( x ); }
   static I myAsin( const I& x ) { return mc::asin( x ); }
   static I myAcos( const I& x ) { return mc::acos( x ); }
   static I myAtan( const I& x ) { return mc::atan( x ); }
   static bool myEq( const I& x, const I& y ) { return x==y; }
   static bool myNe( const I& x, const I& y ) { return x!=y; }
   static bool myLt( const I& x, const I& y ) { return x<y; }
   static bool myLe( const I& x, const I& y ) { return x<=y; }
   static bool myGt( const I& x, const I& y ) { return x>y; }
   static bool myGe( const I& x, const I& y ) { return x>=y; }
 };
 }
 #endif 