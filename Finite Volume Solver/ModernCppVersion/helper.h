#ifndef __HELPER_H__
#define __HELPER_H__

#include<iostream>
#include<string>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<iterator>
#include<fstream>
#include<algorithm>
#include<math.h>

typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;

typedef std::vector<int> intVector;
typedef std::vector<intVector> intMatrix;

namespace State {
  enum Value : unsigned int {
    FLUID      = 1 << 0 ,
    NO_SLIP    = 1 << 1,
    FREE_SLIP  = 1 << 2,
    OUTFLOW    = 1 << 3,
    INFLOW     = 1 << 4,
    B_N        = 1 << 5, 
    B_S        = 1 << 7,
    B_W        = 1 << 6,
    B_O        = 1 << 8,
    B_SO       = 1 << 8 && 1 << 6,
    B_NO       = 1 << 8 && 1 << 5,
    B_NW       = 1 << 7 && 1 << 5,
    B_SW       = 1 << 7 && 1 << 6
  };
}


inline State::Value& operator|=(State::Value& a, State::Value b) {
    return a = static_cast<State::Value> (a | b);
}

inline constexpr State::Value operator|(State::Value a, State::Value b) {
    return a = static_cast<State::Value> (a | b);
}

inline constexpr State::Value operator&(State::Value a, State::Value b) {
    return a = static_cast<State::Value> (a & b);
}

template<typename iter>
double norm(iter first, iter last){
  double tmp = 0;
  typename std::iterator_traits<iter>::difference_type n = std::distance(first, last);
  while(n>0){
    tmp += (*first)*(*first);
    ++(*first);
    --n;
  }
  return std::sqrt(tmp);
}
//#############################################################################
//########### Useful print function for vectors and matrices ##################
template<typename T>
void printVector(const T& t){
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, " "));
}

template<typename T>
void printMatrix(const T& t){
    std::for_each(t.cbegin(), t.cend(), printVector<typename T::value_type>);
    std::cout<<std::endl;
}


//#############################################################################
//####################### for INIT_CPP_ File ##################################

void init_matrix(Matrix& M, double a);
intMatrix read_pgm(const char* filename);

//##############################################################################
//##############################################################################

// #define B_O     256  /*100000000*/
// #define B_W     128  /*010000000*/
// #define B_S      64  /*001000000*/
// #define B_N      32  /*000100000*/
// #define INFLOW   16  /*000010000*/
// #define OUTFLOW   8  /*000001000*/
// #define FREE_SLIP 4  /*000000100*/
// #define NO_SLIP   2  /*000000010*/
// #define FLUID     1  /*000000001*/
//
// #define B_SO    320  /*101000000*/
// #define B_NO    288  /*100100000*/
// #define B_SW    192  /*011000000*/
// #define B_NW    160  /*010100000*/

#endif
