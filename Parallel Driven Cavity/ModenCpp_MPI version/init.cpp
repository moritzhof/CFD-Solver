#include"helper.h"
#include"init.h"

#include<utility>
#include<iostream>
#include<iterator>
#include<vector>
#include<math.h>
#include<numeric>

void init_FLAG(const char* geometry, intMatrix& FLAG, intMatrix& PIC){

  PIC = read_pgm(geometry);
  int xsize = FLAG.size();
  int ysize = FLAG[0].size();

  for(std::pair<intMatrix::iterator,intMatrix::iterator >itr(PIC.begin(), FLAG.begin()); itr.first != PIC.end(); ++itr.first, ++itr.second){
    for(std::pair<std::vector<int>::iterator, std::vector<int>::iterator> itc(itr.first->begin(), itr.second->begin()); itc.first != itr.first->end(); ++itc.first, ++itc.second){
      switch (*itc.first) {
        case 0: *itc.second = 1<<1; break;
        case 1: *itc.second = 1<<2; break;
    		case 2: *itc.second = 1<<3; break;
    		case 3: *itc.second = 1<<4; break;
    		case 4: *itc.second = 1<<0; break;
      }
    }
  }

  for(int i = 0; i<xsize; ++i){
    for(int j = 0; j<ysize; ++j){
      if(FLAG[i][j] & State::NO_SLIP){
          if(i != 0       && (FLAG[i-1][j] & State::FLUID)) FLAG[i][j] |= State::B_W;
          if(i != xsize-1 && (FLAG[i+1][j] & State::FLUID)) FLAG[i][j] |= State::B_O;
          if(j != 0       && (FLAG[i][j-1] & State::FLUID)) FLAG[i][j] |= State::B_S;
          if(j != ysize-1 && (FLAG[i][j+1] & State::FLUID)) FLAG[i][j] |= State::B_N;
      }
    }
  }
}
