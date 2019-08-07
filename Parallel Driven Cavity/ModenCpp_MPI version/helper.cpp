#include"helper.h"
#include<iostream>
#include<string>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<iterator>
#include<algorithm>
#include<fstream>


void init_matrix(Matrix& M, double a){
    int nrow = M.size();
    int ncol = M[0].size();
   for(int i = 0; i<nrow; ++i){
      for(int j = 0; j<ncol; ++j)
      M[i][j] = a;
  }
}


intMatrix read_pgm(const char* filename){
  std::FILE* input = nullptr;
  char line[1024];
  int levels;
  int xsize, ysize;

  if ((input = std::fopen(filename, "rb")) == 0) {
      char szBuff[80];
      std::sprintf(szBuff, "Can not read file %s !!!", filename);
  }

  /* check for the right "magic number" */
  if (std::fread(line, 1, 3, input) != 3) {
      std::fclose(input);
  }

  do {
    std::fgets(line, sizeof line, input);
  } while(*line == '#');

  std::sscanf(line, "%d %d\n", &xsize, &ysize);
  std::cout<<" Image size: " << xsize << "x" << ysize <<std::endl;

  std::fgets(line, sizeof line, input);
  std::sscanf(line, "%d\n", &levels);

  intMatrix pic(xsize, std::vector<int>(ysize));

  std::cout<< "Initializing Image ... " << std::endl;
  for (int j = ysize-1; j >= 0; --j) {
    for (int i = 0; i < xsize; ++i) {
        int value;
        std::fscanf(input, "%d", &value);

        if (value == EOF) {
            fclose(input);
            std::cerr<< "read of geometry file failed!" <<std::endl;
        }
        pic[i][j] = value;
    }
  }
  std::cout<< "Initializing of image completed" << std::endl;
  std::fclose(input);
  return pic;
}
