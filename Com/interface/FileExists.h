#ifndef FileExists_h
#define FileExists_h
#include <iostream>
#include <assert.h>
#include <fstream>
  /**
   * Will test whether a file exists
   * if fail is set to true, executable will exit
   */

static bool fexists(const string filename, bool fail)
{
  ifstream ifile(filename.c_str());
  if (fail && !ifile) {
    std::cout << "File does not exist: " << filename << endl;
    assert(ifile);
  }
  return ifile;
}

#endif
