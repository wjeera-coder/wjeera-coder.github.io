//mytest.cc
#include "mylib.h"
int main(){
  XPrint XP;
  XP.print("Printing using XPrint");
  XP.print("Testing XRand ... ");
  XRand XR;
  XR.test(10);
  print("Testing print() ...");
  print("Testing myrand() ...");
  test_myrand(10);
}
