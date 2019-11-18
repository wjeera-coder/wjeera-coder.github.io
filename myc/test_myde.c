//test_myde.c
#include "mylib.h"

int main(){
  print("Testing mylib.h & myde (struct)");
  init_myrand();
  myde DE; init_myde(&DE,100,30,-100.0,100.0);
  set_myde_SFCR(&DE,0.5,0.9);
  set_myde_VTR(&DE,1.0e-10);
  feval_myde(&DE);
  run_myde(&DE,100000);
  return 0;
}
