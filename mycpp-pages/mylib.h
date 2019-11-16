//mylib.h
#include<iostream>
#include<string>
#include<random>
using namespace std;
mt19937 mt_rand(time(0));

class XPrint{
  public:
    void print(string text){
      cout<<text<<endl;
    }
};

class XRand{
  public:
    double rand(){
      return (1.0*mt_rand()/mt_rand.max());
    }
    void test(int n){
      for(int i=1;i<=n;i++){
        cout<<i<<" : "<<rand()<<endl;
      }
    } 
};

void print(string text){
  cout<<text<<endl;
}

double myrand(){
  return (1.0*mt_rand()/mt_rand.max());
}

void test_myrand(int n){
  for(int i=1;i<=n;i++){
    printf("%2d : %f\n",i,myrand());
   }
} 

