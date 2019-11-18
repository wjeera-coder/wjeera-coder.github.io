//mylib.h
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

//--------------------------- myrand()
/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,1)-real-interval */
double genrand_float32_notone(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

//---------------------------------------------------------
void init_myrand(){
  unsigned long s = time(0);
  printf("seed = %ld\n",s);
  init_genrand(s);
}
double myrand(){
  return genrand_float32_notone();
}
void test_myrand(int n){
  for(int i=1;i<=n;i++){
    printf("%2d : %12.10f\n",i,myrand());
  }
}

//--------------------------- for printing
void print(char* text){
  printf("%s\n",text);
}

void print_vec(double* x,int nd){
  for(int i=0;i<nd;i++){
      printf("%12.5f ",x[i]);
  }
  printf("\n");
}

void print_mat(double** A,int nr,int nc){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      printf("%12.5f ",A[i][j]);
    }
    printf("\n");
  }
}

void print_mat_vec(double** A,int nr,int nc,double* VF){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      printf("%12.5f ",A[i][j]);
    }
    printf("%20.5f\n",VF[i]);
  }
}

void rand_mat(double** A,int nr,int nc,double LB,double UB){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      A[i][j] = LB+(UB-LB)*myrand();
    }
  }
}

//------------------------------------- (struct) dvec
typedef struct {
  int size;
  double* valueAt; 
} dvec;

void init_dvec2(dvec* x,int size){
  x->size = size;
  x->valueAt = malloc(size*sizeof(double));
}

void init_dvec(dvec* x,int size,double a){
  x->size = size;
  x->valueAt = malloc(size*sizeof(double));
  for(int i=0;i<size;i++){ x->valueAt[i] = a; }
}

void show_dvec(const dvec* x){
  printf("(struct) dvec : size = %d\n",x->size);
  for(int i=0;i<x->size;i++){
    printf("%12.5f ",x->valueAt[i]);
  }
  printf("\n");
}

void show_dvec2(const dvec* x){
  for(int i=0;i<x->size;i++){
    printf("%12.5f ",x->valueAt[i]);
  }
  printf("\n");
}

void init_dvec_rand(dvec* x,int size,double LB,double UB){
  x->size = size;
  x->valueAt = malloc(size*sizeof(double));
  for(int i=0;i<size;i++){ x->valueAt[i] = LB+(UB-LB)*myrand(); }
}

//------------------------------------- (struct) dmat
typedef struct {
  int nr, nc;
  dvec* dvecAt;
} dmat;

void init_dmat(dmat* A,int nr,int nc,double a){
  A->nr = nr; A->nc = nc;
  A->dvecAt = malloc(nr*sizeof(dvec));
  for(int i=0;i<nr;i++){ init_dvec(&(A->dvecAt[i]),nc,a); }
}

void init_dmat_rand(dmat* A,int nr,int nc,double LB,double UB){
  A->nr = nr; A->nc = nc;
  A->dvecAt = malloc(nr*sizeof(dvec));
  for(int i=0;i<nr;i++){ init_dvec_rand(&(A->dvecAt[i]),nc,LB,UB); }
}

void show_dmat(const dmat* A){
  printf("(struct) dmat : nr = %d, nc = %d\n",A->nr,A->nc);
  for(int i=0;i<A->nr;i++){
    show_dvec2(&(A->dvecAt[i]));
  }
}

//------------------------------------- DE operation
void copy_vec(double* x1,int nd,double* x2){
  for(int i=0;i<nd;i++){
    x1[i] = x2[i];
  }
}
void mutate(double* xm,int nd,double* x1, double* x2, double * x3,double SF){
  for(int i=0;i<nd;i++){
    xm[i] = x1[i]+SF*(x2[i]-x3[i]);
  }
}
void adjust(double* x,int nd,double LB,double UB){
  for(int i=0;i<nd;i++){
    if(x[i]<LB || x[i]>UB){
      x[i] = LB+(UB-LB)*myrand();
    }
  }
}
void crossover(double* xc,int nd,double* xm, double CR){
  int IC = (int) (nd*myrand());
  for(int i=0;i<nd;i++){
    if (i==IC || myrand()<CR){ xc[i] = xm[i]; }
  }
}

//------------------------------------- Objective function
double f2(double* x,int nd){
  double fx = 0.0;
  double temp;
  for(int i=0;i<nd;i++){
    temp = x[i]-(i+1);
    fx += temp*temp;
  }
  return fx;
}
double f3(double* x,int nd){
  double fx = 0.0;
  double temp1, temp2;
  for(int i=0;i<nd-1;i++){
    temp1 = x[i+1]-x[i]*x[i];
    temp2 = x[i] -1.0;
    fx += (100.0)*temp1*temp1+temp2*temp2;
  }
  return fx;
}
double f_dvec(dvec* x){
  double fx = 0.0;
  double temp1, temp2;
  for(int i=0;i<x->size-1;i++){
    temp1 = x->valueAt[i+1]-(x->valueAt[i])*(x->valueAt[i]);
    temp2 = x->valueAt[i] -1.0;
    fx += (100.0)*temp1*temp1+temp2*temp2;
  }
  return fx;
}

void copy_dvec(dvec* x1,dvec* x2){
  for(int i=0;i<x1->size;i++){
    x1->valueAt[i] = x2->valueAt[i];
  }
}

//---------------------------------------------------------

typedef struct {
  int NP, ND;
  double LB, UB, SF, CR, VTR;
  dmat P;
  dvec xm, xc, xb, VF;
  double fx, fc, fb;
  //int tindex, r1, r2, r3;
  int nf, check, check2;
} myde;

void show_pop(const myde* DE){
  show_dmat(&(DE->P));
}

void init_myde(myde* DE,int np, int nd,double LB,double UB){ 
  DE->NP = np; DE->ND = nd;
  DE->LB = LB; DE->UB = UB;
  init_dmat_rand(&(DE->P),np,nd,LB,UB);
  //show_pop(DE);
  init_dvec2(&(DE->xm),nd); init_dvec2(&(DE->xc),nd);
  init_dvec2(&(DE->xb),nd); init_dvec2(&(DE->VF),np);
}

void set_myde_SFCR(myde* DE,double SF,double CR){
  DE->SF = SF; DE->CR = CR;
}

void set_myde_VTR(myde* DE,double VTR){ DE->VTR = VTR; }

void feval_myde(myde* DE){
  int bindex; double fb, fx;
  for(int i=0;i<DE->NP;i++){
    fx = f_dvec(&(DE->P.dvecAt[i])); DE->VF.valueAt[i] = fx;
    if (i==0){ bindex = i; fb = fx; }
    else {
      if (fx<fb){ bindex = i; fb = fx; }
    }
  }
  DE->fb = fb;
  copy_dvec(&(DE->xb),&(DE->P.dvecAt[bindex]));
  //show_dvec(&(DE->VF));
  print("After initialization ...");
  printf("fb = %f\n",DE->fb);
  show_dvec(&(DE->xb));
  printf("\n\n");
}



//------------------------------------- DE operation

void gen_indices3(int* r1,int* r2,int* r3,int i,int n){
  while(1){
    *r1 = (int) (n*myrand());
    *r2 = (int) (n*myrand());
    *r3 = (int) (n*myrand());
    if(*r1!=i && *r2!=i && *r3!=i && *r1!=*r2 && *r1!=*r3 && *r2!=*r3){ break; }
  }
}

void mutate_myde_xm(myde* DE,int tindex){
  int r1, r2, r3;
  gen_indices3(&r1,&r2,&r3,tindex,DE->NP);
  //double SF = 0.5+(0.2)*myrand();
  //printf("%d %d %d %d\n",tindex,r1,r2,r3);
  for(int i=0;i<DE->ND;i++){
    DE->xm.valueAt[i] = DE->P.dvecAt[r1].valueAt[i]+(DE->SF)*(DE->P.dvecAt[r2].valueAt[i]-DE->P.dvecAt[r3].valueAt[i]);
    //DE->xm.valueAt[i] = DE->P.dvecAt[r1].valueAt[i]+(SF)*(DE->P.dvecAt[r2].valueAt[i]-DE->P.dvecAt[r3].valueAt[i]);
  }
}

void adjust_myde_xm(myde* DE){
  double LB = DE->LB;
  double UB = DE->UB;
  for(int i=0;i<DE->ND;i++){
    if(DE->xm.valueAt[i]<LB || DE->xm.valueAt[i]>UB){
      DE->xm.valueAt[i] = LB+(UB-LB)*myrand();
    }
  }
}

void crossover_myde(myde* DE,int tindex){
  copy_dvec(&(DE->xc),&(DE->P.dvecAt[tindex]));
  int IC = (int) (DE->ND)*myrand();
  for(int i=0;i<DE->ND;i++){
    if (i==IC || myrand()<DE->CR){ DE->xc.valueAt[i] = DE->xm.valueAt[i]; }
  }
  //printf("in crossover : %d\n",tindex);
  //show_dvec(&(DE->xc));
}

void selection_myde(myde* DE,int tindex){
  //print("In selection ...");
  double fc = f_dvec(&(DE->xc)); DE->nf = (DE->nf)+1;
  if (fc<DE->VF.valueAt[tindex]){
    DE->VF.valueAt[tindex] = fc; copy_dvec(&(DE->P.dvecAt[tindex]),&(DE->xc));
    //compare with fb
    if (fc<DE->fb){
      DE->fb = fc; copy_dvec(&(DE->xb),&(DE->xc)); DE->check2 = 1;
      if (fc<DE->VTR){ DE->check = 1; }
    }
  }
}

void run_myde(myde* DE,int maxgen){
  DE->nf = 0; DE->check = 0;
  for(int ig=1;ig<=maxgen;ig++){
    for(int i=0;i<DE->NP;i++){
      mutate_myde_xm(DE,i);
      adjust_myde_xm(DE);
      crossover_myde(DE,i);
      DE->check2 = 0;
      selection_myde(DE,i);
      //if (DE->check2){
      if ((DE->nf)%10000 == 0){
        printf("ig = %d, nf = %d, fb = %e\n",ig,DE->nf,DE->fb);
        //show_dvec(&(DE->xb));
      }
      if (DE->check){ break; }
    }
    if (DE->check){
      print("Reach VTR ..."); break;
    }
  }
  //Done ...
  printf("nf = %d, fb = %e\n",DE->nf,DE->fb);
  show_dvec(&(DE->xb));
}

