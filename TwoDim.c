/*   TwoDim.c 
 *   Written Spring 2012 -- Patrick Malsom
 *   HMC algorithm written in C and OpenMP
 *   Reference: C HMC code
 */

// =========================================
// Library Definitions
// =========================================
//STD Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//OpenMP libraries
#include <omp.h>
//GNU Scientific Libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

//defines the potentials in calcPotentials function
#include "potentials.h"
//constants header defines path length and temp etc.
#include "constants.h"

const int NUMu=NUMBEAD-1;
const int NUMl=NUMBEAD-2;

const double SIGMA=sqrt(2.0*TEMP);
const double SIGMA2=2.0*TEMP;

// =========================================
// Structure Definitions
// =========================================

//Define the "position" struct
//only stores the positions
typedef struct _position
{
  double pos[NUMDIM];
} position;

//Define the "averages" struct
//stores average of the mean and covariance matrix
typedef struct _averages
{
  double mean[NUMDIM];
  double xx;
  double yy;
  double xy;
} averages;

//Define the "config" struct
//stores positions and all potentials
typedef struct _config
{
  double pos[NUMDIM];
  double Energy;
  double G;
  double gradG[NUMDIM];
  double LinvG[NUMDIM];
} config;

// =========================================
// Function Prototypes
// =========================================
// =========================================

void renormBB(double bb[NUMBEAD], double du, double doubleNUMu);

void renorm(position *currentpos, double du, double doubleNUMu);

void calcPotentials(config *currentConfig, int beadIndex);
//Takes the config 'currentConfig' and calculates
// G, gradG, LinvG, Energy for currentConfig[beadIndex].
//Note: calculates potentials for a single bead (beadIndex).

void generateBB(double bb[NUMBEAD], double du, double dt, double GaussRandArray[NUMu], gsl_rng *RanNumPointer);
// generates a standard browinan bridge (starts at zero ends at zero)
// does not renormalize the bridge to correct quadratic variation (see renormBB)

void GenGaussRand(double GaussRand[NUMu], gsl_rng *RanNumPointer, double StdDev);
//generate NUMBEAD of Gaussian random nums; save to GaussRandArray[]

void rotateConfig(config **a, config **b, config **c);
//rotation of pointers.
// b       -> a,    c   -> b,        a   -> c
//(current -> old,  new -> current,  old -> new)
//can now write over new to calculate a new state (discarding old-old state)

void LInverse(config* currentConfig, double doubleNUMl, double du, double vecdg[NUMl], double veci1[NUMl], double veci0[NUMl]);
//L^(-1) y = x. finds the vector y
// function solves for all currentConfig.LinvG of the current config

void MolecularDynamics(config *oldConfig, config *currentConfig, config *newConfig, double du, double dt);
//performs molecular dynamics. uses oldConfig and currentConfig to calculate newConfig

void preconditionSPDE(config* currentConfig, config* newConfig, double du, double dt, double doubleNUMu, double bb[NUMBEAD], double GaussRandArray[NUMu], gsl_rng *RanNumPointer);

void saveConfig(config *currentConfig, config *saveConfig);
//copies ALL elements of currentto saveConfig
//leaves currentConfig unmodified.

void savePosition(position *currentpos, position *savepos);
//copies ALL elements of currentpos to savepos
//leaves currentpos unmodified.

void saveConfigtoPos(config *currentConfig, position *savepos);
//copies ONLY pos elements of currentConfig to savepos
//leaves currentConfig unmodified.

void savePostoConfig(position *currentpos, config *saveConfig);
//copies ONLY pos elements of currentpos to saveConfig
//leaves currentpos unmodified.

void writeConfig(config *newConfig,averages *tubeAve, int MCloopi);
//write the config to file with name position******.dat where ****** is MCloopi
//only writes the modulo-th iteration.

void printPositionPos(position* currentpos, int beadIndex);
//Print the positions of the position struct 

void printConfigPos(config* currentConfig, int beadIndex);
//Print the positions of the config struct 

void printConfigPot(config* currentConfig, int beadIndex);
//Print all information of the config struct calculated in calcPotentials
// (Postitions, grad G, Energy, G)

void printConfigAll(config* currentConfig, int beadIndex);
//Print all information of the config struct

void printDistance(config *newConfig, position *savePos);

void ProbAccRatio(config *currentConfig, config *newConfig, double dt, double du, double *ratio);
// calculate the probability accaptance ratio for the step
//used to determine the acc/rej of the Metropolis-Hastings MC test

void zeroAverages(averages *tubeAve, int *tau);
//zero all elements in the averages struct 

void accumulateAverages(averages *tubeAve, config *newConfig, int *tau);
//sum the averages. This includes the mean (sum of the positions) and the 
//square of the positions. This is enough information to calculate the covariance
//matrix in post processing

void normalizeAverages(averages *tubeAve, int *tau);
//normalize the average struct before writing to file
//simply dividing by the total number of accumulate average calls

void accumulateArrayPlot(int arrayPlot[300][200], config *currentConfig);
//accumulate the averages for the array plot
// if in the bound of the array
//this is not all that general (size of matrix is hardcoded)

void writeArrayPlot(int arrayPlot[300][200], int MCloopi);
//print the arrayPlot to file

//===============================================================
//               MAIN Program
//===============================================================
int main(int argc, char *argv[])
{ 
printf("numbead %d \n",NUMBEAD);

  setbuf(stdout,NULL);

  //===============================================================
  // Declare variables and print to std output for reference
  //===============================================================
  //define the Config sturcts. Example: configOld[n].pos[i][j]
  //where n->Bead i->Particle j->dimension
  //configOld and configCurrent are switched between when doing MD
  config *configOld = calloc(NUMBEAD,sizeof(config));
  config *configCurrent = calloc(NUMBEAD,sizeof(config));
  config *configNew = calloc(NUMBEAD,sizeof(config));

  //Used to save positions for MHMC rejection
  position *savePos = calloc(NUMBEAD,sizeof(position));

  //averages
  averages *tubeAve = calloc(NUMBEAD,sizeof(averages));

  double doubleNUMu=(double)NUMu;
  double doubleNUMl=(double)NUMl;

  //Parameters
  double du=DU;
  double dt=PREDT*DU*DU;
  double h=sqrt(2.0l*dt);

  //Incrimenter Declarations
  int i,j,acc,rej;
  int MDloopi,MCloopi;
  int tau=0;

  //Vectors for doing the L Inverse
  double *vecdg = calloc(NUMl,sizeof(double));;
  double *veci0 = calloc(NUMl,sizeof(double));;
  double *veci1 = calloc(NUMl,sizeof(double));;
  //double veci1[NUMBEAD-2];
  //double veci0[NUMBEAD-2];

  // array to store rand nums in
  double *GaussRandArray = calloc(NUMu,sizeof(double));
  //double GaussRandArray[NUMu];

  // ratio for incrimenting the MH-MC test
  double ratio;

  // storage for brownian bridge
  double *bb = calloc(NUMu,sizeof(double));
  //double bb[NUMBEAD];

  // array plot of the average path
  //int xBinMax=300;
  //int yBinMax=200;
  int arrayPlot[300][200];
  for(i=0;i<300;i++){
    for(j=0;j<200;j++){
      arrayPlot[i][j]=0;
  } }

  //Print parameters for the run in stdout
  printf("=======================================================\n");
  printf("HMC method for 2D potentials \n");
  printf("=======================================================\n");
  printf("TEMPERATURE = %f \n",TEMP);
  printf("=======================================================\n");
  printf("Number of Metropolis Hastings steps: %i\n",NUMMC);
  printf("Number of MD steps: %i \n",NUMMD);
  printf("=======================================================\n");
  printf("Number of Dimensions: %i \n",NUMDIM);
  printf("Number of Beads: %i \n",NUMBEAD);
  printf("Path grid: du = %+.8e \n", DU);
  printf("Sampling Parameters: dt=%f \n",dt);
  printf("=======================================================\n");
  printf("MD step: h=%+.8e \n",h);
  printf("MD time (n*h): %+.8e \n",NUMMD*h);
  printf("=======================================================\n");

  //===============================================================
  // Reading the input configuration file into savepos
  //===============================================================
  //Input file to be read as first command line argument
  if(argv[1]==NULL) { 
    printf("No input file. Exiting!\n");
    exit(-1);
  }
  else {
    printf("Input Configuratrion File: %s\n",argv[1]);
  }
  int lineNum = 0;
  FILE *fptr = fopen(argv[1],"r");
  switch(NUMDIM){
    case 3:  //For 3 Dimensions
      while( EOF != fscanf(fptr,"%lf %lf %lf",
      &(savePos[lineNum].pos[0]),
      &(savePos[lineNum].pos[1]),
      &(savePos[lineNum].pos[2])) ) {
        lineNum++;
      }
      break;
    case 2:   //For 2 Dimensions 
      while( EOF != fscanf(fptr,"%lf %lf",
      &(savePos[lineNum].pos[0]),
      &(savePos[lineNum].pos[1])) ) {
        lineNum++;
      }
      break;
    case 1:  //For 1 Dimension
      while( EOF != fscanf(fptr,"%lf",
      &(savePos[lineNum].pos[0])) ) {
        lineNum++;
      }
      break;
    default:
      printf("ERROR: NUMDIM incorrectly defined. Exiting!\n");
      exit(-1);
  }


  //===============================================================
  // GNU Scientific Library Random Number Setup
  //===============================================================
  // Example shell command$ GSL_RNG_SEED=123 ./a.out
  printf("=======================================================\n");
  const gsl_rng_type * RanNumType;
  gsl_rng *RanNumPointer; 
  gsl_rng_env_setup();
  RanNumType = gsl_rng_default;
  RanNumPointer= gsl_rng_alloc (RanNumType);
  printf("Random Number Generator Type: %s \n", gsl_rng_name(RanNumPointer));
  printf("RNG Seed: %li \n", gsl_rng_default_seed);
  printf("=======================================================\n");

  double randUniform;

  renorm(savePos, DU, doubleNUMu);



  //===============================================================
  //     Start of HMC Loop (loops over Metropolis Hastings - MC steps)
  //===============================================================

  printf("START Hybrid Monte Carlo MAIN LOOP\n");
  printf("=======================================================\n");
  acc=0;
  rej=0;
  zeroAverages(tubeAve,&tau);

  for(MCloopi=1; MCloopi<=NUMMC; MCloopi++)
  {
    //zero ratio for MH MC test
    ratio=0.0l;
    //===============================================================
    //     Perform one SPDE step
    //===============================================================
  
    //store savePos.pos values to configCurrent.pos
    // savePos.pos stores the positions in case of rejection of the MHMC
    savePostoConfig(savePos, configCurrent);
 
    //(calculates potentials in config given the positions)
    #pragma omp parallel for
    for(i=0;i<NUMBEAD;i++) {calcPotentials(configCurrent,i);}
 
    //calculate LinvG for the config
    LInverse(configCurrent, doubleNUMl, du, vecdg, veci1, veci0);

    //do the preconditioned form of the SPDE
    preconditionSPDE(configCurrent, configNew, du, dt, doubleNUMu, bb, GaussRandArray, RanNumPointer);

    //(calculates potentials in config given the positions)
    #pragma omp parallel for
    for(i=0;i<NUMBEAD;i++) {calcPotentials(configNew,i);}

    //calculate LinvG for the config
    LInverse(configNew, doubleNUMl, du, vecdg, veci1, veci0);

    //acc ratio of newconfig
    ProbAccRatio(configCurrent, configNew, dt, du, &ratio);

    //calculate the averages for the tubes estimator
    accumulateAverages(tubeAve,configNew,&tau);
    accumulateArrayPlot(arrayPlot, configNew);

    printf("SPDE ratio: %+0.10f \n",ratio);
    //===============================================================
    //     Start of MD Loop
    //     This loop needs to be focused on for parallelization
    //===============================================================

    for(MDloopi=1;MDloopi<=NUMMD; MDloopi++)
    {
      //rotate the configuration        
      rotateConfig(&configOld, &configCurrent, &configNew);

      //do the MD position update
      MolecularDynamics(configOld, configCurrent, configNew, du, dt);
 
      //(calculates potentials in config given the positions)
      #pragma omp parallel for
      for(i=0;i<NUMBEAD;i++) {calcPotentials(configNew,i);}

      //calculate LinvG for the config
      LInverse(configNew, doubleNUMl, du, vecdg, veci1, veci0);

      //calculate the average distance moved in the step and print to std out
      if(MDloopi%WRITESTDOUT==0){
        printf("MDi: %.5d | MDi*h: %0.5f | MD ratio: %+0.5f | distance: ",MDloopi,MDloopi*sqrt(2*dt),ratio);
        printDistance(configNew, savePos);
      }

      //acc ratio of newconfig
      ProbAccRatio(configCurrent, configNew, dt, du, &ratio);
      //printf("%i  ProbAcc= %+.15e  QV Vel= %0.15e \n", MDloopi, ratio, qvvel);

      //calculate the averages for the tubes estimator
      accumulateAverages(tubeAve,configNew,&tau);
      accumulateArrayPlot(arrayPlot, configNew);
    }
    //===============================================================
    //Metropolis Hastings Monte-Carlo test
    //===============================================================
    randUniform = gsl_rng_uniform(RanNumPointer);
    if( exp(ratio/SIGMA2) > randUniform ){
      acc++;
      saveConfigtoPos(configNew, savePos);
    }
    else{
      rej++;
    }
    printf("rand=%+0.6f  Exp[ratio]=%+0.6f   dt= %+0.5e     acc= %i      rej= %i  \n",randUniform,exp(ratio/SIGMA2),dt,acc,rej);

    // Write the configuration to file
    if(MCloopi % WRITECONFIGS==0){
      normalizeAverages(tubeAve,&tau);
      writeConfig(configNew,tubeAve,MCloopi);
      zeroAverages(tubeAve,&tau);

      writeArrayPlot(arrayPlot, MCloopi);
      for(i=0;i<300;i++){
        for(j=0;j<200;j++){
          arrayPlot[i][j]=0;
      } }
    }

  }
  // GSL random number generator release memory
  gsl_rng_free (RanNumPointer);

  return(0);
}




// ============================================
//     Function Declarations
// ============================================

void rotateConfig(config **a, config **b, config **c)
// rotates the pointers
//configOld     ->  configNew (a->c)
//configCurrent ->  configOld (b->a)
//configNew     ->  configCurrent (c->b)
//the new configNew is then ready to be overwritten with new positions
{
  config *temp;

  temp=*c;
  *c=*a;
  *a=*b;
  *b=temp;
}

// ============================================
void saveConfig(config *currentConfig, config *saveConfig)
//copies ALL elements of currentConfig to saveConfig
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    saveConfig[n].Energy=currentConfig[n].Energy;
    saveConfig[n].G=currentConfig[n].G;
      for(i=0;i<NUMDIM;i++){
        saveConfig[n].pos[i]=currentConfig[n].pos[i];
        saveConfig[n].gradG[i]=currentConfig[n].gradG[i];
        saveConfig[n].LinvG[i]=currentConfig[n].LinvG[i];
} } }

// ============================================
void savePosition(position *currentpos, position *savepos)
//copies ALL elements of currentpos to save pos
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
      savepos[n].pos[i]=currentpos[n].pos[i];
} } }

// ============================================
void saveConfigtoPos(config *currentConfig, position *savepos)
//copies ONLY pos elements of currentConfig to savepos
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
      savepos[n].pos[i]=currentConfig[n].pos[i];
} } }

// ============================================
void savePostoConfig(position *currentpos, config *saveConfig)
//copies ONLY pos elements of currentpos to saveConfig
{
  int i,n;
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
      saveConfig[n].pos[i]=currentpos[n].pos[i];
} } }

// ============================================
void renorm(position *currentpos, double du, double doubleNUMu)
{
  long double alpha;
  double sum, term, term0; 
  double endPtCorr;
  int i,n;

  for(i=0;i<NUMDIM;i++)
  {
    endPtCorr=gsl_pow_int(currentpos[0].pos[i]-currentpos[NUMBEAD-1].pos[i],2)/doubleNUMu;
    sum=0.0l;
    for(n=0;n<NUMu;n++)
    {
      sum+=gsl_pow_int(currentpos[n].pos[i]-currentpos[n+1].pos[i],2);
    }
    alpha=sqrt((doubleNUMu*du*SIGMA2-endPtCorr)/(sum-endPtCorr));
    term=(1.0L-alpha)*(currentpos[NUMBEAD-1].pos[i]-currentpos[0].pos[i])/doubleNUMu;
    term0=(1.0L-alpha)*currentpos[0].pos[i];
    //term and term0 have a subtraction of roughly equal numbers and thus is not very accurate
    // alpha is ~1 with an error of 10^-4 or 5 for sample configs. This makes the routine 
    //nondeterministic between Fortran and C
    for(n=1;n<NUMu;n++)
    {
      currentpos[n].pos[i]=alpha*currentpos[n].pos[i]+term0+(((double)(n))-1.0l)*term;
    }
  }
}

//============================================
void printPositionPos(position* currentpos, int beadIndex)
//Print the positions of the position struct 
{
  printf("position x             position y\n");
  printf("%+.17e %+.17e \n",currentpos[beadIndex].pos[0],currentpos[beadIndex].pos[1]);
}

//============================================
void printConfigPos(config* currentConfig, int beadIndex)
//Print the positions of the config struct 
{
  printf("position x             position y\n");
  printf("%+.17e %+.17e \n",currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
}

//============================================
void printConfigPot(config* currentConfig, int beadIndex)
//Print all information of the config struct calculated in calcPotentials
// (Postitions, grad G, Energy, G)
{
  printf("position x            position y               gradG x               gradG y               Energy                      G \n");
  printf("%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e \n",currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1],currentConfig[beadIndex].gradG[0],currentConfig[beadIndex].gradG[1],currentConfig[beadIndex].Energy,currentConfig[beadIndex].G);
}

//============================================
void printConfigAll(config* currentConfig, int beadIndex)
//Print all information of the config struct
{
  printf("position x          position y           gradG x          gradG y         Energy               G              LinvG x            LinvG y\n");
  printf("%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n",currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1],currentConfig[beadIndex].gradG[0],currentConfig[beadIndex].gradG[1],currentConfig[beadIndex].Energy,currentConfig[beadIndex].G,currentConfig[beadIndex].LinvG[0],currentConfig[beadIndex].LinvG[1]);
}

//============================================
void printDistance(config *newConfig, position *savePos)
{
  int i,n;
  double tempSum;

  tempSum=0.0l;
  #pragma omp parallel for private(i) reduction(+:tempSum)
  for(n=1;n<NUMBEAD-1;n++){
    for(i=0;i<NUMDIM;i++){
      tempSum+=gsl_pow_int(newConfig[n].pos[i]-savePos[n].pos[i],2);
    }
  }
  printf("%+.10f \n",tempSum/(NUMBEAD-2));

}

//============================================
void calcPotentials(config *currentConfig, int beadIndex)
{
  double V    =    VFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vx   =   VxFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vy   =   VyFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vxx  =  VxxFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vxy  =  VxyFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vyy  =  VyyFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vxxx = VxxxFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vxxy = VxxyFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vxyy = VxyyFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  double Vyyy = VyyyFunc(currentConfig[beadIndex].pos[0],currentConfig[beadIndex].pos[1]);
  
  currentConfig[beadIndex].Energy = V;
  currentConfig[beadIndex].G = 0.5*(Vx*Vx +Vy*Vy) - TEMP*(Vxx +Vyy);
  currentConfig[beadIndex].gradG[0] = Vx*Vxx + Vy*Vxy - TEMP*(Vxxx + Vxyy);
  currentConfig[beadIndex].gradG[1] = Vx*Vxy + Vy*Vyy - TEMP*(Vxxy + Vyyy);
}

//============================================
void LInverse(config* currentConfig, double doubleNUMl, double du, double vecdg[NUMl], double veci1[NUMl], double veci0[NUMl])
{
  double lasti;
  int i,n;
  double du2 = du*du;


  for(i=0;i<NUMDIM;i++)
  {
    #pragma omp parallel for
    for(n=0;n<NUMl;n++){
      vecdg[n]=currentConfig[n+1].gradG[i]*du2;
    }

    veci0[0]=vecdg[0];
    veci1[0]=vecdg[0];
    //this for loop is recursive!!!!
    for(n=1;n<NUMl;n++){
      veci0[n]=veci0[n-1]+vecdg[n];
      veci1[n]=veci1[n-1]+((double)(n+1))*vecdg[n];
    }

    lasti=veci0[NUMl-1]-(currentConfig[NUMBEAD-1].pos[i]-currentConfig[0].pos[i]+veci1[NUMl-1])/((double)(NUMl+1));
    #pragma omp parallel for
    for(n=0;n<NUMl;n++){
      currentConfig[n+1].LinvG[i]=currentConfig[0].pos[i]+((double)(n+1))*(veci0[n]-lasti) - veci1[n];
    }
    currentConfig[0].LinvG[i]=currentConfig[0].pos[i];
    currentConfig[NUMBEAD-1].LinvG[i]=currentConfig[NUMBEAD-1].pos[i];

  }
}

//============================================
void preconditionSPDE(config* currentConfig, config* newConfig, double du, double dt, double doubleNUMu, double bb[NUMBEAD], double GaussRandArray[NUMu], gsl_rng *RanNumPointer)
//generates a preconditioned step using the Stochastic Partial Differential Equation
//currentConfig is the incoming configuration that has LinvG and GradG calculated
//newConfig is a temp array that is used to make all of the calsulations without touching currentConfig.pos[*][*]
//newConfig.pos is saved to currentConfig.pos before exiting the function
{

  double qvvel = 0.0l;
  double qvpos = 0.0l;
  int i,n;

  double h=sqrt(2.0l * dt);
  double co=(4.0l-h*h)/(4.0l+h*h);
  double si=(4.0l*h)/(4.0l+h*h); //(h/2)*si
  double hOverTwoSi=(2.0l*h*h)/(4.0l+h*h); //(h/2)*si

  //for grahm shmidt
  double alpha, alphaNum, alphaDenom;

  for(i=0;i<NUMDIM;i++)
  {
    generateBB(bb, du, dt, GaussRandArray, RanNumPointer);
    //need to make the bb orthogonal to pos without the linear term
    //store pos w/o linear term in newconfig.pos temporarily
    #pragma omp parallel for
    for(n=0;n<NUMBEAD;n++){
      newConfig[n].pos[i]=currentConfig[n].pos[i]-currentConfig[0].pos[i]-(((double)(n))*(currentConfig[NUMBEAD-1].pos[i]-currentConfig[0].pos[i]))/((double)(NUMBEAD -1));
    }
  
    //Gram Schmidt orthogonalization
    alphaNum = 0.0l;
    alphaDenom = 0.0l;
    #pragma omp parallel for reduction(+:alphaNum,alphaDenom)
    for(n=1;n<NUMBEAD;n++)
    {
      alphaNum+=(bb[n]-bb[n-1])*(newConfig[n].pos[i]-newConfig[n-1].pos[i]);
      alphaDenom+=(newConfig[n].pos[i]-newConfig[n-1].pos[i])*(newConfig[n].pos[i]-newConfig[n-1].pos[i]);
    }
    alpha=alphaNum/alphaDenom;

    #pragma omp parallel for
    for(n=0;n<NUMBEAD;n++){ bb[n]=bb[n]-alpha*newConfig[n].pos[i];}

    renormBB(bb, du, doubleNUMu);

    #pragma omp parallel for
    for(n=0;n<NUMBEAD;n++){
      newConfig[n].pos[i]=hOverTwoSi*currentConfig[n].LinvG[i] + si*bb[n] + co*currentConfig[n].pos[i];
    }

    //calculate the quadratic variation 
    #pragma omp parallel for reduction(+:qvpos,qvvel)
    for(n=1;n<NUMBEAD;n++){
    qvpos+=gsl_pow_int((newConfig[n].pos[i]-newConfig[n-1].pos[i]),2);
    qvvel+=gsl_pow_int((bb[n]-bb[n-1]),2);
    }
  }

  //print the quadratic variation 
  qvvel *= 0.5/(2.0l*du*((double)(NUMBEAD-1)));
  qvpos *= 0.5/(2.0l*du*((double)(NUMBEAD-1)));
  printf("qvvel=%0.10f      qvpos=%0.10f \n",qvvel,qvpos);

}
  

//============================================
void generateBB(double bb[NUMBEAD], double du, double dt, double GaussRandArray[NUMu], gsl_rng *RanNumPointer)
//generate a Brownian bridge and store to bb
{
  int n;
  double sqdu, xn;

  //generate NUMBEAD of Gaussian random nums; save to GaussRandArray[]
  GenGaussRand(GaussRandArray, RanNumPointer, 1.0l);

  //****************************
  // read in random numbers for testing precondition function
  // remove for real random numbers to be used
  //int linenum=0;
  //FILE *fptr = fopen("Random_Numbers.dat","r");
  //while( EOF != fscanf(fptr,"%lf",
  //&(GaussRandArray[linenum])) ) {
  //  linenum++;}
  //****************************

  sqdu=SIGMA*sqrt(du);
  bb[0]=0.0l;
  for(n=1;n<NUMBEAD;n++)
  {
    bb[n]=bb[n-1]+sqdu*GaussRandArray[n-1];
  }

  xn=bb[NUMu]/((double)(NUMu));
  #pragma omp parallel for
  for(n=1;n<NUMu;n++)
  {
    bb[n]-=((double)(n))*xn;
  }
  bb[NUMBEAD-1]=0.0l;
}

// ============================================
void renormBB(double bb[NUMBEAD], double du, double doubleNUMu)
{
  int n;
  double endPtCorr;
  double sum, term, term0;
  double alpha;

  endPtCorr=gsl_pow_int(bb[0]-bb[NUMBEAD-1],2)/doubleNUMu;

  sum=0.0l;
  #pragma omp parallel for reduction(+:sum)
  for(n=0;n<NUMu;n++)
  {
    sum+=gsl_pow_int(bb[n]-bb[n+1],2);
  }

  alpha=sqrt((doubleNUMu*du*SIGMA2-endPtCorr)/(sum-endPtCorr));
  term=(1.0l-alpha)*(bb[NUMBEAD-1]-bb[0])/doubleNUMu;
  term0=(1.0l-alpha)*bb[0];
  //term and term0 have a subtraction of roughly equal numbers and thus is not very accurate
  // alpha is ~1 with an error of 10^-4 or 5 for sample configs. This makes the routine 
  //nondeterministic between Fortran and C

  #pragma omp parallel for
  for(n=1;n<NUMu;n++)
  {
    bb[n]=alpha*bb[n]+term0+((double)(n-1))*term;
  }
}

//============================================
void GenGaussRand(double GaussRand[NUMu], gsl_rng *RanNumPointer, double StdDev)
//generate NUMu of Gaussian random nums; save to GaussRandArray[]
{
  int i;
  for(i=0;i<NUMu; i++)
  {
    GaussRand[i]= gsl_ran_gaussian(RanNumPointer,StdDev);
  } 
}

//============================================
void ProbAccRatio(config *currentConfig, config *newConfig, double dt, double du, double *ratio)
// calculate the probability accaptance ratio for the step
{

  //there are 3 terms in the error (see notes)
  //lambda1 with h/4 in front (l1h4)
  //lambda1 with h^/16 in front (l1hh16)
  //lambda2 (l2)
  double l1h4, l1hh16, l2;
  int n;

  double h=sqrt(2.0l*dt);
  //double co=(4.0l-(h*h))/(4.0l+(h*h));
  //double si=sqrt(1.0l-(co*co));
  double cot=(4.0l-h*h)/(4.0l*h);
  double csc=(4.0l+h*h)/(4.0l*h);
  
  l1h4=0.0l;
  l1hh16=0.0l;
  l2=0.0l;

  #pragma omp parallel for reduction(+:l1h4,l1hh16,l2)
  for(n=0;n<NUMBEAD;n++)
  {
    //this only works for the two dimensional case
    l1h4+=(cot*currentConfig[n].pos[0]-csc*newConfig[n].pos[0])*currentConfig[n].gradG[0] + (csc*currentConfig[n].pos[0]-cot*newConfig[n].pos[0])*newConfig[n].gradG[0]+(cot*currentConfig[n].pos[1]-csc*newConfig[n].pos[1])*currentConfig[n].gradG[1] + (csc*currentConfig[n].pos[1]-cot*newConfig[n].pos[1])*newConfig[n].gradG[1];
    l1hh16+=currentConfig[n].gradG[0]*currentConfig[n].LinvG[0] - newConfig[n].gradG[0]*newConfig[n].LinvG[0]+currentConfig[n].gradG[1]*currentConfig[n].LinvG[1] - newConfig[n].gradG[1]*newConfig[n].LinvG[1];
    l2+=newConfig[n].G-currentConfig[n].G;
  }
  l1h4*=h*0.25l;
  l1hh16*=h*h*0.0625l;
  l2*=0.5l;

  *ratio+=du*(l1h4+l1hh16+l2);
}

//============================================
void MolecularDynamics(config *oldConfig, config *currentConfig, config *newConfig, double du, double dt)
// perform the MD step
{

  int i,n;
  double h=sqrt(2.0l*dt);
  double co=(4.0-h*h)/(4.0l+h*h);
  double si=sqrt(1.0l-co*co);
  double twoCosMinusOne = 2.0l*co-1.0l;
  double sinH = si*h;

  #pragma omp parallel for private(i)
  for(n=0;n<NUMBEAD;n++) {
    for(i=0;i<NUMDIM;i++)  {
      newConfig[n].pos[i]= (currentConfig[n].pos[i]-oldConfig[n].pos[i]) + twoCosMinusOne*currentConfig[n].pos[i] + sinH*currentConfig[n].LinvG[i];
    }
  }
}

//============================================
void writeConfig(config *newConfig, averages *tubeAve, int MCloopi)
{
//print the configs to file
//File order: posx   posy   meanx   meany   posx^2   posx*posy   posy^2
  char filename[50];
  int i,n;

  sprintf(filename,"%s-T%.2f-pos%07d.dat",PotentialString,TEMP,MCloopi);
  FILE * pWritePos;
  pWritePos = fopen(filename,"w");
  for(n=0;n<NUMBEAD;n++){
    for(i=0;i<NUMDIM;i++){
    fprintf(pWritePos, "%+.15e \t",newConfig[n].pos[i]);
    }
    for(i=0;i<NUMDIM;i++){
    fprintf(pWritePos, "%+.15e \t",tubeAve[n].mean[i]);
    }
    fprintf(pWritePos, "%+.15e \t",tubeAve[n].xx);
    fprintf(pWritePos, "%+.15e \t",tubeAve[n].xy);
    fprintf(pWritePos, "%+.15e",tubeAve[n].yy);
    fprintf(pWritePos, "\n");
  }
  fclose(pWritePos);

  //call the python script to plot the paths
  char * path;
  char pyCall[300];
  path=getenv("PWD");
  if(path != NULL){
    strcpy(pyCall,"2dHMC-plot ");
    strcat(pyCall,path);
    strcat(pyCall,"/");
    strcat(pyCall,filename);
    strcat(pyCall," 1000");
    system(pyCall); 
  }

}

//============================================
void zeroAverages(averages *tubeAve, int *tau)
{
  int n;
  for(n=0;n<NUMBEAD;n++)
  {
    tubeAve[n].mean[0]=0.0l;
    tubeAve[n].mean[1]=0.0l;
    tubeAve[n].xx=0.0l;
    tubeAve[n].yy=0.0l;
    tubeAve[n].xy=0.0l;
  }
  *tau=0;
}

//============================================
void accumulateAverages(averages *tubeAve, config *newConfig, int *tau)
{
  int n;
  *tau=*tau+1;
  #pragma omp parallel for
  for(n=0;n<NUMBEAD;n++)
  {
    tubeAve[n].mean[0]+=newConfig[n].pos[0];
    tubeAve[n].mean[1]+=newConfig[n].pos[1];
    tubeAve[n].xx+=newConfig[n].pos[0]*newConfig[n].pos[0];
    tubeAve[n].yy+=newConfig[n].pos[1]*newConfig[n].pos[1];
    tubeAve[n].xy+=newConfig[n].pos[0]*newConfig[n].pos[1];
  }
}
    
//============================================
void normalizeAverages(averages *tubeAve, int *tau)
{
  int n;
  double oneOverTau=1.0l/((double)(*tau));

  #pragma omp parallel for
  for(n=0;n<NUMBEAD;n++)
  {
    tubeAve[n].mean[0]*=oneOverTau;
    tubeAve[n].mean[1]*=oneOverTau;
    tubeAve[n].xx*=oneOverTau;
    tubeAve[n].yy*=oneOverTau;
    tubeAve[n].xy*=oneOverTau;
  }
}

//============================================
void accumulateArrayPlot(int arrayPlot[300][200], config *currentConfig)
//accumulate the averages for the array plot
// if in the bound of the array
//this is not all that general (size of matrix is hardcoded)
{
  int xbin;
  int ybin;
  int n;
  
  for(n=0;n<NUMBEAD;n++)
  {
    xbin=(int)(floor( (currentConfig[n].pos[0]+1.5)*100.));
    ybin=(int)(floor( (currentConfig[n].pos[1]+0.5)*100.));
    if((xbin>=0) && (xbin<300) && (ybin>=0) && (ybin<200)){
      arrayPlot[xbin][ybin]++;
    }
  }
}

//============================================
void writeArrayPlot(int arrayPlot[300][200], int MCloopi)
{
//print the arrayPlot to file
  char filename[50];
  int i,j;

  sprintf(filename,"arrayPlot%s-T%.2f-pos%07d.dat",PotentialString,TEMP,MCloopi);
  FILE * pWritePos;
  pWritePos = fopen(filename,"w");
  for(i=0;i<300;i++){
    for(j=0;j<200;j++){
    fprintf(pWritePos, "%d \t",arrayPlot[i][j]);
    }
    fprintf(pWritePos, "\n");
  }
  fclose(pWritePos);
}

