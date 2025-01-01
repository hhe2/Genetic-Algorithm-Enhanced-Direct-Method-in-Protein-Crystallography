#define PDB_CODE_REFE_HIST "1rb6"
#define PDB_CODE "1no4"

//cell parameters
#define CRYST_A 114.240
#define CRYST_B 114.164 
#define CRYST_C 60.297
#define CRYST_ALPHA ( 90.00 / 180.00 * PI )
#define CRYST_BETA  ( 90.00 / 180.00 * PI )
#define CRYST_GAMMA ( 90.00 / 180.00 * PI )

//grid size
#define ANUM 104 //mod(ANUM,4)==0
#define BNUM 104 //mod(BNUM,4)==0
#define CNUM 56 //mod(CNUM,4)==0 for P4132

using namespace std;

const float resoCutoff = 2.2;
const float resoCutoffLow = 15.0;

const float solvContIni = 0.70; //sfcheck 70% pdb matthews 77.7%
const float solvContFina = 0.70;
float solvCont = solvContIni;

//weighted observed data
const float sigmWeigObsIni = 1.5;
float sigmWeigObs = sigmWeigObsIni;

//weighted average density
const float sigmWeigAvgIni = 4.0;
const float sigmWeigAvgFina = 3.0;
float sigmWeigAvg = sigmWeigAvgIni;
char updateSpeedWeigAvgDens[] = "normally"; 
//"normally" 100%; "slowly" 50% + 50%; "very slowly" 5% + 95%; "extremely slowly" 1% + 99%

const int numIter = 10000;
const bool shareMeanMissReflFlag = false;

//shrink protMask for HPR and CHIO
float shrinkProtMaskPercentIni = 0.05;
float shrinkProtMaskPercentFina = 0;

//MPI
const int numCrossPosi = 11;
const float crossGenePercent = 0.05; //percent of ASU for each crossover position
const float mutateGenePercent = 0.01; //percent of ASU for mutation

const float goodProtMaskThreshold = 0.90; //protMaskMatch > threshold
const float convThreshold = 70.0; //meanPhaseError < threshold

#define NCSANUM 180
#define NCSBNUM 180
#define NCSCNUM 180

//sigma of weighting function
const float sigmWeigNcs = 5.0;

//sigma of weighting function for updating centMassOrth
const float sigmWeigCentMass = 15.0;

const int numNcsOper = 2;
float ncsOper[numNcsOper][3][4] = 
{  	{ { 1.0, 0.0, 0.0, 0.0 }, 
	  { 0.0, 1.0, 0.0, 0.0 }, 
	  { 0.0, 0.0, 1.0, 0.0 } }, 
	{ { -0.2094,  0.8364, -0.5066,  51.1571},
	  { -0.6932, -0.4923, -0.5263, 141.5947},
	  { -0.6896,  0.2410,  0.6829,  42.6693} } };

const int numCentMass = 1;
const float centMassOrthTrue[numCentMass][3] = {22.242, 0.011, -1.513}; //center of molecule
float centMassOrth[numCentMass][3] = {22.242, 0.011, -1.513};

const float pointAOnNcsAxisTrue[3] = {22.242, 0.011, -1.513};
const float vectOnNcsAxisTrue[3] = {0,0,1}; //{-0.016, 0.036, 27.957};
float pointAOnNcsAxis[3] = {22.242, 0.011, -1.513};
float pointBOnNcsAxis[3] = {};
float vectOnNcsAxis[3] = {0,0,1};;

//self-rotation function at section kappa: lateral angle psi; azimuthal angle phi.
const float kappa = 180.0 / 180.0 * 3.1415926;
const float lateAnglPsi = 26.0 / 180.0 * 3.1415926;
const float azimAnglPhi = -175.0 / 180.0 * 3.1415926; //+-5 deg, and +-175 deg
//const float vectOnNcsAxisIni[3] = {-0.43670, -0.03821, 0.89879};

const int numNcsAxisCand = 101;//candNcsAxis + initNcsAxis

const float ncsAxisDeviDistMax = 10.0;	//10 angstroms
const float ncsAxisDeviAnglMax = 10.0 / 180.0 * 3.1415926; //10 degrees

const float rotaNcsAngl = 180.0 / 180.0 * 3.1415926;
const bool rotaNcsFlag = true;
const bool tranNcsFlag = false;
const float tranNcsDist = 0;

#include "../header/commPara.hpp"
#include "../header/P21212.hpp"

