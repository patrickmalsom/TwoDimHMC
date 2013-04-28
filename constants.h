// =========================================
// Preprocessor Definitions
// =========================================

// Path Definitions
#define NUMDIM    2           // Spacial Dimensions
#define NUMBEAD   262144     // Path Points
#define DU        0.001       // du (see potentials.h for estimate of DU)
#define PREDT     2.0         // dt=PreDT*du^2 (path time)

// Temperature Definition
#define TEMP      0.15

// Incrimenter Definitions
#define NUMMD     1000         // Number of MD steps 
//      NUMMD     ~3/(2*sqrt(2*PreDT*DU^2)) <- Approx optimal value of NUMMD
#define NUMMC     1000        // Number of Metropolis Hastings MC steps

// Constants for writing to stdout and config
#define WRITESTDOUT  50       // How often to print to stdout (# of MD loops)
#define WRITECONFIGS 100       // How often to save config to file (# of MHMC steps)
const char PotentialString[]="2WellSym";// Potential Description 
