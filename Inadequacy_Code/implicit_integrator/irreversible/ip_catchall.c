 /*******************************************************************
 *
 * ip_catchall.c
 *
 * This is the wrapper for the model inadequacy code.  It calls the
 * driver compute.c.
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *******************************************************************/

#include "compute.h"
#include <yaml-cpp/yaml.h>

int main(int argc, char* argv[])
{

  // Initialize QUESO environment
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env =
    new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
  
  // Call application
  computeAllParams(*env);

  // Finalize QUESO environment
  delete env;
  MPI_Finalize();

  return 0;
}
