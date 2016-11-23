/* example file demonstrating the creation of a custom halo reader. 
 * The reader will be called in halo_io.cpp under "my_group_format".
 
 * To use this reader, set
 
    GroupFileFormat my_group_format
 
 * in config file.
 */

#include "../mymath.h"
#include "custom_io.h"

HBTInt MyGroupReader(int SnapshotId, vector <Halo_t> Halos)
{
  HBTInt NumberOfGroups;
  HBTInt *NumberOfParticlesInGroup;
  HBTInt **ParticlesInGroup;
  
/* !!!
 * add some code or function to read your catalogue here.
 * suppose after this you have NumberOfGroups, NumberOfParticlesInGroup[], ParticlesInGroup[][];
 * !!!
*/

  HBTInt ntotal=0;
  Halos.resize(NumberOfGroups);
  for(HBTInt i=0;i<NumberOfGroups;i++)
  {
    Halos[i].Particles.resize(NumberOfParticlesInGroup[i]);
    for(HBTInt j=0;j<NumberOfParticlesInGroup[i];j++)
    {
    Halos[i].Particles[j]=ParticlesInGroup[i][j];
    }
    ntotal+=NumberOfParticlesInGroup[i];
  }
  return ntotal;
}