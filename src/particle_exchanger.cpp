#include "snapshot.h"
#include "particle_exchanger.h"

void create_Mpi_RemoteParticleType(MPI_Datatype& dtype)
{
  /*to create the struct data type for communication*/	
  RemoteParticle_t p;
  #define NumAttr 10
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;
  
  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;
  
  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(Id, MPI_HBT_INT, 1)
  RegisterAttr(ComovingPosition, MPI_HBT_REAL, 3)
  RegisterAttr(PhysicalVelocity, MPI_HBT_REAL, 3)
  RegisterAttr(Mass, MPI_HBT_REAL, 1)
  RegisterAttr(Order, MPI_HBT_INT, 1)
  #undef RegisterAttr
  assert(i<=NumAttr);
  
  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}