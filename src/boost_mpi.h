#ifndef HBT_BOOST_MPI_H_INCLUDED
#define HBT_BOOST_MPI_H_INCLUDED

#include <array>
#include <boost/version.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
// #include <mpi.h>
namespace mpi = boost::mpi;

namespace boost {
namespace serialization {
#if BOOST_VERSION<105600
template<class Archive, class T, size_t N>
void serialize(Archive & ar, std::array<T,N> & a, const unsigned int version)
{
  ar & boost::serialization::make_array(a.data(), a.size());
}
#endif
} // namespace serialization
} // namespace boost


template <class T>
void VectorAllToAll(mpi::communicator &world, vector < vector<T> > &SendVecs, vector < vector <T> > &ReceiveVecs, MPI_Datatype dtype)
{
  vector <int> SendSizes(world.size()), ReceiveSizes(world.size());
  for(int i=0;i<world.size();i++)
	SendSizes[i]=SendVecs[i].size();
  MPI_Alltoall(SendSizes.data(), 1, MPI_INT, ReceiveSizes.data(), 1, MPI_INT, world);
  
  ReceiveVecs.resize(world.size());
  for(int i=0;i<world.size();i++)
	ReceiveVecs[i].resize(ReceiveSizes[i]);
  
  vector <MPI_Datatype> SendTypes(world.size()), ReceiveTypes(world.size());
  for(int i=0;i<world.size();i++)
  {
	MPI_Aint p;
	MPI_Address(SendVecs[i].data(), &p);
	MPI_Type_create_hindexed(1, &SendSizes[i], &p, dtype, &SendTypes[i]);
	MPI_Type_commit(&SendTypes[i]);
	
	MPI_Address(ReceiveVecs[i].data(), &p);
	MPI_Type_create_hindexed(1, &ReceiveSizes[i], &p, dtype, &ReceiveTypes[i]);
	MPI_Type_commit(&ReceiveTypes[i]);
  }
  vector <int> Counts(world.size(),1), Disps(world.size(),0);
  MPI_Alltoallw(MPI_BOTTOM, Counts.data(), Disps.data(), SendTypes.data(),
				MPI_BOTTOM, Counts.data(), Disps.data(), ReceiveTypes.data(), world);
  
  for(int i=0;i<world.size();i++)
  {
	MPI_Type_free(&SendTypes[i]);
	MPI_Type_free(&ReceiveTypes[i]);
  }
  
}

#endif