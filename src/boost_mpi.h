#ifndef HBT_BOOST_MPI_H_INCLUDED
#define HBT_BOOST_MPI_H_INCLUDED

#include <array>
#include <boost/version.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
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

#endif