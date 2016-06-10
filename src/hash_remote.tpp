#include "hash.h"
#include "snapshot.h"

template <class Pair_t, class Val_t>
inline int CompPairWithValue(const Pair_t a, const Val_t b)
{
  return (a.Key<b);
};
template <class Key_t, class Index_t>
template <class ParticleIdList_T>
void MappedIndexTable_t<Key_t, Index_t>::GetIndices(ParticleIdList_T &particles) const
{ 
#define ALWAYS_BATCH_BINARY_SEARCH
  
#ifdef ALWAYS_BATCH_BINARY_SEARCH
  GetIndicesRecursive(particles, 0, particles.size(), Map.begin(), Map.end());//batch-binary-search: is this always faster?
#else  
  if(particles.size()<NumQueryCrit)//do individual binary search
  {
	for(auto &&p: particles)
	  p.Id=GetIndex(p.Id);
	return;
  }
  //otherwise do batch sequential search
  auto &null=BaseClass_t::NullIndex;
  
  auto it_p=particles.begin();
  auto it_map=Map.begin();
  while(true)
  {
	if(it_p==particles.end()) return;
	
	if(it_map==Map.end()) break;
	
	if(it_p->Id<it_map->Key)
	{
	  it_p->Id=null;
	  ++it_p;
	}
	else if(it_p->Id==it_map->Key)
	{
	  it_p->Id=it_map->Index;
	  ++it_p;
	}
	else
	  ++it_map;
  }
  
  while(true)
  {
	  it_p->Id=null;
	  ++it_p;
	  if(it_p==particles.end()) return;
  }
#endif  
}

template <class Key_t, class Index_t>
template <class ParticleIdList_T>
void MappedIndexTable_t<Key_t, Index_t>::GetIndicesRecursive(ParticleIdList_T &particles, HBTInt imin, HBTInt imax, MapIter_t MapBegin,  MapIter_t MapEnd) const
{
  //GetIndices of particles in storage range [imin, imax) from map [MapBegin, MapEnd).
  auto &null=BaseClass_t::NullIndex;
 
  if(MapBegin==MapEnd)
  {
	for(HBTInt i=imin;i<imax;i++)
	  particles[i].Id=null;
	return;
  }
  
  if(imin>=imax) return;
  
  HBTInt imid;
  if(imax-imin==1) 
	imid=imin;
  else
	imid=(imin+imax)/2;
  Key_t key=particles[imid].Id;
  MapIter_t MapMid=lower_bound(MapBegin, MapEnd, key, CompPairWithValue<Pair_t, Key_t>);
  MapIter_t MapEndLeft=MapMid, MapBeginRight=MapMid;
  if(MapMid==MapEnd||MapMid->Key>key)
	particles[imid].Id=null;
  else
  {
	particles[imid].Id=MapMid->Index;
	++MapEndLeft;
  }
	
  GetIndicesRecursive(particles, imin, imid, MapBegin, MapEndLeft);
  GetIndicesRecursive(particles, imid+1, imax, MapBeginRight, MapEnd);
}

template <class Key_t, class Index_t>
template <class ParticleIdList_T>
void FlatIndexTable_t<Key_t, Index_t>::GetIndices(ParticleIdList_T &particles) const
{
  for(auto &&p: particles)
	p.Id=GetIndex(p.Id);
}

template <class ParticleIdList_t>
void ParticleSnapshot_t::GetIndices(ParticleIdList_t& particles) const
{//ParticleIdList_t is a list of particle structs containing at least an Id field
  if(HBTConfig.ParticleIdNeedHash)
	MappedHash.GetIndices(particles);
  else
	FlatHash.GetIndices(particles);
}
