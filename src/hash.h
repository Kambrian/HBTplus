/*hash table data for ID2Index*/
//TODO: rename Index to Value, to avoid confusion when used for (Key,Val) pairs rather than ordered Keys.
#ifndef HASH_HEADER_INCLUDED
#define HASH_HEADER_INCLUDED

#include "datatypes.h"
#include <exception>
#include <string>

class RemoteParticle_t;

template <class Key_t, class Index_t>
class KeyList_t
{
public:
  virtual Key_t GetKey(const Index_t i) const=0;
  virtual Index_t GetIndex(const Index_t i) const=0;
  virtual Index_t size() const=0;
};

class InvalidPIdException_t : public exception
{
private:
  HBTInt PId;
public:
  InvalidPIdException_t(HBTInt pid)
  {
	PId=pid;
  };
  const char * what () const throw ()
  {
	stringstream msg;
	msg<<"Invalid Particle Id "<<PId<<" for index lookup\n";
	return msg.str().c_str();
  };
    virtual ~InvalidPIdException_t()throw()
  {};
};

template <class Key_t, class Index_t>
struct IndexedKey_t
{
  Key_t Key;
  Index_t Index;  
};

template <class Key_t, class Index_t>
class IndexTable_t
{
public: 
//   typedef HBTInt Key_t;
//   typedef HBTInt Index_t;
  Index_t NullIndex;
  
  virtual void Fill(const KeyList_t<Key_t, Index_t> &Keys, Index_t null_index=SpecialConst::NullParticleId)=0;
  virtual void Clear()=0;
  virtual Index_t GetIndex(const Key_t key) const =0;
  typedef vector <RemoteParticle_t> ParticleIdList_T; 
  virtual void GetIndices(ParticleIdList_T &particles) const =0;
};

template <class Key_t, class Index_t>
class FlatIndexTable_t: public IndexTable_t <Key_t, Index_t>
{
private:
  typedef IndexTable_t<Key_t, Index_t> BaseClass_t;
  typedef typename BaseClass_t::ParticleIdList_T ParticleIdList_T;
  Index_t * Index;
  Index_t Offset;
  Key_t KeySpan, KeyMax, KeyMin;
public:
  FlatIndexTable_t(): Index(), Offset(0), KeySpan(0)
  {
  }
  void Fill(const KeyList_t <Key_t, Index_t> &Keys, Index_t null_index=SpecialConst::NullParticleId);
  void Clear();
  Index_t GetIndex(const Key_t key) const;
  void GetIndices(ParticleIdList_T &particles) const;
  ~FlatIndexTable_t()
  {
	Clear();
  }
};

template <class Key_t, class Index_t>
class MappedIndexTable_t: public IndexTable_t <Key_t, Index_t>
{
public:  
  typedef IndexedKey_t<Key_t, Index_t> Pair_t;  
private:
  HBTInt NumQueryCrit;
  typedef IndexTable_t<Key_t, Index_t> BaseClass_t;
  typedef typename BaseClass_t::ParticleIdList_T ParticleIdList_T;
  vector <Pair_t> Map; 
  typedef typename vector <Pair_t>::const_iterator MapIter_t;
  void GetIndicesRecursive(ParticleIdList_T &particles, HBTInt imin, HBTInt imax, MapIter_t MapBegin,  MapIter_t MapEnd) const;
public:
  MappedIndexTable_t(): Map(), NumQueryCrit()
  {
  }
  void Fill(const KeyList_t <Key_t, Index_t> &Keys, Index_t null_index=SpecialConst::NullParticleId);
  void Clear();
  Index_t GetIndex(const Key_t key) const;
  void GetIndices(ParticleIdList_T &particles) const;
  ~MappedIndexTable_t()
  {
	Clear();
  }
};

#include "hash.tpp"

#endif