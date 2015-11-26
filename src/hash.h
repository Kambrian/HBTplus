/*hash table data for ID2Index*/
#ifndef HASH_HEADER_INCLUDED
#define HASH_HEADER_INCLUDED

#include "datatypes.h"
#include <exception>
#include <string>

template <class Key_t, class Index_t>
class KeyList_t
{
public:
  virtual Key_t GetKey(const Index_t i) const=0;
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
  static const Index_t NullIndex=SpecialConst::NullParticleId;
  
  virtual void Fill(const KeyList_t<Key_t, Index_t> &Keys)=0;
  virtual void Clear()=0;
  virtual Index_t GetIndex(const Key_t key) const =0;
};

template <class Key_t, class Index_t>
class FlatIndexTable_t: public IndexTable_t <Key_t, Index_t>
{
private:
  typedef IndexTable_t<Key_t, Index_t> BaseClass_t;
  Index_t * Index;
  Index_t Offset;
  Key_t KeySpan;
public:
  FlatIndexTable_t(): Index(), Offset(0), KeySpan(0)
  {
  }
  void Fill(const KeyList_t <Key_t, Index_t> &Keys);
  void Clear();
  Index_t GetIndex(const Key_t key) const;
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
  typedef IndexTable_t<Key_t, Index_t> BaseClass_t;
  vector <Pair_t> Map; 
public:
  MappedIndexTable_t(): Map()
  {
  }
  void Fill(const KeyList_t <Key_t, Index_t> &Keys);
  void Clear();
  Index_t GetIndex(const Key_t key) const;
  ~MappedIndexTable_t()
  {
	Clear();
  }
};

#include "hash.tpp"

#endif