/* if the particles are ordered in the snapshot and labelled from 1 to N in the group file, then define PID_ORDERED to avoid complex IndexTables*/
#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>
#include <cstdlib>

#include "hash.h"
#include "mysort.h"

//=====general ID2Index table======//
/* the hash-table implementation here is by sorting Id and use binsearch to locate keys*/
/* more general functions hcreate(),hsearch()... exists in glibc; but binsearch should be
 * more efficient than the general hsearch() I guess?*/

template <class Key_t, class Index_t>
inline bool CompPair(const IndexedKey_t<Key_t, Index_t> & a, const IndexedKey_t<Key_t, Index_t> & b)
{
  return (a.Key<b.Key);
};
template <class Key_t, class Index_t>
void MappedIndexTable_t<Key_t, Index_t>::Fill(const KeyList_t <Key_t, Index_t> &Keys, Index_t null_index)
{
  BaseClass_t::NullIndex=null_index;
  Index_t n=Keys.size();
  Map.resize(n);
  #pragma omp parallel for
  for(Index_t i=0;i<n;i++)
  {
    Map[i].Key=Keys.GetKey(i);
    Map[i].Index=Keys.GetIndex(i);
  }
  MYSORT(Map.begin(), Map.end(), CompPair<Key_t, Index_t>);
}
template <class Key_t, class Index_t>
void MappedIndexTable_t<Key_t, Index_t>::Clear()
{
  vector <Pair_t>().swap(Map);
}
template <class Key_t, class Index_t>
inline int CompKeyWithPair(const void *a, const void *b)//used to sort Id in ascending order; 
{
  Key_t va=* static_cast<const Key_t *>(a);
  Key_t vb=static_cast<const IndexedKey_t<Key_t, Index_t> *> (b)->Key;
  if(va>vb) return 1;
  if(va<vb) return -1;
  return 0;
};
template <class Key_t, class Index_t>
Index_t MappedIndexTable_t<Key_t, Index_t>::GetIndex(const Key_t key) const
{//maybe implement the exception here? could be slow... test it first.
  if(key<0) return BaseClass_t::NullIndex;
  Pair_t *p=(Pair_t *) bsearch(&key,Map.data(),Map.size(),sizeof(Pair_t),CompKeyWithPair<Key_t, Index_t>);
  if(NULL==p) return BaseClass_t::NullIndex;  //no match
  return p->Index;  
}
template <class Key_t, class Index_t>
void FlatIndexTable_t<Key_t, Index_t>::Fill(const KeyList_t<Key_t, Index_t> &Keys, Index_t null_index)
{
  BaseClass_t::NullIndex=null_index;
  Clear();
  Index_t n=Keys.size();
  if(0==n) return;
  
  Key_t keymin, keymax;
  keymin=keymax=Keys.GetKey(0);
  #pragma omp parallel for reduction(max: keymax) reduction(min: keymin)
  for(Index_t i=1;i<n;i++)
  {
    Key_t key=Keys.GetKey(i);
    if(key>keymax)
      keymax=key;
    if(key<keymin)	//cannot do elseif due to the way omp initalizes the private copy of reductions
      keymin=key;
  }
  KeyMax=keymax; KeyMin=keymin;
  KeySpan=KeyMax-KeyMin+1;
  Index=new Index_t[KeySpan];
  Offset=KeyMin;
  Index-=Offset;/*shift PIndex by KeyMin element so that it is accessed through PIndex[KeyMin~KeyMax],
				  i.e.,its subscript ranges the same as particle ID range. PID_all[PIndex[Id]]=Id;*/
  
  #pragma omp parallel 
  {
  #pragma omp for
  for(Key_t i=KeyMin;i<=KeyMax;i++)
	  Index[i]=BaseClass_t::NullIndex; //initialize with -1, although not necessary
  /*====make ID index for query====*/
  #pragma omp for
  for(Index_t i=0;i<n;i++)
    Index[Keys.GetKey(i)]=Keys.GetIndex(i);		
  }
}
template <class Key_t, class Index_t>
void FlatIndexTable_t<Key_t, Index_t>::Clear()
{
	if(KeySpan)
	{
		KeySpan=0;
		Index+=Offset;
		delete [] Index;
	}
}
template <class Key_t, class Index_t>
Index_t FlatIndexTable_t<Key_t, Index_t>::GetIndex(const Key_t key) const
{
	if(KeySpan==0||key<KeyMin||key>KeyMax) return BaseClass_t::NullIndex;//no match
	return Index[key];
}
