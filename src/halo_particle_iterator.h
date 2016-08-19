#ifndef HALO_PARTICLE_ITERATOR_INCLUDED
#define HALO_PARTICLE_ITERATOR_INCLUDED

template <class HaloIterator>
class HaloParticleIterator_t
{
  typedef vector<Particle_t>::iterator particle_iterator;
  HaloIterator FirstHalo, EndHalo, CurrHalo;
  particle_iterator CurrPart;
public:
  HaloParticleIterator_t()=default;
  HaloParticleIterator_t(const HaloIterator &begin, const HaloIterator &end)
  {
	init(begin, end);
  }
  void init(HaloIterator begin, HaloIterator end)
  {
	while((begin!=end)&&(begin->Particles.size()==0))//skip empty ones, though not necessary for current HBT2
	  ++begin;
	FirstHalo=begin;
	EndHalo=end;
	reset();
  }
  void reset()
  {
	CurrHalo=FirstHalo;
	if(CurrHalo!=EndHalo)
	  CurrPart=FirstHalo->Particles.begin();
  }
  particle_iterator begin()
  {
	return FirstHalo->Particles.begin();
  }
  HaloParticleIterator_t<HaloIterator> & operator ++()//left operator
  {
	++CurrPart;
	while(CurrPart==CurrHalo->Particles.end())//increment halo and skip empty haloes
	{
	  ++CurrHalo;
	  if(CurrHalo==EndHalo) break;
	  CurrPart=CurrHalo->Particles.begin();
	}
	return *this;
  }
  Particle_t & operator *()
  {
	return *CurrPart;
  }
  bool is_end()
  {
	return CurrHalo==EndHalo;
  }
};


template <class HaloIterator>
class HaloNestIterator_t
{
  typedef HBTInt NestMember_t;
  typedef vector<NestMember_t>::iterator nest_iterator;
  HaloIterator FirstHalo, EndHalo, CurrHalo;
  nest_iterator CurrPart;
public:
  HaloNestIterator_t()=default;
  HaloNestIterator_t(const HaloIterator &begin, const HaloIterator &end)
  {
	init(begin, end);
  }
  void init(HaloIterator begin, HaloIterator end)
  {
	while((begin!=end)&&(begin->NestedSubhalos.size()==0))//skip empty ones, though not necessary for current HBT2
	  ++begin;
	FirstHalo=begin;
	EndHalo=end;
	reset();
  }
  void reset()
  {
	CurrHalo=FirstHalo;
	if(CurrHalo!=EndHalo)
	  CurrPart=FirstHalo->NestedSubhalos.begin();
  }
  nest_iterator begin()
  {
	return FirstHalo->NestedSubhalos.begin();
  }
  HaloNestIterator_t<HaloIterator> & operator ++()//left operator
  {
	++CurrPart;
	while(CurrPart==CurrHalo->NestedSubhalos.end())//increment halo and skip empty haloes
	{
	  ++CurrHalo;
	  if(CurrHalo==EndHalo) break;
	  CurrPart=CurrHalo->NestedSubhalos.begin();
	}
	return *this;
  }
  NestMember_t & operator *()
  {
	return *CurrPart;
  }
  bool is_end()
  {
	return CurrHalo==EndHalo;
  }
};

#endif