#include "mymath.h"

//TODO:discard the fortran-style ll; use struct or indexed table to parallelize the linklist!

class PositionData_t
{
public:
  virtual const HBTxyz & operator [](HBTInt i) const=0;
  /*virtual const HBTReal GetPos(HBTInt i, int j) const
  {
    return (*this)[i][j];
  }*/
  virtual size_t size() const=0;
};
class Linkedlist_t
{
private:
  int NDiv, NDiv2;
  HBTInt NumPart;
  bool PeriodicBoundary;
  HBTReal BoxSize, BoxHalf;
  vector <HBTInt> HOC;
  vector <HBTInt> List;
  HBTReal Range[3][2];
  HBTReal Step[3];
  PositionData_t *Particles;
  int RoundGridId(int i)
  //to correct for rounding error near boundary
  {
    return i<0?0:(i>=NDiv?NDiv-1:i);
  }
  int ShiftGridId(int i)
  /*to correct for periodic conditions; 
  only applicable when def PERIODIC_BDR and ll.UseFullBox=1 */
  {
	i=i%NDiv;
	if(i<0) i+=NDiv;
	return i;
  }
  int FixGridId(int i)
  {
    if(PeriodicBoundary)
      return ShiftGridId(i);
    return RoundGridId(i);
  }
  HBTInt Sub2Ind(int i, int j, int k)
  {
    return i+j*NDiv+k*NDiv2;
  }
  HBTInt GetHOC(int i, int j, int k)
  {
    return HOC[Sub2Ind(i,j,k)];
  }
  HBTInt GetHOCSafe(int i, int j, int k)
  {
    return HOC[Sub2Ind(FixGridId(i), FixGridId(j), FixGridId(k))];
  }
  HBTReal Distance(const HBTxyz &x, const HBTxyz &y)
  {
	  HBTxyz dx;
	  dx[0]=x[0]-y[0];
	  dx[1]=x[1]-y[1];
	  dx[2]=x[2]-y[2];
	  if(PeriodicBoundary)
	  {
	    #define _NEAREST(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))
	    dx[0]=_NEAREST(dx[0]);
	    dx[1]=_NEAREST(dx[1]);
	    dx[2]=_NEAREST(dx[2]);
	    #undef _NEAREST
	  }
	  return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  }
public:
  Linkedlist_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false)
  {
    init(ndiv, data, boxsize, periodic);
  }
  init(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false) 
  {
    NumPart=data->size();
    NDiv=ndiv;
    Particles=data;
    PositionData_t &particles=*Particles;
    HOC.resize(ndiv*ndiv*ndiv);
    List.resize(data->size());
    BoxSize=boxsize;
    PeriodicBoundary=periodic;
    if(BoxSize==0.) PeriodicBoundary=false; //only effective when boxsize is specified
    BoxHalf=BoxSize/2.;
    
    NDiv2=NDiv*NDiv;
    HBTInt i,j,grid[3];
    HBTInt ind;
    //~ float range[3][2],step[3];
    cout<<"creating linked list..\n";
    /*determining enclosing cube*/
    if(BoxSize)
    {
      for(i=0;i<3;i++)
      {
	Range[i][0]=0.;
	Range[i][1]=BoxSize;
      }
      for(j=0;j<3;j++)
	Step[j]=BoxSize/NDiv;	
    }
    else
    {
      for(i=0;i<3;i++)
	for(j=0;j<2;j++)
	  Range[i][j]=particles[0][j];
      for(i=1;i<np;i++)
	for(j=0;j<3;j++)
	{
	  auto x=particles[i][j];
	  if(x<Range[j][0])
	    Range[j][0]=x;
	  else if(x>Range[j][1])
	    Range[j][1]=x;
	}
      for(j=0;j<3;j++)
	Step[j]=(Range[j][1]-Range[j][0])/NDiv;
    }
    /*initialize hoc*/
    HBTInt *phoc=HOC;
    for(i=0;i<NDiv*NDiv*NDiv;i++,phoc++)
	    *phoc=-1;
	    
    for(i=0;i<np;i++)
    {
	    for(j=0;j<3;j++)
	    {
		    grid[j]=floor((particles[i][j]-Range[j][0])/Step[j]);
		    grid[j]=FixGridId(grid[j]);
	    }
	    ind=Sub2Ind(grid[0],grid[1],grid[2]);
	    List[i]=HOC[ind];
	    HOC[ind]=i; /*use hoc[ind] as swap varible to temporarily 
								    store last ll index, and finally the head*/
    }
  }
  void SearchSphere(HBTReal radius, const HBTxyz &searchcenter, vector <HBTInt> &found_ids, int nmax_guess=8)
  {//nmax_guess: initial guess for the max number of particles to be found. for memory allocation optimization purpose.
	PositionData_t &particles=*Particles;
	HBTReal dr;
	int i,j,k,subbox_grid[3][2];
	
	found_ids.clear();
	found_ids.reserve(nmax_guess);
		
	for(i=0;i<3;i++)
	{
	  subbox_grid[i][0]=floor((searchcenter[i]-radius-Range[i][0])/Step[i]);
	  subbox_grid[i][1]=floor((searchcenter[i]+radius-Range[i][0])/Step[i]);
	  if(!PeriodicBoundary)
	  {//do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
	    subbox_grid[i][0]=FixGridId(subbox_grid[i][0],ll);
	    subbox_grid[i][1]=FixGridId(subbox_grid[i][1],ll);
	  }	
	}
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
	  for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
	    for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
	    {
	      HBTInt pid=GetHOCSafe(i,j,k); //in case the grid-id is out of box, in the periodic case
	      while(pid>=0)
	      {
		dr=Distance(particles[pid],searchcenter);
		if(dr<radius)  found_ids.push_back(pid);
		pid=List[pid];
	      }
	    }
  }
};


