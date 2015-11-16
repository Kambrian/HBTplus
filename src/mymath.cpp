#include <iostream>
#include <glob.h>

#include "mymath.h"

int count_pattern_files(char *filename_pattern)
{
   glob_t globbuf;

   globbuf.gl_offs = 0;
   glob(filename_pattern, GLOB_ERR, NULL, &globbuf);
//    printf("couting %s...\n found %d files\n", filename_pattern, globbuf.gl_pathc);
   globfree(&globbuf);
   return globbuf.gl_pathc;  
}
/*
std::ostream& operator << (std::ostream& o, HBTxyz &a)
{
   o << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
   return o;
};
*/
HBTInt compile_offset(HBTInt Len[], HBTInt Offset[], HBTInt n)
{//fill offset info, and return total length.
  HBTInt i,offset;
  for(i=0,offset=0;i<n;i++)
  {
    Offset[i]=offset;
    offset+=Len[i];
  }
  return offset;
}

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((unsigned) (x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
//bit shift operation is invalid for 8byte+ data

void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte)
/*This function is used to switch endian, for data2swap[nel] with element size mbyte*/
{
  size_t i,j;
  char *data, *old_data;//by definition, sizeof(char)=1, one byte
  
  data=(char *)data2swap;
  
  switch(mbyte)
  {
	  case 1 :break;
	  case 2 :
	  		for(j=0;j<nel;j++)
			FIX_SHORT(data[j*2]);
			break;
	  case 4 :
	  		for(j=0;j<nel;j++)
			FIX_LONG(data[j*4]);
			break;
	  default :
			old_data=new char[mbyte];
			for(j=0;j<nel;j++)
			{
			  memcpy(&old_data[0],&data[j*mbyte],mbyte);
			  for(i=0;i<mbyte;i++)
				{
				  data[j*mbyte+i]=old_data[mbyte-i-1];
				}
			}
			delete old_data;
	}
}

void spherical_basisD(double dx[3],double er[3],double et[3],double ef[3])
/*er: radial;
 *et: azimuthal,(0~2*pi)
 *ef: elevation (0~pi)
 * */ 
{
	int i;
	double dr,dxy2,modet,modef;//unit vectors of spherical coord.
	dxy2=dx[0]*dx[0]+dx[1]*dx[1];
	dr=sqrt(dxy2+dx[2]*dx[2]);
	for(i=0;i<3;i++) er[i]=dx[i]/dr;
	modet=sqrt(dxy2);
	if(0.==modet)
	{
		et[0]=1.;
		et[1]=et[2]=0.;
		ef[1]=1.;
		ef[0]=ef[2]=0.;
	}
	else
	{
		et[0]=-dx[1]/modet;
		et[1]=dx[0]/modet;
		et[2]=0;
		if(0.==dx[2])
		{
		ef[0]=ef[1]=0.;
		ef[2]=1.;	
		}
		else
		{
		modef=sqrt(dxy2+dxy2*dxy2/dx[2]/dx[2]);
		ef[0]=-dx[0]/modef;
		ef[1]=-dx[1]/modef;
		ef[2]=dxy2/dx[2]/modef;
		}
	}
}

HBTReal vec_prod(HBTReal *a,HBTReal *b,HBTInt dim)
{
	HBTInt i;
	HBTReal c;
	c=0.;
	for(i=0;i<dim;i++)
		c+=a[i]*b[i];
	return c;
}
void vec_cross(HBTxyz a,HBTxyz b,HBTxyz c)
{
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[0]*b[2]-a[2]*b[0];
	c[2]=a[0]*b[1]-a[1]*b[0];
}

/* From Numerical Recipes*/
HBTReal psort(HBTInt k, HBTInt numel, HBTReal arr[])
{/* partition sort: sort arr into two parts: k small elements and n-k large elements,
  * with arr[k-1] as the k-th smallest value , also the return value.  
  * slightly modified from the original version to have arr[0~n-1] rather than arr[1~n]*/
	HBTInt i,ir,j,l,mid;
	HBTReal a,temp;
	#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
	
	arr--;//so that it's accessed with arr[1~n]
	l=1;
	ir=numel;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
	#undef SWAP
}

int LargestRootFactor(int N, int dim)
/*find an integer factor of N that is the largest subject to x<=N**(1./dim)*/
{
  int x=floor(pow(N, 1./dim));
  for(;x>0;x--)
	if(N%x==0) return x;
}
vector <int> ClosestFactors(int N, int dim)
/* return a factorization of `N` into `dim` factors that are as close as possible to each other*/
{
  vector <int> factors;
  for(;dim>0;dim--)
  {
	int x=LargestRootFactor(N, dim);
	factors.push_back(x);
	N/=x;
  }
//   sort(factors.begin(), factors.end());
  return factors;
}
