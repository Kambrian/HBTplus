#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define BIGENDIAN

#ifdef BIGENDIAN
#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | (((x) >> 8) & 0x0000ff00) | ((unsigned) (x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))

void swap_Nbyte(char *data,int n,int m)
/*This function is used to switch endian*/
{
  int i,j;
  char old_data[16];
  
  switch(m)
  {
	  case 1 :break;
	  case 2 :
	  		for(j=0;j<n;j++)
			FIX_SHORT(data[j*2]);
			break;
	  case 4 :
	  		for(j=0;j<n;j++)
			FIX_LONG(data[j*4]);
			break;
	  default :
			for(j=0;j<n;j++)
			{
			  memcpy(&old_data[0],&data[j*m],m);
			  for(i=0;i<m;i++)
				{
				  data[j*m+i]=old_data[m-i-1];
				}
			}
	}
}
size_t fread_BE(void *buf,size_t Nsize,size_t Nbuf,FILE *fp)
{
	size_t Nread;
	Nread=fread(buf,Nsize,Nbuf,fp);
	swap_Nbyte((char *)buf,Nbuf,Nsize);
	return Nread;
}
#else
#define fread_BE fread
#endif

typedef struct
{
   int Np;   //number of particles in the simulation
   int ips;  //time step number, or snapshot num
   float ztp; // current redshift
   float Omegat;  // current omega_m, mass density
   float Lambdat; // current Omega_lambda, dark energy density
   float rLbox;  // boxsize in unit of Mpc/h
   float xscale;
   float vscale;
   float Hz; //Hubble Param at ztp
   float vunit; //velocity unit to get physical peculiar velocity in km/s
   /*==extension needed for calculating binding energy:==*/
   float time;//current reduced scale factor 
   float mass[2];//Gas (mass[0]) and DM (mass[1]) particle masses, in units of 10^10Msun/h
}IO_HEADER; //header of jing's data structure
	
	#define G 43007.1
	#define HUBBLE0 0.1    //H_0 in internal units
	
int main(int argc,char **argv)
{	
	FILE *fp;
	int i,j,dummy,dummy2;
	IO_HEADER header;
	float PMass;
	char snapfile[1024];//="/data/A4700r3d1/ypjing/pos6121.5000";
	
	sprintf(snapfile,"%s",argv[1]);
	printf("reading %s\n",snapfile);
	#ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread_BE(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread_BE(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){printf("error!record brackets not match for file %s\t%d,%d\n",snapfile,dummy,dummy2);fflush(stdout);exit(1);} 

	if(!(fp=fopen(snapfile,"r")))
	{
		printf("Error:cannot open file %s\n",snapfile);
		exit(1);
	}
	SKIP;printf("%d\t",dummy);
	fread_BE(&header.Np,sizeof(int),1,fp);
	fread_BE(&header.ips,sizeof(int),1,fp);
	fread_BE(&header.ztp,sizeof(float),1,fp);
	fread_BE(&header.Omegat,sizeof(float),1,fp);
	fread_BE(&header.Lambdat,sizeof(float),1,fp);
	fread_BE(&header.rLbox,sizeof(float),1,fp);
	fread_BE(&header.xscale,sizeof(float),1,fp);
	fread_BE(&header.vscale,sizeof(float),1,fp);
	SKIP2;printf("%d\n",dummy2);
	CHECK;
	PMass = header.Omegat*3.*HUBBLE0*HUBBLE0/8./3.1415926/G*header.rLbox*header.rLbox*header.rLbox*1.0e9/header.Np;
	//particle mass in units of 10^10Msun/h; only valid when using params at z=0;
	printf("Np=%d,ips=%d,Omegat=%g,Lambdat=%g,Lbox=%g\nztp=%g,xscale=%g,vscale=%g\nPartMass=%g (only valid for ztp=0!)\n",
		header.Np,header.ips,header.Omegat,header.Lambdat,header.rLbox,header.ztp,header.xscale,header.vscale,PMass);
	SKIP;printf("%d\t",dummy);
	fseek(fp,-4L,SEEK_END);
	SKIP2;
	CHECK;printf("%d\n",dummy2);
	return 0;
}
