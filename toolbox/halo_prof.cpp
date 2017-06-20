/*to compute the density profile of each halo (or halo-matter correlation) out to RMAX in logrithmic bins. 
 * the bins are generated as logspace(RMIN, RMAX, NBIN+1), with the innermost bin edge replaced by 0.
 * output the count in each bin for each halo.
 */

#include <cmath>
#include <iostream>
#include <string>

#include "../src/datatypes.h"
#include "../src/config_parser.h"
#include "../src/snapshot.h"
#include "../src/halo.h"
#include "../src/subhalo.h"
#include "../src/mymath.h"
#include "../src/linkedlist_parallel.h"
#include "../src/geometric_tree.h"

#define RMIN 1e-2 //only for deciding binwidth; not used to set the innermost bin edge.
#define RMAX 20.
#define NBIN 20
//#define USE_LL //algorithm: whether to use linkedlist or geotree for spatial search.

struct HaloSize_t
{
  HBTInt TrackId;
  HBTInt HaloId;
  HBTInt n[NBIN];
  HaloSize_t(): TrackId(-1), HaloId(-1), n{0}
  {
  }
  void Compute(HBTxyz &cen, LinkedlistPara_t &ll);
  void Compute(HBTxyz &cen, GeoTree_t &tree);
};
void BuildHDFHaloSize(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk);

vector <float> rbin, rbin2;
void logspace(double xmin,double xmax,int N, vector <float> &x);
void save(vector <HaloSize_t> HaloSize, int isnap, int ifile, int nfiles);
int main(int argc, char **argv)
{
  if(argc!=3)
  {
    cerr<<"Usage: "<<endl;
    cerr<<" "<<argv[0]<<" [config_file] [snapshot_number]"<<endl;
    cerr<<"    If snapshot_number<0, then it's counted from final snapshot in reverse order"<<endl; 
    cerr<<"    (i.e., FinalSnapshot=-1,... FirstSnapshot=-N)"<<endl;
    return 1;
  }
  HBTConfig.ParseConfigFile(argv[1]);
  int isnap=atoi(argv[2]);
  if(isnap<0) isnap=HBTConfig.MaxSnapshotIndex+isnap+1;
  
  SubhaloSnapshot_t subsnap(isnap, SubReaderDepth_t::SubTable);;
  ParticleSnapshot_t partsnap(isnap);
  auto &Cosmology=partsnap.Cosmology;
 
#ifdef USE_LL
  SnapshotPos_t PartPos(partsnap);
  LinkedlistPara_t searcher(16, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);//memory consumption by ll is much lighter, even with 256**3 grids and 24 threads (abount 3GB of HoC, plus np ints)
  cout<<"linked list compiled\n";
#else
  GeoTree_t searcher;//memory consumption can be heavy, HBTInt[8]*Np, bigger than snapshot!
  searcher.Build(partsnap);
  cout<<"tree built\n";
#endif
  
  logspace(RMIN, RMAX, NBIN+1, rbin);
  rbin2.resize(rbin.size());
  for(int i=0;i<rbin.size();i++)
    rbin2[i]=rbin[i]*rbin[i];
  
  
  vector <HaloSize_t> HaloSize; 
  const int nfiles=1; //do it in blocks to save memory
  HBTInt nmax=subsnap.MemberTable.SubGroups.size();
  HBTInt blocksize=nmax/nfiles+1;
  for(int ifile=0;ifile<nfiles;ifile++)
  {
    HBTInt grpidmin=blocksize*ifile;
    HBTInt grpidmax=min(grpidmin+blocksize, nmax);
    HaloSize.resize(blocksize);
    #pragma omp parallel for schedule(dynamic,1)
    for(HBTInt grpid=grpidmin;grpid<grpidmax;grpid++)
    {
      HaloSize[grpid].HaloId=grpid;//need to adjust this in MPI version..
      auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
      if(subgroup.size()==0||subsnap.Subhalos[subgroup[0]].Nbound<1000) 
	continue;
      HaloSize[grpid].TrackId=subgroup[0];
      HaloSize[grpid].Compute(subsnap.Subhalos[subgroup[0]].ComovingMostBoundPosition, searcher); 
    }
    auto it=HaloSize.begin(), it_save=HaloSize.begin();
    for(;it!=HaloSize.end();++it)
    {
      if(it->TrackId>=0)
      {
	if(it!=it_save)
	  *it_save=move(*it);
	++it_save;
      }
    }
    HaloSize.resize(it_save-HaloSize.begin());
    
    save(HaloSize, isnap, ifile, nfiles);
  }

}

void HaloSize_t::Compute(HBTxyz &cen, LinkedlistPara_t &ll)
{
    vector <LocatedParticle_t> founds;
    founds.reserve(1024);
    ll.SearchSphereSerial(RMAX, cen, founds);
    for(auto &&p: founds)
    {
      for(int i=NBIN-1;i>=0;i--)
      {
	if(p.d>rbin2[i])
	{
	  n[i]++;
	  break;
	}
      }
    }
}

void HaloSize_t::Compute(HBTxyz &cen, GeoTree_t &tree)
{
    vector <LocatedParticle_t> founds;
    founds.reserve(1024);
    tree.Search(cen, RMAX, founds);
    for(auto &&p: founds)
    {
      for(int i=NBIN-1;i>=0;i--)
      {
	if(p.d>rbin2[i])
	{
	  n[i]++;
	  break;
	}
      }
    }
}

void save(vector <HaloSize_t> HaloSize, int isnap, int ifile, int nfiles)
{
  string filename=HBTConfig.SubhaloPath+"/HaloSize";
  mkdir(filename.c_str(), 0755);
  if(ifile==0&&nfiles==0)
    filename=filename+"/HaloProf_"+to_string(isnap)+".hdf5";
  else
    filename=filename+"/HaloProf_"+to_string(isnap)+"."+to_string(ifile)+".hdf5";
  hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t dim_atom[]={1}, dims[]={HaloSize.size()};
  hid_t H5T_HaloSizeInMem, H5T_HaloSizeInDisk;
  BuildHDFHaloSize(H5T_HaloSizeInMem, H5T_HaloSizeInDisk);
  writeHDFmatrix(file, &nfiles, "NumberOfFiles", 1, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, HaloSize.data(), "HostHalos", 1, dims, H5T_HaloSizeInMem, H5T_HaloSizeInDisk);
  dims[0]=rbin.size();
  writeHDFmatrix(file, rbin.data(), "RadialBins", 1, dims, H5T_NATIVE_FLOAT);
  H5Tclose(H5T_HaloSizeInDisk);
  H5Tclose(H5T_HaloSizeInMem);
  H5Fclose(file);
}

void BuildHDFHaloSize(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk)
{
  H5T_dtypeInMem=H5Tcreate(H5T_COMPOUND, sizeof (HaloSize_t));
  hsize_t dims[]={NBIN};
  hid_t H5T_HBTIntArray=H5Tarray_create2(H5T_HBTInt, 1, dims);

  #define InsertMember(x,t) H5Tinsert(H5T_dtypeInMem, #x, HOFFSET(HaloSize_t, x), t)//;cout<<#x<<": "<<HOFFSET(HaloSize_t, x)<<endl
  InsertMember(TrackId, H5T_HBTInt);
  InsertMember(HaloId, H5T_HBTInt);
  InsertMember(n, H5T_HBTIntArray);
  #undef InsertMember	

  H5T_dtypeInDisk=H5Tcopy(H5T_dtypeInMem);
  H5Tpack(H5T_dtypeInDisk); //clear fields not added to save disk space

  H5Tclose(H5T_HBTIntArray);
}
void logspace(double xmin,double xmax,int N, vector <float> &x)
{
  x.resize(N);
  int i;
  double dx;
  x[0]=xmin;x[N-1]=xmax;
  xmin=log(xmin);
  xmax=log(xmax);
  dx=exp((xmax-xmin)/(N-1));
  for(i=1;i<N-1;i++)
  {
    x[i]=x[i-1]*dx;
  }
}