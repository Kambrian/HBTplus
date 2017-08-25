/*to compute the density profile of each halo (or halo-matter correlation) out to RMAX in logrithmic bins. 
 * the bin edges are [0, r1, r2, ...rn) where (r1, ...rn) are generated as logspace(RMIN, RMAX, NBIN).
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

#define NBOUND_MIN 1000 //min number of host bound particles to consider.
#define RMIN 1e-2 //only for deciding binwidth; not used to set the innermost bin edge.
#define RMAX 20.
#define NBIN 20
// #define USE_LL //algorithm: whether to use linkedlist or geotree for spatial search.

class SnapshotOffset_t: public Snapshot_t
{
public:
  HBTInt Offset;
  HBTInt N;
  Snapshot_t & Snapshot;
  SnapshotOffset_t(HBTInt offset, Snapshot_t & fullsnapshot): Offset(offset), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot)
  {
    N=Snapshot.size()-Offset;
  }
  HBTInt ReSize(HBTInt n)
  {
    N=n;
  }
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetMemberId(const HBTInt i) const
  {
	return i+Offset;
  }
  HBTReal GetMass(const HBTInt i) const
  {
	return Snapshot.GetMass(GetMemberId(i));
  }
  const HBTxyz & GetPhysicalVelocity(const HBTInt i) const
  {
	return Snapshot.GetPhysicalVelocity(GetMemberId(i));
  }
  const HBTxyz & GetComovingPosition(const HBTInt i) const
  {
	return Snapshot.GetComovingPosition(GetMemberId(i));
  }
};

struct HaloSize_t
{
  HBTInt TrackId;
  HBTInt HaloId;
  HBTInt n[NBIN];
  HaloSize_t(){};
  void Compute(HBTxyz &cen, LinkedlistPara_t &ll);
  void Compute(HBTxyz &cen, GeoTree_t &tree);
};
void BuildHDFHaloSize(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk);

float DlnX=logf(RMAX/RMIN)/(NBIN-1)*2; //spacing in ln-space for r**2
void save(vector <HaloSize_t> &HaloSize, int isnap, int ifile=0, int nfiles=0);
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
    
  vector <HaloSize_t> HaloSize(subsnap.MemberTable.SubGroups.size()); 
  #pragma omp parallel for schedule(dynamic)
  for(HBTInt grpid=0;grpid<HaloSize.size();grpid++)
  {
    HaloSize[grpid].HaloId=grpid;//need to adjust this in MPI version..
    auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
    if(subgroup.size()==0||subsnap.Subhalos[subgroup[0]].Nbound<NBOUND_MIN) 
    {
      HaloSize[grpid].TrackId=-1;
      continue;
    }
    HaloSize[grpid].TrackId=subgroup[0];
    fill(begin(HaloSize[grpid].n), end(HaloSize[grpid].n), 0);
  }
  
  const int nloop=256; //do it in blocks to save memory
  HBTInt nmax=partsnap.size();
  HBTInt blocksize=nmax/nloop+1;
  for(HBTInt offset=0,iloop=0;offset<nmax;offset+=blocksize, iloop++)
  {
    HBTInt np=min(blocksize, nmax-offset);
    SnapshotOffset_t snapslice(offset, partsnap);//split the snapshot and do it piece by piece
    snapslice.ReSize(np);
    cout<<"i="<<iloop<<endl;
    
  #ifdef USE_LL
    SnapshotPos_t PartPos(snapslice);
    LinkedlistPara_t searcher(16, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    cout<<"linked list compiled\n";
  #else
    GeoTree_t searcher;
    searcher.Build(snapslice);
    cout<<"tree built\n";
  #endif
    
    #pragma omp parallel for schedule(dynamic,1)
    for(HBTInt grpid=0;grpid<HaloSize.size();grpid++)
    {
      if(HaloSize[grpid].TrackId<0) continue;
      auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
      HaloSize[grpid].Compute(subsnap.Subhalos[subgroup[0]].ComovingMostBoundPosition, searcher); 
    }
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
    
    save(HaloSize, isnap);
    
    return 0;
}

void HaloSize_t::Compute(HBTxyz &cen, LinkedlistPara_t &ll)
{
    vector <LocatedParticle_t> founds;
    founds.reserve(1024);
    ll.SearchSphereSerial(RMAX, cen, founds);
    for(auto &&p: founds)
    {
	int ibin=ceilf(logf(p.d2/RMIN/RMIN)/DlnX);
	if(ibin<0) ibin=0;
	else if(ibin>=NBIN) ibin=NBIN-1;
	n[ibin]++;
    }
}

void HaloSize_t::Compute(HBTxyz &cen, GeoTree_t &tree)
{
    vector <LocatedParticle_t> founds;
    founds.reserve(1024);
    tree.Search(cen, RMAX, founds);
    for(auto &&p: founds)
    {
	int ibin=ceilf(logf(p.d2/RMIN/RMIN)/DlnX);
	if(ibin<0) ibin=0;
	else if(ibin>=NBIN) ibin=NBIN-1;
	n[ibin]++;
    }
}

void save(vector <HaloSize_t> &HaloSize, int isnap, int ifile, int nfiles)
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
  vector <float> rbin;
  logspace(RMIN, RMAX, NBIN, rbin);
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