//outputs subhalo lists 
#include <cmath>
#include <iostream>
#include <string>

#include "../../src/datatypes.h"
#include "../../src/config_parser.h"
#include "../../src/snapshot.h"
#include "../../src/halo.h"
#include "../../src/subhalo.h"
#include "../../src/mymath.h"
#include "../../src/linkedlist_parallel.h"

#define RMAX 1.5	 //statistics done in r<RMAX*rvir

#define MHOST M200Crit
#define RHOST R200CritComoving

struct Satellite_t
{
    float DBound2HostBound;
    float DCoM2HostCoM;
    float DBound2HostCoM;
    float VBound2HostBound;
    float VCoM2HostCoM;
    float VBound2HostCoM;
    float DBoundCoMHost;
    float VBoundCoMHost;
    float MHost;
    float RHost;
    float Mbound, LastMaxMass, LastMaxVmax;
    float Vmax, VmaxHost;
    HBTInt HostId;
    HBTInt TrackId;
    HBTInt CentralTrackId;
};

class SubhaloPos_t: public PositionData_t
{
    vector <Subhalo_t> &Subhalos;
public:
    SubhaloPos_t(vector <Subhalo_t> &subhalos):Subhalos(subhalos)
    {}
    const HBTxyz & operator [](HBTInt i) const
    {   return Subhalos[i].ComovingMostBoundPosition;
    }
    size_t size() const
    {   return Subhalos.size();
    }
};

struct HaloSize_t
{
    HBTInt HaloId;
    HBTxyz CenterComoving;
    float MVir, M200Crit, M200Mean, RVirComoving, R200CritComoving, R200MeanComoving;
    float RmaxComoving, VmaxPhysical;
};
vector <HaloSize_t> HaloSize;
void LoadHaloSize(vector<HaloSize_t> &HaloSize, int isnap);
void save(vector <Satellite_t> &Satellites, int isnap);
void collect_submass(int grpid, const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll, vector <Satellite_t> &satellites);

int main(int argc,char **argv)
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
    SubhaloSnapshot_t subsnap(isnap, SubReaderDepth_t::SubTable);
    cout<<"data loaded\n";

    LoadHaloSize(HaloSize,isnap);
    if(HaloSize.size()!=subsnap.MemberTable.SubGroups.size())
    {
        cerr<<HaloSize.size()<<"!="<<subsnap.MemberTable.SubGroups.size()<<endl;
        throw(runtime_error("number of halos mismatch\n"));
    }

    SubhaloPos_t SubPos(subsnap.Subhalos);
    LinkedlistPara_t ll(200, &SubPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    cout<<"linked list compiled\n";

    vector <Satellite_t> Satellites;
    Satellites.reserve(subsnap.Subhalos.size());

    auto &subgroups=subsnap.MemberTable.SubGroups;
    #pragma omp parallel for
    for(HBTInt grpid=0; grpid<subgroups.size(); grpid++)
    {
        if(subgroups[grpid].size()==0) continue;
        vector <Satellite_t> satlist;
        collect_submass(grpid, subsnap, ll, satlist);
        #pragma omp critical
        Satellites.insert(Satellites.end(), satlist.begin(), satlist.end());
    }

    save(Satellites, isnap);

    return 0;
}

void collect_submass(int grpid, const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll, vector <Satellite_t> &satellites)
{
    satellites.clear();
    auto &host=HaloSize[grpid];
    float rhost=host.RHOST;
    float mhost=host.MHOST;
    if(mhost<=0) return;
    satellites.reserve(subsnap.MemberTable.SubGroups[grpid].size());

    HBTInt cenid=subsnap.MemberTable.SubGroups[grpid][0];
    auto &central=subsnap.Subhalos[cenid];
    vector <LocatedParticle_t> sublist;
    ll.SearchSphereSerial(RMAX*rhost, central.ComovingMostBoundPosition, sublist, 8);
    float d0=PeriodicDistance(central.ComovingMostBoundPosition, central.ComovingAveragePosition);
    float v0=Distance(central.PhysicalMostBoundVelocity, central.PhysicalAverageVelocity);
    for(auto &&subid: sublist)
    {
//         if(subid.id!=cenid)
        {
            auto &sub=subsnap.Subhalos[subid.id];
            Satellite_t sat;
	    sat.CentralTrackId=central.TrackId;
	    sat.TrackId=sub.TrackId;
            sat.Mbound=sub.Mbound;
            sat.DBound2HostBound=PeriodicDistance(sub.ComovingMostBoundPosition, central.ComovingMostBoundPosition);
            sat.DBound2HostCoM=PeriodicDistance(sub.ComovingMostBoundPosition, central.ComovingAveragePosition);
            sat.DCoM2HostCoM=PeriodicDistance(sub.ComovingAveragePosition, central.ComovingAveragePosition);
            sat.VBound2HostBound=Distance(sub.PhysicalMostBoundVelocity, central.PhysicalMostBoundVelocity);
            sat.VBound2HostCoM=Distance(sub.PhysicalMostBoundVelocity, central.PhysicalAverageVelocity);
            sat.VCoM2HostCoM=Distance(sub.PhysicalAverageVelocity, central.PhysicalAverageVelocity);
            sat.MHost=mhost;
            sat.RHost=rhost;
            sat.DBoundCoMHost=d0;
            sat.VBoundCoMHost=v0;
            sat.HostId=grpid;
	    sat.LastMaxMass=sub.LastMaxMass;
	    sat.LastMaxVmax=sub.LastMaxVmaxPhysical;
	    sat.Vmax=sub.VmaxPhysical;
// 	    sat.VmaxHost=host.VmaxPhysical;
	    sat.VmaxHost=central.VmaxPhysical;
            satellites.push_back(sat);
        }
    }
}

void BuildHDFSatellite(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk)
{
    H5T_dtypeInMem=H5Tcreate(H5T_COMPOUND, sizeof (Satellite_t));
    hsize_t dims[2]= {3,3};
    hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);

#define InsertMember(x,t) H5Tinsert(H5T_dtypeInMem, #x, HOFFSET(Satellite_t, x), t)
    InsertMember(DBound2HostBound, H5T_NATIVE_FLOAT);
    InsertMember(DCoM2HostCoM, H5T_NATIVE_FLOAT);
    InsertMember(DBound2HostCoM, H5T_NATIVE_FLOAT);
    InsertMember(VBound2HostBound, H5T_NATIVE_FLOAT);
    InsertMember(VCoM2HostCoM, H5T_NATIVE_FLOAT);
    InsertMember(VBound2HostCoM, H5T_NATIVE_FLOAT);
    InsertMember(DBoundCoMHost, H5T_NATIVE_FLOAT);
    InsertMember(VBoundCoMHost, H5T_NATIVE_FLOAT);
    InsertMember(MHost, H5T_NATIVE_FLOAT);
    InsertMember(RHost, H5T_NATIVE_FLOAT);
    InsertMember(Mbound, H5T_NATIVE_FLOAT);
    InsertMember(HostId, H5T_HBTInt);
    InsertMember(TrackId, H5T_HBTInt);
    InsertMember(CentralTrackId, H5T_HBTInt);
    InsertMember(LastMaxMass, H5T_NATIVE_FLOAT);
    InsertMember(LastMaxVmax, H5T_NATIVE_FLOAT);
    InsertMember(Vmax, H5T_NATIVE_FLOAT);
    InsertMember(VmaxHost, H5T_NATIVE_FLOAT);
#undef InsertMember

    H5T_dtypeInDisk=H5Tcopy(H5T_dtypeInMem);
    H5Tpack(H5T_dtypeInDisk); //clear fields not added to save disk space

    H5Tclose(H5T_HBTxyz);
}

void save(vector <Satellite_t> &Satellites, int isnap)
{
    string outdir=HBTConfig.SubhaloPath+"/analysis/";
    mkdir(outdir.c_str(),0755);
#define xstr(s) str(s)
#define str(s) #s
    string filename=outdir+"massList"+to_string(isnap)+"." xstr(MHOST) ".hdf5";
#undef xstr
#undef str
    hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dtype, dtype2;
    BuildHDFSatellite(dtype, dtype2);
    hsize_t dims[]= {Satellites.size()};
    writeHDFmatrix(file, Satellites.data(), "Satellites", 1, dims, dtype, dtype2);
    H5Tclose(dtype);
    H5Tclose(dtype2);
    H5Fclose(file);
}

void BuildHDFHaloSize(hid_t &H5T_dtype)
{
    H5T_dtype=H5Tcreate(H5T_COMPOUND, sizeof (HaloSize_t));
    hsize_t dims[2]= {3,3};
    hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);

#define InsertMember(x,t) H5Tinsert(H5T_dtype, #x, HOFFSET(HaloSize_t, x), t)//;cout<<#x<<": "<<HOFFSET(HaloSize_t, x)<<endl
    InsertMember(HaloId, H5T_HBTInt);
    InsertMember(CenterComoving, H5T_HBTxyz);
    InsertMember(R200CritComoving, H5T_NATIVE_FLOAT);
    InsertMember(R200MeanComoving, H5T_NATIVE_FLOAT);
    InsertMember(RVirComoving, H5T_NATIVE_FLOAT);
    InsertMember(M200Crit, H5T_NATIVE_FLOAT);
    InsertMember(M200Mean, H5T_NATIVE_FLOAT);
    InsertMember(MVir, H5T_NATIVE_FLOAT);
    InsertMember(RmaxComoving, H5T_NATIVE_FLOAT);
    InsertMember(VmaxPhysical, H5T_NATIVE_FLOAT);
#undef InsertMember

    H5Tclose(H5T_HBTxyz);
}

void LoadHaloSize(vector <HaloSize_t> &HaloSize, int isnap)
{
    string filename=HBTConfig.SubhaloPath+"/HaloSize/HaloSize_"+to_string(isnap)+".hdf5";
    hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t H5T_HaloSize;
    BuildHDFHaloSize(H5T_HaloSize);

    hid_t dset=H5Dopen2(file, "HostHalos", H5P_DEFAULT);
    hsize_t dims[1];
    GetDatasetDims(dset, dims);
    HaloSize.resize(dims[0]);
    if(dims[0])	H5Dread(dset, H5T_HaloSize, H5S_ALL, H5S_ALL, H5P_DEFAULT, HaloSize.data());
    H5Dclose(dset);

    H5Tclose(H5T_HaloSize);
    H5Fclose(file);
}