//to test the entry_ptr->size of hdf5

#include <iostream>
#include <string>
#include <vector>

using namespace std;
#include "../hdf_wrapper.h"

#define N 2

int main(void)
{

    cout<<"long size="<<sizeof(long)<<endl;

    long lens[]={2, 4386789};
    long totlen=lens[0]+lens[1];

    vector <long> data;
    data.resize(totlen);

    vector <hvl_t> vl(N);
    vl[0].len=lens[0];
    vl[0].p=data.data();
    vl[1].len=lens[1];
    vl[1].p=&data[lens[0]];

    hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_NATIVE_LONG);
    hid_t file=H5Fcreate("test_h5size.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dim[1];
    dim[0]=totlen;
    //this always work
    writeHDFmatrix(file, data.data(), "AllParticles", 1, dim, H5T_NATIVE_LONG);
    cout<<"single array written"<<endl<<flush;
    dim[0]=N;
    //this leads to error with hdf5-1.12.0 when size of vl entry is larger than 32MB (32*1k*1k bytes), but works in other versions including 1.12.1 and 1.8
    writeHDFmatrix(file, vl.data(), "SrchaloParticles", 1, dim, H5T_HBTIntArr);
    cout<<"vl written"<<endl;
    H5Fclose(file);
    H5Tclose(H5T_HBTIntArr);

    return 0;
}
