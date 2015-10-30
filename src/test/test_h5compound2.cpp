#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>
#include <vector>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
using std::vector;
#endif  // H5_NO_STD
#endif
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

const H5std_string FILE_NAME( "test_compound2.hdf5" );
const H5std_string DATASET_NAME( "data" );
const int   LENGTH = 5;
const int   RANK = 1;

#define ShowField(s,f){\
	cout << endl<<"Field "<<#f<<" : " << endl; \
	for(int i = 0; i < LENGTH; i++)\
	  cout<<s[i].f<<" ";\
	  cout<<endl;\
	}
	
int main(void)
{

    struct s_t
    {
        int    a;
        float  b;
        int c;
    };
    CompType mtype( sizeof(s_t) );
    /*only insert a,b, do not insert c*/
    mtype.insertMember( "a", HOFFSET(s_t, a), PredType::NATIVE_INT);
    mtype.insertMember( "b", HOFFSET(s_t, b), PredType::NATIVE_FLOAT);
	/*note field c is not inserted!*/
	CompType mtype2;
	mtype2.copy(mtype);
	mtype2.pack();
	cout<<mtype.getSize()<<","<<mtype2.getSize()<<endl;

    hsize_t dim[] = {LENGTH};
    vector <s_t> datain(LENGTH);
    for(int i=0; i<LENGTH; i++)/* init data*/
    {
        datain[i].a=i;
        datain[i].b=i*i;
        datain[i].c=-i;
    }
    cout<<"==========Data initialized=============\n";
	ShowField(datain, a);
	ShowField(datain, b);
	ShowField(datain, c);
	
    /*write to file*/
    {
        DataSpace space( RANK, dim );
        H5File file( FILE_NAME, H5F_ACC_TRUNC );
        DataSet dset(file.createDataSet(DATASET_NAME, mtype2, space)); //remove padding here!
        dset.write( datain.data(), mtype );
    }

    /*read back*/
    H5File file( FILE_NAME, H5F_ACC_RDONLY );
    DataSet dset(file.openDataSet( DATASET_NAME ));
    vector <s_t> dataout(LENGTH);

	dset.read( dataout.data(), mtype );

    cout<<"\n===========Data Read==========\n";
    ShowField(dataout,a);
    ShowField(dataout,b);
    ShowField(dataout,c); 
		
	CompType btype(sizeof(float));
	btype.insertMember("b", 0, PredType::NATIVE_FLOAT);
	float b[LENGTH];
	dset.read(b, btype);
	cout<<"b \n";
	for(int i=0;i<LENGTH;i++) cout<<b[i]<<" ";
	cout<<endl;
/*Output:
==========Data initialized=============

Field a : 
0 1 2 3 4 

Field b : 
0 1 4 9 16 

Field c : 
0 -1 -2 -3 -4 

===========Data Read==========

Field a : 
0 1 2 3 4 

Field b : 
0 1 4 9 16 

Field c : 
0 -1 -2 -3 -4 
*/
    return 0;
}
