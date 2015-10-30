/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * This example shows how to create a compound datatype,
 * write an array which has the compound datatype to the file,
 * and read back fields' subsets.
 */

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
#include "hdf5.h"
#include "hdf5_hl.h"	
#include "H5Cpp.h"
// #include "H5VarLenType.h"	

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const H5std_string FILE_NAME( "SDScompound.hdf5" );
const H5std_string DATASET_NAME( "ArrayOfStructures" );
const H5std_string MEMBER1( "a" );
const H5std_string MEMBER2( "b" );
const H5std_string MEMBER3( "c" );
const int   LENGTH = 10;
const int   RANK = 1;

int main(void)
{
   /* First structure  and dataset*/
   typedef struct s1_t {
	int    a;
	float  b;
	double c;
   } s1_t;

   /* Second structure (subset of s1_t)  and dataset*/
   struct s2_t {
	double c;
	int    a;
   };
   
   // Try block to detect exceptions raised by any of the calls inside it
   try
   {
      /*
       * Initialize the data
       */
      int  i;
      s1_t s1[LENGTH];
      for (i = 0; i< LENGTH; i++)
      {
         s1[i].a = i;
         s1[i].b = i*i;
         s1[i].c = 1./(i+1);
      }

      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
//       Exception::dontPrint();

      /*
       * Create the data space.
       */
      hsize_t dim[] = {LENGTH};   /* Dataspace dimensions */
      DataSpace space( RANK, dim );

      /*
       * Create the file.
       */
      H5File file( FILE_NAME, H5F_ACC_TRUNC );

      /*
       * Create the memory datatype.
       */
      CompType mtype1( sizeof(s1_t) );
      mtype1.insertMember( MEMBER1, HOFFSET(s1_t, a), PredType::NATIVE_INT);
      mtype1.insertMember( MEMBER3, HOFFSET(s1_t, c), PredType::NATIVE_DOUBLE);
      mtype1.insertMember( MEMBER2, HOFFSET(s1_t, b), PredType::NATIVE_FLOAT); //only field exposed by CompType will get written

      /*
       * Create the dataset.
       */
      DataSet* dataset;
      dataset = new DataSet(file.createDataSet(DATASET_NAME, mtype1, space));

      /*
       * Write data to the dataset;
       */
      dataset->write( s1, mtype1 );

      /*
       * Release resources
       */
      delete dataset;

      /*
       * Open the file and the dataset.
       */
      H5File file2( FILE_NAME, H5F_ACC_RDONLY );
      dataset = new DataSet (file2.openDataSet( DATASET_NAME ));

      /*
       * Create a datatype for s2
       */
      CompType mtype2( sizeof(s2_t) );

      mtype2.insertMember( MEMBER3, HOFFSET(s2_t, c), PredType::NATIVE_DOUBLE);
      mtype2.insertMember( MEMBER1, HOFFSET(s2_t, a), PredType::NATIVE_INT);

      /*
       * Read two fields c and a from s1 dataset. Fields in the file
       * are found by their names "c_name" and "a_name".
       */
      s2_t s2[LENGTH];
      dataset->read( s2, mtype2 );

      /*
       * Display the fields
       */
#define ShowField(s,x){\
      cout << endl << "Field "<<#x<<" : " << endl;\
      for( i = 0; i < LENGTH; i++)\
	 cout << s[i].x << " ";\
      cout << endl;\
   }
   ShowField(s2, a);
   ShowField(s2, c);

      /*
       * Create a datatype for s3.
       */
      CompType mtypef( sizeof(float) );

      mtypef.insertMember( MEMBER2, 0, PredType::NATIVE_FLOAT);

      /*
       * Read field b from s1 dataset. Field in the file is found by its name.
       */
      float sf[LENGTH];  // Third "structure" - used to read float field of s1
      dataset->read( sf, mtypef );

      /*
       * Display the field
       */
      cout << endl << "Field b : " << endl;
      for( i = 0; i < LENGTH; i++)
	 cout << sf[i] << " ";
      cout << endl;

	  struct s3_t{
		float b;
		double c;
		int a;
		int d;
		float x[3];
		vector <int> i; //this invalidates POD requirement of s3_t
// 		float *f;
// 		s2_t *s2;
//  		CompType mtype2;
	  };//you just cannot be derived class??
	  struct ct_t{
		CompType ct;
		ct_t(): ct(sizeof(s3_t))
		{
		}
		size_t getsize()
		{
		  return ct.getSize();
		}
	  };
	  ct_t mtype3container;
	  cout<<mtype3container.getsize()<<endl;
	  VarLenType vlt_i(&PredType::NATIVE_INT);
	  cout<<"s3_t size="<<sizeof(s3_t)<<endl;
	  cout<<"vecotor i size="<<sizeof(vector <int>)<<endl;
	  vector <int> ii(10);
	  cout<<"vector ii size="<<sizeof(ii)<<endl;
      CompType mtype3(sizeof(s3_t));
      mtype3.insertMember("a", HOFFSET(s3_t, a), PredType::NATIVE_INT);
//       mtype3.insertMember(MEMBER2, HOFFSET(s3_t, b), PredType::NATIVE_FLOAT);
      mtype3.insertMember(MEMBER3, HOFFSET(s3_t, c), PredType::NATIVE_DOUBLE);
      mtype3.insertMember("d", HOFFSET(s3_t, d), PredType::NATIVE_INT);
	  hsize_t dims=3;
	  mtype3.insertMember("x", HOFFSET(s3_t, x), ArrayType(PredType::NATIVE_INT, 1, &dims));
// 	  mtype3.insertMember("i", HOFFSET(s3_t, i), vlt_i); //this is not useful, since i is the vector object, not the memory storage.
      
      vector <s3_t> s3(LENGTH);

	  s3[0].i.assign(3, 1);
	  s3[3].i.assign(2, -1);
      dataset->read(s3.data(), mtype3);
      
      cout << endl << "Field d : " << endl;
      for( i = 0; i < LENGTH; i++)
	 cout << s3[i].d << " ";
      cout << endl;

      cout << endl << "Field a : " << endl;
      for( i = 0; i < LENGTH; i++)
	 cout << s3[i].a << " ";
      cout << endl;
	  
	  cout << endl << "Field b : " << endl;
      for( i = 0; i < LENGTH; i++)
	 cout << s3[i].b << " ";
      cout << endl;
	  
	  cout << endl << "Field x : " << endl;
      for( i = 0; i < LENGTH; i++)
		{
		  cout <<i<<":";
		  for(auto x: s3[i].x)
			cout << x << " ";
		  cout << endl;
		}
		
	  cout << endl << "Field i : " << endl;
      for( i = 0; i < LENGTH; i++)
		if(s3[i].i.size())
		{
		  cout <<i<<":";
		  for(auto x: s3[i].i)
			cout << x << " ";
		  cout << endl;
		}
      
      /*
       * Release resources
       */
      delete dataset;
	  
	  {
	  H5File file3( "VL.hdf5", H5F_ACC_TRUNC );
      DataSet dset(file3.createDataSet("S3", mtype3, space));
	  dset.write( s3.data(), mtype3 );
	  //better avoid variable length array since it is requires newest h5py/pytables etc and may be different for Fortraners.
	  hvl_t vl[LENGTH];
	  for(i=0;i<LENGTH;i++)
	  {
		vl[i].p=s3[i].i.data();
		vl[i].len=s3[i].i.size();
	  }
	  file3.createGroup("/data");
	  DataSet dset2(file3.createDataSet("/data/S3Data", vlt_i, space));
	  dset2.write( vl, vlt_i);
// 	  dset2.setComment("comment","variable length");
	  H5LTset_attribute_string(file3.getId(),"/data/S3Data","Unit","(km/s)^2, -GM/R_physical");
	  }
	  {
		H5File file3("VL.hdf5", H5F_ACC_RDONLY);
		DataSet dset(file3.openDataSet("S3"));
		dset.read(s3.data(), mtype3);
		// now how to read S3Data back??...
		DataSet dset2(file3.openDataSet("/data/S3Data"));
		DataSpace dspace=dset2.getSpace();
		dspace.getSimpleExtentDims(dim);
		hsize_t count[]={1}, offset[]={0}, stride[]={1}, block[]={1};
		for(offset[0]=0;offset[0]<dim[0];offset[0]++)
		{
		  dspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
		  cout<<dset2.getVlenBufSize(vlt_i, dspace)/vlt_i.getSuper().getSize()<<" ";
		}
		cout<<endl;
		cout<<dset2.getStorageSize()<<','<<dset2.getInMemDataSize()<<endl;
		cout<<vlt_i.getSize()<<endl;
		cout<<dim[0]<<endl;
		hvl_t vl[dim[0]];
		dset2.read(vl, VarLenType(&PredType::NATIVE_FLOAT)); //this does not have to be the same as the storage type!
		for(i=0;i<LENGTH;i++)
		{
		  cout<<"data "<<i<<":";
		  for(int j=0;j<vl[i].len;j++)
			cout<<((float *)vl[i].p)[j]<<" ";
		  cout<<endl;
		}
	  }
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
      error.printError();
      return -1;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printError();
      return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printError();
      return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printError();
      return -1;
   }

   return 0;
}