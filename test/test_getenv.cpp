using namespace std;
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>

int main(int argc, char **argv)
{
 
  string dir=".";
  char * versionstr=getenv("USER");
  if(versionstr) 
  {
	cout<<dir+"/VER"+versionstr;
	ofstream version_file(dir+"/VER"+versionstr, fstream::trunc);
// 	version_file<<versionstr;
// 	version_file.close();
  }
  else
  {
	cout<<"HBT_VERSION not specified"<<endl;
  }
  
  return 0;
}
