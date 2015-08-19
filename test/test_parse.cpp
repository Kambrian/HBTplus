#include "../config_parser.h"

int main()
{
	Parameter_t config;
	config.ParseConfigFile("test.conf");
	cout<<config.SnapshotPath<<endl;
	cout<<config.SubhaloPath<<endl;
	cout<<config.MinNumPartOfSub<<endl;
}