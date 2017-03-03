// task4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "data.h"

void main()
{
	data mydata;
	//mydata.print_data(cout);
	ifstream f1("E:\\talk\\biology\\task4\\chi1.txt",ios::in);
	ifstream f2("E:\\talk\\biology\\task4\\alphabet.txt",ios::in);
	mydata.readfile(f1,f2);
	mydata.print_data(cout);
}

/*void printa(int **dataa,int length){
	for (int i=0;i<1296;i++)
	{
		for (int j=0;j<length;j++)
		{
			cout<<*((int *)dataa+i*length+j)<<endl;
		}
	}
}*/

void main0(){

	cout<<true;

	//int mydata[1296][20];
	/*vector<vector<double>> mydata;
	mydata.resize(1296);
	for (int i=0;i<1296;i++)
	{
		mydata[i].resize(20);
		for (int j=0;j<20;j++)
		{
			mydata[i][j] = 1;
		}
	}
	for (int i=0;i<1296;i++)
	{
		for (int j=0;j<20;j++)
		{
			++mydata[i][j];
			cout<<mydata[i][j]<<endl;
		}
	}*/
}