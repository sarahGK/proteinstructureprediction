#include "stdafx.h"
#include "data.h"

data::data(){
	sd.resize(LA);
	se.resize(LA);
	for (int i=0;i<LA;i++)
	{
		for (int j=0;j<20;j++)
		{
			sa[i][j] = 1;
			sb[i][j] = 1;
			sc[i][j] = 1;
		}
		sd[i].resize(LA);
		se[i].resize(LA);
		for (int j=0;j<LA;j++)
		{
			sd[i][j] = 1;
			se[i][j] = 1;
		}
		for (int j=0;j<3;j++)
		{
			sf[i][j] = 1;
			sg[i][j] = 1;
			sh[i][j] = 1;
		}
	}
	for(int i=0;i<8;i++)
		mark[i] = false;
}

void data::set_data(int s, int a, int b, int c, int d, int e, int f, int g, int h){
	++sa[s][a];
	++sb[s][b];
	++sc[s][c];
	++sd[s][d];
	++se[s][e];
	++sf[s][f];
	++sg[s][g];
	++sh[s][h];

}

void data::set_data(double **dataarray,int length,int no){
	double sum = 0.0;
	double sums[LA];//
	double sumac[LA];
	//double joint[LA][LA];
	double pre,post;
	pre= post = 0.0;
	for (int i=0;i<LA;i++)
	{
		sums[i]=0.0;
		for (int j=0;j<length;j++)
		{
			sums[i] = sums[i] + *((double *)dataarray + i*length +j);
		}
		sum = sum + sums[i];
	}
	for (int i=0;i<length;i++)
	{	
		
		sumac[i] = 0.0;
		for (int j=0;j<LA;j++)
		{
			sumac[i] = sumac[i] + *((double *)dataarray + i*length +j);
			*((double *)dataarray + j*length +i) = *((double *)dataarray + j*length +i) / sum;
		}
	}
	for (int i=0;i<LA;i++)
	{
		sums[i] = sums[i]/sum;
		pre = pre -sums[i] * log(sums[i]);
		for (int j=0;j<length;j++)
		{
			//*(*(dataarray+i)+j) = joint[i][j] /sumac[j] ;
			post = post-*((double *)dataarray + i*length +j)*log(*((double *)dataarray + i*length +j)/sumac[j]);
		}
	}

	if(post>pre) mark[no] = true;
	
}

void data::set_data(std::vector<vector<double> > &de,int no){

	double sum = 0.0;
	double sums[LA];
	double sumac[LA];
	//double joint[LA][LA];
	double pre,post;
	pre= post = 0.0;
	for (int i=0;i<LA;i++)
	{
		sums[i]=0.0;
		for (int j=0;j<LA;j++)
		{
			sums[i] = sums[i] + de[i][j];
		}
		sum = sum + sums[i];
	}
	for (int i=0;i<LA;i++)
	{	

		sumac[i] = 0.0;
		for (int j=0;j<LA;j++)
		{
			sumac[i] = sumac[i] + de[i][j];
			de[j][i] = de[j][i] / sum;
		}
	}
	for (int i=0;i<LA;i++)
	{
		sums[i] = sums[i]/sum;
		pre = pre -sums[i] * log(sums[i]);
		for (int j=0;j<LA;j++)
		{
			//de[i][j] = joint[i][j] /sumac[j] ;
			post = post-de[i][j]*log(de[i][j]/sumac[j]);
		}
	}

	if(post>pre) mark[no] = true;
}

/*void data::set_data(int s, int f, int g, int h){
	++sf[s][f];
	++sg[s][g];
	++sh[s][h];
}*/

/*void data::readin(ifstream &f1,ifstream &f2){

}*/

void data::print_data(ostream & of){
	/*for (int i=0;i<LA;i++)
	{
		for (int j=0;j<20;j++)
		{
			of<<sa[i][j] ;
			of<<sb[i][j] ;
			of<<sc[i][j] ;
		}

		for (int j=0;j<LA;j++)
		{
			of<<sd[i][j];
			of<<se[i][j];
		}
		for (int j=0;j<3;j++)
		{
			of<<sf[i][j] ;
			of<<sg[i][j] ;
			of<<sh[i][j] ;
		}
	}
	of<<endl;
	of<<endl;*/
	for (int i=0;i<8;i++)
	{
		of<<mark[i]<<endl;
	}
}

int BinarySearchChar(char c)
{
	int left,right,middle;
	left=0;
	right=20;
	if(c<rname_one_rank20[0])
		return -1;
	if(c>rname_one_rank20[19])
		return -1;
	while(right-left>1)
	{
		middle=(left+right)/2;
		if(c==rname_one_rank20[middle])
			return middle;
		if(c<rname_one_rank20[middle])
			right=middle;
		else
			left=middle;
	}
	if(c==rname_one_rank20[left])
		return left;
	else if(c==rname_one_rank20[right])
		return right;
	else 
		return -1;
}

void data::readfile(ifstream & file1,ifstream & file2){
	bool flag1,flag2;
	flag1 = flag2 = true;
	char af,am,ab;
	int pf,pm,pb,chif,chim,chib,sf1,sm,sb1;
	int i=0;
	while (!file1.eof())
	{
		i++;
	
		string line1,line2 ;
		getline(file1,line1);
		getline(file2,line2);
		if(line1.length()<=10){
			for(;!file2.eof()&&strcmp(line2.c_str(),line1.c_str())!=0;getline(file2,line2));			
			continue;
		}
		if(flag1){
			af = line1.substr(9,1).c_str()[0];
			pf = BinarySearchChar(af);
			chif = atoi(line1.substr(25,4).c_str());
			sf1 = atoi(line2.substr(4,4).c_str());
		}
		if(-1==pf||-999==chif||-999==sf1){
			flag1 = flag2=true;
			continue;
		}
		if(flag1){
		flag1 = false;
		getline(file1,line1);
		getline(file2,line2);
		if(line1.length()<=10){
			for(;!file2.eof()&&strcmp(line2.c_str(),line1.c_str())!=0;getline(file2,line2));			
			continue;
		}
		am = line1.substr(9,1).c_str()[0];
		pm = BinarySearchChar(am);
		chim = atoi(line1.substr(25,4).c_str());
		sm = atoi(line2.substr(4,4).c_str());
		}
		if(-1==pm||-999==chim||-999==sm){
			flag1 =flag2= true;
			continue;
		}
		if(flag2){
		flag2 = false;
		getline(file1,line1);
		getline(file2,line2);
		if(line1.length()<=10){
			for(;!file2.eof()&&strcmp(line2.c_str(),line1.c_str())!=0;getline(file2,line2));			
			continue;
		}
		}
		ab = line1.substr(9,1).c_str()[0];
		pb = BinarySearchChar(ab);
		chib = atoi(line1.substr(25,4).c_str());
		sb1 = atoi(line2.substr(4,4).c_str());		
		if(-1==pb||-999==chib||-999==sb1){
			flag1=flag2=true;
			continue;
		}
		set_data(sm,pf,pm,pb,sf1,sb1,chif,chim,chib);
		sf1=sm;
		pf=pm;
		chif=chim;
		sm=sb1;
		pm=pb;
		chim=chib;
		cout<<"now we are reading data now"<<i<<endl;
	}

	set_data((double **)sa,20,0);
	set_data((double **)sb,20,1);
	set_data((double **)sc,20,2);
	set_data(sd,3);
	set_data(se,4);
	set_data((double **)sf,3,5);
	set_data((double **)sg,3,6);
	set_data((double **)sh,3,7);
}