#include "amino_acid.h"

char residulename321(const char *three)
{
	if(strlen(three)!=3) return 'U';
	char tempthree[4];
	for(int i=0;i<3;i++)
		tempthree[i]=toupper(three[i]);
	tempthree[3]='\0';
	int po;
	po=BinarySearchStr(rname_three_rank,three,0,NUM_RESIDULE_TYPE-1);
	if(po!=-1) return three2one[po];
	else return 'U';
}
void residulename123(const char c,char *three)
{
	int po;
	po=BinarySearchChar(rname_one_rank,c,0,NUM_RESIDULE_TYPE-1);
	if(po!=-1) strcpy_s(three,4,one2three[po]);
	else strcpy_s(three,4,"UNK");
}
int Rname2ID(const char c)
{
	int po;
	char tempc;
	tempc=toupper(c);
	po=BinarySearchChar(rname_one_rank,tempc,0,NUM_RESIDULE_TYPE-1);
	if(po==-1) po=18;
	return po;
}
int Rname2ID(const char* three)
{
	char c;
	c=residulename321(three);
	return Rname2ID(c);
}
int Rname2ID20(const char c)
{
	int po;
	po=Rname2ID(c);
	if(po==0) return po;
	//ASX to ASN
	else if(po==1) po=11;
	//GLX to GLN
	else if(po==22) po=13;
	else if((po>1)&&(po<18)) po-=1;
	else if((po>18)&&(po<22)) po-=2;
	else po=-1;
	return po;
}
int Rname2ID20(const char* three)
{
	char c;
	c=residulename321(three);
	return Rname2ID20(c);
}
int Atomname2ID(const char* atom)
{
	int po;
	po=BinarySearchStr(atom_name_rank,atom,0,NUM_ATOM_TYPE-1);
	if(po==-1)
	{
#ifdef ERROR_INFO
		cout<<"did not find the atom"<<atom<<" the invalid rank number"<<endl;
#endif
		;
	}
	return po;
}
void Upper(string& str)
{
	int i,len;
	len=(int)str.length();
	for(i=0;i<len;i++)
		if((str[i]>='a')&&(str[i]<='z')) str[i]=toupper(str[i]);
}

int BinarySearchStr(const char threename[][4],const char *str,int l,int r)
{
	int left,right,middle;
	left=l;
	right=r;
	if(strcmp(str,threename[l])<0) return -1;
	if(strcmp(str,threename[r])>0) return -1;
	while(right-left>1)
	{
		middle=(left+right)/2;
		if(strcmp(str,threename[middle])==0) return middle;
		if(strcmp(str,threename[middle])<0) right=middle;
		else left=middle;
	}
	if(strcmp(str,threename[left])==0) return left;
	else if(strcmp(str,threename[right])==0) return right;
	else return -1;
}
int BinarySearchChar(const char *str,const char c,int l,int r)
{
	int left,right,middle;
	left=l;
	right=r;
	if(c<str[l]) return -1;
	if(c>str[r]) return -1;
	while(right-left>1)
	{
		middle=(left+right)/2;
		if(c==str[middle]) return middle;
		if(c<str[middle]) right=middle;
		else left=middle;
	}
	if(c==str[left]) return left;
	else if(c==str[right]) return right;
	else return -1;
}
