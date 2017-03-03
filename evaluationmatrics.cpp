/*************************************************
 * 说明：此程序是按照计算的能量进行评估
 * 如果按照概率评估则需要修改部分代码，已注释
 ***********************************************/
//定义结构：代表一个蛋白质的天然结构和对应的非天然结构
typedef struct
{
	string pdbid;//ID
	int num_aa;	
	int num_decoys;//扰动结构个数
	vector<int> decoys;	//扰动结构标识
	vector<double> potvalue_decoys;//每个扰动结构计算的能量(概率)
	vector<double> rmsd;	//每个扰动结构的RMSD
	vector<int> potrank;	//能量的排序
	double potvalue_native;//天然结构对应的能量
	
	//评价标准
	int rank;
	double z;
	double drmsd;
	double ef;
	double cc;
} HR_genel;

{
	vector<HR_genel> allhr;

	//evaluated and output
	string pdbid;
	cout<<"evaluate..."<<endl;
	double mean,mean2,dix,rmsd;
	double accuracy,averank,avez,avedrmsd,aveef,avecc;
	double ef=0.1;			//10% enrichment fraction
	int num_ef;
	string tfname1,tfname2;
	ofstream fileout(outname.c_str());
	if(!fileout.is_open())
	{
		cout<<"cannot open"<<outname.c_str()<<endl;
		return 0;
	}
	accuracy=averank=avez=avedrmsd=aveef=avecc=0;
	fileout<<"pdbid\trank\tz\tdrmsd\tfe\tcc"<<endl;
	for(i=0;i<total;i++)
	{
		pdbid=allhr[i].pdbid;
		allhr[i].rank=1;
		num_decoys=allhr[i].num_decoys;
		//rank and Z-score
		mean=mean2=0;
		for(j=0;j<num_decoys;j++)
		{
			if(allhr[i].potvalue_native>allhr[i].potvalue_decoys[j])//概率将>改为<
				allhr[i].rank++;
			mean+=allhr[i].potvalue_decoys[j];
			mean2+=allhr[i].potvalue_decoys[j]*allhr[i].potvalue_decoys[j];
		}
		mean/=num_decoys;
		mean2/=num_decoys;
		dix=sqrt(mean2-mean*mean);
		allhr[i].z=(mean-allhr[i].potvalue_native)/dix;
		
		if(allhr[i].rank==1)
			accuracy+=1.0;
			
		Index_Sort_Reverse(allhr[i].potvalue_decoys,allhr[i].potrank);////概率按照从大到小排序
		
		if(num_decoys>0)
		{
			//drmsd
			sprintf(str,"%d",allhr[i].decoys[0]);
			tfname1=hrdir+pdbid+"/"+pdbid+"."+str+".coord.pdb";
			sprintf(str,"%d",allhr[i].decoys[allhr[i].potrank[0]]);
			tfname2=hrdir+pdbid+"/"+pdbid+"."+str+".coord.pdb";
			if(!Ev_Rmsd(tfname1,tfname2,1,rmsd))
			{
				cout<<"waring: cannot comput rmsd"<<endl;
			}
			allhr[i].drmsd=rmsd;
			//ef
			num_ef=int(num_decoys*ef);
			int ef1=0;
			for(j=0;j<num_ef;j++)
			{
				if(allhr[i].potrank[j]<num_ef)
					ef1++;
			}
			allhr[i].ef=(ef1/(double)num_ef)/ef;	
			//cc
			double meanx,meany,cx,cy,ta,tb,tc;
			meanx=meany=0;
			for(j=0;j<num_decoys;j++)
			{
				meanx+=allhr[i].rmsd[j];
				meany+=allhr[i].potvalue_decoys[j];
			}		
			meanx/=num_decoys;meany/=num_decoys;
			ta=tb=tc=0;
			for(j=0;j<num_decoys;j++)
			{
				cx=allhr[i].rmsd[j]-meanx;
				cy=allhr[i].potvalue_decoys[j]-meany;
				ta+=cx*cy;tb+=cx*cx;tc+=cy*cy;
			}
			allhr[i].cc=ta/sqrt(tb*tc);				
		}
		averank+=allhr[i].rank;
		avez+=allhr[i].z;
		avedrmsd+=allhr[i].drmsd;
		aveef+=allhr[i].ef;
		avecc+=allhr[i].cc;
		fileout<<allhr[i].pdbid<<'\t'<<allhr[i].rank<<'\t'<<allhr[i].z<<'\t'<<allhr[i].drmsd<<'\t'<<allhr[i].ef<<'\t'<<allhr[i].cc<<endl;
	}
	accuracy/=total;
	averank/=total;
	avez/=total;
	avedrmsd/=total;
	aveef/=total;
	avecc/=total;
	fileout<<"ave accuracy: "<<accuracy<<endl;
	fileout<<"ave rank: "<<averank<<endl;
	fileout<<"ave z: "<<avez<<endl;
	fileout<<"ave drmsd: "<<avedrmsd<<endl;
	fileout<<"ave FE: "<<aveef<<endl;	
	fileout<<"ave CC: "<<avecc<<endl;
}
/***************************************************************************
 * 索引排序：对va中的数从小到大排序，排序索引放于index中，不改变va中的数
 ***************************************************************************/
void Index_Sort_Reverse(const vector<double>& va,vector<int>& index)
{
	int i,j,len;
	int temp;
	len=(int)va.size();

	if(len==0)
	{
		index.clear();
		return;
	}
	index.resize(len);
	for(i=0;i<len;i++)
		index[i]=i;
	//起泡排序
	for(i=len-1;i>0;i--)
	{
		for(j=0;j<i;j++)
		{
			if(va[index[j]]>va[index[j+1]])//从大到小排序只需将〉改为<
			{
				temp=index[j];
				index[j]=index[j+1];
				index[j+1]=temp;
			}
		}
	}
}