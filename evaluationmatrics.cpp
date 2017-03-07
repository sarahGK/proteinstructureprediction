/*************************************************
 * This program is to calculate the performance metrics 
 ************************************************/

// define a structure to represent one native protein structure and its decoys
typedef struct
{
	string pdbid;                   //ID
	int num_aa;	
	int num_decoys;                //number of decoys
	vector<int> decoys;	           // all the decoy ids
	vector<double> potvalue_decoys;// all the potential values of decoys
	vector<double> rmsd;	         // all the RMSD values of decoys
	vector<int> potrank;	         // the rank of potential values
	double potvalue_native;        // the potential value of the native structure
	
// all the performance metrics
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
  //10% enrichment fraction
	double ef=0.1;			
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
			if(allhr[i].potvalue_native>allhr[i].potvalue_decoys[j]) // if compare the probability, it's '<' instead of '>'
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
			
		Index_Sort_Reverse(allhr[i].potvalue_decoys,allhr[i].potrank); // if compare the probability,
		
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
 * sort va by bubble sorting and save the values into index
 * ***************************************************************************/
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
	//bubble sorting
	for(i=len-1;i>0;i--)
	{
		for(j=0;j<i;j++)
		{
			if(va[index[j]]>va[index[j+1]])// change into '<' if sort from the biggest 
			{
				temp=index[j];
				index[j]=index[j+1];
				index[j+1]=temp;
			}
		}
	}
}
