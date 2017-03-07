#include "rd_pdb.h"
#include "direct.h"

//main1 将一个蛋白质以Fasta格式输出
void main1()
{
	Protein protein;
	cout<<"正在读入"<<endl;
	protein.LoadFromPDB("pdb1ngr.ent");
	cout<<"正在写入"<<endl;
	protein.WriteSeqFasta("pdb1ngr.txt");
}
//main11 将一个蛋白质以Fasta格式输出，带参数
void main11(int argv,char **argc)
{
	Protein protein;
	string input,output;
	if(argv!=3)
	{
		cout<<"用法：程序名 蛋白质结构文件 输出序列文件"<<endl;
		return;
	}
	input=argc[1];
	output=argc[2];
	cout<<"正在读入"<<input.c_str()<<"...";
	if(protein.LoadFromPDB(input))
		cout<<"成功"<<endl;
	else
	{
		cout<<"失败"<<endl;
		return;
	}
	cout<<"正在写入"<<output.c_str()<<"...";
	if(protein.WriteSeqFasta(output))
	{
		cout<<"成功"<<endl;
	}
	else
		cout<<"失败"<<endl;
}
//main2 将一个目录下的所有蛋白质以Fasta格式输出到另一个目录下
void main2()
{
	intptr_t ptr;
	_finddata_t data;
	string sourcedir="11";//pdb文件所在的目录
	string desdir="22";  //fasta格式的文件存放的目录
	string temp1,temp2,findstr;

	_mkdir(desdir.c_str());
	findstr=sourcedir+"\\*.*";
	ptr=_findfirst(findstr.c_str(),&data);
	if(ptr==-1)
	{
		cout<<"目录不存在"<<sourcedir.c_str()<<endl;
		return;
	}
	do
	{
		if(!(data.attrib&_A_SUBDIR))
		{
			temp1=sourcedir+"\\"+data.name;
			temp2=desdir+"\\"+data.name+".faa";
			Protein protein;	
			protein.LoadFromPDB(temp1);	//读入PDB蛋白质
			//protein.WriteSeqFasta(temp2);//以fasta格式输出每个链的序列信息
			protein.WriteSeqFasta_total(temp2);//以fasta格式输出所有链的序列信息，链之间没有间隔
		}
	}while(_findnext(ptr,&data)==0);
}

void main3()
{
	Protein protein;
	cout<<"正在读入"<<endl;
	protein.LoadFromPDB("d1a04a2.ent");
	cout<<"正在写入"<<endl;
	//protein.SaveChains("1.txt");
	protein.SaveToFile("1.txt");

	/*cout<<"正在读入"<<endl;
	protein.LoadFromPDB("save1");

	cout<<"正在写入"<<endl;
	protein.SaveToFile("save2");*/
	/*Protein protein;
	cout<<"正在读入"<<endl;
	protein.LoadFromPDB("pdb1lba.ent");
	protein.WriteSeqFasta("1lbaeq.txt");*/
	/*cout<<"正在计算contact"<<endl;
	vector<Contact> all_contact;
	protein.ObtainContact(0,all_contact);
	protein.SaveContact(0,"contact.txt",all_contact);*/

	cout<<"完毕"<<endl;
}
//4 计算SAS,经过验证与DSSP计算的结果一致
void main4()
{
	string fname="T0295.pdb";
	Protein protein;

	if(!protein.LoadFromPDB(fname))
	{
		cout<<"无法读入蛋白质"<<endl;
		return;
	}
	if(!protein.CalSAA())
	{
		cout<<"计算过程出错"<<endl;
		return;
	}
	ofstream fileout("samplepdb_saresidue.txt");
	int num_chain=protein.num_chain;
	int i,j,num_aa;
	for(i=0;i<num_chain;i++)
	{
		num_aa=protein.m_chains[i].num_residue ;
		for(j=0;j<num_aa;j++)
		{
			fileout<<j<<"\t"<<protein.m_chains[i].m_residues[j].c_name<<"\t"<<protein.m_chains[i].m_residues[j].saa<<endl;
		}
	}
	int index,num_atom,k;
	int num_sa;
	bool find;
	for(i=0;i<num_chain;i++)
	{
		num_aa=protein.m_chains[i].num_residue ;
		num_sa=0;
		for(j=0;j<num_aa;j++)
		{
			num_atom=protein.m_chains[i].m_residues[j].num_atom;
			find=false;
			for(k=0;k<num_atom;k++)
			{
				index=protein.m_chains[i].m_residues[j].m_atoms[k];
				if(protein.m_chains[i].m_atoms[index].saa >2)
				{
					find=true;
					break;
				}
			}
			if(find)
				num_sa++;
		}
		cout<<"链"<<i<<"氨基酸总数"<<num_aa<<" 表面氨基酸总数"<<num_sa<<endl;
	}
}
//5 计算两面角 包括主链和侧链 经验证，与MODELLER软件包的计算结果一致
void main()
{
	Protein p;
	if(!p.LoadFromPDB("1A1x0.ent"))
	{
		cout<<"cannot load"<<endl;
		return;
	}
	if(!p.m_chains[0].CalculateDihedral())
	{
		cout<<"cannot calculate "<<endl;
		return;
	}
	int num_aa,i;
	num_aa=p.m_chains[0].num_residue;
	ofstream fileout("logdih.txt");
	fileout<<"索引\t氨基酸\tphi\tpsi\tchi1\tchi2\tchi3\tchi4"<<endl;
	for(i=0;i<num_aa;i++)
	{
		fileout<<i<<"\t"<<p.m_chains[0].m_residues[i].c_name<<"\t"<<p.m_chains[0].m_residues[i].phi<<"\t"<<p.m_chains[0].m_residues[i].psi<<"\t"<<p.m_chains[0].m_residues[i].chi_1<<"\t"<<p.m_chains[0].m_residues[i].chi_2<<"\t"<<p.m_chains[0].m_residues[i].chi_3<<"\t"<<p.m_chains[0].m_residues[i].chi_4<<endl;
	}
}