#include "rd_pdb.h"
#include "direct.h"

//main1 ��һ����������Fasta��ʽ���
void main1()
{
	Protein protein;
	cout<<"���ڶ���"<<endl;
	protein.LoadFromPDB("pdb1ngr.ent");
	cout<<"����д��"<<endl;
	protein.WriteSeqFasta("pdb1ngr.txt");
}
//main11 ��һ����������Fasta��ʽ�����������
void main11(int argv,char **argc)
{
	Protein protein;
	string input,output;
	if(argv!=3)
	{
		cout<<"�÷��������� �����ʽṹ�ļ� ��������ļ�"<<endl;
		return;
	}
	input=argc[1];
	output=argc[2];
	cout<<"���ڶ���"<<input.c_str()<<"...";
	if(protein.LoadFromPDB(input))
		cout<<"�ɹ�"<<endl;
	else
	{
		cout<<"ʧ��"<<endl;
		return;
	}
	cout<<"����д��"<<output.c_str()<<"...";
	if(protein.WriteSeqFasta(output))
	{
		cout<<"�ɹ�"<<endl;
	}
	else
		cout<<"ʧ��"<<endl;
}
//main2 ��һ��Ŀ¼�µ����е�������Fasta��ʽ�������һ��Ŀ¼��
void main2()
{
	intptr_t ptr;
	_finddata_t data;
	string sourcedir="11";//pdb�ļ����ڵ�Ŀ¼
	string desdir="22";  //fasta��ʽ���ļ���ŵ�Ŀ¼
	string temp1,temp2,findstr;

	_mkdir(desdir.c_str());
	findstr=sourcedir+"\\*.*";
	ptr=_findfirst(findstr.c_str(),&data);
	if(ptr==-1)
	{
		cout<<"Ŀ¼������"<<sourcedir.c_str()<<endl;
		return;
	}
	do
	{
		if(!(data.attrib&_A_SUBDIR))
		{
			temp1=sourcedir+"\\"+data.name;
			temp2=desdir+"\\"+data.name+".faa";
			Protein protein;	
			protein.LoadFromPDB(temp1);	//����PDB������
			//protein.WriteSeqFasta(temp2);//��fasta��ʽ���ÿ������������Ϣ
			protein.WriteSeqFasta_total(temp2);//��fasta��ʽ�����������������Ϣ����֮��û�м��
		}
	}while(_findnext(ptr,&data)==0);
}

void main3()
{
	Protein protein;
	cout<<"���ڶ���"<<endl;
	protein.LoadFromPDB("d1a04a2.ent");
	cout<<"����д��"<<endl;
	//protein.SaveChains("1.txt");
	protein.SaveToFile("1.txt");

	/*cout<<"���ڶ���"<<endl;
	protein.LoadFromPDB("save1");

	cout<<"����д��"<<endl;
	protein.SaveToFile("save2");*/
	/*Protein protein;
	cout<<"���ڶ���"<<endl;
	protein.LoadFromPDB("pdb1lba.ent");
	protein.WriteSeqFasta("1lbaeq.txt");*/
	/*cout<<"���ڼ���contact"<<endl;
	vector<Contact> all_contact;
	protein.ObtainContact(0,all_contact);
	protein.SaveContact(0,"contact.txt",all_contact);*/

	cout<<"���"<<endl;
}
//4 ����SAS,������֤��DSSP����Ľ��һ��
void main4()
{
	string fname="T0295.pdb";
	Protein protein;

	if(!protein.LoadFromPDB(fname))
	{
		cout<<"�޷����뵰����"<<endl;
		return;
	}
	if(!protein.CalSAA())
	{
		cout<<"������̳���"<<endl;
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
		cout<<"��"<<i<<"����������"<<num_aa<<" ���氱��������"<<num_sa<<endl;
	}
}
//5 ��������� ���������Ͳ��� ����֤����MODELLER������ļ�����һ��
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
	fileout<<"����\t������\tphi\tpsi\tchi1\tchi2\tchi3\tchi4"<<endl;
	for(i=0;i<num_aa;i++)
	{
		fileout<<i<<"\t"<<p.m_chains[0].m_residues[i].c_name<<"\t"<<p.m_chains[0].m_residues[i].phi<<"\t"<<p.m_chains[0].m_residues[i].psi<<"\t"<<p.m_chains[0].m_residues[i].chi_1<<"\t"<<p.m_chains[0].m_residues[i].chi_2<<"\t"<<p.m_chains[0].m_residues[i].chi_3<<"\t"<<p.m_chains[0].m_residues[i].chi_4<<endl;
	}
}