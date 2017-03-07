#include "rd_pdb.h"

//write protein information in the format with Fasta
void main1()
{
	Protein protein;
	cout<<"Reading data ..."<<endl;
	protein.LoadFromPDB("pdb1ngr.ent");
	cout<<"Writing into file now"<<endl;
	protein.WriteSeqFasta("pdb1ngr.txt");
}

void main2()
{
	Protein protein;
	cout<<"Reading"<<endl;
	protein.LoadFromPDB("d1a04a2.ent");
	cout<<"Writing"<<endl;
	protein.SaveToFile("1.txt");
	cout<<"Finish"<<endl;
}
//calculate SAS, and results are the same as DSSP
void main3()
{
	string fname="T0295.pdb";
	Protein protein;

	if(!protein.LoadFromPDB(fname))
	{
		cout<<"Can not read Protein"<<endl;
		return;
	}
	if(!protein.CalSAA())
	{
		cout<<"Error in Calculation"<<endl;
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
		cout<<"Chain"<<i<<"Number of amino acids"<<num_aa<<"the number of SA"<<num_sa<<endl;
	}
}
// calculate backbone and side chain dihedral angles and the results are the same as MODEllER
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
	fileout<<"index\t°amino acid\tphi\tpsi\tchi1\tchi2\tchi3\tchi4"<<endl;
	for(i=0;i<num_aa;i++)
	{
		fileout<<i<<"\t"<<p.m_chains[0].m_residues[i].c_name<<"\t"<<p.m_chains[0].m_residues[i].phi<<"\t"<<p.m_chains[0].m_residues[i].psi<<"\t"<<p.m_chains[0].m_residues[i].chi_1<<"\t"<<p.m_chains[0].m_residues[i].chi_2<<"\t"<<p.m_chains[0].m_residues[i].chi_3<<"\t"<<p.m_chains[0].m_residues[i].chi_4<<endl;
	}
}
