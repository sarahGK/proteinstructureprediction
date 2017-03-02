#include "rd_pdb.h"
#include "amino_acid.h"
#include "vector.h"
#include "rdfasta.h"
#include "calcsurf.h"
#include "dssp\\calcaccsurf.h"
Point::Point()
{
	x=0;y=0;z=0;
}
Atom::Atom()
{
	chain_seqm=0;
	c_resname='x';
	residule_chain_seqm=0;
	chain_identifier=-1;
	radium=0;
}
Residue::Residue()
{
	c_name='X';
	chain_seqm=0;
	num_atom=0;	
	alpha=phi=psi=omega=chi_1=chi_2=chi_3=chi_4=chi_5=-999.0;
	for(int i=0;i<NUM_ATOM_TYPE;i++)
		atom_index[i]=-1;
	num_contact=0;
}
Chain::Chain()
{
	chain_identifier='X';
	num_residue=0;
	num_atom=0;
}
Protein::Protein()
{
	num_chain=0;
}
void Residue::Clear()
{
	name.clear();
	c_name='X';
	chain_seqm=0;
	num_atom=0;
	m_atoms.clear();
	alpha=phi=psi=omega=chi_1=chi_2=chi_3=chi_4=chi_5=-999.0;
	for(int i=0;i<NUM_ATOM_TYPE;i++)
		atom_index[i]=-1;
	num_contact=0;
	saa=0;
	side_center.x=side_center.y=side_center.z=0;
	sidecenter_te13.x=sidecenter_te13.y=sidecenter_te13.z=0;;
	sidecenter_HRSC.x=sidecenter_HRSC.y=sidecenter_HRSC.z=0;
}
void Chain::Clear()
{
	chain_identifier=-1;
	m_residues.clear();
	num_residue=0;
	m_atoms.clear();
	num_atom=0;
	seq.clear();
	
}
void Protein::Clear()
{
	m_chains.clear();
	num_chain=0;
	PDB_id.clear();
	method.clear();
	resolution=0;
}
/*************************************************
 * read data from PDB files 
 * return true if sucess, otherwise return false
*************************************************/
bool  Protein::LoadFromPDB(string filename)
{
	ifstream filein;
	string line;
	bool model=false;//if multiple models,only deal with the first one
	vector<Atom> all_atom;
	if(filename.empty())
		return false;
	filein.open(filename.c_str(),ios::in);
	if(!filein.is_open())
		return false;

	Clear();
	while(!filein.eof())
	{
		getline(filein,line);		
		if(strcmp(line.substr(0,6).c_str(),"HEADER")==0) 
			process_header(line);
		else if(strcmp(line.substr(0,6).c_str(),"SOURCE")==0) 
			process_source(line);
		else if(strcmp(line.substr(0,6).c_str(),"EXPDTA")==0)
			process_expdta(line);
		else if(strcmp(line.substr(0,6).c_str(),"REMARK")==0)
			process_remark(line);
		else if(strcmp(line.substr(0,6).c_str(),"ATOM  ")==0) 
		{
			if(!model)
				process_atom(line,all_atom);
		}
		else if(strcmp(line.substr(0,6).c_str(),"ENDMDL")==0)
			model=true;
		else process_field(line);
	}

	AddAtom(all_atom);
	return true;
}
/*********************************************************
 * save information into files and return true if success
 ********************************************************/
bool Protein::SaveToFile(string filename)
{
	ofstream file;
	string line;
	Atom atom;
	int i,j;
	if(filename.empty())
		return false;
	file.open(filename.c_str(),ios::out);
	if(!file.is_open())
		return false;
	OutputPreField(file);
	int num_atom;
	for(i=0;i<num_chain;i++)
	{
		num_atom=m_chains[i].num_atom;
		for(j=0;j<num_atom;j++)
		{
			atom=m_chains[i].m_atoms[j];
			Atom2Line(atom,line);
			file<<line<<endl;
		}
		Atom2Ter(atom,line);
		file<<line<<endl;
	}
	OutputSucField(file);
	return true;
}

void Protein::process_header(const string& line)
{
	if((int)line.size()>=66)
		PDB_id=line.substr(62,4);
	process_field(line);
}
void Protein::process_expdta(const string& line)
{
	if(line.find(NMR)!=-1)
		method=NMR;
	else if(line.find(X_RAY)!=-1)
		method=X_RAY;
}
void Protein::process_source (const string& line)
{
	process_field(line);
}
void Protein::process_remark(const string& line)
{
	string id=line.substr(6,4);
	int r=atoi(id.c_str());
	double re;
	if(r==2)//处for resolution field
	{
		if(resolution<=0.0001)// haven't assigned value yet
		{
			id=line.substr(23,4);
			re=-1;
			re=atof(id.c_str());
			if(re<=0)
				resolution=0;
			else
				resolution=re;
		}
	}
	process_field(line);
}
void Protein::process_atom (const string& line, vector<Atom>& all_atom)
{
	int i,len;
	len=(int)line.size();
	Atom atom;
	string subline;
	if(len<7)
		return;
	subline=line.substr(6,5);
	atom.chain_seqm=atoi(subline.c_str());
	subline=line.substr(12,4);
	for(i=0;i<4;i++)
	{
		if(subline[i]!=' ')
			atom.name.push_back(subline[i]);
		Upper(atom.name);
	}
	// determine radium based on atom type and get the value from DSSP
	if(atom.name=="N")
		atom.radium=RN;
	else if(atom.name=="CA")
		atom.radium=RCA;
	else if(atom.name=="C")
		atom.radium=RC;
	else if(atom.name=="O")
		atom.radium=RO;
	else
		atom.radium=RSIDEATOM;

	subline=line.substr(17,3);
	atom.residule_name=subline;
	Upper(atom.residule_name);
	atom.c_resname=residulename321(subline.c_str());
	atom.chain_identifier=line.length()>=22?line[21]:' ';
	subline=line.substr(22,5);
	if(IsDigit(subline))
		atom.residule_chain_seqm=atoi(subline.c_str());
	else
		//atom.residule_chain_seqm=-1;//ignore the overlap 
		return;
	subline=line.substr(30,8);
	atom.pt.x=atof(subline.c_str());
	subline=line.substr(38,8);
	atom.pt.y=atof(subline.c_str());
	subline=line.substr(46,8);
	atom.pt.z=atof(subline.c_str());
	if(len>=60)
	{
		subline=line.substr(54,6);
		atom.occupancy=atof(subline.c_str());
	}
	if(len>=66)
	{
		subline=line.substr(60,6);
		atom.tempFactor=atof(subline.c_str());
	}
	if(len>=73)
		atom.segID=line.substr(72,4);
	if(len>=77)
		atom.element=line.substr(76,2);
	if(len>=79)
		atom.charge=line.substr(78,2);
	all_atom.push_back(atom);
}
/************************************************************
 * Add the atoms into protein by the order of reading sequence
 ************************************************************/
bool Protein::AddAtom(const vector<Atom>& all_atom)
{
	int len,i,index;
	Atom atom;
	Chain chain;
	Residue residue;

	chain.chain_identifier=NOCHAIN;
	residue.chain_seqm=NOSEQNUM;
	len=(int)all_atom.size();
	if(len<=0)//if no atom
		return true;
	for(i=0;i<len;i++)
	{
		atom=all_atom[i];
		if(atom.residule_chain_seqm!=residue.chain_seqm)
		{
			if(residue.chain_seqm!=NOSEQNUM)
			{
				chain.m_residues.push_back(residue);
				chain.num_residue++;			
				chain.seq.push_back(residue.c_name);
				residue.Clear();
			}
				residue.chain_seqm=atom.residule_chain_seqm;
				residue.name=atom.residule_name;
				residue.c_name=atom.c_resname;
				residue.num_atom=0;			
		}
		if(atom.chain_identifier!=chain.chain_identifier)
		{
			if(chain.chain_identifier!=NOCHAIN)
			{
				m_chains.push_back(chain);
				chain.Clear();			
				num_chain++;
			}
				chain.chain_identifier=atom.chain_identifier;
				chain.num_residue=0;	
				chain.num_atom=0;			
		}		
		atom.residule_seqm=chain.num_residue;
		chain.m_atoms.push_back(atom);
		residue.m_atoms.push_back(chain.num_atom);
		//residue.chain_seqm=atom.residule_chain_seqm;
		index=Atomname2ID(atom.name.c_str());
		if(index!=-1)
		{
			residue.atom_index[index]=chain.num_atom;
		}
		chain.num_atom++;		
		residue.num_atom++;
	}
	chain.m_residues.push_back(residue);
	chain.num_residue++;
	chain.seq.push_back(residue.c_name);
	m_chains.push_back(chain);
	num_chain++;
	return true;
}

void Protein::Atom2Line(const Atom& atom,string& line)
{
	char str[81];
	char name[5]="    ";
	if((int)atom.name.length()==4)//if the length of the name of atom is 4
	{
		for(int i=0;i<(int)atom.name.length();i++)
			name[i]=atom.name[i];
	}
	else//if the length is not 4,then the fisrt letter will be blank
	{
		for(int i=0;i<(int)atom.name.length()&&i<3;i++)
			name[i+1]=atom.name[i];
	}
	memset(str,' ',81);
	str[80]='\0';
	sprintf_s(str,"%s%5d %s %s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %s%s%s","ATOM  ",atom.chain_seqm,name,atom.residule_name.c_str(),atom.chain_identifier,atom.residule_chain_seqm,atom.pt.x,atom.pt.y,atom.pt.z,atom.occupancy,atom.tempFactor,atom.segID.c_str(),atom.element.c_str(),atom.charge.c_str());
	line=str;
}
void Protein::Atom2Ter(const Atom& atom,string& line)
{
	char str[81];	
	memset(str,' ',81);
	str[80]='\0';
	sprintf_s(str,"%s%5d      %s %c%4d                                                      ","TER   ",atom.chain_seqm+1,atom.residule_name.c_str(),atom.chain_identifier,atom.residule_chain_seqm);
	line=str;
}
/********************************************
  * get all the contacts and return the total number
 ********************************************/
int Protein::ObtainContact(int chainnumber,vector<Contact>& all_contact)
{
	if(chainnumber<num_chain)
	{
		return m_chains[chainnumber].ObtainContact(all_contact);
	}
	else
		return 0;
}
/**************************************************
 * save contacts into files
 * return 0 if sucess, otherwise return -1
 ***********************************************/
int Protein::SaveContact(int chainnumber,const string& filename,const vector<Contact>& all_contact)
{
	if(chainnumber<num_chain)
		return m_chains[chainnumber].SaveContact(filename,all_contact);
	else
		return -1;
}
int Chain::ObtainContact(vector<Contact >& all_contact)
{
	int i,j,len;
	int id;
	double sider,d,t;
	bool findca,findcb;
	Point capt,cbpt,sidept;
	len=num_residue;
	//
	for(i=0;i<len;i++)
	{
		id=Rname2ID(m_residues[i].c_name);
		sider=AA_SideR[id];
		m_residues[i].side_radium=sider;
		findca=findcb=false;
		for(j=0;j<m_residues[i].num_atom;j++)
		{
			if(m_atoms[m_residues[i].m_atoms[j]].name=="CA")
			{
				capt=m_atoms[m_residues[i].m_atoms[j]].pt;
				findca=true;
			}
			if(m_atoms[m_residues[i].m_atoms[j]].name=="CB")
			{
				cbpt=m_atoms[m_residues[i].m_atoms[j]].pt;
				findcb=true;
			}
		}
		if(findca&&findcb)
		{
			d=Distance(capt,cbpt);
			t=sider/d;
			sidept.x=capt.x+(cbpt.x-capt.x)*t;
			sidept.y=capt.y+(cbpt.y-capt.y)*t;
			sidept.z=capt.z+(cbpt.z-capt.z)*t;
			m_residues[i].side_center=sidept;
		}
		else if(findca)
		{
			m_residues[i].side_center=capt;
		}				
	}
	Contact contact;
	double r1,r2;
	for(i=0;i<len;i++)
	{
		r1=m_residues[i].side_radium;
		for(j=i+1;j<len;j++)
		{
			r2=m_residues[j].side_radium;
			d=Distance(m_residues[i].side_center,m_residues[j].side_center);
			if(d<r1+r2+2.8)
			{
				contact.i1=i;
				contact.i2=j;
				all_contact.push_back(contact);
			}
		}
	}
	return (int)all_contact.size();
}
/*******************************************
 * 将一个链的contact信息存入文件中
 * 存储格式为：
 * 氨基酸序列
 * contact总数
 * 所有的contact
 * 成功时返回0，否则，返回非0，表示
 * 1 打不开文件
 * 2 无法写入文件
 ******************************************/
int Chain::SaveContact(const string& filename,const vector<Contact>& all_contact)
{
	ofstream fileout;
	int len;
	if(filename.empty())
		return 1;
	fileout.open(filename.c_str(),ios::out);
	if(!fileout.is_open())
		return 1;
	fileout<<seq.c_str()<<endl;
	len=(int)all_contact.size();
	fileout<<len<<endl;
	for(int i=0;i<len;i++)
		fileout<<all_contact[i].i1<<' '<<all_contact[i].i2<<' ';	
	return 0;
}
int Protein::Get_Num_chain()
{
	return num_chain;
}
void Protein::Get_PDBID(string& id)
{
	id=PDB_id;
}
bool Protein::GetSeq(int chain,string& seq)
{
	if(chain>num_chain)
		return false;
	seq=m_chains[chain].seq;
	return true;
}
bool Chain::CalculateDihedral()
{
	int len,i,index,atom_index,rindex;
	Point p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11;
	len=(int)m_residues.size();
	for(i=0;i<len;i++)
	{
		//计算alpha
		if((i>=1)&&(i<len-2))
		{
			index=Atomname2ID("CA");
			atom_index=m_residues[i-1].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i-1<<"中找不到CA原子"<<endl;
				return false;
			}
			p1=m_atoms[atom_index].pt;
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i<<"中找不到CA原子"<<endl;
				return false;
			}
			p2=m_atoms[atom_index].pt;
			atom_index=m_residues[i+1].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i+1<<"中找不到CA原子"<<endl;
				return false;
			}
			p3=m_atoms[atom_index].pt;
			atom_index=m_residues[i+2].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i+2<<"中找不到CA原子"<<endl;
				return false;
			}
			p4=m_atoms[atom_index].pt;
			m_residues[i].alpha=Dihedralangle(p1,p2,p3,p4);
		}		

		//计算phi,psi,omega
		//if(i>0)
		{
			if(i>0)
			{
				index=Atomname2ID("C");
				atom_index=m_residues[i-1].atom_index[index];
				if(atom_index==-1)
				{
					cout<<"计算两面角时出错：在氨基酸"<<i-1<<"中找不到C原子"<<endl;
					return false;
				}
				p1=m_atoms[atom_index].pt;
			}
			index=Atomname2ID("N");
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i<<"中找不到N原子"<<endl;
				return false;
			}
			p2=m_atoms[atom_index].pt;
			index=Atomname2ID("CA");
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i<<"中找不到CA原子"<<endl;
				return false;
			}
			p3=m_atoms[atom_index].pt;
			index=Atomname2ID("C");
			atom_index=m_residues[i].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算两面角时出错：在氨基酸"<<i<<"中找不到C原子"<<endl;
				return false;
			}
			p4=m_atoms[atom_index].pt;
			if(i>0)
				m_residues[i].phi =Dihedralangle(p1,p2,p3,p4);
			if(i<len-1)
			{
				index=Atomname2ID("N");
				atom_index=m_residues[i+1].atom_index[index];
				if(atom_index==-1)
				{
					cout<<"计算两面角时出错：在氨基酸"<<i+1<<"中找不到C原子"<<endl;
					return false;
				}
				p5=m_atoms[atom_index].pt;
				m_residues[i].psi=Dihedralangle(p2,p3,p4,p5);
			
				index=Atomname2ID("CA");
				atom_index=m_residues[i+1].atom_index[index];
				if(atom_index==-1)
				{
					cout<<"计算两面角时出错：在氨基酸"<<i+1<<"中找不到CA原子"<<endl;
					return false;
				}
				p6=m_atoms[atom_index].pt;
				m_residues[i].omega=Dihedralangle(p3,p4,p5,p6);
			}
			//计算侧链角度
			rindex=Rname2ID20(m_residues[i].c_name);
			if(rindex!=-1)
			{
				if(haschi1[rindex])//计算X1
				{
					index=Atomname2ID("CB");
					atom_index=m_residues[i].atom_index[index];
					if(atom_index==-1)
					{
						cout<<"计算侧链两面角时出错：在氨基酸"<<i<<"中找不到CB原子"<<endl;
						return false;
					}
					p7=m_atoms[atom_index].pt;

					index=atomchi1[rindex];
					atom_index=m_residues[i].atom_index[index];
					if(atom_index==-1)
					{
						cout<<"计算侧链X1两面角时出错：在氨基酸"<<i<<"中找不到指定的原子"<<endl;
						return false;
					}
					p8=m_atoms[atom_index].pt;
					m_residues[i].chi_1=Dihedralangle(p2,p3,p7,p8);
				}
				if(haschi2[rindex])//计算X2
				{
					index=atomchi2[rindex];
					atom_index=m_residues[i].atom_index[index];
					if(atom_index==-1)
					{
						cout<<"计算侧链X1两面角时出错：在氨基酸"<<i<<"中找不到指定的原子"<<endl;
						return false;
					}
					p9=m_atoms[atom_index].pt;
					m_residues[i].chi_2=Dihedralangle(p3,p7,p8,p9);
				}
				if(haschi3[rindex])//计算X3
				{
					index=atomchi3[rindex];
					atom_index=m_residues[i].atom_index[index];
					if(atom_index==-1)
					{
						cout<<"计算侧链X1两面角时出错：在氨基酸"<<i<<"中找不到指定的原子"<<endl;
						return false;
					}
					p10=m_atoms[atom_index].pt;
					m_residues[i].chi_3=Dihedralangle(p7,p8,p9,p10);
				}
				if(haschi4[rindex])//计算X4
				{
					index=atomchi4[rindex];
					atom_index=m_residues[i].atom_index[index];
					if(atom_index==-1)
					{
						cout<<"计算侧链X1两面角时出错：在氨基酸"<<i<<"中找不到指定的原子"<<endl;
						return false;
					}
					p11=m_atoms[atom_index].pt;
					m_residues[i].chi_4=Dihedralangle(p8,p9,p10,p11);
				}
			}
			
		}
		
	
	}
		return true;
}
bool Chain::WriteDihedral(const string filename)
{
	ofstream file(filename.c_str());
	if(!file.is_open())
		return false;
	int i,len;
	len=(int)m_residues.size();
	file<<"两面角"<<endl;
	
	file<<"序号"<<'\t'<<"氨基酸名称"<<'\t'<<"alpha"<<'\t'<<"phi"<<'\t'<<"psi\t"<<"omega\t"<<endl;
	for(i=0;i<len;i++)
	{
		file<<i<<'\t'<<m_residues[i].name<<'\t';
		file.precision(7);
		file<<m_residues[i].alpha<<'\t';
		file.precision(7);
		file<<m_residues[i].phi<<'\t';
		file.precision(7);
		file<<m_residues[i].psi<<'\t';
		file.precision(7);
		file<<m_residues[i].omega <<'\t';

		file<<endl;
	}
	return true;
}
/*********************************************************
 * 根据两面角进行蛋白质三维结构的坐标重建
 * 具体方法请参考生物信息学学习笔记2中的专题8.8
 * 需要说明的是,由于在重建过程中使用的键长键角
 * 都是在能量最低时的经验值，重建后的坐标并不能
 * 和原来坐标完全重合
 * 校准的键长键角取自AMBER势能函数，见蛋白质结构
 * 准则一文。
 * 另外，即使是按照PDB文件中给出的坐标来计算键长
 * 和键角，也不能和标准状态下的值，完全一样
 * 目前仅完成对主链原子(除了O原子)进行重建
 * 实验结果显示：这种重建方法和实验值基本吻合
 * 但是会出现累计误差，后面的原子的坐标偏差较大
 * 大约为2埃
 *********************************************************/
bool Chain::Rebuildxyz_Dihedral()
{
	int len=(int)m_residues.size();
	int i;
	int index,atom_index;
	Point qN,qCA,qC;//前面氨基酸的Ca,C,N原子坐标
	Point N,CA,C,O;//当前氨基酸的N,Ca,C,O原子坐标
	double theta,tao,bond_len;

	//对第一个氨基酸的相应原子赋值，这里拷贝原有坐标

	for(i=1;i<len-1;i++)//依次对随后的氨基酸中的原子进行坐标重建
	{
		//计算当前氨基酸的N
		index=Atomname2ID("N");
		atom_index=m_residues[i-1].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"坐标重建时时出错：在氨基酸"<<i-1<<"中找不到N原子"<<endl;
			return false;
		}
		qN=m_atoms[atom_index].pt;

		index=Atomname2ID("CA");
		atom_index=m_residues[i-1].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"坐标重建时时出错：在氨基酸"<<i-1<<"中找不到CA原子"<<endl;
			return false;
		}
		qCA=m_atoms[atom_index].pt;
		index=Atomname2ID("C");
		atom_index=m_residues[i-1].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"坐标重建时时出错：在氨基酸"<<i-1<<"中找不到C原子"<<endl;
			return false;
		}
		qC=m_atoms[atom_index].pt;
		
		//计算当前的N原子坐标
		theta=m_residues[i-1].psi;//两面角
		tao=116.6;//Ca-C-N键角
		bond_len=1.335;//C-N键长
		N=LastXYZ(qN,qCA,qC,theta,tao,bond_len);
		index=Atomname2ID("N");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"坐标重建时时出错：在氨基酸"<<i<<"中找不到N原子"<<endl;
			return false;
		}
		m_atoms[atom_index].pt=N;

		//对当前的CA原子进行坐标重建
		theta=m_residues[i-1].omega;
		tao=121.9;//C-N-CA键角
		bond_len=1.449;//N-CA键长
		CA=LastXYZ(qCA,qC,N,theta,tao,bond_len);
		index=Atomname2ID("CA");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"坐标重建时时出错：在氨基酸"<<i<<"中找不到CA原子"<<endl;
			return false;
		}
		m_atoms[atom_index].pt=CA;

		//对当前的C原子进行坐标重建
		theta=m_residues[i].phi;
		tao=110.3;//N-CA-C键角
		bond_len=1.522;//CA-C键长
		C=LastXYZ(qC,N,CA,theta,tao,bond_len);
		index=Atomname2ID("C");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"坐标重建时时出错：在氨基酸"<<i<<"中找不到C原子"<<endl;
			return false;
		}
		m_atoms[atom_index].pt=C;


	}
	return true;
}
bool Protein::WriteSeqFasta(string filename)
{
	if(filename.empty())
		return false;
	ofstream file;
	file.open(filename.c_str());
	if(!file.is_open())
		return false;
	string header,header1;
	header1=this->PDB_id;
	if(header1.empty())
		header1=filename;
	int i;
	for(i=0;i<(int)m_chains.size();i++)
	{
		header=header1+m_chains[i].chain_identifier;
		if(!WriteFastaSeq(file,m_chains[i].seq,header))
			return false;
	}
	return true;
}
bool Protein::WriteSeqFasta_total(string filename)
{
	if(filename.empty())
		return false;
	string seq;
	ofstream file;
	file.open(filename.c_str());
	if(!file.is_open())
		return false;
	string header,header1;
	header1=this->PDB_id;
	if(header1.empty())
		header1=filename;
	int i;
	for(i=0;i<(int)m_chains.size();i++)
	{		
		seq+=m_chains[i].seq;			
	}
	if(!WriteFastaSeq(file,seq,header1))
		return false;
	return true;
}
//处理零散的字段，仅仅记录数据
void Protein::process_field(const string& line)
{
	string subline=line.substr(0,6);
	if(subline=="ANISOU")
		anisou.push_back(line);
	else if(subline=="AUTHOR")
		author.push_back(line);
	else if(subline=="CISPEP")
		cispep.push_back(line);
	else if(subline=="COMPND")
		compnd.push_back(line);
	else if(subline=="CONECT")
		conect.push_back(line);
	else if(subline=="CRYST1")
		cryst1.push_back(line);
	else if(subline=="DBREF ")
		dbref.push_back(line);
	else if(subline=="EXPDTA")
		expdta.push_back(line);
	else if(subline=="FORMUL")
		formul.push_back(line);
	else if(subline=="FTNOTE")
		ftnote.push_back(line);
	else if(subline=="HEADER")
		header.push_back(line);
	else if(subline=="HELIX ")
		helix.push_back(line);
	else if(subline=="HET   ")
		het.push_back(line);
	else if(subline=="HETATM")
		hetatm.push_back(line);
	else if(subline=="HETNAM")
		hetnam.push_back(line);
	else if(subline=="HETSYN")
		hetsyn.push_back(line);
	else if(subline=="HYDBND")
		hydbnd.push_back(line);
	else if(subline=="JNRL  ")
		jnrl.push_back(line);
	else if(subline=="JRNL  ")
		jrnl.push_back(line);
	else if(subline=="KEYWDS")
		keywds.push_back(line);
	else if(subline=="LINK  ")
		link.push_back(line);
	else if(subline=="MASTER")
		master.push_back(line);
	else if(subline=="MODRES")
		modres.push_back(line);
	else if(subline=="MTRIX1")
		mtrix1.push_back(line);
	else if(subline=="MTRIX2")
		mtrix2.push_back(line);
	else if(subline=="MTRIX3")
		mtrix3.push_back(line);
	else if(subline=="ORIGX1")
		origx1.push_back(line);
	else if(subline=="ORIGX2")
		origx2.push_back(line);
	else if(subline=="ORIGX3")
		origx3.push_back(line);
	else if(subline=="REMARK")
		remark.push_back(line);
	else if(subline=="REVDAT")
		revdat.push_back(line);
	else if(subline=="SCALE1")
		scale1.push_back(line);
	else if(subline=="SCALE2")
		scale2.push_back(line);
	else if(subline=="SCALE3")
		scale3.push_back(line);
	else if(subline=="SEQADV")
		seqadv.push_back(line);
	else if(subline=="SEQRES")
		seqres.push_back(line);
	else if(subline=="SHEET ")
		sheet.push_back(line);
	else if(subline=="SITE  ")
		site.push_back(line);
	else if(subline=="SLTBRG")
		sltbrg.push_back(line);
	else if(subline=="SOURCE")
		source.push_back(line);
	else if(subline=="SPRSDE")
		sprsde.push_back(line);
	else if(subline=="SSBOND")
		ssbond.push_back(line);
	else if(subline=="TITLE ")
		title.push_back(line);
	else if(subline=="TURN  ")
		turn.push_back(line);
	else if(subline=="TVECT")
		tvect.push_back(line);
	else if(subline=="END   ")
		end.push_back(line);
}
void Protein::OutputPreField(ofstream& file)
{
	int i;
	if(!header.empty())
	{
		for(i=0;i<(int)header.size();i++)
			file<<header[i]<<endl;
	}
	if(!title.empty())
	{
		for(i=0;i<(int)title.size();i++)
			file<<title[i]<<endl;
	}
	if(!compnd.empty())
	{
		for(i=0;i<(int)compnd.size();i++)
			file<<compnd[i]<<endl;
	}
	if(!source.empty())
	{
		for(i=0;i<(int)source.size();i++)
			file<<source[i]<<endl;
	}
	if(!keywds.empty())
	{
		for(i=0;i<(int)keywds.size();i++)
			file<<keywds[i]<<endl;
	}
	if(!expdta.empty())
	{
		for(i=0;i<(int)expdta.size();i++)
			file<<expdta[i]<<endl;
	}
	if(!author.empty())
	{
		for(i=0;i<(int)author.size();i++)
			file<<author[i]<<endl;
	}
	if(!revdat.empty())
	{
		for(i=0;i<(int)revdat.size();i++)
			file<<revdat[i]<<endl;
	}
	if(!sprsde.empty())
	{
		for(i=0;i<(int)sprsde.size();i++)
			file<<sprsde[i]<<endl;
	}
	if(!jnrl.empty())
	{
		for(i=0;i<(int)jnrl.size();i++)
			file<<jnrl[i]<<endl;
	}
	if(!jrnl.empty())
	{
		for(i=0;i<(int)jrnl.size();i++)
			file<<jrnl[i]<<endl;
	}
	if(!remark.empty())
	{
		for(i=0;i<(int)remark.size();i++)
			file<<remark[i]<<endl;
	}
	if(!dbref.empty())
	{
		for(i=0;i<(int)dbref.size();i++)
			file<<dbref[i]<<endl;
	}
	if(!seqadv.empty())
	{
		for(i=0;i<(int)seqadv.size();i++)
			file<<seqadv[i]<<endl;
	}
	if(!seqres.empty())
	{
		for(i=0;i<(int)seqres.size();i++)
			file<<seqres[i]<<endl;
	}
	if(!modres.empty())
	{
		for(i=0;i<(int)modres.size();i++)
			file<<modres[i]<<endl;
	}
	if(!ftnote.empty())
	{
		for(i=0;i<(int)ftnote.size();i++)
			file<<ftnote[i]<<endl;
	}
	if(!het.empty())
	{
		for(i=0;i<(int)het.size();i++)
			file<<het[i]<<endl;
	}
	if(!hetnam.empty())
	{
		for(i=0;i<(int)hetnam.size();i++)
			file<<hetnam[i]<<endl;
	}
	if(!hetsyn.empty())
	{
		for(i=0;i<(int)hetsyn.size();i++)
			file<<hetsyn[i]<<endl;
	}
	if(!formul.empty())
	{
		for(i=0;i<(int)formul.size();i++)
			file<<formul[i]<<endl;
	}
	if(!helix.empty())
	{
		for(i=0;i<(int)helix.size();i++)
			file<<helix[i]<<endl;
	}
	if(!sheet.empty())
	{
		for(i=0;i<(int)sheet.size();i++)
			file<<sheet[i]<<endl;
	}
	if(!turn.empty())
	{
		for(i=0;i<(int)turn.size();i++)
			file<<turn[i]<<endl;
	}
	if(!ssbond.empty())
	{
		for(i=0;i<(int)ssbond.size();i++)
			file<<ssbond[i]<<endl;
	}
	if(!hydbnd.empty())
	{
		for(i=0;i<(int)hydbnd.size();i++)
			file<<hydbnd[i]<<endl;
	}
	if(!sltbrg.empty())
		{
		for(i=0;i<(int)sltbrg.size();i++)
			file<<sltbrg[i]<<endl;
	}	
	if(!link.empty())
	{
		for(i=0;i<(int)link.size();i++)
			file<<link[i]<<endl;
	}
	if(!cispep.empty())
	{
		for(i=0;i<(int)cispep.size();i++)
			file<<cispep[i]<<endl;
	}
	if(!site.empty())
	{
		for(i=0;i<(int)site.size();i++)
			file<<site[i]<<endl;
	}
	if(!cryst1.empty())
	{
		for(i=0;i<(int)cryst1.size();i++)
			file<<cryst1[i]<<endl;
	}
	if(!origx1.empty())
	{
		for(i=0;i<(int)origx1.size();i++)
			file<<origx1[i]<<endl;
	}
	if(!origx2.empty())
	{
		for(i=0;i<(int)origx2.size();i++)
			file<<origx2[i]<<endl;
	}
	if(!origx3.empty())
	{
		for(i=0;i<(int)origx3.size();i++)
			file<<origx3[i]<<endl;
	}
	if(!scale1.empty())
	{
		for(i=0;i<(int)scale1.size();i++)
			file<<scale1[i]<<endl;
	}
	if(!scale2.empty())
	{
		for(i=0;i<(int)scale2.size();i++)
			file<<scale2[i]<<endl;
	}
	if(!scale3.empty())
	{
		for(i=0;i<(int)scale3.size();i++)
			file<<scale3[i]<<endl;
	}
	if(!mtrix1.empty())
	{
		for(i=0;i<(int)mtrix1.size();i++)
			file<<mtrix1[i]<<endl;
	}
	if(!mtrix2.empty())
	{
		for(i=0;i<(int)mtrix2.size();i++)
			file<<mtrix2[i]<<endl;
	}
	if(!mtrix3.empty())
	{
		for(i=0;i<(int)mtrix3.size();i++)
			file<<mtrix3[i]<<endl;
	}
	if(!tvect.empty())
	{
		for(i=0;i<(int)tvect.size();i++)
			file<<tvect[i]<<endl;
	}
}
//anisou暂时未处理
void Protein::OutputSucField(ofstream& file)
{
	int i;
	if(!hetatm.empty())
	{
			for(i=0;i<(int)hetatm.size();i++)
				file<<hetatm[i]<<endl;
	}
	if(!conect.empty())
	{
			for(i=0;i<(int)conect.size();i++)
				file<<conect[i]<<endl;
	}
	if(!master.empty())
	{
			for(i=0;i<(int)master.size();i++)
				file<<master[i]<<endl;
	}
	if(!end.empty())
	{
		for(i=0;i<(int)end.size();i++)
				file<<end[i]<<endl;
	}
	else
		file<<"END                                                                             ";

			
}
/**************************************************
 * 将每个链的信息分别存入文件中
 * 如果filename为默认值，则文件名为pdb标识+链标识.ent
 * 如果文件名不为默认值，则将链标识插入到最后一个.之前
 * 当当前的蛋白质中仅有一个链，则文件名中不加链标识
 * 如果链标识为空格，则在需要表示时，用A表示
 ***************************************************/
bool Protein::SaveChains(string dirname,string filename)
{
	
	string line,fname,tempstr;
	Atom atom;
	int j,po;
	int num_atom;
	int chain_num;
	if(num_chain<1)
		return false;

	for(chain_num=0;chain_num<num_chain;chain_num++)
	{
		//生成文件名
		if(filename!="AUTO")
		{
			fname=filename;
			if(chain_num>1)
			{
				if(m_chains[chain_num].chain_identifier!=' ')
				{
					po=(int)fname.find('.');
					tempstr=m_chains[chain_num].chain_identifier;
					fname.insert(po,tempstr);
				}
				else
				{
					po=(int)fname.find('.');					
					fname.insert(po,"A");
				}
			}
		}
		else
		{
			fname=PDB_id;
			if(num_chain>1)
			{
				if(m_chains[chain_num].chain_identifier!=' ')
				{					
					fname+=m_chains[chain_num].chain_identifier;
				}
				else
					fname+='A';
			}
			fname+=".ent";
		}
		//添加目录标识
		if(dirname!="AUTO")
		{
			_mkdir(dirname.c_str());
			fname=dirname+"\\"+fname;
		}
		ofstream file;		
		file.open(fname.c_str(),ios::out);
		if(!file.is_open())
			return false;
		/*line="HEADER";
		line.resize(80);
		for(i=0;i<(int)PDB_id.length();i++)
			line[62+i]=PDB_id[i];
		file<<line<<endl;*/
		
		OutputPreField(file);		
		num_atom=m_chains[chain_num].num_atom;
		for(j=0;j<num_atom;j++)
		{
			atom=m_chains[chain_num].m_atoms[j];
			Atom2Line(atom,line);
			file<<line<<endl;
		}
		Atom2Ter(atom,line);
		file<<line<<endl;
		OutputSucField(file);
	}
	return true;
}
/**********************************************
 * 统计氨基酸对骨架扭转角的倾向性
 * 骨架扭转角必须以计算
 * 这里骨架扭转角分成36*36份
 * 成功时返回true，否则返回false
 *********************************************/
bool Chain::Pre_Dihedral(vector<vector<double> >& preference)
{
	int len,i,bin,bpsi,bphi,aa;
	double phi,psi;
	len=(int)m_residues.size();
	for(i=0;i<len;i++)
	{
		//计算骨架扭转角类型
		phi=m_residues[i].phi;
		psi=m_residues[i].psi;
		if(fabs(phi+999)<1e-5)//两面角不存在则不统计
			continue;
		if(fabs(psi+999)<1e-5)
			continue;
		bphi=(int)((phi+180)/10);
		if(bphi>=36)
			bphi=36;
		bpsi=(int)((psi+180)/10);
		if(bpsi>=36)
			bpsi=36;
		bin=bphi*36+bpsi;
		if(bin>1295)
			return false;
		//计算氨基酸序号
		aa=Rname2ID(m_residues[i].c_name);		
		preference[bin][aa]+=1.0;
	}
	return true;
}
//计算每个氨基酸的接触个数
bool Chain::Cal_Num_Contact()
{
	int i,j,len;
	int index,atom_index;
	Point pi,pj;
	double distance;
	len=num_residue;

	for(i=0;i<len;i++)
	{
		index=Atomname2ID("CA");
		atom_index=m_residues[i].atom_index[index];
		if(atom_index==-1)
		{
			cout<<"计算接触个数时出错：在氨基酸"<<i+1<<"中找不到CA原子"<<endl;
			return false;
		}
		pi=m_atoms[atom_index].pt;
		for(j=i+1;j<len;j++)
		{
			index=Atomname2ID("CA");
			atom_index=m_residues[j].atom_index[index];
			if(atom_index==-1)
			{
				cout<<"计算接触个数时出错：在氨基酸"<<j+1<<"中找不到CA原子"<<endl;
				return false;
			}
			pj=m_atoms[atom_index].pt;
			distance=Distance(pi,pj);
			if(distance<=8)//距离小于8埃，定义为存在contact
			{
				m_residues[i].num_contact++;
				m_residues[j].num_contact++;
			}
		}
	}
	return true;
}
//统计氨基酸对接触个数类型的出现次数，接触个数必须已经计算
bool Chain::Pre_Num_Contact(vector<vector<double> >& preference)
{
	int i,len;
	int bin,aa;
	len=num_residue;

	for(i=0;i<len;i++)
	{
		bin=m_residues[i].num_contact;
		if(bin<1)
			bin=1;
		else if(bin>25)
			bin=25;
		bin-=1;//转变成序号
		aa=Rname2ID(m_residues[i].c_name);	
		preference[bin][aa]+=1.0;
	}
	return true;
}
bool IsDigit(string str)
{
	int i,len;
	len=(int)str.size();
	for(i=0;i<len;i++)
	{
		if((!isdigit(str[i]))&&(str[i]!=' ')&&(str[i]!='-'))
			return false;
	}
	return true;
}
/*********************************************************************
 * 累加原子类型1的势能
 * 原子类型1的势能指：原子在序列上的间隔距离k分三中，以11和22为分界点
 * 原子在空间的距离类别从0变化到10，以0.5为间隔，20个
 * 仅考虑主链上的重原子对25对
 * 对应的势能向量的长度为1500
 *********************************************************************/
bool Chain::Pre_Num_Atomtype1(vector<double>& preference)
{
	int i1,i2,len;
	int k,s,daa,po,kk,type1,type2;
	Point p1,p2;
	double distance;
	string atomname1,atomname2;

	len=this->num_atom;
	for(i1=0;i1<len;i1++)
	{
		atomname1=m_atoms[i1].name;
		if((atomname1=="CB")||(atomname1=="N")||(atomname1=="O")||(atomname1=="CA")||(atomname1=="C"))//原子名称在读入时已经转换为大写
		{
			p1=m_atoms[i1].pt;
			for(i2=i1+1;i2<len;i2++)
			{
				atomname2=m_atoms[i2].name;
				if((atomname2=="CB")||(atomname2=="N")||(atomname2=="O")||(atomname2=="CA")||(atomname2=="C"))
				{
					p2=m_atoms[i2].pt;
					distance=Distance(p1,p2);
					if(distance>=10.0)//超过该距离的原子对忽略掉
						continue;
					//确定空间距离编号s
					s=int(distance*2);
					if(s>=20)
						s=19;
					//确定氨基酸间隔编号k
					kk=m_atoms[i1].residule_chain_seqm-m_atoms[i2].residule_chain_seqm;
					if(kk<0)
						kk=-kk;
					if(kk<11)
						k=0;
					else if(kk>=11||kk<=22)
						k=1;
					else
						k=2;
					//确定原子对编号
					if(atomname1=="CA")
						type1=0;
					else if(atomname1=="CB")
						type1=1;
					else if(atomname1=="N")
						type1=2;
					else if(atomname1=="C")
						type1=3;
					else //"O"
						type1=4;
					if(atomname2=="CA")
						type2=0;
					else if(atomname2=="CB")
						type2=1;
					else if(atomname2=="N")
						type2=2;
					else if(atomname2=="C")
						type2=3;
					else //"O"
						type2=4;
					daa=type1*5+type2;
					po=k*500+s*25+daa;
					preference[po]+=1.0;
				}
			}
		}
	}
	return true;
}
/*********************************************************************
 * 累加原子类型2的势能
 * 原子类型1的势能指：原子在序列上的间隔距离k分三中，以11和22为分界点
 * 原子在空间的距离类别从1.5变化到11.5，以1为间隔，10个
 * 原子类型共40个,原子对1600对
 * 对应的势能向量的长度为48000
 *********************************************************************/
bool Chain::Pre_Num_Atomtype2(vector<double>& preference)
{
	int i1,i2,len;
	int k,s,daa,po,kk,type1,type2;
	Point p1,p2;
	double distance;
	string atomname1,atomname2;	
	char c1,c2;
	int aaindex,atomindex;

	len=this->num_atom;
	for(i1=0;i1<len;i1++)
	{
		atomname1=m_atoms[i1].name;
		c1=m_atoms[i1].c_resname;
		aaindex=Rname2ID20(c1);//将原子所在的氨基酸名称转换成20个标准氨基酸序号
		if(aaindex==-1)
		{
			cout<<"警告:出现了无效的氨基酸序号第"<<i1<<"个原子所在的氨基酸名称为"<<c1<<endl;
			cout<<"该原子将被忽略"<<endl;
			continue;
		}
		atomindex=Atomname2ID(atomname1.c_str());//原子名称转换成序号
		if(atomindex==36)//忽略OXT原子
			continue;
		if(atomindex==-1)
		{
			cout<<"警告:出现了无效的原子序号第"<<i1<<"个原子所名称为"<<atomname1.c_str()<<endl;
			cout<<"该原子将被忽略"<<endl;
			continue;
		}
		type1=Atomtype40[aaindex][atomindex];
		if(type1==-1)
		{
			cout<<"警告:未能获得原子的正确序号氨基酸名称"<<c1<<" 序号"<<aaindex<<" 原子名称"<<atomname1.c_str()<<" 序号"<<atomindex<<endl;
			cout<<"该原子将被忽略"<<endl;
			continue;
		}
		//if((atomname1=="CB")||(atomname1=="N")||(atomname1=="O")||(atomname1=="CA")||(atomname1=="C"))//原子名称在读入时已经转换为大写
		{
			p1=m_atoms[i1].pt;
			for(i2=i1+1;i2<len;i2++)
			{
				atomname2=m_atoms[i2].name;
				c2=m_atoms[i2].c_resname;
				aaindex=Rname2ID20(c2);//将原子所在的氨基酸名称转换成20个标准氨基酸序号
				if(aaindex==-1)
				{
					cout<<"警告:出现了无效的氨基酸序号第"<<i2<<"个原子所在的氨基酸名称为"<<c2<<endl;
					cout<<"该原子将被忽略"<<endl;
					continue;
				}
				atomindex=Atomname2ID(atomname2.c_str());//原子名称转换成序号
				if(atomindex==36)//忽略OXT原子
					continue;
				if(atomindex==-1)
				{
					cout<<"警告:出现了无效的原子序号第"<<i2<<"个原子所名称为"<<atomname2.c_str()<<endl;
					cout<<"该原子将被忽略"<<endl;
					continue;
				}
				type2=Atomtype40[aaindex][atomindex];
				if(type2==-1)
				{
					cout<<"警告:未能获得原子的正确序号氨基酸名称"<<c2<<" 序号"<<aaindex<<" 原子名称"<<atomname2.c_str()<<" 序号"<<atomindex<<endl;
					cout<<"该原子将被忽略"<<endl;
					continue;
				}
				//if((atomname2=="CB")||(atomname2=="N")||(atomname2=="O")||(atomname2=="CA")||(atomname2=="C"))
				{
					p2=m_atoms[i2].pt;
					distance=Distance(p1,p2);
					if((distance<1.5)||(distance>=11.5))//超过该距离的原子对忽略掉
						continue;
					//确定空间距离编号s
					s=(int)(distance-1.5);
					if(s>=10)
						s=9;
					//确定氨基酸间隔编号k
					kk=m_atoms[i1].residule_chain_seqm-m_atoms[i2].residule_chain_seqm;					
					if(kk<0)
						kk=-kk;
					if(kk<=1) //相邻或者同一个氨基酸之内的不计算
						continue;
					if(kk<11)
						k=0;
					else if(kk>=11||kk<=22)
						k=1;
					else
						k=2;					
					daa=(type1-1)*40+(type2-1);
					po=k*16000+s*1600+daa;
					preference[po]+=1.0;
				}
			}
		}
	}
	return true;
}
bool Chain::Pre_FS(vector<vector<vector< double> > >& preference)
{
	int i1,i2,len;
	int s,type1,type2;
	Point p1,p2;
	double distance,r1,r2;
	string atomname1,atomname2;	
	char c1,c2;
	int aaindex,atomindex;
	int residule_chain_seqm1,residule_chain_seqm2;

	len=this->num_atom;
	for(i1=0;i1<len;i1++)
	{
		atomname1=m_atoms[i1].name;
		c1=m_atoms[i1].c_resname;
		aaindex=Rname2ID20(c1);//将原子所在的氨基酸名称转换成20个标准氨基酸序号
		if(aaindex==-1)
		{
			cout<<"警告:出现了无效的氨基酸序号第"<<i1<<"个原子所在的氨基酸名称为"<<c1<<endl;
			cout<<"该原子将被忽略"<<endl;
			continue;
		}
		atomindex=Atomname2ID(atomname1.c_str());//原子名称转换成序号	
		if(atomindex==36)
			continue;
		if(atomindex==-1)
		{
			//cout<<"警告:出现了无效的原子序号第"<<i1<<"个原子所名称为"<<atomname1.c_str()<<endl;
			//cout<<"该原子将被忽略"<<endl;
			continue;
		}
		type1=Atomtype_FS[aaindex][atomindex];
		
		if(type1==-1)
		{
			cout<<"警告:未能获得原子的正确序号氨基酸名称"<<c1<<" 序号"<<aaindex<<" 原子名称"<<atomname1.c_str()<<" 序号"<<atomindex<<endl;
			cout<<"该原子将被忽略"<<endl;
			continue;
		}
		residule_chain_seqm1=m_atoms[i1].residule_chain_seqm;
		if(atomindex<17)//atom C
			r1=1.4;
		else if(atomindex<26)//atom N
			r1=1.2;
		else if(atomindex<37)//atom O
			r1=1.1;
		else //atom S
			r1=1.4;
		//if((atomname1=="CB")||(atomname1=="N")||(atomname1=="O")||(atomname1=="CA")||(atomname1=="C"))//原子名称在读入时已经转换为大写
		{
			p1=m_atoms[i1].pt;
			for(i2=i1+1;i2<len;i2++)
			{
				atomname2=m_atoms[i2].name;
				c2=m_atoms[i2].c_resname;
				aaindex=Rname2ID20(c2);//将原子所在的氨基酸名称转换成20个标准氨基酸序号
				if(aaindex==-1)
				{
					cout<<"警告:出现了无效的氨基酸序号第"<<i2<<"个原子所在的氨基酸名称为"<<c2<<endl;
					cout<<"该原子将被忽略"<<endl;
					continue;
				}
				atomindex=Atomname2ID(atomname2.c_str());//原子名称转换成序号
				if(atomindex==36)
					continue;
				if(atomindex==-1)
				{
					//cout<<"警告:出现了无效的原子序号第"<<i2<<"个原子所名称为"<<atomname2.c_str()<<endl;
					//cout<<"该原子将被忽略"<<endl;
					continue;
				}
				type2=Atomtype_FS[aaindex][atomindex];				
				if(type2==-1)
				{
					cout<<"警告:未能获得原子的正确序号氨基酸名称"<<c2<<" 序号"<<aaindex<<" 原子名称"<<atomname2.c_str()<<" 序号"<<atomindex<<endl;
					cout<<"该原子将被忽略"<<endl;
					continue;
				}
				residule_chain_seqm2=m_atoms[i2].residule_chain_seqm;
				if(atomindex<17)//atom C
					r2=1.4;
				else if(atomindex<26)//atom N
					r2=1.2;
				else if(atomindex<37)//atom O
					r2=1.1;
				else //atom S
					r2=1.4;
                //if((atomname2=="CB")||(atomname2=="N")||(atomname2=="O")||(atomname2=="CA")||(atomname2=="C"))
				{
					if((int)fabs((double)residule_chain_seqm1-residule_chain_seqm2)<4)//原子在序列上的间隔
						continue;
					p2=m_atoms[i2].pt;
					distance=Distance(p1,p2);
					distance-=(r1+r2);
					if((distance<0)||(distance>=4))//超过该距离的原子对忽略掉
						continue;
					//确定空间距离编号s
					s=int(distance/0.2);
					if(s>19)
						s=19;					
					preference[type1][type2][s]+=1.0;
					preference[type2][type1][s]+=1.0;
				}
			}
		}
	}
	return true;
}
bool Protein::CalSAA()
{
    long i,j,k,index;
    double f,fatom;
    long FORLIM;
    Residue *WITH;
    long FORLIM1,FORLIM2;
	double boxmin[3],boxmax[3],v4[4];
    CAS *cas= CASCreate(2);


    AddToAllAtomRadii(RWATER); /* add the water radius :-) */ //to be added
    CalcResidueCenters();//to be added
    FORLIM = num_chain;
    for (i = 0; i < FORLIM; i++) 
	{
		FORLIM1=m_chains[i].num_residue;
		for(j=0;j<FORLIM1;j++)
		{	
			WITH = &m_chains[i].m_residues[j];
			boxmax[0]=WITH->boxmax.x; boxmax[1]=WITH->boxmax.y; boxmax[2]=WITH->boxmax.z;
			boxmin[0]=WITH->boxmin.x; boxmin[1]=WITH->boxmin.y; boxmin[2]=WITH->boxmin.z;
			FindNeighbourRes(boxmin, boxmax);
			f=0;
			FORLIM2=(int)WITH->m_atoms.size();
			for(k=0;k<FORLIM2;k++)
			{
				index=WITH->m_atoms[k];
				Atom2Pointer4(m_chains[i].m_atoms[index],v4);
				fatom=Surface(cas,v4);
				m_chains[i].m_atoms[index].saa=fatom;
				f+=fatom;
			}
			WITH->saa = (long)floor(f + 0.5);
	    }
	    
	}    
    AddToAllAtomRadii(-RWATER);   /* remove the water :-) */
    CASDelete(cas);
	return true;
}
/************************************
 * 所有链的所有原子半径增加r
 * r为负值时减小
 ***********************************/
void Protein::AddToAllAtomRadii(double r)
{
  long i, j, FORLIM;  
  long FORLIM1;

  FORLIM = num_chain;
  for (i = 0; i < FORLIM; i++) 
  {
	  FORLIM1=(int)m_chains[i].m_atoms.size();
	  for(j=0;j<FORLIM1;j++)
	  {
		  m_chains[i].m_atoms[j].radium+=r;		  
	  }
  }    
}  /* AddToAllAtomRadii */
void Protein::CalcResidueCenters()
{
  long i, j, k,index,FORLIM;
  long FORLIM1,FORLIM2;
  double boxmin[3],boxmax[3],v4[4];
  Residue * re;

  FORLIM = num_chain;
  for(i=0;i<FORLIM;i++)
  {
	  FORLIM1=m_chains[i].num_residue;
	  for(j=0;j<FORLIM1;j++)
	  {
		  boxmax[0] = -1e10;
		  boxmax[1] = -1e10;
		  boxmax[2] = -1e10;
		  boxmin[0] = 1e10;
		  boxmin[1] = 1e10;
		  boxmin[2] = 1e10;
		  FORLIM2=(int)m_chains[i].m_residues[j].m_atoms.size();
		  re=&m_chains[i].m_residues[j];
		  for(k=0;k<FORLIM2;k++)
		  {
			  index=re->m_atoms[k];
			  Atom2Pointer4(m_chains[i].m_atoms[index],v4);
			  MinMax(v4,boxmin,boxmax);
		  }
		  re->boxmax.x=boxmax[0];re->boxmax.y=boxmax[1];re->boxmax.z=boxmax[2];
		  re->boxmin.x=boxmin[0];re->boxmin.y=boxmin[1];re->boxmin.z=boxmin[2];		  
	  }
  }  
}  /* CalcResidueCenters */
void Protein::FindNeighbourRes(double *vmin, double *vmax)//PDB
{
  long i,j,FORLIM,FORLIM1;

  LastNeighbourRes = 0;
  FORLIM = num_chain;
  /* gcc (2.5.7 and 2.4.5) have problems on some machines with */
  /* -funroll-loops here.*/
  /* e.g on SGI 1tim produces wrong results!!!! */
  for (i = 0; i < FORLIM; i++) 
  {
	  FORLIM1=m_chains[i].num_residue;
	  for(j=0;j<FORLIM1;j++)
	  {
		if(
		vmin[0] < m_chains[i].m_residues[j].boxmax.x && vmin[1] < m_chains[i].m_residues[j].boxmax.y &&
		vmin[2] < m_chains[i].m_residues[j].boxmax.z && vmax[0] > m_chains[i].m_residues[j].boxmin.x &&
		vmax[1] > m_chains[i].m_residues[j].boxmin.y && vmax[2] > m_chains[i].m_residues[j].boxmin.z) 
		{
			LastNeighbourRes++;	  
			NeighbourRes[LastNeighbourRes].chain = i;
			NeighbourRes[LastNeighbourRes].res  = j;
		}
	  }
  }
}  /*  FindNeighbourRes*/
  void Protein::Liste(CAS* cas, double *atom)//PDB
{
    long i, j,FORLIM;
    Residue *WITH;
    long FORLIM1;
	double boxmin[3],boxmax[3],v4[4];

    FORLIM = LastNeighbourRes;
    for (i = 1; i <= FORLIM; i++) 
	{
		WITH = &m_chains[NeighbourRes[i].chain].m_residues[NeighbourRes[i].res];
		boxmax[0]=WITH->boxmax.x;boxmax[1]=WITH->boxmax.y;boxmax[2]=WITH->boxmax.z;
		boxmin[0]=WITH->boxmin.x;boxmin[1]=WITH->boxmin.y;boxmin[2]=WITH->boxmin.z;
		if (InBox(atom, boxmin, boxmax)) 
		{
			FORLIM1=(int)WITH->m_atoms.size();
			for(j=0;j<FORLIM1;j++)
			{
				Atom2Pointer4(m_chains[NeighbourRes[i].chain].m_atoms[WITH->m_atoms[j]],v4);
				Listentry(cas, atom, v4);
			}			
		}
	}    
}  /* Liste */
double Protein::Surface(CAS* cas,double *xatom)
{
    double f;
    CASResetAtoms(cas);
    CASSetCenterAtom(cas,xatom[0],xatom[1],xatom[2],xatom[3]);
    Liste(cas, xatom);
    f= CASSurface(cas);
    return f;
}  /* Surface */
bool Chain::ComputSidecenter_Te13()
{
	int i,j,indexaa,num_sideatom;
	int atomindex;
	for(i=0;i<num_residue;i++)
	{
		m_residues[i].sidecenter_te13.x=0;
		m_residues[i].sidecenter_te13.y=0;
		m_residues[i].sidecenter_te13.z=0;
		num_sideatom=0;
		indexaa=Rname2ID20(m_residues[i].c_name);
		if(indexaa==-1)
			return false;
		for(j=0;j<NUM_ATOM_TYPE;j++)
		{
			if(SideAtom[indexaa][j])
			{
				atomindex=m_residues[i].atom_index[j];
				if(atomindex!=-1)
				{
					m_residues[i].sidecenter_te13=Add(m_residues[i].sidecenter_te13,m_atoms[atomindex].pt);
					num_sideatom++;
				}
			}
		}		
		if(indexaa==5)//GLY氨基酸用Ca原子代替
		{
			m_residues[i].sidecenter_te13=m_atoms[m_residues[i].atom_index[1]].pt;
		}
		else
		{
			if(num_sideatom==0)
				return false;
			m_residues[i].sidecenter_te13.x/=num_sideatom;
			m_residues[i].sidecenter_te13.y/=num_sideatom;
			m_residues[i].sidecenter_te13.z/=num_sideatom;
		}
	}
	return true;
}
bool Chain::ComputSidecenter_HRSC()
{
	int i,j,indexaa,num_sideatom;
	int atomindex;
	for(i=0;i<num_residue;i++)
	{
		m_residues[i].sidecenter_HRSC.x=0;
		m_residues[i].sidecenter_HRSC.y=0;
		m_residues[i].sidecenter_HRSC.z=0;
		num_sideatom=0;
		indexaa=Rname2ID20(m_residues[i].c_name);
		if(indexaa==-1)
			return false;
		for(j=0;j<NUM_ATOM_TYPE;j++)
		{
			if(SideAtom_HRSC[indexaa][j])
			{
				atomindex=m_residues[i].atom_index[j];
				if(atomindex!=-1)
				{
					m_residues[i].sidecenter_HRSC=Add(m_residues[i].sidecenter_HRSC,m_atoms[atomindex].pt);
					num_sideatom++;
				}
			}
		}
		if(num_sideatom==0)
		return false;
		m_residues[i].sidecenter_HRSC.x/=num_sideatom;
		m_residues[i].sidecenter_HRSC.y/=num_sideatom;
		m_residues[i].sidecenter_HRSC.z/=num_sideatom;		
	}
	return true;		
}
bool Chain::Pre_DFIRESCM(vector<vector<vector<double > > >& preference)
{
	char aa1,aa2;
	int i,j,aaindex1,aaindex2,disindex;
	double distance;
	Point p1,p2;
	
	int len=this->num_residue;	
	for(i=0;i<len;i++)
	{
		aa1=m_residues[i].c_name;
		aaindex1=Rname2ID20(aa1);
		if(aaindex1==-1)
			return false;
		p1=m_residues[i].sidecenter_te13;
		for(j=i+1;j<len;j++)
		{
			aa2=m_residues[j].c_name;
			aaindex2=Rname2ID20(aa2);
			if(aaindex2==-1)
				return false;
			p2=m_residues[j].sidecenter_te13;
			distance=Distance(p1,p2);
			if((distance>=0)&&(distance<=20))
			{
				disindex=Dis2DFIRESCMndex(distance);
				preference[aaindex1][aaindex2][disindex]+=1.0;
			}
		}
	}		
	return true;
}
int Dis2DFIRESCMndex(double dis)
{
	int i,index;
	double center[20]={1,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.5,9.5,10.5,11.5,12.5,13.5,14.5};
	double radius[20]={1,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
	index=0;
	for(i=0;i<20;i++)
	{
		if((dis>=(center[i]-radius[i]))&&(dis<=(center[i]+radius[i])))
			{
				index=i;
				break;
			}
	}
	if(index<0)
		index=0;
	if(index>=20)
		index=19;
	return index;
}
