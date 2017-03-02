#ifndef _RD_PDB_H
#define _RD_PDB_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <conio.h>
#include <math.h>
#include "amino_acid.h"
//#include "direct.h"
#include "io.h"
#include "calcsurf.h"

using namespace std;
// invalid sequence number
#define NOSEQNUM -9999
// invalid chain label
#define NOCHAIN -1
const string NMR="NMR";
const string X_RAY="X-RAY";

// define a contact
typedef struct S_Contact
{
	int i1;
	int i2;
} Contact;
//define a point in 3D
class Point
{
public:
	Point();
public:
	double x,y,z;
};
//define an atom
class Atom
{
public:
	Atom();
	//basic information for an atom
	Point pt;			//coordinate
	string name;			//name
	int chain_seqm;			//sequence number in the chain
	string residule_name;		//name of the residue
	char c_resname;			//one letter code for the residue
	int residule_chain_seqm;	//the sequence number of the residue
	char chain_identifier;		//the chain identifier
	double occupancy;		//occupancy
	double tempFactor;		//tempFactor
	string segID;			//segID
	string element;			//element
	string charge;			//charge

	int residule_seqm;		// the index of the residue starting from 0
	double saa;			// the area of surface accessible 

	double radium;			// the radium of atom
};
//define the residue class
class Residue
{
	friend class Protein;
public:
	Residue();
	void Clear();
public:
	//basic information
	string name;				//name
	char c_name;				//one letter code
	int chain_seqm;				//sequence number of the chain
	int num_atom;				//the number of atoms
	vector<int> m_atoms;		//all the atom IDs
	int atom_index[NUM_ATOM_TYPE];	//the index of the atoms,starting at 0, -1 means atom is not included
	//dihedral
	double alpha,phi,psi,omega,chi_1,chi_2,chi_3,chi_4,chi_5;//alpha is the dihedral formed by 4 consecutive C-alpha atoms

	//contact
	int num_contact;

	//surface accessible area
	double saa;

	double main_radium;			//radium of the main chain
	Point main_center;			//the center of the main chain
	double side_radium;			//readium of side chain
	Point side_center;			//the center of side chain
	Point sidecenter_te13;	//the centriod for the side chain
	Point sidecenter_HRSC;	//the centriod for the side chain used in HRSC potential 
protected:
	
	//to calculate the saa
	Point boxmin,boxmax;//the box containing the amino acid

};
class Chain
{
public:
	Chain();
	void Clear();
	int AtomIndexInAA(int aa_num,string atom_name);//

	//contact
	// get all the contacts
	int ObtainContact(vector<Contact>& all_contact);
	//save all the contacts into a file
	int SaveContact(const string& filename, const vector<Contact>& all_contact);
	//dihedral
	//calculate all the dihedral for each residue
	bool CalculateDihedral();
	//write all the dihedral into a file
	bool WriteDihedral(const string filename);
	// statistically calculate the backbone torsion angle/residue--frequences
	bool Pre_Dihedral(vector<vector<double> >& preference);//Í³¼Æ°±»ùËá¶Ô¹Ç¼ÜÅ¤×ª½ÇµÄ³öÏÖ´ÎÊý

	//nubmer of contacts
	//calculate the number of contacts for each residue
	bool Cal_Num_Contact();
	//calculate the frequences of contacts/residue
	bool Pre_Num_Contact(vector<vector<double> >& preference);
	//potential
	//the total number of appearance of atom type 1 in the terms of potential
	bool Pre_Num_Atomtype1(vector<double >& preference);
	// the total number of appearance of atom type 2  in the terms of potential
	bool Pre_Num_Atomtype2(vector<double >& preference);
	bool Pre_FS(vector<vector<vector< double> > >& preference);
	bool Pre_DFIRESCM(vector<vector<vector<double > > >& preference);

	// rebuild the coordinates based on dihedral
	bool Rebuildxyz_Dihedral();¨
public:
	//basic information
	char chain_identifier;		//the chain identifier
	vector<Residue> m_residues;	// all the residues
	vector<Atom> m_atoms;		// all the atoms
	int num_atom;				//number of atoms
	int num_residue;			//number of residues
	string seq;					//the sequence

	//¸½¼ÓÐÅÏ¢
	bool ComputSidecenter_Te13();
	bool ComputSidecenter_HRSC();
};
class Protein
{
public:
	Protein();
	// clear all the chains
	void Clear();
	// read from PDB files
	bool LoadFromPDB(string filename);
	//save coordinates into files
	bool SaveToFile(string filename);
	//save different chains into different files
	bool SaveChains(string dirname="AUTO",string filename="AUTO");
	//save sequence into files with format of Fasta
	bool WriteSeqFasta(string filename);
	//save all the sequence into the same file
	bool WriteSeqFasta_total(string filename);
	// get the total number of chains
	int Get_Num_chain();
	// get the PDB ID
	void Get_PDBID(string& id);
	// get the sequence of the chain
	bool GetSeq(int chain,string& seq);

	//contact
	//get the contact for the chain
	int ObtainContact(int chainnumber,vector<Contact>& all_contact);
	//save contact into files
	int SaveContact(int chainnumber,const string& filename, const vector<Contact >& all_contact);	//½«µÚi¸öÁ´µÄcontactÐÅÏ¢Ð´ÈëÎÄ¼þ

	//calculate solvent-accessible surface area
	bool CalSAA();

private:
	
	void process_header(const string& line);
	void process_source(const string& line);
	void process_expdta(const string& line);
	void process_remark(const string& line);
	void process_atom(const string& line,vector<Atom>& all_atom);
	void process_field(const string& line);
	bool AddAtom(const vector<Atom>& all_atom);
	void Atom2Line(const Atom& atom,string& line);
	void Atom2Ter(const Atom& atom,string& line);

	void OutputPreField(ofstream& file);
	void OutputSucField(ofstream& file);

	void AddToAllAtomRadii(double r);
	void CalcResidueCenters();
	void FindNeighbourRes(double *vmin, double *vmax);
	void Liste(CAS* cas, double *atom);
	double Surface(CAS* cas,double *xatom);
	


public:
	string PDB_id;
	vector<Chain> m_chains;
	int num_chain;
	string method;
	double resolution;
// to store all different fields of data in PDB files
	vector<string> header;
	vector<string> title;
	vector<string> compnd;
	vector<string> source;
	vector<string> author;
	vector<string> revdat;
	vector<string> jrnl;
	vector<string> remark;
	vector<string> dbref;
	vector<string> seqadv;
	vector<string> seqres;
	vector<string> het;
	vector<string> hetnam;
	vector<string> formul;
	vector<string> helix;
	vector<string> sheet;
	vector<string> link;
	vector<string> cispep;
	vector<string> origx1;
	vector<string> origx2;
	vector<string> origx3;
	vector<string> scale1;
	vector<string> scale2;
	vector<string> scale3;
	vector<string> hetatm;
	vector<string> concet;
	vector<string> master;
	vector<string> anisou;
	vector<string> cryst1;
	vector<string> expdta;
	vector<string> ftnote;
	vector<string> hetsyn;
	vector<string> hydbnd;
	vector<string> jnrl;
	vector<string> keywds;
	vector<string> modres;
	vector<string> mtrix1;
	vector<string> mtrix2;
	vector<string> mtrix3;	
	vector<string> site;
	vector<string> sltbrg;
	vector<string> sprsde;
	vector<string> ssbond;
	vector<string> turn;
	vector<string> tvect;
	vector<string> conect;
	vector<string> end;
};

bool IsDigit(string str);
void Atom2Pointer4(const Atom& atom,double *v);
int Dis2DFIRESCMndex(double dis);
#endif
