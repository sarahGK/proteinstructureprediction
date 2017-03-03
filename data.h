#define LA 1296
const char rname_one_rank20[21]="ACDEFGHIKLMNPQRSTVWY";

class data{
private:
	double sa[LA][20];//p(si|ai-1)
	double sb[LA][20];//p(si|ai)
	double sc[LA][20];//p(si|ai+1)
	//double sd[LA][LA];//p(si|si-1)
	//double se[LA][LA];//p(si|si+1)
	vector<vector<double>> sd;//p(si|si-1)
	vector<vector<double>> se;//p(si|si+1)
	double sf[LA][3];//p(si|chi-1)
	double sg[LA][3];//p(si|chi)
	double sh[LA][3];//p(si|chi+1)
	bool mark[8];//false为条件不可取，true为可取
public:
	data();
	void set_data(int s,int a,int b,int c,int d,int e, int f, int g, int h);//统计出现的个数
	//void set_data(int s, int f, int g, int h);
	void set_data(double **dataarray,int length,int no);//计算概率值
	void set_data(vector<vector<double>> & de,int no);//计算概率值
	void readfile(ifstream &f1,ifstream &f2);//文件中读取数据
	void print_data(ostream & of);//输出mark

};