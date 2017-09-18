#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

vector < vector <double> > parseMatrix(string str)
{
	istringstream iss(str);

	unsigned m;
	unsigned n;
	iss >> m;
	iss >> n;

	vector < vector <double> > matrix;
	matrix.resize(m);
	for(unsigned i = 0; i < m; ++i)
	{
		matrix[i].resize(n);
	}

	unsigned i = 0;
	double tmp;
    iss >> tmp;
	while (iss)
    {
        matrix[i/n][i % n] = tmp;
        ++i;
        iss >> tmp;
    };

    return matrix;
}

template <typename T>
void showMat(vector < vector <T> > mat)
{
	for(unsigned i = 0; i < mat.size(); ++i)
	{
		for(unsigned j = 0; j < mat[i].size(); ++j)
		{
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
}

template <typename T>
void showMatLine(vector < vector <T> > mat)
{
	cout << mat.size() << " " << mat[0].size();
	for(unsigned i = 0; i < mat.size(); ++i)
	{
		for(unsigned j = 0; j < mat[i].size(); ++j)
		{
			cout << " " << mat[i][j];
		}
	}
	cout << endl;
}

template <typename T>
vector < vector <T> > multMatrix(vector < vector <T> > pA, vector < vector <T> > pB)
{
    vector < vector <T> > mult;
    mult.resize(pA.size());
    for(unsigned i = 0; i < mult.size(); ++i)
    {
        mult[i] = vector<T> (pB[0].size(), 0);
    }

    for(unsigned i = 0; i < pA.size(); ++i)
    {
        for(unsigned j = 0; j < pB[0].size(); ++j)
        {
            for(unsigned k = 0; k < pB.size(); ++k)
            {
                mult[i][j] += pA[i][k] * pB[k][j];
            }
        }
    }
    return mult;
}

vector < vector <double> > emissionProba(vector < vector <double> > pA, vector < vector <double> > pB, vector < vector <double> > pPi)
{
	return multMatrix(multMatrix(pPi, pA), pB);
}

int main(int argc,char **argv)
{
	vector < vector <double> > A;
	vector < vector <double> > B;
	vector < vector <double> > Pi;

	string line;

	getline (cin,line);
	A = parseMatrix(line);
	/*showMat(A);
	cout << endl;*/

	getline (cin,line);
	B = parseMatrix(line);
	/*showMat(B);
	cout << endl;*/

	getline (cin,line);
	Pi = parseMatrix(line);
	/*showMat(Pi);
	cout << endl;*/

	vector < vector <double> > result = emissionProba(A,B,Pi);

	//showMat(result);
	showMatLine(result);

	return 0;
}