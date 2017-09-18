#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <numeric>

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



vector <unsigned> parseVector(string str)
{
	istringstream iss(str);

	unsigned m;
	iss >> m;

	vector <unsigned> vect;
	vect.resize(m);

	for(unsigned i = 0; i < vect.size(); ++i)
	{
		unsigned tmp;
    	iss >> tmp;
    	vect[i] = tmp;
	}

    return vect;
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
void showVect(vector <T> vect)
{
	for(unsigned i = 0; i < vect.size(); ++i)
	{
		cout << vect[i] << " ";
	}
	cout << endl;
}

template <typename T>
vector < vector <T> > transposeMat(vector < vector <T> > mat)
{
	vector < vector <T> > result;
	result.resize(mat[0].size());
	for(unsigned i = 0; i < result.size(); ++i)
	{
		result[i].resize(mat.size());
	}
	for(unsigned i = 0; i < mat.size(); ++i)
	{
		for(unsigned j = 0; j < mat[i].size(); ++j)
		{
			result[j][i] =  mat[i][j];
		}
	}
	return result;
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

template <typename T>
vector<T> multVector(vector<T> pA, vector<T> pB)
{
    vector<T> mult;
    mult.resize(pA.size());
    for(unsigned i = 0; i < mult.size(); ++i)
    {
        mult[i] = pA[i]*pB[i];
    }

    return mult;
}

template <typename T>
vector <T> multVectMatrix(vector <T> pA, vector < vector <T> > pB)
{
    vector <T> mult(pB[0].size(), 0);

    for(unsigned i = 0; i < pB[0].size(); ++i)
    {
        for(unsigned j = 0; j < pB.size(); ++j)
        {
            mult[i] += pA[j] * pB[j][i];
        }
    }

    return mult;
}

vector < vector <double> > emissionProba(vector < vector <double> > pA, vector < vector <double> > pB, vector < vector <double> > pPi)
{
	return multMatrix(multMatrix(pPi, pA), pB);
}

vector< vector<double> > alpha(vector < vector <double> > pA, vector < vector <double> > pBT, vector < vector <double> > pPi, vector<unsigned> pO)
{
	vector < vector<double> > result;
	result.resize(pO.size()); 
	result[0] =  multVector(pPi[0], pBT[pO[0]]);
	for(unsigned i = 1; i < pO.size(); ++i)
	{
		result[i] = multVector(multVectMatrix(result[i-1], pA), pBT[pO[i]]);
	}
	double sum = 0;
	for(unsigned i = 0; i < result[0].size(); ++i)
	{
		sum += result[result.size()-1][i];
	}
	//cout << sum << endl;
	return result;
}

int main(int argc,char **argv)
{
	vector < vector <double> > A;
	vector < vector <double> > B;
	vector < vector <double> > Pi;
	vector <unsigned> O;

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

	getline (cin,line);
	O = parseVector(line);
	//showVect(O);

	vector< vector<double> > tmp = alpha(A, transposeMat(B), Pi, O);

	cout << accumulate(tmp[tmp.size()-1].begin(), tmp[tmp.size()-1].end(), 0.0) << endl;


	return 0;
}