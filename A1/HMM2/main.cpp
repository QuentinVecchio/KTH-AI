#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

vector < vector <long double> > parseMatrix(string str)
{
	istringstream iss(str);

	unsigned m;
	unsigned n;
	iss >> m;
	iss >> n;

	vector < vector <long double> > matrix;
	matrix.resize(m);
	for(unsigned i = 0; i < m; ++i)
	{
		matrix[i].resize(n);
	}

	unsigned i = 0;
	long double tmp;
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

vector < vector <long double> > emissionProba(vector < vector <long double> > pA, vector < vector <long double> > pB, vector < vector <long double> > pPi)
{
	return multMatrix(multMatrix(pPi, pA), pB);
}

double alpha(vector < vector <long double> > pA, vector < vector <long double> > pBT, vector < vector <long double> > pPi, vector<unsigned> pO)
{
	vector < vector<long double> > result;
	result.resize(1); 
	result[0] =  multVector(pPi[0], pBT[pO[0]]);
	for(unsigned i = 1; i < pO.size(); ++i)
	{
		result[0] = multVector(multMatrix(result, pA)[0], pBT[pO[i]]);
	}
	long double sum = 0;
	for(unsigned i = 0; i < result[0].size(); ++i)
	{
		sum += result[0][i];
	}
	return sum;
}

void viterbi(vector < vector <long double> > pA, vector < vector <long double> > pB, vector < vector <long double> > pPi, vector<unsigned> pO)
{
	vector < vector <long double> > pAT = transposeMat(pA);
	vector < vector <long double> > pBT = transposeMat(pB);
	vector<unsigned> statesSequence;
	vector< vector<long double> > delta;
	vector< vector<unsigned> > deltaState;
	delta.resize(pO.size());
	deltaState.resize(pO.size());
	for(unsigned i = 0; i < delta.size(); ++i)
	{
		delta[i].resize(pA.size());
		deltaState[i].resize(pA.size());
	}

	delta[0] = multVector(pPi[0], pBT[pO[0]]);

	//showVect(delta[0]);

	unsigned max = distance(delta[0].begin(), max_element(delta[0].begin(), delta[0].end()));

	for(unsigned t = 1; t < delta.size(); ++t)
	{
		for(unsigned i = 0; i < pA.size(); ++i)
		{
			vector<long double> tmp = vector<long double>(pA.size());
			for(unsigned j = 0; j < pA.size(); ++j)
			{
				tmp[j] = delta[t-1][j] * pA[j][i] * pB[i][pO[t]];
			}

			delta[t][i] =  *max_element(tmp.begin(), tmp.end());
			deltaState[t][i] = distance(tmp.begin(), max_element(tmp.begin(), tmp.end()));
		}
		//showVect(delta[t]);
	}

	max = distance(delta[delta.size()-1].begin(), max_element(delta[delta.size()-1].begin(), delta[delta.size()-1].end()));
	statesSequence.resize(delta.size());
	statesSequence[delta.size()-1] = max;

	for(unsigned i = delta.size()-1; i > 0; --i)
	{
		vector<long double> tmp;
		tmp = multVector(delta[i-1], pAT[statesSequence[i]]);
		statesSequence[i-1] = distance(tmp.begin(), max_element(tmp.begin(), tmp.end()));
	}

	for(unsigned i = 0; i < statesSequence.size(); ++i)
	{
		cout << statesSequence[i] << " ";
	}
	cout << endl;

}

int main(int argc,char **argv)
{
	vector < vector <long double> > A;
	vector < vector <long double> > B;
	vector < vector <long double> > Pi;
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

	viterbi(A, B, Pi, O);

	return 0;
}