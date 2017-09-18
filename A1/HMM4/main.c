#include <math.h>
#include <stdio.h>
#include <stdlib.h>

short unsigned T, M, N;

double** parseMatrix(short unsigned m, short unsigned n)
{
	double **matrix = (double**)malloc(m*sizeof(double));
	for(short unsigned i = 0; i < m; ++i)
	{
		matrix[i] = (double*)malloc(n*sizeof(double));
		for(short unsigned j = 0; j < n; ++j)
		{
			scanf("%lf", &matrix[i][j]);
		}
	}

    return matrix;
}



short unsigned* parseVector(short unsigned m)
{
	short unsigned* vect = (short unsigned*)malloc(m*sizeof(short unsigned));

	for(short unsigned i = 0; i < m; ++i)
	{
    	scanf("%hu",&vect[i]);
	}

    return vect;
}

void showMatLine(double **mat, short unsigned m, short unsigned n)
{
	printf("%u %u", m, n);
	for(short unsigned i = 0; i < m; ++i)
	{
		for(short unsigned j = 0; j < n; ++j)
		{
			printf(" %f", mat[i][j]);
		}
	}
	printf("\n");
}


void alpha(double **result, double **pA, double **pB, double **pPi, short unsigned* pO, double*pC)
{
	for(short unsigned i = 0; i < N; ++i)
	{
		result[0][i] = pPi[0][i]*pB[i][pO[0]];
	}

	pC[0] = 0;
	for(short unsigned i = 0; i < N; ++i)
	{
		pC[0] += result[0][i];
	}
	pC[0] = 1/pC[0];

	for(short unsigned i = 0; i < N; ++i)
	{
		result[0][i] *= pC[0];
	}


	for(short unsigned t = 1; t < T; ++t)
	{
		pC[t] = 0;
		for(short unsigned i = 0; i < N; ++i)
		{
			result[t][i] = 0;
			for(short unsigned j = 0; j < N; ++j)
			{
				result[t][i] += result[t-1][j]*pA[j][i];
			}
			result[t][i] *= pB[i][pO[t]];
			pC[t] += result[t][i];
		}
		pC[t] = 1/pC[t];
		for(short unsigned i = 0; i < N; ++i)
		{
			result[t][i] *= pC[t];
		}
	}

}

void beta(double **result, double **pA, double **pB, double **pPi, short unsigned* pO, double* pC)
{
	for(short unsigned i = 0; i < N; ++i)
	{
		result[T-1][i] = pC[T-1];
	}

	for(int t = T-2; t >= 0; --t)
	{
		for(short unsigned i = 0; i < N; ++i)
		{
			result[t][i] = 0;
			for(short unsigned j = 0; j < N; ++j)
			{
				result[t][i] += result[t+1][j]*pB[j][pO[t+1]]*pA[i][j];
			}
			result[t][i] *= pC[t];
		}
	}
}


void baumWelch(double **pA, double **pB, double **pPi, short unsigned* pO)
{
	short unsigned i, j, t;

	double ***mDiGamma;
	double **mAlpha = (double**)malloc(T*sizeof(double));
	for(t = 0; t < T; ++t)
	{
		mAlpha[t] = (double*)malloc(N*sizeof(double));
	}
	double **mBeta = (double**)malloc(T*sizeof(double));
	for(t = 0; t < T; ++t)
	{
		mBeta[t] = (double*)malloc(N*sizeof(double));
	}
	double **mGamma;
	double *mC = (double*)malloc(T*sizeof(double));

	mDiGamma = (double***)malloc(T*sizeof(double));
	for(t = 0; t < T; ++t)
	{
		mDiGamma[t] = (double**)malloc(N*sizeof(double));
		for(j = 0; j < N; ++j)
		{
			mDiGamma[t][j] = (double*)malloc(N*sizeof(double));
		} 
	}

	mGamma = (double**)malloc(T*sizeof(double));
	for(t = 0; t < T; ++t)
	{
		mGamma[t] = (double*)malloc(N*sizeof(double));
		for(i = 0; i < N; ++i)
		{
			mGamma[t][i] = 0;
		}
	}

	
	double logProb;
	double oldLogProg = -100000;

	short unsigned idx = 0;
	for(;;)
	{
		alpha(mAlpha, pA, pB, pPi, pO, mC);

		beta(mBeta, pA, pB, pPi, pO, mC);

		double denom;
		for(t = 0; t < T-1; ++t)
		{
			denom = 0;
			for(i = 0; i < N; ++i)
			{
				for(j = 0; j < N; ++j)	
				{
					denom += mAlpha[t][i] * pA[i][j] * pB[j][pO[t+1]] * mBeta[t+1][j];
				}
			}
			for(i = 0; i < N; ++i)
			{
				mGamma[t][i] = 0;
				for(j = 0; j < N; ++j)
				{
					mDiGamma[t][i][j] = mAlpha[t][i]*pA[i][j]*pB[j][pO[t+1]]*mBeta[t+1][j]/denom;
					mGamma[t][i] += mDiGamma[t][i][j];
				}
			}
		}


		denom = 0;
		for(i = 0; i < N; ++i)
		{
			denom += mAlpha[N-1][i];
		}
		for(i = 0; i < N; ++i)
		{
			mGamma[T-1][i] = mAlpha[T-1][i]/denom;
		}


		//Re-estimate Pi
		pPi[0] = mGamma[0];

		double sumDiGamma = 0;
		double sumGamma = 0;
		//Re-estimate A
		for(i = 0; i < N; ++i)
		{
			for(j = 0; j < N; ++j)
			{
				sumDiGamma = 0;
				sumGamma = 0;
				for(t = 0; t < T-1; ++t)
				{
					sumDiGamma += mDiGamma[t][i][j];
					sumGamma += mGamma[t][i];
				}
				pA[i][j] = sumDiGamma/sumGamma;
			} 
		}


		double sum1;
		double sum2;
		//Re-estimate B
		for(i = 0; i < N; ++i)
		{
			for(j = 0; j < M; ++j)
			{
				sum1 = 0;
				sum2 = 0;
				for(t = 0; t < T; ++t)
				{
					if(pO[t] == j)
						sum1 += mGamma[t][i];
					sum2 += mGamma[t][i];
				}
				pB[i][j] = sum1/sum2;
			}
		}
		
		//Compute log[P(O|y)]
		logProb = 0;
		for(i = 0; i < T; ++i)
		{
			logProb += log(mC[i]);
		}
		logProb = -logProb;

		++idx;

		if(idx >= 30 || logProb <= oldLogProg)
			break;

		oldLogProg = logProb;
	}

	showMatLine(pA, N, N);
	showMatLine(pB, N, M);

}

int main(int argc,char **argv)
{
	double **A;
	double **B;
	double **Pi;
	short unsigned *O;
	short unsigned tmp;

	scanf("%hu", &N);
	scanf("%hu", &tmp);
	A = parseMatrix(N, N);
	/*showMat(A);
	cout << endl;*/

	scanf("%hu", &tmp);
	scanf("%hu", &M);
	B = parseMatrix(N, M);
	/*showMat(B);
	cout << endl;*/

	scanf("%hu", &tmp);
	scanf("%hu", &tmp);
	Pi = parseMatrix(1, N);
	/*showMat(Pi);
	cout << endl;*/

	scanf("%hu", &T);
	O = parseVector(T);
	//showVect(O);

	baumWelch(A, B, Pi, O);

	return 0;
}