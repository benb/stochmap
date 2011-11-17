#undef complex
struct complex
	{
	MrBFlt re;
	MrBFlt im;
	};

typedef struct complex complex;

void    AddTwoMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
int     AreDoublesEqual (MrBFlt x, MrBFlt y, MrBFlt tol);
MrBFlt  **AllocateDoubleMatrix (int r, int c);
MrBFlt  *AllocateDoubleVector (int dim);
int **AllocateIntegerMatrix (int r, int c);
int     *AllocateIntegerVector (int dim);
complex **AllocateSquareComplexMatrix (int dim);
MrBFlt  **AllocateSquareDoubleMatrix (int dim);
int     **AllocateSquareIntegerMatrix (int dim);
int     AutodGamma (MrBFlt *M, MrBFlt rho, int K);
void    BetaBreaks (MrBFlt alpha, MrBFlt beta, MrBFlt *values, int K);
void    CalcCijk (int dim, MrBFlt *c_ijk, MrBFlt **u, MrBFlt **v);
void    CopyComplexMatrices (int dim, complex **from, complex **to);
void    CopyDoubleMatrices (int dim, MrBFlt **from, MrBFlt **to);
void    DirichletRandomVariable (MrBFlt *alp, MrBFlt *z, int n, long int *seed);
int     DiscreteGamma (MrBFlt *rK, MrBFlt alfa, MrBFlt beta, int K, int median);
void    FreeDoubleMatrix (MrBFlt **m);
void    FreeIntegerMatrix (int **m);
void    FreeSquareComplexMatrix (complex **m);
void    FreeSquareDoubleMatrix (MrBFlt **m);
void    FreeSquareIntegerMatrix (int **m);
int     GetEigens (int dim, MrBFlt **q, MrBFlt *eigenValues, MrBFlt *eigvalsImag, MrBFlt **eigvecs, MrBFlt **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs);
MrBFlt  LnGamma (MrBFlt alp);
void    MrBayesPrint (char *format, ...);
void    MultiplyMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
void    MultiplyMatrixByScalar (int dim, MrBFlt **a, MrBFlt scalar, MrBFlt **result);
int     MultiplyMatrixNTimes (int dim, MrBFlt **Mat, int power, MrBFlt **Result);
void    PrintDoubleVector (int dim, MrBFlt *vec);
void    PrintSquareDoubleMatrix (int dim, MrBFlt **matrix);
void    PrintSquareIntegerMatrix (int dim, int **matrix);
MrBFlt  QuantileGamma (MrBFlt x, MrBFlt alfa, MrBFlt beta);
MrBFlt  RandomNumber (long int *seed);

