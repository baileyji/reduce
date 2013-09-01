#define MIN(a,b) ((a<b)?a:b)
#define MAX(a,b) ((a>b)?a:b)

int bandsol(int Argc, void *Argv[])
{
  double *a, *r, aa;
  int n, nd, i, j, k;

/* Double arrays are passed by reference. */
  a = (double *)Argv[0];
  r = (double *)Argv[1];
/* The size of the system and the number of diagonals are passed by value */
  n = *(int *)Argv[2];
  nd = *(int *)Argv[3];

/*
   bandsol solve a sparse system of linear equations with band-diagonal matrix.
   Band is assumed to be symmetrix relative to the main diaginal. Usage:
   CALL_EXTERNAL('bandsol.so', 'bandsol', a, r, n, nd)
   where a is 2D array [n,m] where n - is the number of equations and nd
           is the width of the band (3 for tri-diagonal system),
           nd is always an odd number. The main diagonal should be in a(*,nd/2)
           The first lower subdiagonal should be in a(1:n-1,nd-2-1), the first
           upper subdiagonal is in a(0:n-2,nd/2+1) etc. For example:
                  / 0 0 X X X \
                  | 0 X X X X |
                  | X X X X X |
                  | X X X X X |
              A = | X X X X X |
                  | X X X X X |
                  | X X X X X |
                  | X X X X 0 |
                  \ X X X 0 0 /
         r is the array of RHS of size n.
*/

/* Forward sweep */
  for(i=0; i<n-1; i++)
  {
    aa=a[i+n*(nd/2)];
    r[i]/=aa;
    for(j=0; j<nd; j++) a[i+j*n]/=aa;
    for(j=1; j<MIN(nd/2+1,n-i); j++)
    {
      aa=a[i+j+n*(nd/2-j)];
      r[i+j]-=r[i]*aa;
      for(k=0; k<n*(nd-j); k+=n) a[i+j+k]-=a[i+k+n*j]*aa;
    }
  }

/* Backward sweep */
  r[n-1]/=a[n-1+n*(nd/2)];
  for(i=n-1; i>0; i--)
  {
    for(j=1; j<=MIN(nd/2,i); j++) r[i-j]-=r[i]*a[i-j+n*(nd/2+j)];
    r[i-1]/=a[i-1+n*(nd/2)];
  }

  r[0]/=a[n*(nd/2)];

  return 0;
}
