#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef double (*sistema)();
#define N 4

double f0(double *x)
{
  return 4*x[0] - x[1] + x[2] - x[0]*x[3];
  
}

double f1(double *x)
{
  return -x[0] + 3*x[1] - 2*x[2] - x[1]*x[3];

}

double f2(double *x)
{
  return x[0] - 2*x[1] + 3*x[2] - x[2]*x[3];
  
}

double f3(double *x)
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] -1;
}

double* LeVetor(char *nome,int *m)
{
  FILE *fp;
  int i;
  double *v;
  
  fp=fopen(nome,"r");
  
  fscanf(fp,"%d",m);

  v=malloc(sizeof(double) * (*m));
  
  for(i=0;i<*m; i++)
  {
    fscanf(fp,"%lf",&v[i]);
  }
  
  return v;
}

void ImprimeVetor(double *V, int m)
{
  int i;
  for(i=0; i<m; i++) printf("%g\t",V[i]);
  //puts("");
}

double NormaVetor(double *v,int m,double p)
{
  double norma=0;
  int i;
  
  if(p==0) 
  {
    norma=fabs(v[0]);
    for(i=1;i<m;i++)
    {
      if(fabs(v[i])>norma) norma=fabs(v[i]);
    }
  }
  else
  {
    for(i=0; i<m; i++)
    {
      norma+=pow(fabs(v[i]), p);
    }
    norma=pow(norma, 1/p);
  }
  
  return norma;
}

void TrocaLinhaMestre(double **k, double **l)
{
  double *temp;
  
  temp=*k;
  *k=*l;
  *l=temp;
}

void LU(double **M, double ***L, double ***U, double **B, int m,int n)
{
  int i,j,k,p;
  double l,pivot;

  *L=calloc(m,sizeof(double *));
  for(i=0; i<m; i++) (*L)[i]=calloc(m, sizeof( double *));

  *U=malloc(m*sizeof(double *));
  for(i=0; i<m; i++) (*U)[i]=malloc(m* sizeof( double *));

  *B=malloc(m*sizeof(double *));
  for(i=0; i<m; i++) (*B)[i]=M[i][n-1];

  for(i=0; i<m; i++)
  {
    (*L)[i][i]=1;
    for(j=0; j<n; j++)
      {
        (*U)[i][j]=M[i][j];
      }
   }

  for(j=0; j<m-1; j++)
  {
    pivot=(*U)[j][j];
    p=0;

    for(i=j+1; i<m; i++)
    {
      if(fabs((*U)[i][j])>fabs(pivot))
      {
        pivot= fabs((*U)[i][j]);
        p=i;
      }    
    }

    if(p!=0)
    {
      TrocaLinhaMestre(&(*U)[p], &(*U)[j]);
      l=(*B)[p];
      (*B)[p]=(*B)[j];
      (*B)[j]=l; 
    }

    for(i=j+1; i<m; i++)
    {
      l= (*L)[i][j]= (*U)[i][j]/(*U)[j][j];
      for(k=j; k<m; k++) (*U)[i][k]-= l*(*U)[j][k]; 
    }
  }
}

double *SubsRev(double **M, int m, double *B)
{
  double *S,t;
  int i,j;
  
  S= malloc(m*sizeof(double));
  
  for(i=m-1;i>=0;i--)
  {
    t=0;
    for(j=i+1;j<=m-1;j++) 
    {
      t+=M[i][j]*S[j];
    }
    S[i]=(B[i]-t)/M[i][i];
  }
  
  return S;
}

double *SubsDir(double **M, int m, double *B)
{
  double *S,t;
  int i,j;

  S= malloc(m*sizeof(double));
  
  for(i=0; i<m; i++)
  {
    t=0;
    for(j=(i-1); j>=0; j--) 
    {
      t+=M[i][j]*S[j];
    }
    S[i]=(B[i]-t)/M[i][i];
  }
 
  return S;
}

void Jacobiana(double **J, double *x)
{
  J[0][0]=4-x[3];
  J[0][1]=-1;
  J[0][2]=1;
  J[0][3]=-x[0];
    J[1][0]=-1;
    J[1][1]=3-x[3];
    J[1][2]=-2;
    J[1][3]=-x[1];
      J[2][0]=1;
      J[2][1]=-2;
      J[2][2]=3-x[3];
      J[2][3]=-x[2];
        J[3][0]=2*x[0];
        J[3][1]=2*x[1];
        J[3][2]=2*x[2];
        J[3][3]=0;
}

void Newton(double *x, double tol, sistema f[N])
{
  double **J, **L, **U, *b, norma;
  int i, itera=0;
    
  J = (double **) malloc(N*sizeof(double *));
  for (i=0; i<N; i++) J[i] = (double *) malloc((N+1)*sizeof(double));
  
  do{
    Jacobiana(J,x);
    
    for(i=0; i<N; i++) J[i][N] = -f[i](x);
    
    LU(J, &L, &U, &b, N, N+1);
    b = SubsDir(L, N, b);
    b = SubsRev(U, N, b);
    
    itera++;
    
    norma = NormaVetor(b, N, 2);
    
    for (i=0; i<N; i++) x[i] += b[i];
    
    printf("\n%2d |b|= %g\n\t x= ",itera, norma);
    ImprimeVetor(x, N);
    puts("\n");
  }while (norma>tol);
}
  
int main(int argc, char **argv) 
{
  int it=1,m;
  double *x, tol=1e-10,norma; 
  sistema S[N]={f0,f1,f2,f3};
  
  x=LeVetor(argv[1],&m);
  
  Newton(x, tol, S);
  return 0;
}