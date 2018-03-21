#include "Solvers.h"


void SolveByScalarRun(int N,double *u_temp,double *A,double *C,double *B,double *F)
{
    double alfa[N+1];
    double beta[N+1];
    alfa[1]=-B[0]/C[0];
    beta[1]=F[0]/C[0];
    for(int i = 1; i<=N; i++)
    {
        alfa[i+1]= -B[i]/(A[i]*alfa[i]+C[i]);
        beta[i+1]=(F[i]-A[i]*beta[i])/(A[i]*alfa[i]+C[i]);
    }
    u_temp[N]=(F[N]-A[N]*beta[N])/(C[N]+A[N]*alfa[N]);
    for(int i = N-1;i>-1; i--)
    {
        u_temp[i] = alfa[i+1]*u_temp[i+1]+beta[i+1];
    }
}
