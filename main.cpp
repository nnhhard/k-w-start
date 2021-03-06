﻿#include <iostream>
#include <cmath>
#include "Solvers.h"

#include <map>
#include <string>
#include <fstream>


using namespace std;

char name_out[100];


double f_t(double p1, double p2, double a_x, double b_x) {
    return (p2 - p1) / (b_x - a_x);
}
double d2dx2(double **vector, int i, int j, double *hx) {
    return ( ( ( ( vector[i+1][j] - vector[i][j] ) / (hx[i]) ) - ( ( vector[i][j] - vector[i-1][j] ) / (hx[i-1] ) ) ) / ( ( hx[i] + hx[i-1] ) / 2.0 ) );
}
double d2dy2(double **vector, int i, int j, double *hy) {
    return ( ( ( ( vector[i][j+1] - vector[i][j] ) / (hy[j]) ) - ( ( vector[i][j] - vector[i][j-1] ) / (hy[j-1] ) ) ) / ( ( hy[j] + hy[j-1] ) / 2.0 ) );
}
double ddxParabola(double **vector,int i,int j,double *hx) {
    double lamda = hx[i] / hx[i-1];
    return (vector[i+1][j] - vector[i][j] * (1 - lamda * lamda) - vector[i-1][j] * lamda * lamda) / (lamda * (hx[i] + hx[i-1]));
}
double ddyParabola(double **vector,int i,int j,double *hy) {
    double lamda = hy[j] / hy[j-1];
    return (vector[i][j+1] - vector[i][j] * (1 - lamda * lamda) - vector[i][j-1] * lamda * lamda) / (lamda * (hy[j] + hy[j-1]));
}
void null(double **vector, int lN, int lM) {
    for(int i = 0; i <= lN; i++)
        for(int j = 0; j <= lM; j++)
            vector[i][j]=0;
}
double norm(double **vector,int lN,int lM) {
    double max=fabs(vector[0][0]);

    for(int i=0;i<=lN;i++)
    {
        for(int j=0;j<=lM;j++)
        {
            if(fabs(vector[i][j])>max)
            {
                max=fabs(vector[i][j]);
            }
        }
    }
    return max;
}
double equationLine(double x1, double x2, double y1, double y2, double x) {
    double y = ( (x - x1) / (x2 - x1) * (y2 - y1) ) + y1;
    return y;
}
double operatorPressure(int lN, int lM, double **vector, int i, int j, double *hx, double *hy) {
    if(i == 0) {
        return (vector[i+1][j]+vector[i][j])/2.0;
    } else if(i == lN) {
        return (vector[i][j]+vector[i-1][j])/2.0;
    } else if(j == 0) {
        return -1.0 * (vector[i][j+1]-vector[i][j])/(hy[j]);
    } else if(j == lM) {
        return (vector[i][j]-vector[i][j-1])/(hy[j-1]);
    } else {
        return -1.0 * ( ( ( ( ( vector[i+1][j] - vector[i][j] ) / (hx[i]) ) - ( ( vector[i][j] - vector[i-1][j] ) / (hx[i-1] ) ) ) / ( ( hx[i] + hx[i-1] ) / 2.0 ) )
         + ( ( ( ( vector[i][j+1] - vector[i][j] ) / (hy[j]) ) - ( ( vector[i][j] - vector[i][j-1] ) / (hy[j-1] ) ) ) / ( ( hy[j] + hy[j-1] ) / 2.0 ) ) );
    }
}
double rightPartPressure(double p1, double p2, double TAU, int lN, int lM, double **u, double **v, int i, int j, double *hx, double *hy) {
    if(i == 0) {
        return p1;
    } else if(i == lN) {
        return p2;
    } else if(j == 0) {
        return 0;
    } else if(j == lM) {
        return 0;
    } else {
        return -1 * ( ( (u[i][j]-u[i-1][j])/(hx[i-1]) + (v[i][j]-v[i][j-1])/(hy[j-1]) ) / TAU );
    }
}
double sc(double **x, double **y,int lN,int lM)
{
    double s=0;
    for(int i=0;i<=lN;i++)
    {
        for(int j=0;j<=lM;j++)
        {
            s+=x[i][j]*y[i][j];
        }
    }
    return s;
}
void solveByBCGstab(double p1, double p2, double TAU, int lN, int lM, double **p,double **u,double **v,double *hx,double *hy,double *u_hx,double *v_hy)
{
    double **Rn = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        Rn[i] = new double [lM + 1];
    double **_Rn = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        _Rn[i] = new double [lM + 1];
    double **Pn = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        Pn[i] = new double [lM + 1];
    double **Pn1 = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        Pn1[i] = new double [lM + 1];
    double **Vn = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        Vn[i] = new double [lM + 1];
    double **Sn = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        Sn[i] = new double [lM + 1];
    double **Tn = new double *[lN+1];
    for(int i = 0; i <= lN; i++)
        Tn[i] = new double [lM + 1];
    double pn,an,wn,pn1,bn;

    null(Rn, lN, lM);
    null(_Rn, lN, lM);
    null(Pn, lN, lM);
    null(Pn1, lN, lM);
    null(Vn, lN, lM);
    null(Sn, lN, lM);
    null(Tn, lN, lM);


    for(int i = 0; i <= lN; i++)
        for(int j=0; j <= lM;j++)
            if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                Rn[i][j] = rightPartPressure(p1, p2, TAU, lN, lM, u, v, i, j, u_hx, v_hy) - operatorPressure(lN, lM, p, i, j, hx, hy);
    for(int i = 0; i <= lN; i++)
        for(int j=0; j <= lM;j++)
            if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                _Rn[i][j] = Rn[i][j];

    pn = 1.;
    an = 1.;
    bn = 1.;
    wn = 1.;
    int n1 = 0;
    double eps = 1e-12;

    while(norm(Rn, lN, lM) > eps)
    {
        n1++;
        pn1 = sc(_Rn, Rn, lN, lM);
        bn = (pn1/pn)*(an/wn);

        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    Pn1[i][j] = Rn[i][j] + bn*(Pn[i][j]-wn*Vn[i][j]);
        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    Vn[i][j] = operatorPressure(lN, lM, Pn1, i, j, hx, hy);
        an = pn1/sc(_Rn, Vn, lN, lM);
        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    Sn[i][j] = Rn[i][j] - an*Vn[i][j];
        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    Tn[i][j] = operatorPressure(lN, lM, Sn, i, j, hx, hy);
        wn = sc(Tn, Sn, lN, lM)/sc(Tn, Tn, lN, lM);
        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    p[i][j] += an*Pn1[i][j] + wn*Sn[i][j];
        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    Rn[i][j] = Sn[i][j] - wn*Tn[i][j];
        pn = pn1;
        for(int i = 0; i <= lN; i++)
            for(int j=0; j <= lM;j++)
                if(!(   (i==0 && j==0)   ||   (i == 0 && j == lM)   ||    (i == lN && j == 0)   ||   (i == lN && j == lM)     ))
                    Pn[i][j] = Pn1[i][j];
        //if(n1 % 10000 == 0)
            //printf("%d %.20lf\n", n1, norm(Rn, lN, lM));
    }
    for(int i = 0; i <= lN; i++)
    {
        delete Rn[i];
        delete _Rn[i];
        delete Pn[i];
        delete Pn1[i];
        delete Vn[i];
        delete Sn[i];
        delete Tn[i];

    }
    delete []Rn;
    delete []_Rn;
    delete []Pn;
    delete []Pn1;
    delete []Vn;
    delete []Sn;
    delete []Tn;
    Rn=NULL;
    _Rn=NULL;
    Pn=NULL;
    Pn1=NULL;
    Vn=NULL;
    Sn=NULL;
    Tn=NULL;
}
void solveTransportU(double p1, double p2, double TAU, double a_x, double b_x, double **u,double **u_n,double **v_n, double **nuTurbulent, int lN,int lM,double *hx,double *hy)
{
    double *u_temp_x=new double[lN+1];
    double *u_temp_y=new double[lM+1];

    double **u_1_2=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        u_1_2[i]=new double[lM+1];


    for(int i=0;i<=lN;i++)
        for(int j=0;j<=lM;j++)
        {
            u_1_2[i][j] = 0;
        }

    double *A=new double[lM+1];
    double *B=new double[lM+1];
    double *C=new double[lM+1];
    double *F=new double[lM+1];

    for(int j=0;j<=lM;j++)
    {
        A[j]=0;
        B[j]=0;
        C[j]=0;
        F[j]=0;
    }

    for(int i=1;i<lN;i++)
    {
        for(int j = 0; j<=lM; j++)
        {
            if(j == 0)
            {
                A[j]=0;
                B[j]=1.0/2.0;
                C[j]=1.0/2.0;
                F[j]=0 + TAU * f_t(p1, p2, a_x, b_x);
            }
            else if(j == lM)
            {
                A[j]=1.0/2.0;
                B[j]=1.0/2.0;
                C[j]=0;
                F[j]=0 + TAU * f_t(p1, p2, a_x, b_x);
            }
            else
            {
                double lamda = hy[j] / hy[j-1];
                A[j] =  - v_n[i][j] * lamda / (hy[j] + hy[j-1]) - 2.0 / ( (hy[j] + hy[j-1]) * hy[j-1] ) * nuTurbulent[i][j] + ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] )  / ( (hy[j] + hy[j-1]) * hy[j-1] );
                B[j] =  2.0/TAU - nuTurbulent[i][j] * (- 1.0 / hy[j] - 1.0 / hy[j-1] ) / ( (hy[j] + hy[j-1])/2.0 ) + v_n[i][j] * ( - (1 - lamda * lamda) / (lamda * (hy[j] + hy[j-1])) ) + /*ddyParabola(nuTurbulent, i, j, hy)*/ ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] ) * ( - (1 - lamda * lamda) / (lamda * (hy[j] + hy[j-1])) );
                C[j] =  v_n[i][j] / (lamda * (hy[j] + hy[j-1])) - 2.0 / ( (hy[j] + hy[j-1]) * hy[j] ) * nuTurbulent[i][j] + /*ddyParabola(nuTurbulent, i, j, hy)*/ ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] ) / ( (hy[j] + hy[j-1]) * hy[j] );
                F[j] =  - u_n[i][j] * ddxParabola(u_n,i,j,hx) + nuTurbulent[i][j] * d2dx2(u_n,i,j,hx) + 2.0*u_n[i][j]/TAU + /*ddxParabola(nuTurbulent, i, j, hx)*/( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) * ddxParabola(u_n, i, j, hx);
            }
        }
        SolveByScalarRun(lM,u_temp_y,A,B,C,F);

        for(int j=0;j<=lM;j++)
            u_1_2[i][j]=u_temp_y[j];

        }




        for(int i=0;i<=lN;i++)
        {
            for(int j=0;j<=lM;j++)
            {
                u[i][j]=u_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;


        A=new double[lN+1];
        B=new double[lN+1];
        C=new double[lN+1];
        F=new double[lN+1];


        for(int j=1;j<lM;j++)
        {
            for(int i = 0; i<=lN; i++)
            {
                if(i == 0)
                {
                    A[i] = 0;
                    B[i] = -1.0/hx[i];
                    C[i] = 1.0/hx[i];
                    F[i] = 0;
                }
                else if(i == lN)
                {
                    A[i] = -1.0/hx[i-1];
                    B[i] = 1.0/hx[i-1];
                    C[i] = 0;
                    F[i] = 0;
                }
                else
                {
                    double lamda = hx[i] / hx[i-1];
                    A[i] =  - u_n[i][j] * lamda / (hx[i] + hx[i-1]) - 2.0 / ( (hx[i] + hx[i-1]) * hx[i-1] ) * nuTurbulent[i][j] + /*ddxParabola(nuTurbulent, i, j, hx)*/( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) / ( (hx[i] + hx[i-1]) * hx[i-1] );
                    B[i] =  2.0/TAU - nuTurbulent[i][j] * (- 1.0 / hx[i] - 1.0 / hx[i-1] ) / ( (hx[i] + hx[i-1])/2.0 ) + u_n[i][j] * ( - (1 - lamda * lamda) / (lamda * (hx[i] + hx[i-1])) ) + /*ddxParabola(nuTurbulent, i, j, hx)*/ ( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) * ( - (1 - lamda * lamda) / (lamda * (hx[i] + hx[i-1])) );
                    C[i] =  u_n[i][j] / (lamda * (hx[i] + hx[i-1])) - 2.0 / ( (hx[i] + hx[i-1]) * hx[i] ) * nuTurbulent[i][j] + /*ddxParabola(nuTurbulent, i, j, hx)*/ ( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) / ( (hx[i] + hx[i-1]) * hx[i] );
                    F[i] =  - v_n[i][j] * ddyParabola(u_1_2,i,j,hy) + nuTurbulent[i][j] * d2dy2(u_1_2,i,j,hy) +2.0*u_1_2[i][j]/TAU + ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] ) * ddyParabola(u_1_2, i, j, hy);

                }
            }
        SolveByScalarRun(lN,u_temp_x,A,B,C,F);

        for(int i=0;i<=lN;i++)
            u[i][j]=u_temp_x[i];
        }

        for(int i=1;i<lN;i++)
        {
            u[i][0]= ( 0 + TAU * f_t(p1, p2, a_x, b_x) ) * 2  - u[i][1];
            u[i][lM]= ( 0 + TAU * f_t(p1, p2, a_x, b_x) ) * 2  - u[i][lM-1];
        }


    for(int i=0;i<=lN;i++)
    {
        delete []u_1_2[i];
    }
    delete []u_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp_x;
    delete []u_temp_y;
    u_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    u_temp_x=NULL;
    u_temp_y=NULL;
}
void solveTransportV(double TAU, double **u,double **u_n,double **v_n, double **nuTurbulent, int lN,int lM,double *hx,double *hy) {
    double *u_temp_x=new double[lN+1];
    double *u_temp_y=new double[lM+1];

    double **u_1_2=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        u_1_2[i]=new double[lM+1];


    for(int i=0;i<=lN;i++)
        for(int j=0;j<=lM;j++)
        {
            u_1_2[i][j] = 0;
        }

    double *A=new double[lM+1];
    double *B=new double[lM+1];
    double *C=new double[lM+1];
    double *F=new double[lM+1];

    for(int j=0;j<=lM;j++)
    {
        A[j]=0;
        B[j]=0;
        C[j]=0;
        F[j]=0;
    }

    for(int i=1;i<lN;i++)
    {
        for(int j = 0; j<=lM; j++)
        {
            if(j == 0) {
                A[j]=0;
                B[j]=1;
                C[j]=0;
                F[j]=0;
             } else if(j == lM) {
                A[j]=0;
                B[j]=1;
                C[j]=0;
                F[j]=0;
             } else {
                double lamda = hy[j] / hy[j-1];
                A[j] =  - v_n[i][j] * lamda / (hy[j] + hy[j-1]) - 2.0 / ( (hy[j] + hy[j-1]) * hy[j-1] ) * nuTurbulent[i][j] + /*ddyParabola(nuTurbulent, i, j, hy)*/ ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] )  / ( (hy[j] + hy[j-1]) * hy[j-1] );
                B[j] =  2.0/TAU - nuTurbulent[i][j] * (- 1.0 / hy[j] - 1.0 / hy[j-1] ) / ( (hy[j] + hy[j-1])/2.0 ) + v_n[i][j] * ( - (1 - lamda * lamda) / (lamda * (hy[j] + hy[j-1])) ) + /*ddyParabola(nuTurbulent, i, j, hy)*/ ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] ) * ( - (1 - lamda * lamda) / (lamda * (hy[j] + hy[j-1])) );
                C[j] =  v_n[i][j] / (lamda * (hy[j] + hy[j-1])) - 2.0 / ( (hy[j] + hy[j-1]) * hy[j] ) * nuTurbulent[i][j] + /*ddyParabola(nuTurbulent, i, j, hy)*/ ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] ) / ( (hy[j] + hy[j-1]) * hy[j] );
                F[j] =  - u_n[i][j] * ddxParabola(u_n,i,j,hx) + nuTurbulent[i][j] * d2dx2(u_n,i,j,hx) + 2.0*u_n[i][j]/TAU + /*ddxParabola(nuTurbulent, i, j, hx)*/( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) * ddxParabola(u_n, i, j, hx);
            }
        }
        SolveByScalarRun(lM,u_temp_y,A,B,C,F);

        for(int j=0;j<=lM;j++)
            u_1_2[i][j]=u_temp_y[j];

        }

        for(int i=0;i<=lN;i++)
        {
            for(int j=0;j<=lM;j++)
            {
                u[i][j]=u_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;


        A=new double[lN+1];
        B=new double[lN+1];
        C=new double[lN+1];
        F=new double[lN+1];


        for(int j=1;j<lM;j++)
        {
            for(int i = 0; i<=lN; i++)
            {
                if(i == 0) {
                    A[i]=0;
                    B[i]=1.0/2.0;
                    C[i]=1.0/2.0;
                    F[i]=0;
            } else if(i == lN) {
                A[i]=1.0/2.0;
                B[i]=1.0/2.0;
                C[i]=0;
                F[i]=0;
            } else {
                double lamda = hx[i] / hx[i-1];
                A[i] =  - u_n[i][j] * lamda / (hx[i] + hx[i-1]) - 2.0 / ( (hx[i] + hx[i-1]) * hx[i-1] ) * nuTurbulent[i][j] + /*ddxParabola(nuTurbulent, i, j, hx)*/( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) / ( (hx[i] + hx[i-1]) * hx[i-1] );
                B[i] =  2.0/TAU - nuTurbulent[i][j] * (- 1.0 / hx[i] - 1.0 / hx[i-1] ) / ( (hx[i] + hx[i-1])/2.0 ) + u_n[i][j] * ( - (1 - lamda * lamda) / (lamda * (hx[i] + hx[i-1])) ) + /*ddxParabola(nuTurbulent, i, j, hx)*/ ( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) * ( - (1 - lamda * lamda) / (lamda * (hx[i] + hx[i-1])) );
                C[i] =  u_n[i][j] / (lamda * (hx[i] + hx[i-1])) - 2.0 / ( (hx[i] + hx[i-1]) * hx[i] ) * nuTurbulent[i][j] + /*ddxParabola(nuTurbulent, i, j, hx)*/ ( (nuTurbulent[i][j] - nuTurbulent[i-1][j]) / hx[i-1] ) / ( (hx[i] + hx[i-1]) * hx[i] );
                F[i] =  - v_n[i][j] * ddyParabola(u_1_2,i,j,hy) + nuTurbulent[i][j] * d2dy2(u_1_2,i,j,hy) +2.0*u_1_2[i][j]/TAU + ( (nuTurbulent[i][j] - nuTurbulent[i][j-1]) / hy[j-1] ) * ddyParabola(u_1_2, i, j, hy);
            }
        }
        SolveByScalarRun(lN,u_temp_x,A,B,C,F);

        for(int i=0;i<=lN;i++)
            u[i][j]=u_temp_x[i];
        }

    for(int i=0;i<=lN;i++)
    {
        delete []u_1_2[i];
    }
    delete []u_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp_x;
    delete []u_temp_y;
    u_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    u_temp_x=NULL;
    u_temp_y=NULL;
}
void solveTransportK(double BETTA, double TAU, double NU, double SIGMA, double **k, double **k_n, double **omega, double **u, double **uClassic, double **v, double **vClassic, double **tau11, double **tau12, double **tau21, double **tau22, int lN,int lM, double *hx, double *hy, double *u_hx, double *v_hy)
{
    double *k_temp_x=new double[lN+1];
    double *k_temp_y=new double[lM+1];

    double **k_1_2=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        k_1_2[i]=new double[lM+1];
    double **func=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        func[i]=new double[lM+1];


    for(int i=0;i<=lN;i++)
        for(int j=0;j<=lM;j++)
        {
            k_1_2[i][j] = 0;
            func[i][j] = NU + SIGMA * k_n[i][j] / omega[i][j];
        }


    double *A=new double[lM+1];
    double *B=new double[lM+1];
    double *C=new double[lM+1];
    double *F=new double[lM+1];

    for(int j=0;j<=lM;j++)
    {
        A[j]=0;
        B[j]=0;
        C[j]=0;
        F[j]=0;
    }

    for(int i=1;i<lN;i++)
    {
        for(int j = 0; j<=lM; j++)
        {
            if(j == 0)
            {
                A[j]=0;
                B[j]=1.0/2.0;
                C[j]=1.0/2.0;
                F[j]=0;
            }
            else if(j == lM)
            {
                A[j]=1.0/2.0;
                B[j]=1.0/2.0;
                C[j]=0;
                F[j]=0;
            }
            else
            {
                A[j] = -v[i][j] / (hy[j-1]) - ( (func[i][j] - func[i][j-1]) / hy[j-1] * ( - 1.0 / hy[j-1] ) - func[i][j] * 2.0 / ( (hy[j] + hy[j-1]) * hy[j-1] ) );
                B[j] = 2.0/TAU + v[i][j] / (hy[j-1]) - (func[i][j] - func[i][j-1]) / hy[j-1] *     1.0 / hy[j-1] - func[i][j] * (- 1.0 / hy[j] - 1.0 / hy[j-1] ) / ( (hy[j] + hy[j-1])/2.0 );
                C[j] = - func[i][j] * 2.0 / ( (hy[j] + hy[j-1]) * hy[j] ) ;
                F[j] =  - u[i][j] * ( k_n[i][j] - k_n[i-1][j] ) / hx[i-1] + 2.0*k_n[i][j]/TAU +
                        (
                            tau11[i][j] * ( uClassic[i][j] - uClassic[i-1][j] ) / u_hx[i-1] +
                            tau12[i][j] * ddyParabola(u, i, j, hy)                        +
                            tau21[i][j] * ddxParabola(v, i, j, hx)                        +
                            tau22[i][j] * ( vClassic[i][j] - vClassic[i][j-1] ) / v_hy[j-1]
                        )
                - BETTA * k_n[i][j] * omega[i][j] + ( (func[i][j] - func[i-1][j] ) / hx[i-1] * (k_n[i][j] - k_n[i-1][j] ) / hx[i-1] + func[i][j] * d2dx2(k, i, j, hx) );
            }
        }
        SolveByScalarRun(lM,k_temp_y,A,B,C,F);

        for(int j=0;j<=lM;j++)
            k_1_2[i][j] = k_temp_y[j];

        }

        for(int i=0;i<=lN;i++)
        {
            for(int j=0;j<=lM;j++)
            {
                k[i][j]=k_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;


        A=new double[lN+1];
        B=new double[lN+1];
        C=new double[lN+1];
        F=new double[lN+1];


        for(int j=1;j<lM;j++)
        {
            for(int i = 0; i<=lN; i++)
            {
                if(i == 0)
                {
                    A[i] = 0;
                    B[i] = 1.0 / 2.0;
                    C[i] = 1.0 / 2.0;
                    F[i] = 0;
                }
                else if(i == lN)
                {
                    A[i] = -1.0 / hx[i-1];
                    B[i] = 1.0 / hx[i-1];
                    C[i] = 0;
                    F[i] = 0;
                }
                else
                {
                    A[i] = -u[i][j] / (hx[i-1]) - ( (func[i][j] - func[i-1][j]) / hx[i-1] * ( - 1.0 / hx[i-1] ) - func[i][j] * 2.0 / ( (hx[i] + hx[i-1]) * hx[i-1] ) );
                    B[i] = 2.0/TAU - func[i][j] * (- 1.0 / hx[i] - 1.0 / hx[i-1] ) / ( (hx[i] + hx[i-1])/2.0 ) + u[i][j] / (hx[i-1]) -  (func[i][j] - func[i-1][j]) / hx[i-1] * 1.0 / hx[i-1] ;
                    C[i] =  - func[i][j] * 2.0 / ( (hx[i] + hx[i-1]) * hx[i]   ) ;
                    F[i] =  - v[i][j] * ( k_1_2[i][j] - k_1_2[i][j-1] ) / hy[j-1] + 2.0*k_1_2[i][j]/TAU +
                            (
                                tau11[i][j] * ( uClassic[i][j] - uClassic[i-1][j] ) / u_hx[i-1] +
                                tau12[i][j] * ddyParabola(u, i, j, hy)                        +
                                tau21[i][j] * ddxParabola(v, i, j, hx)                        +
                                tau22[i][j] * ( vClassic[i][j] - vClassic[i][j-1] ) / v_hy[j-1]
                            )
                    - BETTA * k_1_2[i][j] * omega[i][j] + ( (func[i][j] - func[i][j-1] ) / hy[j-1] * (k_n[i][j] - k_n[i][j-1] ) / hy[j-1] + func[i][j] * d2dy2(k, i, j, hy) );
                }
            }
        SolveByScalarRun(lN,k_temp_x,A,B,C,F);

        for(int i=0;i<=lN;i++)
            k[i][j] = k_temp_x[i];
        }

        for(int i=1;i<lN;i++)
        {
            k[i][0]= - k[i][1];
            k[i][lM]= - k[i][lM-1];
        }


    for(int i=0;i<=lN;i++)
    {
        delete []k_1_2[i];
    }
    delete []k_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []k_temp_x;
    delete []k_temp_y;
    k_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    k_temp_x=NULL;
    k_temp_y=NULL;
}
void solveTransportW(double ALPHA, double BETTA, double TAU, double NU, double SIGMA, double **w, double **w_n, double **k, double **u, double **uClassic, double **v, double **vClassic, double **tau11, double **tau12, double **tau21, double **tau22, int lN,int lM, double *hx, double *hy, double *u_hx, double *v_hy)
{
    double *w_temp_x=new double[lN+1];
    double *w_temp_y=new double[lM+1];

    double **SIGMAd=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        SIGMAd[i]=new double[lM+1];
    double **w_1_2=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        w_1_2[i]=new double[lM+1];
    double **func=new double*[lN+1];
    for(int i=0;i<=lN;i++)
        func[i]=new double[lM+1];


    for(int i=0;i<=lN;i++)
        for(int j=0;j<=lM;j++)
        {
            w_1_2[i][j] = 0;
            func[i][j] = NU + SIGMA * k[i][j] / w_n[i][j];
        }
    null(SIGMAd, lN, lM);
    double temp = 0;
    for(int i=1;i<lN;i++)
        for(int j=1;j<lM;j++)
        {
            temp = ( ( k[i][j] - k[i-1][j] ) / hx[i-1] * ( w_n[i][j] - w_n[i-1][j] ) / hx[i-1] )
                 + ( ( k[i][j] - k[i][j-1] ) / hy[j-1] * ( w_n[i][j] - w_n[i][j-1] ) / hy[j-1] );
            if(temp < 0) {
               SIGMAd[i][j] = 0;
            } else {
               SIGMAd[i][j] = 0.3;
            }
        }

    double *A=new double[lM+1];
    double *B=new double[lM+1];
    double *C=new double[lM+1];
    double *F=new double[lM+1];

    for(int j=0;j<=lM;j++)
    {
        A[j]=0;
        B[j]=0;
        C[j]=0;
        F[j]=0;
    }

    for(int i=1;i<lN;i++)
    {
        for(int j = 0; j<=lM; j++)
        {
            if(j == 0)
            {
                A[j]=0;
                B[j]=1.0/2.0;
                C[j]=1.0/2.0;
                F[j]= 10 * 6 * NU / (BETTA * hy[j] * hy[j]);
            }
            else if(j == lM)
            {
                A[j]=1.0/2.0;
                B[j]=1.0/2.0;
                C[j]=0;
                F[j]= 10 * 6 * NU / (BETTA * hy[j-1] * hy[j-1]);
            }
            else
            {
                A[j] = -v[i][j] / (hy[j-1]) + SIGMAd[i][j] / w_n[i][j] * ( k[i][j] - k[i][j-1] ) / hy[j-1] * 1.0 / hy[j-1] + ( func[i][j] - func[i][j-1] ) / hy[j-1] * 1.0 / hy[j-1] - func[i][j] * 2.0 / (hy[j-1] * (hy[j] + hy[j-1]));
                B[j] = 2.0/TAU + v[i][j] / (hy[j-1]) - SIGMAd[i][j] / w_n[i][j] * ( k[i][j] - k[i][j-1] ) / hy[j-1] * 1.0 / hy[j-1] - ( func[i][j] - func[i][j-1] ) / hy[j-1] * 1.0 / hy[j-1] - func[i][j] * (- 1.0 / hy[j] - 1.0 / hy[j-1] ) / ( (hy[j] + hy[j-1])/2.0 );
                C[j] = - func[i][j] * 2.0 / ( (hy[j] + hy[j-1]) * hy[j]  );
                F[j] =  - u[i][j] * ( w_n[i][j] - w_n[i-1][j] ) / hx[i-1] + 2.0*w_n[i][j]/TAU + ALPHA * w_n[i][j] / k[i][j] *
                        (
                            tau11[i][j] * ( uClassic[i][j] - uClassic[i-1][j] ) / u_hx[i-1] +
                            tau12[i][j] * ddyParabola(u, i, j, hy)                        +
                            tau21[i][j] * ddxParabola(v, i, j, hx)                        +
                            tau22[i][j] * ( vClassic[i][j] - vClassic[i][j-1] ) / v_hy[j-1]
                        )
                + ( (func[i][j] - func[i-1][j] ) / hx[i-1] * (w_n[i][j] - w_n[i-1][j] ) / hx[i-1] + func[i][j] * d2dx2(w_n, i, j, hx) )
                + (  SIGMAd[i][j] / w_n[i][j] * (k[i][j] - k[i-1][j] ) / hx[i-1] * (w_n[i][j] - w_n[i-1][j] ) / hx[i-1] );
            }
        }
        SolveByScalarRun(lM,w_temp_y,A,B,C,F);

        for(int j=0;j<=lM;j++)
            w_1_2[i][j] = w_temp_y[j];

        }

        for(int i=0;i<=lN;i++)
        {
            for(int j=0;j<=lM;j++)
            {
                w[i][j]=w_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;


        A=new double[lN+1];
        B=new double[lN+1];
        C=new double[lN+1];
        F=new double[lN+1];


        for(int j=1;j<lM;j++)
        {
            for(int i = 0; i<=lN; i++)
            {
                if(i == 0)
                {
                    A[i] = 0;
                    B[i] = 1.0 / 2.0;
                    C[i] = 1.0 / 2.0;
                    F[i] = 1.0 / (10 * 5.0 / 1.0);
                }
                else if(i == lN)
                {
                    A[i] = -1.0 / hx[i-1];
                    B[i] = 1.0 / hx[i-1];
                    C[i] = 0;
                    F[i] = 0;
                }
                else
                {
                    A[i] = -u[i][j] / (hx[i-1]) + SIGMAd[i][j] / w_n[i][j] * ( k[i][j] - k[i-1][j] ) / hx[i-1] * 1.0 / hx[i-1] + ( func[i][j] - func[i-1][j] ) / hx[i-1] * 1.0 / hx[i-1] - func[i][j] * 2.0 / (hx[i-1] * (hx[i] + hx[i-1]));
                    B[i] = 2.0/TAU + u[i][j] / (hx[i-1]) - SIGMAd[i][j] / w_n[i][j] * ( k[i][j] - k[i-1][j] ) / hx[i-1] * 1.0 / hx[i-1] - ( func[i][j] - func[i-1][j] ) / hx[i-1] * 1.0 / hx[i-1] - func[i][j] * (- 1.0 / hx[i] - 1.0 / hx[i-1] ) / ( (hx[i] + hx[i-1])/2.0 );
                    C[i] = - func[i][j] * 2.0 / ( (hx[i] + hx[i-1]) * hx[i]  );
                    F[i] =  - v[i][j] * ( w_1_2[i][j] - w_1_2[i][j-1] ) / hy[j-1] + 2.0*w_1_2[i][j]/TAU + ALPHA * w_1_2[i][j] / k[i][j] *
                            (
                                tau11[i][j] * ( uClassic[i][j] - uClassic[i-1][j] ) / u_hx[i-1] +
                                tau12[i][j] * ddyParabola(u, i, j, hy)                        +
                                tau21[i][j] * ddxParabola(v, i, j, hx)                        +
                                tau22[i][j] * ( vClassic[i][j] - vClassic[i][j-1] ) / v_hy[j-1]
                            )
                    + ( (func[i][j] - func[i][j-1] ) / hy[j-1] * (w_1_2[i][j] - w_n[i][j-1] ) / hy[j-1] + func[i][j] * d2dy2(w_1_2, i, j, hy) )
                    + (  SIGMAd[i][j] / w_n[i][j] * (k[i][j] - k[i][j-1] ) / hy[j-1] * (w_1_2[i][j] - w_1_2[i][j-1] ) / hy[j-1] );
                }
            }
        SolveByScalarRun(lN,w_temp_x,A,B,C,F);

        for(int i=0;i<=lN;i++)
            w[i][j] = w_temp_x[i];
        }

        for(int i=1;i<lN;i++)
        {
            w[i][0]  = (10 * 6 * NU / (BETTA * hy[0] * hy[0]) ) * 2 - w[i][1];
            w[i][lM] = (10 * 6 * NU / (BETTA * hy[lM-1] * hy[lM-1]) ) * 2 - w[i][lM-1];
        }


    for(int i=0;i<=lN;i++)
    {
        delete []w_1_2[i];
    }
    delete []w_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []w_temp_x;
    delete []w_temp_y;
    w_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    w_temp_x=NULL;
    w_temp_y=NULL;
}


void Print(FILE *f,int t, int lN, int lM, double *x, double *y, double **u,double **v, double **p) {
    fprintf(f, "TITLE = \"PROTEKANIE\"\n");
    fprintf(f, "VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"P\"\n");
    fprintf(f, "ZONE T=\"Zone%d\",I=%d,J=%d,F=POINT\n",10000 + t, lM + 1, lN + 1);

    for (int i = 0; i <= lN; i++)
        {
        for (int j = 0; j <= lM ; j++)
        {
            fprintf(f, "%lf,%lf,%lf,%lf,%lf\n", x[i], y[j], u[i][j], v[i][j], p[i][j]);
        }
    }
}
void out(int it, int N, int M, int input, double *x, double *y, double **p_out, double **u_out, double **v_out, double *hx, double *hy, double **p, double **u, double **v) {
    for(int i=0;i<=N;i++) {
        for(int j=0;j<=M;j++) {
            if(   (i==0 && j==0)   ||   (i==0 && j==M)   ||    (i == N && j == 0)   ||   (i == N && j == M)     ) {
                p_out[i][j] = 0;
                v_out[i][j]=0;
                u_out[i][j]=0;
            } else if(i == 0 || i == N) {
                u_out[i][j] = equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, u[i][j], u[i][j+1], 0);
                v_out[i][j] = (v[i][j] + v[i+1][j]) / 2.0;
                p_out[i][j] =
                                (
                                    equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, p[i][j], p[i][j+1], 0)
                                  + equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, p[i+1][j], p[i+1][j+1], 0)
                                ) / 2.0;
            } else if(j == 0 || j == M) {
                u_out[i][j] = (u[i][j] + u[i][j+1]) / 2.0;
                v_out[i][j] = equationLine( -hx[i-1] / 2.0, hx[i] / 2.0, v[i][j], v[i+1][j], 0);
                p_out[i][j] =
                                (
                                    equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, p[i][j], p[i][j+1], 0)
                                  + equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, p[i+1][j], p[i+1][j+1], 0)
                                ) / 2.0;
            } else {
                u_out[i][j] = equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, u[i][j], u[i][j+1], 0);
                v_out[i][j] = equationLine( -hx[i-1] / 2.0, hx[i] / 2.0, v[i][j], v[i+1][j], 0);
                p_out[i][j] = equationLine
                        (
                            -hx[i-1] / 2.0,
                            hx[i] / 2.0,
                            equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, p[i][j], p[i][j+1], 0),
                            equationLine( -hy[j-1] / 2.0, hy[j] / 2.0, p[i+1][j], p[i+1][j+1], 0),
                            0
                        );
            }
        }
    }
    p_out[0][0] = p_out[0][1];
    p_out[0][M] = p_out[0][M-1];
    p_out[N][0] = p_out[N][1];
    p_out[N][M] = p_out[N][M-1];

    sprintf(name_out,"./Solve/Zone_%d.dat",it);
    if( it % input == 0 ) {
        FILE *f = fopen(name_out,"w");
        Print(f, it, N, M, x, y, u_out, v_out, p_out);
        fclose(f);
    }
}

string extractLeft(const string& buf, size_t i) {
    size_t left = 0;
    while(buf[left] == ' ')
        ++left;
    size_t right = i-1;
    while(buf[right] == ' ')
        --right;
    return buf.substr(left, right-left+1);
}
string extractRight(const string& buf, size_t i) {
    size_t left = i+1;
    while(buf[left] == ' ')
        ++left;
    size_t right = buf.size()-1;
    while(buf[right] == ' ')
        --right;
    return buf.substr(left, right-left+1);
}
map<string, string> parseDict(const string& filename) {
    ifstream file(filename);
    if(!file.is_open())
        throw runtime_error("Cant open file: " + filename);
    string buffer;
    map<string, string> tmp;
    while(!file.eof()) {
        getline(file, buffer);
        //cout << buffer << endl;
        size_t i = buffer.find('=', 0);
        if(i != string::npos) {
            tmp.emplace(extractLeft(buffer, i), extractRight(buffer, i));
        }
    }
    return tmp;
}

const int input = 10;
const double T = 10;

const double SIGMA  = 0.5;
const double ALPHA  = 5.0 / 9.0;

const double BETTAS = 0.09;
const double BETTA  = 0.075;

int main() {
    map<string, string> dict = parseDict("input.txt");
    const int N = stoi(dict["Number X"]);
    cout << "N: " << N << endl;
    const int M = stoi(dict["Number Y"]);
    cout << "M: " << M << endl;
    const double p1 = stod(dict["Pressure inlet"]);
    cout << "p1: " << p1 << endl;
    const double p2 = stod(dict["Pressure outlet"]);
    cout << "p2: " << p2 << endl;
    const double NU = stod(dict["NU"]);
    cout << "NU: " << NU << endl;
    const double TAU = stod(dict["TAU"]);
    cout << "TAU: " << TAU << endl;
    const double a_x = stod(dict["AX"]);
    cout << "a_x: " << a_x << endl;
    const double b_x = stod(dict["BX"]);
    cout << "b_x: " << b_x << endl;
    const double a_y = stod(dict["AY"]);
    cout << "a_y: " << a_y << endl;
    const double b_y = stod(dict["BY"]);
    cout << "b_y: " << b_y << endl;


    const int uN = N;
    const int uM = M + 1;

    const int vN = N + 1;
    const int vM = M;

    const int pN = N + 1;
    const int pM = M + 1;

    double *hx = new double[N];
    double *hy = new double[M];

    double *u_hx = new double[uN];
    double *u_hy = new double[uM];

    double *v_hx = new double[vN];
    double *v_hy = new double[vM];

    double *p_hx = new double[pN];
    double *p_hy = new double[pM];

    double h = (b_x - a_x) / (N);
    for(int i = 0; i <= N-1; i++) {
        hx[i] = h;
    }
    h = (0.1 - a_y) / (M / 3 - 1);
    for(int i = 0; i < M / 3 - 1; i++) {
        hy[i] = h;
    }
    h = (b_y - (b_y - 0.1)) / (M / 3 - 1);
    for(int i = (M - 1) - (M / 3 - 1); i <= (M - 1); i++) {
        hy[i] = h;
    }
    h = ( (b_y - 0.1) - (a_y + 0.1) ) / (M - (M / 3 - 1) * 2 );
    for(int i = M / 3 - 1; i <= (M - 1) - (M / 3 - 1); i++) {
        hy[i] = h;
    }

    /// Заполнение сеток hx и hy для вектора u (start)
    for(int i = 0; i <= uN-1; i++) {
        u_hx[i] = hx[i];
    }

    u_hy[0] = hy[0];
    u_hy[uM-1] = hy[M-1];

    for(int i = 1; i <= uM-2; i++) {
        u_hy[i] = (hy[i-1] + hy[i]) / 2.0;
    }
    /// Заполнение сеток hx и hy для вектора u (end)

    /// Заполнение шагов сетки hx и hy для вектора v (start)
    for(int i = 0; i <= vM-1; i++) {
        v_hy[i] = hy[i];
    }

    v_hx[0] = hx[0];
    v_hx[vN-1] = hx[N-1];
    for(int i = 1; i <= vN-2; i++) {
        v_hx[i] = (hx[i-1] + hx[i]) / 2.0;
    }
    /// Заполнение шагов сетки hx и hy для вектора v (end)

    /// Заполнение шагов сетки hx и hy для функции p (start)
    p_hy[0] = hy[0];
    p_hy[pM-1] = hy[M-1];

    for(int i = 1; i <= pM-2; i++) {
        p_hy[i] = (hy[i-1] + hy[i]) / 2.0;
    }

    p_hx[0] = hx[0];
    p_hx[pN-1] = hx[N-1];
    for(int i = 1; i <= pN-2; i++) {
        p_hx[i] = (hx[i-1] + hx[i]) / 2.0;
    }
    /// Заполнение шагов сетки hx и hy для функции p (end)

    double *x = new double[N+1];
    x[0]=a_x;
    for(int i = 1; i <= N; i++)
    {
        x[i] = x[i-1] + hx[i-1];
    }
    double *y = new double[M+1];
    y[0]=a_y;
    for(int i = 1; i <= M; i++)
    {
        y[i] = y[i-1] + hy[i-1];
    }

    double **u = new double*[uN + 1];
    for(int i = 0; i <= uN; i++)
        u[i] = new double[uM + 1];

    double **nutProecU = new double*[uN + 1];
    for(int i = 0; i <= uN; i++)
        nutProecU[i] = new double[uM + 1];

    double **u_n = new double*[uN + 1];
    for(int i = 0; i <= uN; i++)
        u_n[i] = new double[uM + 1];

    double **v = new double*[vN + 1];
    for(int i = 0; i <= vN; i++)
        v[i] = new double[vM + 1];

    double **nutProecV = new double*[vN + 1];
    for(int i = 0; i <= vN; i++)
        nutProecV[i] = new double[vM + 1];

    double **v_n = new double*[vN + 1];
    for(int i = 0; i <= vN; i++)
        v_n[i] = new double[vM + 1];

    double **p = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        p[i] = new double[pM + 1];

    double **k = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        k[i] = new double[pM + 1];

    double **k_n = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        k_n[i] = new double[pM + 1];

    double **S11 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        S11[i] = new double[pM + 1];
    double **S12 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        S12[i] = new double[pM + 1];
    double **S21 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        S21[i] = new double[pM + 1];
    double **S22 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        S22[i] = new double[pM + 1];

    double **tau11 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        tau11[i] = new double[pM + 1];
    double **tau12 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        tau12[i] = new double[pM + 1];
    double **tau21 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        tau21[i] = new double[pM + 1];
    double **tau22 = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        tau22[i] = new double[pM + 1];

    double **w = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        w[i] = new double[pM + 1];

    double **w_n = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        w_n[i] = new double[pM + 1];

    double **nuTurbulent = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        nuTurbulent[i] = new double[pM + 1];

    double **u_proecNut = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        u_proecNut[i] = new double[pM + 1];

    double **v_proecNut = new double*[pN + 1];
    for(int i = 0; i <= pN; i++)
        v_proecNut[i] = new double[pM + 1];

    double **v_proec = new double *[uN + 1];
    for(int i = 0; i <= uN; i++)
        v_proec[i] = new double[uM + 1];

    double **u_proec = new double *[vN + 1];
    for(int i = 0; i <= vN; i++)
        u_proec[i] = new double[vM + 1];

    double **p_out = new double *[N + 1];
    for(int i = 0; i <= N; i++)
        p_out[i] = new double[M+1];

    double **u_out = new double *[N + 1];
    for(int i = 0; i <= N; i++)
        u_out[i] = new double[M+1];
    double **v_out = new double *[N + 1];
    for(int i = 0; i <= N; i++)
        v_out[i] = new double[M+1];

    null(tau11, pN, pM);
    null(tau12, pN, pM);
    null(tau21, pN, pM);
    null(tau22, pN, pM);
    null(S11, pN, pM);
    null(S12, pN, pM);
    null(S21, pN, pM);
    null(S22, pN, pM);
    null(u_proecNut, pN, pM);
    null(v_proecNut, pN, pM);
    null(u_out, N, M);
    null(v_out, N, M);
    null(p_out, N, M);
    null(u_n, uN, uM);
    null(u, uN, uM);
    null(p, pN, pM);
    null(k, pN, pM);
    null(w, pN, pM);
    null(nuTurbulent, pN, pM);
    null(v_n, vN, vM);
    null(v, vN, vM);
    null(u_proec, vN, vM);
    null(v_proec, uN, uM);

    for(int i = 0; i <= pN; i++) {
        for(int j = 0; j <= pM; j++) {
            w[i][j] = 1;
        }
    }

    int it = 0;
    double t =  TAU;
    while(t <= T) {
        it++;


        /// Интерпалировать nut на сетки u и v
        for(int i = 0; i <= uN; i++) {
            for(int j = 0; j <= uM; j++) {
                if( (i == 0) || (i == uN) ) {
                    nutProecU[i][j] = (nuTurbulent[i+1][j] + nuTurbulent[i][j]) / 2.0 + NU;
                } else {
                    nutProecU[i][j] = equationLine( -u_hx[i-1] / 2.0, u_hx[i] / 2.0, nuTurbulent[i][j], nuTurbulent[i+1][j], 0) + NU;
                }
            }
        }

        for(int i = 0; i <= vN; i++) {
            for(int j = 0; j <= vM; j++) {
                if( (j == 0) || (j == vM) ) {
                    nutProecV[i][j] = (nuTurbulent[i][j+1] + nuTurbulent[i][j]) / 2.0 + NU;
                } else {
                    nutProecV[i][j] = equationLine( -v_hy[i-1] / 2.0, v_hy[i] / 2.0, nuTurbulent[i][j], nuTurbulent[i][j+1], 0) + NU;
                }
            }
        }

        for(int i=1;i<uN;i++) {
            for(int j=1;j<uM;j++) {
               v_proec[i][j] = (equationLine( -u_hx[i-1] / 2.0, u_hx[i] / 2.0, v_n[i][j-1], v_n[i+1][j-1], 0)
                              + equationLine( -u_hx[i-1] / 2.0, u_hx[i] / 2.0, v_n[i][j], v_n[i+1][j], 0)) / 2.0;
            }
        }

        for(int i=1;i<vN;i++) {
            for(int j=1;j<vM;j++) {
                u_proec[i][j] = (equationLine( -v_hy[j-1] / 2.0, v_hy[j] / 2.0, u_n[i-1][j],u_n[i-1][j+1], 0)
                              + equationLine( -v_hy[j-1] / 2.0, v_hy[j] / 2.0, u_n[i][j], u_n[i][j+1], 0)) / 2.0;
            }
        }

        solveTransportU(p1, p2, TAU, a_x, b_x, u, u_n, v_proec, nutProecU, uN, uM, u_hx, u_hy);
        solveTransportV(TAU, v, v_n, u_proec, nutProecV, vN, vM, v_hx, v_hy);
        solveByBCGstab(p1, p2, TAU, pN, pM, p, u, v, p_hx, p_hy, u_hx, v_hy);

        for(int i = 0; i <= uN; i++)
            for(int j = 0; j <= uM; j++)
                if(!(   (i==0 && j==0)   ||   (i==0 && j==uM)   ||    (i==uN && j==0)   ||   (i==uN && j==uM)     )) {
                        u[i][j] = u[i][j] - TAU * (p[i+1][j] - p[i][j]) / p_hx[i];
                }

        for(int i = 0; i <= vN; i++)
            for(int j = 0; j <= vM; j++)
                if(!(   (i==0 && j==0)   ||   (i==0 && j==vM)   ||    (i==vN && j==0)   ||   (i==vN && j==vM)     )) {
                    v[i][j] = v[i][j] - TAU * (p[i][j+1] - p[i][j]) / p_hy[j];
                }

        /// Интерполировать u и v в центры ячеек
        for(int i = 1; i < pN; i++) {
            for(int j = 0; j <= pM; j++) {
                u_proecNut[i][j] = (u[i][j] + u[i-1][j]) / 2.0;
            }
        }
        for(int i = 0; i <= pN; i++) {
            for(int j = 1; j < pM; j++) {
                v_proecNut[i][j] = (v[i][j] + v[i][j-1]) / 2.0;
            }
        }
        for(int i = 1; i < pN; i++) {
            for(int j = 1; j < pM; j++) {
                S11[i][j] = 1.0 / 2.0 * ( (u[i][j] - u[i-1][j]) / u_hx[i-1] + (u[i][j] - u[i-1][j]) / u_hx[i-1]);

                S12[i][j] = 1.0 / 2.0 * ( ddyParabola(u_proecNut, i, j, p_hy) + ddxParabola(v_proecNut, i, j, p_hx) );
                S21[i][j] = 1.0 / 2.0 * ( ddyParabola(u_proecNut, i, j, p_hy) + ddxParabola(v_proecNut, i, j, p_hx) );

                S22[i][j] = 1.0 / 2.0 * ( (v[i][j] - v[i][j-1]) / v_hy[j-1] + (v[i][j] - v[i][j-1]) / v_hy[j-1]);
            }
        }
        for(int i = 1; i < pN; i++) {
            for(int j = 1; j < pM; j++) {
                tau11[i][j] = 2 * nuTurbulent[i][j] * S11[i][j] - 2.0/3.0 * k[i][j] * 1;
                tau12[i][j] = 2 * nuTurbulent[i][j] * S12[i][j] - 2.0/3.0 * k[i][j] * 0;
                tau21[i][j] = 2 * nuTurbulent[i][j] * S21[i][j] - 2.0/3.0 * k[i][j] * 0;
                tau22[i][j] = 2 * nuTurbulent[i][j] * S22[i][j] - 2.0/3.0 * k[i][j] * 1;
            }
        }
        /// Найти k
        solveTransportK(BETTAS, TAU, NU, SIGMA, k, k_n, w, u_proecNut, u, v_proecNut, v, tau11, tau12, tau21, tau22, pN, pM, p_hx, p_hy, u_hx, v_hy);
        for(int i = 0; i <= pN; i++) {
            for(int j = 0; j <= pM; j++) {
                printf("%lf ", k[i][j]);
            }
            printf("\n");
        }

        printf("\n");
        printf("\n");
        printf("\n");
        /// Найти w
        solveTransportW(ALPHA, BETTA, TAU, NU, SIGMA, w, w_n, k, u_proecNut, u, v_proecNut, v, tau11, tau12, tau21, tau22, pN, pM, p_hx, p_hy, u_hx, v_hy);
        for(int i = 0; i <= pN; i++) {
            for(int j = 0; j <= pM; j++) {
                printf("%lf ", w[i][j]);
            }
            printf("\n");
        }
   return 0;

        out(it,N,M,input, x,y,p_out,u_out, v_out, hx, hy, p, u, v);
        for(int i = 0; i <= vN; i++) {
            for(int j = 0; j <= vM; j++) {
                v_n[i][j] = v[i][j];
            }
        }

        /// Найти турбулентную вязкость
        for(int i = 0; i <= pN; i++) {
            for(int j = 0; j <= pM; j++) {
                k_n[i][j] = k[i][j];
                w_n[i][j] = w[i][j];
                nuTurbulent[i][j] = k[i][j] / w[i][j];
            }
        }

        for(int i = 0; i <= uN; i++) {
            for(int j = 0; j <= uM; j++) {
                u_n[i][j] = u[i][j];
            }
        }

        double divV = 0;
        for(int i = 1; i < pN; i++) {
            for(int j = 1; j < pM; j++) {
                divV+=((u_n[i][j]-u_n[i-1][j])/(u_hx[i-1]) + (v_n[i][j]-v_n[i][j-1])/(v_hy[j-1]));
            }
        }

        printf("time = %2.5lf div(V*) = %2.10lf\n",t,fabs(divV));

        t = t + TAU;
    }


    for(int i = 0; i <= N; i++) {
        delete []u_out[i];
        delete []v_out[i];
        delete []p_out[i];
    }
    delete []u_out;
    delete []v_out;
    delete []p_out;
    u_out = NULL;
    v_out = NULL;
    p_out = NULL;

    for(int i = 0; i <= uN; i++) {
        delete []u[i];
        delete []u_n[i];
        delete []v_proec[i];
    }
    delete []u;
    delete []u_n;
    delete []v_proec;
    u = NULL;
    u_n = NULL;
    v_proec = NULL;

    for(int i = 0; i <= vN; i++) {
        delete []v[i];
        delete []v_n[i];
        delete []u_proec[i];
    }
    delete []v;
    delete []v_n;
    delete []u_proec;
    v = NULL;
    v_n = NULL;
    u_proec = NULL;

    for(int i = 0; i <= pN; i++) {
        delete []p[i];
    }
    delete []p;
    p = NULL;

    delete []x;
    delete []y;
    x = NULL;
    y = NULL;

    return 0;
}

