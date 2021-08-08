// Program to compute the DFT  of a signal using Decimation in time FFT Algorithm
#include<stdio.h>
#include<math.h>
double PI = acos(-1);
struct complex
{ double real;
  double imag;
};
struct complex exponential(int n, int N )
{ struct complex z;
  z.real = cos(2*PI*n/N);
  z.imag = sin(2*PI*n/N);
  return z;
};
//function to compute bit reversal
void bit_reversal(int b[],int N)
{    int i,a,k;
     double M;
     M = log(N)/log(2);
     for(i=0;i<N;i++)
       { a=i;
         b[i]=0;
         k=0;
         while(a!=0)
        {    b[i]+=(int)(a%2)*pow(2,M-k-1);
             a=(int)(a/2);
             k++;
        }
     }
}
//function to compute the output of a butterfly diagram
void butterfly(struct complex a[], struct complex A[])
{    (*A).real = (*a).real+((*(a+2)).real * (*(a+1)).real-(*(a+2)).imag * (*(a+1)).imag);
     (*A).imag = (*a).imag+(*(a+2)).imag * (*(a+1)).real+(*(a+2)).real * (*(a+1)).imag;
     (*(A+1)).real = (*a).real-((*(a+2)).real * (*(a+1)).real-(*(a+2)).imag * (*(a+1)).imag);
     (*(A+1)).imag = (*a).imag-((*(a+2)).imag * (*(a+1)).real+(*(a+2)).real * (*(a+1)).imag);
}

int main()
{   int i,x[1024],b[1024],N,L,k,m,j;
    //cmplx x[N]={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};
    double M,X_mag[1024];
    char sign;
    struct complex X[1024],temp[1024];
    struct complex A[2],a[3],w;
    printf("Enter the size of the input signal  ");
    scanf("%d",&L);
    //input the signal
    printf("Enter the input signal x(n)= ");
    for(i=0;i<L;i++)
        scanf("%d",&x[i]);
    printf("Enter the size of the DFT signal  ");
    scanf("%d",&N);
    M = log(N)/log(2);// M=number of stages
    printf("M=%lf ",M);
    if(L<N)
       for(i=L;i<N;i++)
            x[i]=0;
    //Display the signal
    printf("x(n)=[ ");
    for(i=0;i<N;i++)
        printf("%d ",x[i]);
    printf("] \n");
    //Bit reversal
    bit_reversal(b,N);
    for(i=0;i<N;i++)
        temp[i].real=x[b[i]];
    //Compute the dft using fft decimation in time algorithm
    for(m=1;m<M+1;m++) // iterating through stages in butterfly
       {for(i=0;i<(int)(pow(2,M-m));i++)  // iterating through sections in a stage
            for(j=0;j<(int)(pow(2,m-1));j++)  // iterating through butterflies in a section
               {a[0]=temp[(int)(i*pow(2,m)+j)];
                a[1]=temp[(int)(i*pow(2,m)+j+pow(2,m-1))];
                w=exponential(-j,(int)pow(2,m));
                a[2]=w;
                butterfly(a,A);
                X[(int)(i*pow(2,m)+j)]=A[0];
                X[(int)(i*pow(2,m)+j+pow(2,m-1))]=A[1];
                }
        for(i=0;i<N;i++)
            temp[i]=X[i];
       }
   //Display the DFT
    printf("X(k)=[ ");
    for(k=0;k<N;k++)
        {sign = X[k].imag<0? '-':'+';
         printf("%f %c %fi   ",X[k].real,sign,fabs(X[k].imag));
        }
    printf("]\n ");
    //Display the magnitude of DFT
    printf("X_mag=[ ");
    for(k=0;k<N;k++)
        printf(" %lf   ",sqrt(pow(X[k].real,2)+pow(X[k].imag,2)));
    printf("] ");
return 0;
}
