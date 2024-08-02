# include <stdio.h>
# include <math.h>
# include "complex.h"


int main()
{
    int i,j,k;
    double x,y,z;
    
    fcomplex a,b,c;
    
    a=Complex(1.0,3.0);
    
    b.r=0.5;
    b.i=0.5;


    
    c=Cmul(a,b);
    
    printf("(%lf+%lf i)*(%lf+%lf i)=%lf+%lf i\n\n",a.r,a.i,b.r,b.i,c.r,c.i);
    
    c=Cadd(a,b);
    
    printf("(%lf+%lf i)+(%lf+%lf i)=%lf+%lf i\n\n",a.r,a.i,b.r,b.i,c.r,c.i);
    
    c=RCmul(3.,b);
    printf("3*(%lf+%lf i)=%lf+%lf i\n\n",b.r,b.i,c.r,c.i);
}
