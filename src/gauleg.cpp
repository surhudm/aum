#include "gauleg.h"

/// Adapted by Surhud More from gsl_integration_glfixed routine in gsl
void gauleg(const double x1, const double x2, double x[], double w[], int n)
{
    double A=0.5*(x2-x1);
    double B=0.5*(x2+x1);

    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc (n);

    int m = (n + 1) >> 1;

    int put=0;

    if(n&1){ /*n - odd*/
	for(int i=1;i<m;i++)
	{
	    x[put]=B+A*t->x[i];
	    w[put]=A*t->w[i];
	    put++;

	    x[put]=B-A*t->x[i];
	    w[put]=A*t->w[i];
	    put++;
	}
    }else{ /*n - even*/
	for(int i=0;i<m;i++)
	{
	    x[put]=B+A*t->x[i];
	    w[put]=A*t->w[i];
	    put++;

	    x[put]=B-A*t->x[i];
	    w[put]=A*t->w[i];
	    put++;
	}
    }

    gsl_integration_glfixed_table_free(t);

}
