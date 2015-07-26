#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


double f (double px, void * par) {

            double ep=-.9, ed=.0, tsp=1.63, tpd=1.13, tpp=.2;

            struct my_f_params {  double x; double y; double z; double v;};
            struct my_f_params * params  = (struct my_f_params *)par;

             double T2 = (params->x);
             double es = (params->y);
             double eF = (params->z);
             double py = (params->v);

            
                    double xi,argument;
                    double sx=2.*sin(px/2.),sy=2.*sin(py/2.);

                    double data[] = {   ed ,  0.   ,  tpd*sx   ,  -tpd*sy,
                                    0., es    ,  tsp*sx   ,   tsp*sy,
                                tpd*sx, tsp*sx,  ep       ,  -tpp*sx*sy,
                               -tpd*sy, tsp*sy, -tpp*sx*sy,   ep };

               gsl_matrix_view m = gsl_matrix_view_array (data, 4, 4);
               gsl_vector *eval = gsl_vector_alloc (4);
               gsl_matrix *evec = gsl_matrix_alloc (4, 4);
               gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (4);
               gsl_eigen_symmv (&m.matrix, eval, evec, w);
               gsl_eigen_symmv_free (w);
               gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

                     double e = gsl_vector_get (eval, 2);
                     gsl_vector_view evec_i   = gsl_matrix_column (evec, 2);
                     xi = gsl_vector_get(&evec_i.vector,0)*gsl_vector_get(&evec_i.vector,1);

                gsl_vector_free (eval);
                gsl_matrix_free (evec);

                double eeF=e-eF;
                if ( eeF == 0 )
                    argument=1./T2;
                else
                    argument=(tanh(eeF/T2))/eeF;
                    argument*=xi*xi;

	return argument;
}


double g(double py, void * p)
{
                struct my_f_params { double a; double b; double c; };
                struct my_f_params * params  = (struct my_f_params *)p;

                double T2 = (params->a);
                double es = (params->b);
                double eF = (params->c);

                //printf("%g %g %g\n",a,b,c);
                struct your_f_params { double x; double y; double z; double v; };

                struct your_f_params parametri = { T2, es , eF , py };

                double result, error;
                gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000000);
                gsl_function F;


                F.function = &f;
                F.params = &parametri;


                gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-7, 5000000, w, &result, &error);
                gsl_integration_workspace_free (w);

return result;
}


double dvoenintegral(double T2, double es, double eF)
{
                double dvoenintegral, error;
                struct my_f_params { double a; double b; double c; };

                gsl_integration_workspace * v = gsl_integration_workspace_alloc (5000000);

                gsl_function F;
                struct my_f_params params = {T2, es , eF};

                F.function = &g;
                F.params = &params;

                gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-7, 5000000, v, &dvoenintegral, &error);
                gsl_integration_workspace_free (v);

return (dvoenintegral)/(2*M_PI*M_PI);
}


double optim(double L, void * pp)
{

                struct my_f_params { double a; double b; double c; };
                struct my_f_params * params  = (struct my_f_params *)pp;

                double Jsd = (params->a);
                double es = (params->b);
                double eF = (params->c);
                
        double TT=2.*exp(-L);
        double tmp=dvoenintegral(TT,es,eF);

return Jsd*tmp-1;
}

double nameri_nula( double Jsd, double es, double eF)
{
       int status;
       double Kelvin=1/11604.;
       struct my_f_params { double a; double b; double c; };
       struct my_f_params params = { Jsd, es , eF};
       int iter = 0, max_iter = 100;
       const gsl_root_fsolver_type *T;
       gsl_root_fsolver *s;
       double r = 0;
       double x_lo = 1E-6, x_hi = 7;
       gsl_function F;

       F.function = &optim;
       F.params = &params;


       T = gsl_root_fsolver_brent;
       s = gsl_root_fsolver_alloc (T);
       gsl_root_fsolver_set (s, &F, x_lo, x_hi);

       //printf ("using %s method\n", gsl_root_fsolver_name (s));

       //printf ("%5s [%9s, %9s] %9s %10s\n", "iter", "lower", "upper", "root", "err");

       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

           //if (status == GSL_SUCCESS)
           //printf ("Converged:\n");

           //printf ("%5d [%.7f, %.7f] %.7f %+.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
         }
       while (status == GSL_CONTINUE && iter < max_iter);

	//printf("T_c=%g\n",exp(-r)/Kelvin);

       gsl_root_fsolver_free (s);

       return (2*exp(-r)/Kelvin);
}

double granichni_tochki(double a, double b,double es)
{

    double sx=2.*sin(a/2.);
    double sy=2.*sin(b/2.);

    double ep=-.9, ed=.0, tsp=1.63, tpd=1.13, tpp=.2;
    
    double data[] = {   ed ,  0.   ,  tpd*sx   ,  -tpd*sy,
                            0., es    ,  tsp*sx   ,   tsp*sy,
                        tpd*sx, tsp*sx,  ep       ,  -tpp*sx*sy,
                       -tpd*sy, tsp*sy, -tpp*sx*sy,   ep };

       gsl_matrix_view m = gsl_matrix_view_array (data, 4, 4);
       gsl_vector *eval = gsl_vector_alloc (4);
       gsl_matrix *evec = gsl_matrix_alloc (4, 4);
       gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (4);
       gsl_eigen_symmv (&m.matrix, eval, evec, w);
       gsl_eigen_symmv_free (w);
       gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

         gsl_vector_free (eval);
         gsl_matrix_free (evec);

         return ( gsl_vector_get (eval, 2) );
}

double pd(double e, double es)
{
    double ep=-.9, ed=.0, tsp=1.63, tpd=1.13, tpp=.2;
    double A,B,C,xd,eep,eed,ees,ppd;

eep=e-ep;
eed=e-ed;
ees=e-es;

A=16.*(4.*tpd*tpd*tsp*tsp + 2.*tsp*tsp*tpp*eed-2.*tpd*tpd*tpp*ees-tpp*tpp*eed*ees);
B=-4.*eep*(tsp*tsp*eed+tpd*tpd*ees);
C=ees*eep*eep*eed;
xd=(-B+sqrt(B*B-A*C))/A;
//printf("xd=%g\n",xd);
ppd=2*asin(sqrt(xd));
//printf("xd=%g ppd=%g e=%g es=%g\n",xd,ppd,e,es);
return ppd;
}



double fermi(double px, void * p) 
{
                struct my_f_params { double a; double b; };
                struct my_f_params * params  = (struct my_f_params *)p;


                double eF = (params->a);
                double es = (params->b);
                


    double A,B,C,eep,eed,ees,e,F;
    double ep=-.9, ed=.0, tsp=1.63, tpd=1.13, tpp=.2;

    e=eF;
eep=e-ep;
eed=e-ed;
ees=e-es;


A=16.*(4.*tpd*tpd*tsp*tsp + 2.*tsp*tsp*tpp*eed-2.*tpd*tpd*tpp*ees-tpp*tpp*eed*ees);
B=-4.*eep*(tsp*tsp*eed+tpd*tpd*ees);
C=ees*eep*eep*eed;

double x=sin(px/2.)*sin(px/2.);
F=(px-2*asin(sqrt(- (C+B*x)/(A*x+B)  )));

return F;
}



double ffactor( double eF , double es )
{
                double res, error;
                struct my_f_params { double a; double b; };

                gsl_integration_workspace * v = gsl_integration_workspace_alloc (1000);

                gsl_function F;
                struct my_f_params params = {eF , es};

                F.function = &fermi;
                F.params = &params;

                gsl_integration_qag (&F, pd(eF,es) , M_PI , 0, 1e-7, 1000, 6 , v, &res, &error);
                gsl_integration_workspace_free (v);


return  (8*res)/(4*M_PI*M_PI) ;
}

double ffactornula( double eF , void * p )
{
                double res, error;
                double es = *(double *) p;

                struct my_f_params { double a; double b; };

                gsl_integration_workspace * v = gsl_integration_workspace_alloc (100000);

                gsl_function F;
                struct my_f_params params = {eF , es};

                F.function = &fermi;
                F.params = &params;

                gsl_integration_qags (&F, pd(eF,es) , M_PI , 0, 1e-7, 100000,  v, &res, &error);
                gsl_integration_workspace_free (v);


return  (8*res)/(4*M_PI*M_PI) - 0.66;
}

double ferminula( double evhs, double etop, double es)
{
    int status;

       int iter = 0, max_iter = 100;
       const gsl_root_fsolver_type *T;
       gsl_root_fsolver *s;
       double r = 0;
       double x_lo = evhs+1E-7, x_hi = etop-1E-7;
       gsl_function F;

       F.function = &ffactornula;
       F.params = &es;

       T = gsl_root_fsolver_brent;
       s = gsl_root_fsolver_alloc (T);
       gsl_root_fsolver_set (s, &F, x_lo, x_hi);

       //printf ("using %s method\n", gsl_root_fsolver_name (s));

       //printf ("%5s [%9s, %9s] %9s %10s\n", "iter", "lower", "upper", "root", "err");

       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

           //if (status == GSL_SUCCESS)
             //printf ("Converged:\n");

           //printf ("%5d [%.7f, %.7f] %.7f %+.7f\n", iter, x_lo, x_hi, r, x_hi - x_lo);
         }
       while (status == GSL_CONTINUE && iter < max_iter);

	
       gsl_root_fsolver_free (s);

       return r;
}

int main(void)
{
    double Kelvin=1/11604.;
    double es=6.0, eF=1.301531002,tsp=1.63,ep=-.9,tpd=1.13;
    double Jsd,Tc,evhs,etop,s,ss;
    //double tmp;
/*
    Jsd=1/dvoenintegral(90*Kelvin,es,eF);
    printf ("Jsd= %g\n",Jsd);
    Tc=nameri_nula(Jsd,es,eF);
    printf ("Tc= %g\n",Tc);
    evhs=granichni_tochki(0,M_PI,es);
    printf("evhs=%g\n",evhs);
    etop=granichni_tochki(M_PI,M_PI,es);
    printf("etop=%g\n",etop);
    ff=ffactor(evhs, es);
    printf("filling factor(evhs)=%g\n",ff);
    printf("%g\n",ferminula(evhs,etop,es));
*/
    Jsd=1/dvoenintegral(90*Kelvin,es,eF);

    for(es=5;es<8.5;es+=0.1)
    {
    evhs=granichni_tochki(0,M_PI,es);
    etop=granichni_tochki(M_PI,M_PI,es);
    eF=ferminula(evhs,etop,es);
    Tc=nameri_nula(Jsd,es,eF);
    //s=((es-eF)*(eF-ep))/(4*tsp*tsp);
    ss=((es-eF)*(eF-ep))/(4*tsp*tpd);
    //printf("Tc=%g es=%g eF=%g\n",Tc,es,eF);
    printf("%g %g\n",ss,log(Tc));
    //printf("%g %g\n",1/(2*(1+s)),Tc );
    }


    return 0;
}
