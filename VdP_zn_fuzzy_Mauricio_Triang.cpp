/*
 *  Descricao	: Programa de um contolador do oscilador de Van der Pol com zona morta e controlador Fuzzy
 *  Arquivo	: VdP_zn_fuzzy.cpp
 *  Autores	: Maurício Gandarela
 *  Atualizacao	: 2011-07-03
 */

#include <stdio.h>
#include <math.h>

void rk4_1(double, double, double, double *, double *, double (*)(double, double, double, double));
double func(double, double, double, double);
double fuzzy(double, double);

main()
{
	FILE *out51, *out52, *out53, *out54, *out55, *out56, *out57;
	
	double x, v, t, tf, mi, b, lam, csr, cst, ssr, sst, xd, vd, ad, e, de, e2, de2, u, du, stf;

	if ((out51 = fopen("out51.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out51.dat!\n");
	if ((out52 = fopen("out52.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out52.dat!\n");
	if ((out53 = fopen("out53.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out53.dat!\n");
	if ((out54 = fopen("out54.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out54.dat!\n");
    if ((out55 = fopen("out55.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out55.dat!\n");
    if ((out56 = fopen("out56.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out56.dat!\n"); 
    if ((out57 = fopen("out57.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo out57.dat!\n");       	
		
	x = 5;		/* Condicao inicial (posiï¿½ï¿½o inicial) */
	v = 0.2;        /* Condicao inicial (velocidade inicial) */
	t = 0.0;		/* Instante inicial */
	tf = 60.0;		/* Instante final */
	mi = 1.5;       /* Coeficiente nï¿½o linear da forï¿½a de amortecimento*/
    b = 2.0;        /* Ganho do controlador*/
    lam = 10;      /* constante estritamente positiva*/
    csr = 500.0;	/* Taxa de amostragem do controlador */
	cst = 1.0/csr;	/* Periodo de amostragem de controlador */
	ssr = 1000.0;	/* Taxa de amostragem do simulador */
	sst = 1.0/ssr;	/* Periodo de amostragem do simulador */
      
	while (t < tf) {
        xd = (t<=15)?t:(t>15&&t<30)?-t+30:(t>=30&&t<=45)?t-30:-t+60;        /* Trajetï¿½ria desejada(m)*/
        vd =(t<=15)?1:(t>15&&t<30)?-1:(t>=30&&t<=45)?1:-1;        /* Velocidade desejada(m/s)*/
        ad = 0;       /* Aceleraï¿½ï¿½o desejada(m/sï¿½)*/
        e = x - xd;         /* Erro*/
  	    de = v - vd;        /* Derivada do erro*/
  	    u = (-(mi*(1.0 - (x*x))*v) + x + (ad - 2.0*lam*de - ((lam*lam)*e))) * (1.0 / b); /* Lei de controle para o OscVDP*/
  	                
  	    
  	  
  	        // Controle Fuzzy
        e2 = xd - x;
        de2 = vd - v;
  	    du = fuzzy(e2, de2);
  	    u += du;               
  	    
  	    stf = t + cst;
        
        while (t < stf) {
            rk4_1(t, t+sst, u, &x, &v, func);
            fprintf(out51, "%.5e %.5e\n", t, x);
            fprintf(out52, "%.5e %.5e\n", t, xd);
		    fprintf(out53, "%.5e %.5e\n", e, de);
		    fprintf(out54, "%.5e %.5e\n", t, u);
		    fprintf(out55, "%.5e %.5e\n", t, du);
		    fprintf(out56, "%.5e %.5e\n", t, e);
		    fprintf(out57, "%.5e %.5e\n", t, de);
            t += sst;
		}
    }
	
	fclose(out51);
	fclose(out52);
	fclose(out53);
	fclose(out54);
	fclose(out55);
	fclose(out56);
	fclose(out57);
	
	return 0;
}
// Runge-Kutta de 4a ordem para EDOs de 2a ordem
void rk4_1(double t, double tf, double u, double *ptx, double *ptv, double (*func)(double, double, double, double))
{
	double x, v, h, k1, k2, k3, k4, q1, q2, q3, q4;
	x = *ptx;
	v = *ptv;
	h = tf - t;
	k1 = h * func(t, x, v, u);
	q1 = h * v;
	k2 = h * func(t + h/2.0, x + q1/2.0, v + k1/2.0, u);
	q2 = h * (v + k1/2.0);
	k3 = h * func(t + h/2.0, x + q2/2.0, v + k2/2.0, u);
	q3 = h * (v + k2/2.0);
	k4 = h * func(t + h, x + q3, v + k3, u);
	q4 = h * (v + k3);
	v += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
	x += (q1 + 2.0 * q2 + 2.0 * q3 + q4) / 6.0;
	*ptx = x;
	*ptv = v;
}
                 
// Funcao v' = (u * b) - x + (mi*(1.0 - (x*x))*v)
double func(double t, double x, double v, double u)
{
		double f, mi, b, dl, dr, ni, E;
	E = 1.0 + 0.2*(sin(3*t*fabs(x)));  /* Incerteza com variaï¿½ï¿½o de 20%*/
    mi = 1.5*E;       /* Coeficiente nï¿½o linear da forï¿½a de amortecimento*/
    b = 2.0;        /* Ganho do controlador*/
    dl = -0.25;      /* Limite inferior da zona morta*/
    dr = 0.25;       /* Limite superior da zn*/
        
    if (u >= dr)
		ni = u-dr;
	else if ( u > dl)
		ni = 0.0;
	else
		ni = u-dl;
        	
	f = (ni * b) - x + (mi*(1.0 - (x*x))*v);

	return f;
}

// Controlador Fuzzy
double fuzzy(double e, double de)
{
   double r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,prod,som,emin,emax,demin,demax;
   double min1,max1,min2,max2,min3,max3,min4,max4,min5,max5,min6,max6,min7,max7,du;
   double pmin1,pmax1,pmin2,pmax2,pmin3,pmax3,pmin4,pmax4,pmin5,pmax5, pmin6,pmax6, pmin7,pmax7;
   int n=7,i;
   int m=7,j;
   int p=49,k;
   double /*ER[n], DR[n],*/ GP[n], GPP[n];
   double W[p];
 /* Tabela do Endrew*/     /*double U[49] = {-20.0,    -20.0,    -20.0,   -20.0,    -20.0,    -20.0,    -5.0,
                                           -20.0,    -20.0,    -5.0,   -5.0,    -5.0,    -5.0,    -2.5,
                                           -20.0,    -5.0,    -5.0,   -2.5,    -2.5,    -2.5,     2.5,
                                           -2.5,    -2.5,    -2.5,    0.0,     2.5,     2.5,     2.0,
                                           -2.5,    2.5,     2.5,     2.5,     5.0,     5.0,     20.0,
                                            2.5,    5.0,     5.0,     5.0,     5.0,     20.0,     20.0,
                                            5.0,    20.0,     20.0,     20.0,     20.0,     20.0,     20.0};*/
 /* Tabela do Prof. Wallace*/  //double U[49] = {-0.5, -0.5, -0.1, -0.05, 0.1, 0.5, 0.5,/**/ -0.5, -0.1, -0.1, -0.05, 0.1, 0.1, 0.5, /**/-0.1, -0.05, -0.05, -0.05, 0.05, 0.05, 0.1, /**/0.05, 0.05, 0.0, 0.0, 0.0, -0.05, -0.05, /**/0.1, 0.05, 0.05, 0.05, -0.05, -0.05, -0.1, /**/0.5, 0.1, 0.1, 0.05, -0.1, -0.1, -0.5, /**/0.5, 0.5, 0.1, 0.05, -0.1, -0.5, -0.5};
 /* Tabela do Livro*/  double U[49] = {-0.20,    -0.20,    -0.20,    -0.005,    -0.005,    -0.025,    0.0,
                                         -0.20,    -0.20,    -0.005,     -0.005,    -0.025,    0.0,     0.025, 
                                         -0.20,    -0.005,     -0.005,     -0.025,    0.0,     0.025,     0.005,
                                         -0.005,     -0.005,     -0.025,     0.0,     0.025,     0.005,     0.005,
                                         -0.005,     -0.025,     0.0,      0.025,     0.005,     0.005,     0.20, 
                                         -0.025,     0.0,      0.025,      0.005,     0.005,     0.20,    0.20,
                                          0.0,     0.025,      0.005,      0.005,     0.20,    0.20,    0.20};
      
      emin = -0.001;             /* Erro mínimo*/
      emax = 0.001;              /* Erro máximo*/
      demin = -0.01;              /* Derivada do erro mínimo*/
      demax = 0.01;               /* Derivada do erro máximo*/     

    double  ER[7] = {-0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1};          /* Vetor erro com n componentes*/
    double  DR[7] = {-0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1};          /* Vetor derivada do erro com n componentes*/

//   for (i = 0; i <= (n-1); i++) {
//       ER[i] = emin + (double)i*((emax-emin)/((double)n-1));          /* Vetor erro com n componentes*/
//       DR[i] = demin + (double)i*((demax-demin)/((double)n-1));        /* Vetor derivada do erro com n componentes*/
//       }         
                                  
    //grau de pertinencia do erro
             r1=(ER[1]-e)/(ER[1]-ER[0]);    
             r2=(e-ER[0])/(ER[1]-ER[0]);    
             r3=(ER[2]-e)/(ER[2]-ER[1]);    
             r4=(e-ER[1])/(ER[2]-ER[1]);    
             r5=(ER[3]-e)/(ER[3]-ER[2]);    
             r6=(e-ER[2])/(ER[3]-ER[2]);    
             r7=(ER[4]-e)/(ER[4]-ER[3]);    
             r8=(e-ER[3])/(ER[4]-ER[3]);
             r9=(ER[5]-e)/(ER[5]-ER[4]);    
             r10=(e-ER[4])/(ER[5]-ER[4]);
             r11=(ER[6]-e)/(ER[6]-ER[5]);    
             r12=(e-ER[5])/(ER[6]-ER[5]);    
               min1=(r1<1?r1:1);
               max1=(min1>0?min1:0);    
                  GP[0]=max1;
               min2=(r2<r3?r2:r3);
               max2=(min2>0?min2:0);
                  GP[1]=max2;
               min3=(r4<r5?r4:r5);
               max3=(min3>0?min3:0);  
                  GP[2]=max3;
               min4=(r6<r7?r6:r7);
               max4=(min4>0?min4:0);  
                  GP[3]=max4;
               min5=(r8<r9?r8:r9);
               max5=(min5>0?min5:0);
                  GP[4]=max5;
               min6=(r10<r11?r10:r11);
               max6=(min6>0?min6:0);
                  GP[5]=max6;
               min7=(r12<1?r12:1);
               max7=(min7>0?min7:0);
                  GP[6]=max7;      
    //grau de pertinencia da derivada do erro
               r1=(DR[1]-de)/(DR[1]-DR[0]);    
               r2=(de-DR[0])/(DR[1]-DR[0]);    
               r3=(DR[2]-de)/(DR[2]-DR[1]);    
               r4=(de-DR[1])/(DR[2]-DR[1]);    
               r5=(DR[3]-de)/(DR[3]-DR[2]);    
               r6=(de-DR[2])/(DR[3]-DR[2]);    
               r7=(DR[4]-de)/(DR[4]-DR[3]);    
               r8=(de-DR[3])/(DR[4]-DR[3]);
               r9=(DR[5]-de)/(DR[5]-DR[4]);    
               r10=(de-DR[4])/(DR[5]-DR[4]);
               r11=(DR[6]-de)/(DR[6]-DR[5]);    
               r12=(de-DR[5])/(DR[6]-DR[5]);    
               pmin1=(r1<1?r1:1);
               pmax1=(pmin1>0?pmin1:0);    
                  GPP[0]=pmax1;
               pmin2=(r2<r3?r2:r3);
               pmax2=(pmin2>0?pmin2:0);
                  GPP[1]=pmax2;
               pmin3=(r4<r5?r4:r5);
               pmax3=(pmin3>0?pmin3:0);  
                  GPP[2]=pmax3;
               pmin4=(r6<r7?r6:r7);
               pmax4=(pmin4>0?pmin4:0);  
                  GPP[3]=pmax4;
               pmin5=(r8<r9?r8:r9);
               pmax5=(pmin5>0?pmin5:0);
                  GPP[4]=pmax5;
               pmin6=(r10<r11?r10:r11);
               pmax6=(pmin6>0?pmin6:0);
                  GPP[5]=pmax6;
               pmin7=(r12<1?r12:1);
               pmax7=(pmin7>0?pmin7:0);
                  GPP[6]=pmax7;      
                  
      //defuzzificaçao
      k=0;
          for (i = 0; i <= (m-1); i++){          
            for (j=0;j<=(n-1);j++){
                W[k]=(GP[i]<GPP[j]?GP[i]:GPP[j]);
                k += 1;
            }          
          }

          for (i = 0; i <= (p-1); i++){          
                prod = prod + (U[i] * W[i]);
                som = som + W[i];
            }          
                du = prod/som;
                som = 0.0;
                prod = 0.0;
      
   return du;
}
         
