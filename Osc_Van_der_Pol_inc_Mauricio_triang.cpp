/*
 *  Descricao	: Programa de um contolador do oscilador de Van der Pol com parametro incerto atrav�s de uma EDO de 2a orden por Runge-Kutta de 4a ordem
 *  Arquivo	: Osc_Van_der_Pol_inc_Nauricio.cpp
 *  Autores	: Maurício Gandarela do Espírito Santo
 *  Atualizacao	: 13-07-2020
 */

#include <stdio.h>
#include <math.h>

void rk4_1(double, double, double, double *, double *, double (*)(double, double, double, double));
double func(double, double, double, double);

main()
{
	FILE *outT21, *outT22, *outT23, *outT24, *outT25, *outT26;
	
	double x, v, t, tf, mi, b, lam, csr, cst, ssr, sst, xd, vd, ad, e, de, u, stf;

	if ((outT21 = fopen("outT21.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outT21.dat!\n");
	if ((outT22 = fopen("outT22.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outT22.dat!\n");
	if ((outT23 = fopen("outT23.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outT23.dat!\n");
	if ((outT24 = fopen("outT24.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outT24.dat!\n");
    if ((outT25 = fopen("outT25.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outT25.dat!\n");
    if ((outT26 = fopen("outT26.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outT26.dat!\n");        	
		
	x = 5;		/* Condicao inicial (posi��o inicial) */
	v = 0.2;        /* Condicao inicial (velocidade inicial) */
	t = 0.0;		/* Instante inicial */
	tf = 60.0;		/* Instante final */
	mi = 1.5;       /* Coeficiente n�o linear da for�a de amortecimento*/
    b = 2.0;        /* Ganho do controlador*/
    lam = 10;      /* constante estritamente positiva*/
    csr = 500.0;	/* Taxa de amostragem do controlador */
	cst = 1.0/csr;	/* Periodo de amostragem de controlador */
	ssr = 1000.0;	/* Taxa de amostragem do simulador */
	sst = 1.0/ssr;	/* Periodo de amostragem do simulador */
      
	while (t < tf) {
        xd = (t<=15)?t:(t>15&&t<30)?-t+30:(t>=30&&t<=45)?t-30:-t+60;        /* Trajet�ria desejada(m)*/
        vd = (t<=15)?1:(t>15&&t<30)?-1:(t>=30&&t<=45)?1:-1;        /* Velocidade desejada(m/s)*/
        ad = 0;       /* Acelera��o desejada(m/s�)*/
        e = x - xd;         /* Erro*/
  	    de = v - vd;        /* Derivada do erro*/
  	    u = (-(mi*(1.0 - (x*x))*v) + x + (ad - 2.0*lam*de - ((lam*lam)*e))) * (1.0 / b); /* Lei de controle para o OscVDP*/
  	    stf = t + cst;
        
        while (t < stf) {
            rk4_1(t, t+sst, u, &x, &v, func);
            fprintf(outT21, "%.5e %.5e\n", t, x);
            fprintf(outT22, "%.5e %.5e\n", t, xd);
		    fprintf(outT23, "%.5e %.5e\n", e, de);
		    fprintf(outT24, "%.5e %.5e\n", t, u);
		    fprintf(outT25
	, "%.5e %.5e\n", x, v);
		    fprintf(outT26, "%.5e %.5e\n", xd, vd);
            t += sst;
		}
    }
	
	fclose(outT21);
	fclose(outT22);
	fclose(outT23);
	fclose(outT24);
	fclose(outT25);
	fclose(outT26);
	
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
	double f, mi, b, E;
	E = 1.0 + 0.25*(sin(3*t*fabs(x)));  /* Incerteza com varia��o de 25%*/
	mi = 1.5*E;                    /* Coeficiente n�o linear da for�a de amortecimento*/
    b = 2.0;                         /* Ganho do controlador*/
   
        	
	f = (u * b) - x + (mi*(1.0 - (x*x))*v);

	return f;
}

         
