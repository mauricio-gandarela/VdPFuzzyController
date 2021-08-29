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
	FILE *outM21, *outM22, *outM23, *outM24, *outM25, *outM26;
	
	double x, v, t, tf, mi, b, lam, csr, cst, ssr, sst, xd, vd, ad, e, de, u, stf;

	if ((outM21 = fopen("outM21.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM21.dat!\n");
	if ((outM22 = fopen("outM22.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM22.dat!\n");
	if ((outM23 = fopen("outM23.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM23.dat!\n");
	if ((outM24 = fopen("outM24.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM24.dat!\n");
    if ((outM25 = fopen("outM25.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM25.dat!\n");
    if ((outM26 = fopen("outM26.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM26.dat!\n");        	
		
	x = 0.2;		/* Condicao inicial (posi��o inicial) */
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
        xd = sin(3*t);        /* Trajet�ria desejada(m)*/
        vd = 3*cos(3*t);        /* Velocidade desejada(m/s)*/
        ad = -9*sin(3*t);       /* Acelera��o desejada(m/s�)*/
        e = x - xd;         /* Erro*/
  	    de = v - vd;        /* Derivada do erro*/
  	    u = (-(mi*(1.0 - (x*x))*v) + x + (ad - 2.0*lam*de - ((lam*lam)*e))) * (1.0 / b); /* Lei de controle para o OscVDP*/
  	    stf = t + cst;
        
        while (t < stf) {
            rk4_1(t, t+sst, u, &x, &v, func);
            fprintf(outM21, "%.5e %.5e\n", t, x);
            fprintf(outM22, "%.5e %.5e\n", t, xd);
		    fprintf(outM23, "%.5e %.5e\n", e, de);
		    fprintf(outM24, "%.5e %.5e\n", t, u);
		    fprintf(outM25, "%.5e %.5e\n", x, v);
		    fprintf(outM26, "%.5e %.5e\n", xd, vd);
            t += sst;
		}
    }
	
	fclose(outM21);
	fclose(outM22);
	fclose(outM23);
	fclose(outM24);
	fclose(outM25);
	fclose(outM26);
	
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

         
