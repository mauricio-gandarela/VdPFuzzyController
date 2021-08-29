/*
 *  Descricao	: Programa de um contolador do oscilador de Van der Pol com zona morta atrav�s de uma EDO de 2a orden por Runge-Kutta de 4a ordem
 *  Arquivo	: Osc_Van_der_Pol_Nauricio.cpp
 *  Autores	: Maurício Gandarela do Espírito Santo
 *  Atualizacao	: 13-07-2020
 */
#include <stdio.h>
#include <math.h>

void rk4_1(double, double, double, double *, double *, double (*)(double, double, double, double));
double func(double, double, double, double);

main()
{
	FILE *outM31, *outM32, *outM33, *outM34, *outM35, *outM36, *outM37;
	
	double x, v, t, tf, mi, b, lam, csr, cst, ssr, sst, xd, vd, ad, e, de, u, stf;

	if ((outM31 = fopen("outM31.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM31.dat!\n");
	if ((outM32 = fopen("outM32.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM32.dat!\n");
	if ((outM33 = fopen("outM33.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM33.dat!\n");
	if ((outM34 = fopen("outM34.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM34.dat!\n");
    if ((outM35 = fopen("outM35.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM35.dat!\n");
    if ((outM36 = fopen("outM36.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM36.dat!\n"); 
    if ((outM37 = fopen("outM37.dat", "w")) == NULL)
		printf("Nao foi possivel abrir o arquivo outM37.dat!\n");             	
		
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
            fprintf(outM31, "%.5e %.5e\n", t, x);
            fprintf(outM32, "%.5e %.5e\n", t, xd);
		    fprintf(outM33, "%.5e %.5e\n", e, de);
		    fprintf(outM34, "%.5e %.5e\n", t, u);
		    fprintf(outM35, "%.5e %.5e\n", x, v);
		    fprintf(outM36, "%.5e %.5e\n", xd, vd);
		    fprintf(outM37, "%.5e %.5e\n", t, e);
            t += sst;
		}
    }
	
	fclose(outM31);
	fclose(outM32);
	fclose(outM33);
	fclose(outM34);
	fclose(outM35);
	fclose(outM36);
	fclose(outM37);
	
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
	E = 1.0 + 0.2*(sin(3*t*fabs(x)));  /* Incerteza com varia��o de 20%*/
    mi = 1.5*E;       /* Coeficiente n�o linear da for�a de amortecimento*/
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

         
