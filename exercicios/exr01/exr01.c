/* *************************************************************
 *  exercicios/exr01/exr01.c
 * 
 *  autor: G. F. Fornel <guilherme.fornel@ufrgs.br>
 * 
 * 
 *  Programa de Pós-Graduação em Matemática Aplicada
 *  Instituto de Matemática / UFRGS
 * 
 *  MAP0202 - Métodos Numéricos para Eq. Diferenciais
 * 
 *  Exercício 01: Solução do P.V.I.
 * 
 *    du/dt = f(t,u(t)) ,
 * 
 *    f(t,u) = (1 + exp(-t)) * u*u / (1 + u*u*u*u)
 * 
 *    u(t0) = u0
 * 
 *    pelo método de Runge-Kutta RK2 com passo uniforme.
 * 
 *  @param
 *    t0: tempo inicial
 *    tf: tempo final
 *    u0: valor inicial
 *    h : passo (opcional; default 1e-2)
 *
 * 
 *  sistema: Linux (recomendável Ubuntu 20.04LTS ou Mint 20)
 *  compilador: gnu gcc
 *  depurador: gnu gdb (geralmente é instalado junto com gcc)
 * 
 * 
 *  dependências:
 *    - biblioteca gráfica gnuplot (através de um pipe)
 * 
 * 
 *  compilar com 'gcc e01.c -o e01 -lm -Wall'
 *  executar com './e01 t0 tf u0' ou com './e01 t0 tf u0 h'
 * 
 * 
 *  obs. 1: arquivos tasks.json e launch.json para compilação e
 *          depuração no vscode inclusos no diretório .vscode
 * 
 *  obs. 2: memória alocada no heap para evitar stack overflow         
 * 
 * ************************************************************* */

#define VERBOSE 1
#define VERB_STEP 10

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double f(double t, double z)
{
  double zz = z*z;
  return (1 + exp(-t)) * zz / (1 + zz*zz);
}


int main(int argc, char *argv[])
{
  printf("\n\tMAP0202 - Métodos numéricos para Equações Diferenciais\n\n");
  printf("\tExemplo 01: Solução do P.V.I.\n\n");
  printf("\t\tdu/dt = f(t,u(t)) ,\n\n");
  printf("\t\tf(t,u) = (1 + exp(-t)) * u*u / (1 + u*u*u*u)\n\n");
  printf("\t\tu(t0) = u0\n\n");
  printf("\tpelo método de Runge-Kutta RK2 com passo uniforme.\n\n\n");
  printf("(pressione ENTER para continuar...)\n\n");
  getchar();

  if (argc < 4)
    {
      printf("Os parâmetros são t0, tf, u0 [, h]:\n\n");
      printf("t0: tempo inicial\n");
      printf("tf: tempo final\n");
      printf("u0: valor inicial\n");
      printf("h : passo (opcional; default 1e-2)\n");
      exit(-1);
    }

  double t0 = atof(argv[1]);
  double tf = atof(argv[2]);
  double u0 = atof(argv[3]);
  double h = 1e-2;
  if (argc > 4) h = atof(argv[4]);

  size_t nt = ceil( (tf - t0) / h ) + 1; /* dimensão dos arranjos */
  h = (tf - t0) / (nt - 1); /* recalcula o passo */

  #if VERBOSE == 1
    printf("\t! tempo inicial    t0 = %e\n", t0);
    printf("\t! tempo final      tf = %e\n", tf);
    printf("\t! valor inicial    u0 = %e\n", u0);
    getchar();
    printf("\t! passo utilizado  h  = %e\n", h);
    printf("\t! numero de pontos nt = %ld\n\n", nt);
    getchar();
  #endif

  double *t = calloc(nt,sizeof(double));
  if (t == NULL)
    {
      printf("erro: não foi possível alocar memória para t[]\n");
      exit(-1);
    }
  t[0] = t0;
  for (size_t n = 1; n < nt; n++) t[n] = h + t[n-1];

  #if VERBOSE == 1
    printf("\t+ Imprimindo arranjo t[]...\n\n");
    for (size_t n = 0; n < nt; n+=nt/VERB_STEP)
      printf("\tt[%ld] = %e\n", n, t[n]);
    printf("\n");
    getchar();
  #endif

  double *u = calloc(nt,sizeof(double));
  if (u == NULL)
    {
      printf("erro: não foi possível alocar memória para u[]\n");
      free(t);
      exit(-1);
    }
  u[0] = u0;

  #if VERBOSE == 1
    printf("\t+ Avançando no tempo...\n\n");
  #endif

  /* Loop do Runge-Kutta RK2 */ 
  for (size_t n = 1; n < nt; n++)
    {
      #if VERBOSE == 1
        if (n % (nt/VERB_STEP) == 0)
          printf("\tn = %ld\n", n);
      #endif
      double k1n = h * f(t[n-1], u[n-1]);
      double k2n = h * f(t[n-1] + 0.5*h, u[n-1] + 0.5*k1n);
      u[n] = u[n-1] + k2n;
    }

  #if VERBOSE == 1
  printf("\n");
  getchar();
  printf("\t+ Imprimindo arranjo u[] (solução)...\n\n");
    for (size_t n = 0; n < nt; n+=nt/10)
      printf("\tu[%ld] = %e\n", n, u[n]);
  printf("\n");
  getchar();
  #endif


  /* Plota o gráfico da solução com gnuplot */
  FILE *outfp;
  FILE *pltpip;

  if ( (outfp = fopen("exr01.dat","w+b") ) == NULL )
    {
      printf("erro: não foi possível criar o arquivo e01.dat\n");
      free(t);
      free(u);
      exit(-1);
    }
  else
    {
      fprintf(outfp, "# Saída de dados\n#\n");
      fprintf(outfp, "#\tExercício 01: Solução do P.V.I.\n#\n");
      fprintf(outfp, "#\t\tdu/dt = f(t,u(t)) ,\n#\n");
      fprintf(outfp, "#\t\tf(t,u) = (1 + exp(-t)) * u*u / (1 + u*u*u*u)\n#\n");
      fprintf(outfp, "#\t\tu(t0) = u0\n#\n");
      fprintf(outfp, "#\tpelo método de Runge-Kutta RK2 com passo uniforme.\n#\n#\n");
      fprintf(outfp, "#\t! tempo inicial    t0 = %e\n", t0);
      fprintf(outfp, "#\t! tempo final      tf = %e\n", tf);
      fprintf(outfp, "#\t! valor inicial    u0 = %e\n#\n", u0);
      fprintf(outfp, "#\t! passo utilizado  h  = %e\n", h);
      fprintf(outfp, "#\t! numero de pontos nt = %ld\n#\n", nt);
      fprintf(outfp, "# t\t\tu\n");
      for (size_t n = 0; n < nt; n++)
        {
          fprintf(outfp, "%e\t%e\n", t[n], u[n]);
        }
      fclose(outfp);
      free(t);
      free(u);

      if ( (pltpip = popen("gnuplot -persist", "w") ) == NULL )
        {
          printf("erro: não foi possível abrir o pipe para gnuplot\n");
        }
      else
        {
          fprintf(pltpip, "set title \'exr01.dat\' font \',10\'\n");
          fprintf(pltpip, "set style data lines\n");
          fprintf(pltpip, "plot 'exr01.dat'\n");
          pclose(pltpip);
        }
    }


  return 0;
}
