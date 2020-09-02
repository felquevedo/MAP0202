/* *************************************************************
 *  e02/e02.c
 * 
 *  autor: G. F. Fornel <guilherme.fornel@ufrgs.br>
 * 
 * 
 *  Programa de Pós-Graduação em Matemática Aplicada
 *  Instituto de Matemática / UFRGS
 * 
 *  MAP0202 - Métodos Numéricos para Eq. Diferenciais
 * 
 *  Exemplo 02: Solução do P.V.I.
 * 
 *    du/dt = f(t,u(t),v(t)) ,
 * 
 *    dv/dt = g(t,u(t),v(t)) ,
 * 
 *    f(t,u,v) = (1 + exp(-t)) * u * v / (1 + u*v*v) ,
 * 
 *    g(t,u,v) =  v / u,
 * 
 *    u(t0) = u0 , v(t0) = v0
 * 
 *    pelo método de Euler com passo uniforme.
 * 
 *  @param
 *    t0: tempo inicial
 *    tf: tempo final
 *    u0: valor inicial para u
 *    v0: valor inicial para v
 *    h : passo (opcional; default 1e-3)
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
 *  compilar com 'gcc -g e02.c -o e02 -lm -Wall'
 *  executar com './e02 t0 tf u0 v0' ou com './e02 t0 tf u0 v0 h'
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


double f(double t, double u, double v)
{
  double vv = v*v;
  return (1 + exp(-t)) * u*v / (1 + u*vv);
}


double g(double t, double u, double v)
{
  if (u == 0) return INFINITY;
  return v/u;
}


int main(int argc, char *argv[])
{
  printf("\n\tMAP0202 - Métodos numéricos para Equações Diferenciais\n\n");
  printf("\tExemplo 02: Solução do P.V.I.\n\n");
  printf("\t\tdu/dt = f(t,u(t),v(t)) ,\n\n");
  printf("\t\tdv/dt = g(t,u(t),v(t)) ,\n\n");
  printf("\t\tf(t,u,v) = (1 + exp(-t)) * u*v / (1 + u*v*v)\n\n");
  printf("\t\tg(t,u,v) = v/u\n\n");
  printf("\t\tu(t0) = u0 ,  v(t0) = v0\n\n");
  printf("\tpelo método de Euler com passo uniforme.\n\n\n");
  printf("(pressione qualquer tecla para continuar...)\n\n");
  getchar();

  if (argc < 5)
    {
      printf("Os parâmetros são t0, tf, u0, v0, [, h]:\n\n");
      printf("t0: tempo inicial\n");
      printf("tf: tempo final\n");
      printf("u0: valor inicial para u\n");
      printf("v0: valor inicial para v\n");
      printf("h : passo (opcional; default 1e-3)\n");
      exit(-1);
    }

  double t0 = atof(argv[1]);
  double tf = atof(argv[2]);
  double u0 = atof(argv[3]);
  double v0 = atof(argv[4]);
  double h = 1e-3;
  if (argc > 5) h = atof(argv[5]);

  size_t nt = ceil( (tf - t0) / h ) + 1; /* dimensão dos arranjos */
  h = (tf - t0) / (nt - 1); /* recalcula o passo */

  #if VERBOSE == 1
    printf("\t! tempo inicial         t0 = %e\n", t0);
    printf("\t! tempo final           tf = %e\n", tf);
    printf("\t! valor inicial para u  u0 = %e\n", u0);
    printf("\t! valor inicial para v  v0 = %e\n", v0);
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

  double *v = calloc(nt,sizeof(double));
  if (v == NULL)
    {
      printf("erro: não foi possível alocar memória para v[]\n");
      free(t);
      free(u);
      exit(-1);
    }
  v[0] = v0;

  #if VERBOSE == 1
    printf("\t+ Avançando no tempo...\n\n");
  #endif

  for (size_t n = 1; n < nt; n++)
    {
      #if VERBOSE == 1
        if (n % (nt/VERB_STEP) == 0) printf("\tn = %ld\n", n);
      #endif
      u[n] = u[n-1] + h * f(t[n-1], u[n-1], v[n-1]);
      v[n] = v[n-1] + h * g(t[n-1], u[n-1], v[n-1]);
    }

  #if VERBOSE == 1
    printf("\n");
    getchar();
    printf("\t+ Imprimindo arranjos u[] e v[] (solução)...\n\n");
    for (size_t n = 0; n < nt; n+=nt/10)
      printf("\tu[%ld] = %e\n", n, u[n]);
    printf("\n");
    for (size_t n = 0; n < nt; n+=nt/10)
      printf("\tv[%ld] = %e\n", n, v[n]);
    printf("\n");
    getchar();
  #endif


  /* Plota o gráfico da solução com gnuplot */
  FILE *outfp;
  FILE *pltpip;

  if ( (outfp = fopen("e02.dat","w+b") ) == NULL )
    {
      printf("erro: não foi possível criar o arquivo e01.dat\n");
      free(t);
      free(u);
      free(v);
      exit(1);
    }
  else
    {
      fprintf(outfp, "# t  u  v\n");
      for (size_t n = 0; n < nt; n++)
        {
          fprintf(outfp, "%e\t%e\t%e\n", t[n], u[n], v[n]);
        }
      fclose(outfp);
      free(t);
      free(u);
      free(v);

      if ( (pltpip = popen("gnuplot -persist", "w") ) == NULL )
        {
          printf("erro: não foi possível abrir o pipe para gnuplot\n");
        }
      else
        {
          fprintf(pltpip, "set title \'e02.dat\' font \',10\'\n");
          fprintf(pltpip, "set style data lines\n");
          fprintf(pltpip, "set xlabel 't'\n");
          fprintf(pltpip, "plot 'e02.dat' using 1:2 title 'u(t)' with lines, "
                          " '' using 1:3 title 'v(t)' with lines \n");
          pclose(pltpip);
        }
    }


  return 0;
}
