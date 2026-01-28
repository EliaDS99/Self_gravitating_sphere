#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EXIT_FAILURE 1
#define R_max 6.96e8
#define G 6.67e-11
#define kb 1.38e-23
#define mH 1.67e-27
#define w 2.36e2
#define zzz 3.68e21
#define c 7.538e-8
#define iii 2.17896e-18

typedef struct {
  double u0, T, x0, rho0, dr, scelta, n, d, r, mi;
  int max_iter, a;
}par;

void rungek4a (double* x, double* rho, double* u, par p, FILE* fp);
void rungek4b (double* x, double* v, double* f, par p, FILE* fp);
void nube (par p);

int main(){
  double* x;
  double* rho;
  double* u;
  double* f;
  double* v;
  int nsteps;
  par pp;
  FILE* fp;
  
  printf("Ciao sono Agata, una volta che mi avrai fornito i valori richiesti sarò in grado di mostrarti come variano densità e massa di un sistema autogravitante in funzione del raggio\n");
  printf("\n");
  printf("Per prima cosa mi devi dire se hai intenzione di visualizzare la sfera isoterma (in tal caso premi 1), o un altro tipo di soluzione (in tal caso premi 2), o se preferisci vedere il potenziale collasso di una nube (premi 3)\n");
  do{
  scanf("%d", &pp.a);
  if(pp.a!=1 && pp.a!=2 && pp.a!=3){
    printf("Il valore inserito non è valido. Reinserire\n");
  }
  }while(pp.a!=1 && pp.a!=2 && pp.a!=3);
  if (pp.a==1){
    printf("Per favore inserisci in quest'ordine i seguenti dati della stella (dovrai usare le unità di misura del SI):\n M(0)\n T\n R(0)\n Rho(0)\n dr\n numero di iterazioni\n mi\n");
    printf("Per una corretta costruzione suggerisco i seguenti valori:\n M(0)=0\n T=15000000\n R(0)=0\n Rho(0)=150000\n dr=1000\n numero di iterazioni 1000000\n 0.5<mi<2\n");
    scanf("%lf %lf %lf %lf %lf %d %lf", &pp.u0, &pp.T, &pp.x0, &pp.rho0, &pp.dr, &pp.max_iter, &pp.mi);
  } else if (pp.a==2){
    printf("Per favore inserisci in quest'ordine i seguenti dati della stella (dovrai usare le unità di misura del SI):\n R(0)\n Rho(0)\n dr\n numero di iterazioni\n n\n");
    printf("Per una corretta costruzione suggerisco i seguenti valori:\n R(0)=0\n Rho(0)=150000\n dr=1000\n numero di iterazioni 1000000\n 0<=n<10\n");
    scanf("%lf %lf %lf %d %lf", &pp.x0, &pp.rho0, &pp.dr, &pp.max_iter, &pp.n);
  } else if (pp.a==3){
    printf("Per favore inserisci in quest'ordine i seguenti dati della nube (dovrai usare i ly, klvin e atomi/m^3):\n R\n T\n n\n");
    printf("Per una corretta costruzione suggerisco i seguenti valori:\n R=1.5ly\n T=10\n 10^6<n<10^10\n");
    scanf("%lf %lf %lf", &pp.r, &pp.T, &pp.d);
  }
  if(pp.a==1 || pp.a==2){
  nsteps = (int)(pp.max_iter*pp.dr+0.5);
  x = (double*)calloc(nsteps+1, sizeof(double));
  rho = (double*)calloc(nsteps+1, sizeof(double));
  u = (double*)calloc(nsteps+1, sizeof(double));
  f = (double*)calloc(nsteps+1, sizeof(double));
  v = (double*)calloc(nsteps+1, sizeof(double));
  if(x==NULL){
    printf("Errore\n");
    exit(EXIT_FAILURE);
  }
  }
  if((fp=fopen("uscita.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
  if (pp.a==1){
    rungek4a(x, rho, u, pp, fp);
    printf("I valori si trovano su un file chiamato uscita.dat\n");
  }
  if (pp.a==2){
    rungek4b(x, f, v, pp, fp);
    printf("I valori si trovano su un file chiamato uscita.dat\n");
  }
  if(pp.a==3){
    nube(pp);
  }
  fclose(fp);
}

void rungek4a(double* x, double* rho, double* u, par pp, FILE* fp){
  double l1, l2, l3, l4, l, g1, g2, g3, g4, g, dmi=(pp.mi-1)/(R_max/pp.dr), mi=pp.mi;
  double z = (1.210e-4)*mi/pp.T, k = 8.38e-10, dr = 0.01, aus = 0, app, aiu; 
  int i;
  x[0] = pp.x0;
  rho[0] = pp.rho0;
  u[0] = pp.u0;
  for (i = 0; i < pp.max_iter; i++){
    if(i == 0){
      app = -(2*M_PI*G*pp.rho0/z)*aus*dr;
      aiu = pp.rho0 + exp(app);
      aus += dr;
      app = -(2*M_PI*G*pp.rho0/z)*aus*dr;
      rho[i] = aiu + exp(app);
    }
    else {
    x[i] = x[i-1] + pp.dr;
    g1 = -z/(((x[i]))*((x[i])))*((rho[i-1]))*((u[i-1]));
    l1 = ((x[i])*(x[i]))*(rho[i-1])*k;
    g2 = -z/(((x[i]))*((x[i])))*((rho[i-1]))*(((u[i-1]))+(l1*pp.dr*0.5));
    l2 = ((x[i])*(x[i]))*((rho[i-1])+g1*pp.dr*0.5)*k;
    g3 = -z/(((x[i]))*((x[i])))*((rho[i-1]))*(((u[i-1]))+(l2*pp.dr*0.5));
    l3 = ((x[i])*(x[i]))*((rho[i-1])+g2*pp.dr*0.5)*k;
    g4 = -z/(((x[i]))*((x[i])))*((rho[i-1]))*(((u[i-1]))+(l3*pp.dr));
    l4 = ((x[i])*(x[i]))*((rho[i-1])+g3*pp.dr)*k;
    g = (g1 + 2*g2 + 2*g3 + g4)/6;
    l = (l1 + 2*l2 + 2*l3 + l4)/6;
    rho[i] = rho[i-1] + g*pp.dr;
    u[i] = u[i-1] + l*pp.dr;
    mi = mi - dmi;
    }
    if(x[i] >= R_max){
      break;
    }
    fprintf(fp,"%lf %lf %lf\n", x[i], rho[i], u[i]/G);
  }
}


void rungek4b(double* x, double* f, double* v, par pp, FILE* fp){
  double l1, l2, l3, l4, l, g1, g2, g3, g4, g, xi, dr, m=0,T=1.5e7,kkkk,e,L=0.;
  double mu,/*a=sqrt((1.4880e17)*(pp.n+1));a = sqrt(((2.23e22)/(pow(pp.rho0,2)))*(pp.n+1));*/a=sqrt((pp.n+1)*(3.1616e25)/(pp.rho0*pp.rho0)); //Direi che questo a è il migliore
  int i;
  f[0] = 1;
  v[0] = 0;
  x[0] = 0;
  dr = pp.dr/a;
  mu = (4*mH*(pow(T,1.5))*exp(-iii/(kb*T)))/((4*mH*(pow(T,1.5))*exp(-iii/(kb*T)))+(pp.rho0*c*c*c));
  printf("%lf\n",mu);
  for (i = 1; i <= pp.max_iter; i++){
    x[i] = x[i-1] + pp.dr;
    xi = (x[i])/a;
    g1 = ((v[i-1]))*dr;
    l1 = (-2/xi)*((v[i-1]))*dr-(pow((f[i-1]),pp.n))*dr;
    g2 = (v[i-1] + l1*dr*0.5)*dr;
    l2 = (-2/xi)*(v[i-1] + l1*dr*0.5)*dr-(pow((f[i-1] + g1*dr*0.5),pp.n))*dr;
    g3 = (v[i-1] + l2*dr*0.5)*dr;
    l3 = (-2/xi)*(v[i-1] + l2*dr*0.5)*dr-(pow((f[i-1] + g2*dr*0.5),pp.n))*dr;
    g4 = (v[i-1] + l3*dr)*dr;
    l4 = (-2/xi)*(v[i-1]  + l3*dr)*dr-(pow((f[i-1] + g3*dr),pp.n))*dr;
    g = (g1 + 2*g2 + 2*g3 + g4)/6;
    l = (l1 + 2*l2 + 2*l3 + l4)/6;
    m += -4*M_PI*a*a*a*pp.rho0*xi*xi*(v[i-1])*dr;
    f[i] = f[i-1] + g;
    v[i] = v[i-1] + l;
    if(x[i] > R_max){
      break;
    }
    if(f[i] < 0){
      printf("%12.11lf\n", fabs(sqrt(6)-xi)); //Serviva per stimare l'errore
      break;
      }
    mu = (4*mH*(pow(T*f[i],1.5))*exp(-iii/(kb*T*f[i])))/((4*mH*(pow(T*f[i],1.5))*exp(-iii/(kb*T*f[i])))+((pp.rho0*pow(f[i],pp.n))*c*c*c));
    if(xi>=6.6){
      mu=0.9995;
    }
    L+=4*M_PI*(pow(f[i],pp.n))*xi*xi*R_max*R_max*pp.dr*pow((T/(10e6))*(f[i]),16);
    fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf\n", xi, f[i], pow(f[i],pp.n), m,T*f[i],mu,L);
  }
}

void nube(par pp){
  double r=pp.r*9.461e15, Mj, t, mi=1, rho=mH*pp.d, m=(4/3)*M_PI*(pow(r,3))*rho;
  Mj = (pow((5.*kb*pp.T)/(G*mi*mH),3./2))*(pow(3./(4*M_PI*rho),1./2));
  if(m > Mj){
    t = pow((3.*M_PI)/(32.*G*rho),1./2);
    t /= 31536000;
    printf("La nube collassa\n La massa  della nube è M = %g kg\n La massa di Jeans è Mj = %g kg\n Il tempo di free fall è t = %g y\n", m, Mj, t);
  }	    
  if(m <= Mj){
    printf("La nube non collassa\n La massa della nube è M = %g mentre la massa di Jeans è Mj = %g kg\n", m, Mj);
  }
}
