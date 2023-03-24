
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <ctype.h>

#include <string.h>

#include "model.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------------*/

void cinet(double  rhot, double    dt, double   *c, double *qa, double *al,
           double *beta, double betat, double alam,     int ngr       )     {

/* ------ vars ------------*/
   double ci[ngr];
   double qai;
   double x;
   double deps;
   double asds;
/* ------------------------*/

   if(fabs(rhot) < 1.e-08) return;  

   double rbl=(rhot-betat)/alam;

   for(int i=1;i<ngr;++i)
      ci[i]=c[i];

   qai=*qa;

/*  calculo de wo(x) */

   double x1 = 0.0;
   double x3 = 0.0;
   int it = 0;
   int nt = 0;
   int nc = 0;

   deps=1.e-07;

   double ref=5.0*betat;

   if((rhot<=ref)||(rhot>-ref) ) {

      if(rhot<0.0) 
         x1 = -al[0]+deps;

      if(rhot>=0.0) 
         x3 = 9.*betat/alam;

      x=0.0;

      while(1){
         double sum  = 0.0;
         double sum2 = 0.0;

         for(int i=0;i<ngr;++i){
            double dsum = beta[i]/(x+al[i]);
            sum  += dsum;
            sum2 += dsum/(x+al[i]);
         }

         double fx  = x*(alam+sum)-rhot;
         double flx = alam+sum-x*sum2;
         double ffl = fx/flx;
         x   = x-ffl;

         nc=nc+1;

         if(fabs(ffl)<10.0*deps) break;

         if(x>x3) x=x3;
         if(x<x1) x=x1;
         if(x<x1) nt=nt+1;

         if(nt>5) break;
      }
   }
   else{ 
      x = -al[0] + deps;
      if(rhot>5.0*betat) x=rbl + deps;
   }

/* calculo de qa e c(i)'s */

   double wo   = x;

   double d1   = 0.01,
          arbl = fabs(rbl),
          dtc11;

   if(arbl<alam)  dtc11 = alam;
   if(arbl>=alam) dtc11 = 0.2/arbl;

   double dtc22 = 0.20/fabs(wo),
          dtc1  = dtc11,
          dtc2  = dtc22;

   if(dtc22<dtc11) dtc1 = dtc22;
   if(dtc22<dtc11) dtc2 = dtc11;

   double tt  = 0.0,
          dtc = 0.5*dt;

   if(dtc<dtc1) dtc=dtc1;

/*----- ---- ------*/
/*----- vars ------*/
/*----- ---- ------*/
   double erbl,ewoh,qfac,s,xj,yj;
/*----------------------------------------------------------------------------*/

point_1:

   if(dtc>dtc2) dtc=dtc2;

point_2:

   tt=tt+dtc;

   if(tt>dt+deps) {
      dtc = dt-(tt-dtc);
      tt  = dt;
   }

   erbl = exp(rbl*dtc);
   ewoh = exp(wo*dtc);
   qfac = (ewoh-erbl)/(wo-rbl);
   s=0.0;

   for(int i=0;i<ngr;++i){
      s=s+al[i]*ci[i];
      double ealh=exp(-al[i]*dtc);
      c[i]=ci[i]*ealh+qai*beta[i]*(ewoh-ealh)/alam/(wo+al[i]);
   }

   *qa = qai*erbl+qfac*s;

   it = it+1;
   xj = fabs(*qa/qai-ewoh);

/*  as vezes precursores nao seguem o nivel de potencia em  dtc
    portanto comparacao de  yj  foi cancelado                                 
        yj=fabs(c[6]/ci[6]-ewoh);
        if(yj>xj) xj = yj;                                                    */

   if(xj<1.0e-10) xj = 1.0e-10;

   if(xj<d1){
      double dif=fabs(1.0-tt/dt);

      if(dif<0.00101) return;

      if( d1/xj > 10.0) xj = 0.1*d1;
      dtc=dtc*d1/xj;

point_3:
      qai=*qa;

      for(int i=0; i<ngr;++i)
         ci[i]=c[i];

      goto point_1;

   }else if(dtc>dtc1) {

      tt=tt-dtc;
      dtc=0.5*dtc;

      goto point_2;
   }
   dtc=dtc1;

   goto point_3;
}

/*----------------------------------------------------------------------------*/


void thinit ( double *tcl0,  // temperatura pastilha (canal médio) (linha central) ------------------
              double *tfr,   // temperatura pastilha (nn)(canal médio) ------------------
              double *tc0,   // temperatura revestimento (canal médio) ------------------
              double *tm0,   // temperatura refrigereante (canal médio) ------------------
              double *ts0,   // temperatura de saida refrigereante (canal médio) ------------------
              double  te,    // temperatura de entrada ------------------
              double *ptcl0, // temperatura pastilha (canal quente) (linha central) ------------------
              double *ptfr,  // temperatura pastilha (nn)(canal quente)------------------
              double *ptc0,  // temperatura revestimento (canal quente)------------------
              double *ptm0,  // temperatura refrigereante (canal quente)------------------
              double *pts0,  // temperatura de saida refrigereante (canal quente)------------------
              double *dr,    // raio eff de transferencia de calor  ---------------------
              double *af,    // area de transferencia combus. (nn) -------------
              double  ac,    // area de transferencia revestimento 
              double  afs,   // area de transferencia pastilha
              double *tf0,   // temperatura média combustivel (canal medio)
              double *ptf0,  // temperatura média combustivel (canal quente) 
              data    dv  ){

   
/* calculo das temperaturas nos canais medio e quente */

   double qa = dv.qa;
      
   int flag  = 0; 

   double ts,
          tm,
          tc,
          sum;

   double xa[dv.nn+1],
          xb[dv.nn+1],
          xc[dv.nn+1],
          xp[dv.nn+1],
          xt[dv.nn+1];

   int nn1,
       nnn;

/* inicializacao das temperaturas em funcao de te */

   ini_point:
       
   ts = qa/(dv.gm*dv.cpm)+te; 
   tm=(ts+te)*0.5;

// -------temperatura inicial do revestimento ------------------------

   tc=0.97*qa/(ac*dv.hm)+tm;

/* -------temperatura inicial do combustivel-------------------------
          como a pastilha e discretizada tem que resolver para todos
          os pontos simultaneamente. utiliza-se o metodo de elimina-
          cao de gauss.               
// -------calculo dos coeficientes------------------------------------*/

   for(int n=0;n<(dv.nn-1);++n){
      xb[n+1] = -af[n+1]*dv.kf/dr[n+1];
      xc[n+1] = -af[n]*dv.kf/dr[n];
      xa[n+1] = -(xb[n+1]+xc[n+1]);
      xp[n+1] =  0.97*qa/dv.nn;
   }

/* ---------temperatura no proximo ponto dentro da pastilha ----------*/

   xp[0] =  0.5*0.97*qa/dv.nn;
   xb[0] =  xc[1];
   xa[0] = -xb[0];
       
   xa[dv.nn] = -xb[dv.nn-1]+afs*dv.hg;
   xc[dv.nn] =  xb[dv.nn-1];
   xp[dv.nn] =  xp[0]+afs*dv.hg*tc;

/* ---------temperatura dos pontos restantes na pastilha---------------
             eliminacao progressiva */

   for(int n=1;n<(dv.nn+1);++n){
      double dum  = xc[n]/xa[n-1];
      xa[n]= xa[n]-xb[n-1]*dum;
      xp[n]= xp[n]-xp[n-1]*dum;
   }

/*        eliminacao regressiva */

   xt[dv.nn]=xp[dv.nn]/xa[dv.nn];      

   for(int n=dv.nn-1;n>=0;--n)
      xt[n]=(xp[n]-xb[n]*xt[n+1])/xa[n];

/* --------xt(1) = temperatura na linha central--------------------
//---------xt(2) = temperatura no ponto 1--------------------------
//           .             .
//           .             .
//           .             .
//---------xt(nn1) = temperatura no ponto nn-2---------------------
//         xt(nn)  =       "          "   nn-1--------------------
//         xt(nnn) =       "          "   nn---------------------- */
//
/*----------calculo da temperatura media do combustivel----------- */

   sum=0.0;
   for(int n=1;n<dv.nn;++n)
      sum+=xt[n];

   sum+=(xt[0]+xt[dv.nn])*0.5;
   sum=sum/dv.nn;

   if(flag == 0) {
 
/* ..........saida para o canal medio ------------------------------*/

      for(int n=0;n<dv.nn;++n)
         tfr[n]=xt[n];

      *tf0  = sum;
      *tc0  = tc;
      *tm0  = tm;
      *ts0  = ts;
      *tcl0 = xt[0];

      flag = 1;
      qa   = dv.qa*dv.pfat;

      goto ini_point;
   }

// -----------saida para o canal quente -------------------------------  

   for(int n=0;n<dv.nn;++n)
      ptfr[n]=xt[n+1];

   *ptc0  = tc;
   *ptm0  = tm;
   *pts0  = ts;
   *ptcl0 = xt[0];
   *ptf0  = sum;
}

void banner(){

   printf("                                                                          \n");
   printf("                                                                          \n");
   printf("     ▄████▄   ██▓ ███▄    █ ▓█████▄▄▄█████▓ ██░ ██  ██▓ ▄████▄   ▄▄▄      \n");
   printf("    ▒██▀ ▀█  ▓██▒ ██ ▀█   █ ▓█   ▀▓  ██▒ ▓▒▓██░ ██▒▓██▒▒██▀ ▀█  ▒████▄    \n");
   printf("    ▒▓█    ▄ ▒██▒▓██  ▀█ ██▒▒███  ▒ ▓██░ ▒░▒██▀▀██░▒██▒▒▓█    ▄ ▒██  ▀█▄  \n");
   printf("    ▒▓▓▄ ▄██▒░██░▓██▒  ▐▌██▒▒▓█  ▄░ ▓██▓ ░ ░▓█ ░██ ░██░▒▓▓▄ ▄██▒░██▄▄▄▄██ \n");
   printf("    ▒ ▓███▀ ░░██░▒██░   ▓██░░▒████▒ ▒██▒ ░ ░▓█▒░██▓░██░▒ ▓███▀ ░ ▓█   ▓██▒\n");
   printf("    ░ ░▒ ▒  ░░▓  ░ ▒░   ▒ ▒ ░░ ▒░ ░ ▒ ░░    ▒ ░░▒░▒░▓  ░ ░▒ ▒  ░ ▒▒   ▓▒█░\n");
   printf("      ░  ▒    ▒ ░░ ░░   ░ ▒░ ░ ░  ░   ░     ▒ ░▒░ ░ ▒ ░  ░  ▒     ▒   ▒▒ ░\n");
   printf("    ░         ▒ ░   ░   ░ ░    ░    ░       ░  ░░ ░ ▒ ░░          ░   ▒   \n");
   printf("    ░ ░       ░           ░    ░  ░         ░  ░  ░ ░  ░ ░            ░  ░\n");
   printf("    ░                     ░                                               \n");

}



