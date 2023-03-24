#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M_PI 3.14159265358979323846

#include "model.h"
#include "helper.h"
#include "inp.h"

/*------------------------------------------------------------------------------

Definições:

  subindices:  f = combustivel
               m = moderador
               c = revestimento

  exemplo:  mf= massa do combustivel
            mm= massa do moderador
            mc= massa de revestimento

    nota:
          variaveis exclusivas do canal quente comecam com o prefixo p

    tmax = tempo maximo do acidente (seg)

    dtprt = intervalo maximo entre prints(seg)
    fqprt = fracao maxima entre potencias para print.
    lftt  = comprimento total do combustivel no reator
    hg    = coeficiente de conveccao do gap
 
------------------------------------------------------------------------------*/


int main (int argc, char **argv) {

    banner();


    data dv; // data values

    if (inp_parse("test.inp", handler, &dv) < 0) {
        printf("Can't load 'test.inp'\n");
        return 1;
    }


   /*-------------------------*/
   /*variaveis usadas internas*/
   /*-------------------------*/
   double    te = 0, // Temperatura REFRIGERANTE ENTRADA (canal medio )
            pte = 0, // Temperatura REFRIGERANTE ENTRADA (canal quente)
          ptmp0 = 0;
   /*-------------------------*/

      // PZERO=QA                          <<<<<<<<<<<<<<<<<<<<<<                       
      // IF(PFAT.EQ.0.0) PFAT=1.0       vai entender isso 
      // IF(IRST.NE.0) READ(5,13) C     pensar aqui

   double beta_t = 0.0;
   for(int i=0;i<dv.ngr;++i)  //obtendo beta total              
      beta_t += dv.beta[i]; 

// RHOT=RHOI + A1*T + A2*T*T + A3*T*T*T + ALFAF*(TF(T)-TF0) + ALFAM*(TM(T)-TM0)


/*                                   
   DT    = INTERVALO DE TEMPO PARA EULER (SEG)                             
   TMAX  = TEMPO MAXIMO DO ACIDENTE (SEG)                                
   DTPRT = INTERVALO MAXIMO ENTRE PRINTS(SEG)                           
   FQPRT = FRACAO MAXIMA ENTRE POTENCIAS PARA PRINT.                    
   LFTT  = COMPRIMENTO TOTAL DO COMBUSTIVEL NO REATOR                    
   HG    = COEFICIENTE DE CONVECCAO DO GAP                                 

*/                                   

//      PTMP0=TMP0                    00011710
//      PTMP1=TMP1                    00011720
//      PTMP2=TMP2                    00011730
//      IF(TCVCC.LE.0.0) TCVCC=1000.  00011800

   ptmp0=dv.tmp0;
                

//      PHM=HM                        00012000


   /*-------------------------*/
   /*variaveis usadas internas*/
   /*-------------------------*/
   double    c [dv.ngr];      // concentação dos precursores

   double    r [dv.nn];      // RAIOS nos INTERVALO nn (VOLUME CONSTANTE)
   double    dr[dv.nn];      // 
   double    af[dv.nn];      // 

   double    afs  = 0.0;     // AREA DE TRANSFERENCIA DE CALOR DA PASTILHA
   double    ac   = 0.0;     // AREA DE TRANSFERENCIA DE CALOR DO REVESTIMENTO
   double    dmf  = 0.0;     // 
   double    dt   = 0.0;     // DT = INTERVALO DE TEMPO PARA EULER (SEG)
   
   /*-------------------------*/

   if(dv.itemp!=0){    // calcula temperatura de entrada do refrigerante considerando variavel
      double tmp = dv.tmp0+dv.tmp1*dv.qa; // TEMP. DE ENTRADA DO MODERADOR 
      te  = tmp-dv.qa/(2.0*dv.gm*dv.cpm);
      pte = te;
   }
   else{
      pte = ptmp0;
      te  = dv.tmp0;  // mantem temperatura de entrada do refrigerante const
   }


   if(dv.iflux == 0)
      printf("*** PROBLEMA SEM FLUXO DE REFRIGERANTE (iflux = %i) *** \n",dv.iflux  );
   else 
      printf("*** PROBLEMA COM FLUXO DE REFRIGERANTE CONSTANTE (iflux = %i) **** \n",dv.iflux  );
   
   


   printf("MASSA DO COMBUSTIVEL.................: %8.4e  [Kg]      \n",dv.mf  );
   printf("MASSA DO REVESTIMENTO................: %8.4e  [Kg]      \n",dv.mc  );
   printf("MASSA DO MODERADOR...................: %8.4e  [Kg]      \n",dv.mm  );
   printf("CALOR ESPECIFICO DO COMBUSTIVEL..... : %8.4e  [J/Kg/oC] \n",dv.cpf );
   printf("CALOR ESPECIFICO DO REVESTIMENTO.....: %8.4e  [J/Kg/oC] \n",dv.cpc );
   printf("CALOR ESPECIFICO DO MODERADOR........: %8.4e  [J/Kg/oC] \n",dv.cpm );
   printf("VAZAO  DO  REFRIGERANTE..............: %8.4e  [Kg/s]    \n",dv.gm  );
   printf("TEMPERATURA DE ENTRADA DO MODERADOR..: %8.4e  [oC]      \n",dv.tmp0);
   printf("TEMP. DE ENT. DO MOD. NO CANAL QUENTE: %8.4e  [oC]      \n",ptmp0  );
   printf("RAIO  DO COMBUSTIVEL.................: %8.4e  [mm]      \n",dv.rf  );
   printf("RAIO DO REVESTIMENTO.................: %8.4e  [mm]      \n",dv.rc  );
   printf("ESPESSURA DO REVESTIMENTO............: %8.4e  [mm]      \n",dv.dc  );
   printf("COMPRIMENTO TOTAL DO COMBUSTIVEL.....: %8.4e  [m]       \n",dv.lftt);
   printf("CONDUTIVIDADE DO COMBUSTIVEL.........: %8.4e  [W/m/oC]  \n",dv.kf  );
   printf("CONDUTIVIDADE DO REVESTIMENTO........: %8.4e  [W/m/oC]  \n",dv.kc  );
   printf("COEFICIENTE DE CONVECCAO DO GAP......: %8.4e  [W/m²/oC] \n",dv.hg  );
   printf("COEFICIENTE DE CONVECCAO DO MODERADOR: %8.4e  [W/m²/oC] \n",dv.hm  );
   printf("COEFICIENTE DE CONVECCAO EM EBULICAO.: %8.4e  [W/m²/oC] \n",dv.hdnb);
   printf("FATOR DE PICO........................: %8.4e  [1]       \n",dv.pfat);

   printf("\nCONVECCAO COM EBULICAO A PARTIR DE %8.4e [oC] NO REVESTIMENTO. \n\n",dv.tcvcc);
 
//  INICIALIZAR  C(I), I=1,NGR PARA PROBLEMAS SEM RESTART                

// retirado ... somente sem restart ...   if (dv.irst == 0){
   double fact = dv.qa/dv.alam;
   for(int j=0;j<dv.ngr;++j)
      c[j]=fact*dv.beta[j]/dv.al[j];

   printf("CONCENTRACAO INICIAL DOS PRECURSORES(WATT). \n\n");

   printf("GRUPO     CONCENTRACAO \n" );
   for(int j=0;j<dv.ngr;++j){
      printf("  %i......%8.4e \n", j+1, c[j]);
   }

   printf("\nAS SAIDAS SERAO EXECUTADAS NO MAXIMO A CADA %.4e [s] OU ANTES DE %.4e [s] VEZES O AUMENTO DE POTENCIA. \n", dv.dtprt, dv.fqprt);


//
//   CALCULO DE AREAS DE TRANSFERENCIA DE CALOR DA PASTILHA, AFS ,       
//   E  DO REVESTIMENTO, AC , ----------------------------               
//          
   afs = 2.0*M_PI*dv.rf*dv.lftt;
   ac  = 2.0*M_PI*dv.rc*dv.lftt; 

//   CALCULO DE RAIOS DO INTERVALO TAL QUE VOLUME CONSTANTE EM NN INTERVALOS    

          
   { 
      double dvol  = M_PI*dv.rf*dv.rf/dv.nn;
      double r[dv.nn];

      dmf = dv.mf/dv.nn;
            
      //   PI*R(I)*R(I)= I*DVOL -----------

      for (int i=0;i<dv.nn;++i)
          r[i]= sqrt((i+1)*dvol/M_PI);

      dr[0] = r[0];
      for (int i=1;i<dv.nn;++i)
         dr[i] = r[i]-r[i-1];

      af[0]= 2.0*M_PI*r[0]*dv.lftt/1.4142;
      for (int i=1;i<dv.nn;++i){
         double r12 = sqrt(r[i]*r[i]+r[i-1]*r[i-1])/1.4142;
         af[i]= 2.0*M_PI*r12*dv.lftt;
      }



      printf("AS TEMPERATURAS SAO CALCULADAS NOS PONTOS: \n\n");

      printf("  #    LOCALIZACAO RADIAL (mm) \n" );
      for(int j=0;j<dv.nn;++j){
         printf("  %i.........%8.4e \n", j+1, r[j]);
      }

//                             
//----- CRITERIO PARA  MAXIMO  DT  -----------------------------   
//                         

      double dtmax1 = dmf*dv.cpf*0.5/((dv.kf/dr[dv.nn-1]+dv.hg)*afs);
      double dtmax2 = dv.mc*dv.cpc/(afs*dv.hg+ac*dv.hm);
      double dtmax=dtmax1;
  
      if(dtmax1>dtmax2)
         dtmax=dtmax2; 

      dt = dv.dt0;

      if( dt > 0.10*dtmax)
         dt = 0.10*dtmax; 


      printf("\nINTERVALO DE TEMPO SOLICITADO PELO USUARIO PARA OS CALCULOS TERMOHIDRAULICOS: %8.4e [s]\n", dv.dt0);
      printf("INTERVALO DE TEMPO MAXIMO PERMITIDO PELO PROGRAMA PARA SE OBTER CONVERGENCIA: %8.4e [s]\n", dtmax );
      printf("INTERVALO DE TEMPO EFETIVAMENTE UTILIZADO PELO PROGRAMA PARA OS CALCULOS    : %8.4e [s]\n", dt    );

   }

 double        tcl0, tfr[dv.nn],  
                tc0,        tm0,  
                ts0,      ptcl0, 
        ptfr[dv.nn],       ptc0,  
               ptm0,       pts0,  
                tf0,       ptf0;  

// ***  thinit ( &tcl0,  tfr,  &tc0, &tm0, &ts0, te, &ptcl0, ptfr, &ptc0, &ptm0, &pts0, 
// ***              qa, pfat,    dr,   kf,   hg, hm,     af,   ac,   afs,    gm,   cpm,
// ***              nn, &tf0, &ptf0              );


   thinit ( &tcl0,  tfr,  &tc0, &tm0, &ts0, te, &ptcl0, ptfr, &ptc0, &ptm0, &pts0, dr, af, ac, afs, &tf0, &ptf0, dv );


   if(dv.iflux==0) {

// ------------INICIALIZAR TFR(I),I=1,NN,  SE TFR0.NE.0.0-----------

      if(tf0>0.0){ 
          for(int i=0;i<(dv.nn);++i){
             tfr[i]  = tf0;             
             ptfr[i] = ptf0;           
          }   
      }

      double sum=0.0;                
      double psum=0.0;               

      for(int i=0;i<(dv.nn);++i){
         sum  = sum+tfr[i];          
         psum = psum+ptfr[i];       
      }
      tf0=sum/dv.nn;             
      ptf0=psum/dv.nn;           
                           
   }

//   INICIALIZANDO AS TEMPERATURAS

   double tcl  = tcl0;
   double tf   = tf0;                 
   double tc   = tc0;                 
   double tm   = tm0;                 
   double ts   = ts0;                 
   double ptcl = ptcl0;             
   double ptf  = ptf0;               
   double ptc  = ptc0;              
   double ptm  = ptm0;              
   double pts  = pts0;              

   double  t; 
      int  it; 
      int  itmax;
   double  ipaf;
   double  tprt;
   double  qprt;
   double  rof; 
   double  rom; 
   double  rhot;
// *** 
   double  qaint; 
   double  dqa; 
   double  pqa; 
   double  pdqa; 
   double  qex; 
   double  pqex; 
   double  qin; 
   double  fr0; 
   double  pqin; 
   double  ptfr0[dv.nn]; 
   double  tfr0[dv.nn]; 
   double  phm=dv.hm; 
// *** 
   double  poterm;  
   double  tmp;  
   double  axm;  
   double  ppterm;  
   double  ptmp;  
// *** 
   double  tfmax; //*
   double  ptfmax;//*
// *** 
   double  roex;
   double  rex;
   double  rhos;
   double  rofs;
   double  roms;
   double  iprt;
   double  pmw;
   double  ptmw;
   double  veloc;
   double  dz;
   double  fccm;
   double  tsc;
   double  tec;
   double  ipag;


   printf("\nACIDENTE OU TRANSIENTE COM DURACAO DE %8.4e [s], \n", dv.tmax);

   printf("\nTEMPERATURAS INICIAIS \n");
   printf("\nMEDIA  PICO\n");

   printf("COMBUSTIVEL ..........  %8.4e  %8.4e   [oC], \n", tf,ptf);
   printf("REVESTIMENTO .........  %8.4e  %8.4e   [oC], \n", tc,ptc);
   printf("REFRIGERANTE MEDIA ...  %8.4e  %8.4e   [oC], \n", tm,ptm);
   printf("REFRIGERANTE SAIDA ...  %8.4e  %8.4e   [oC], \n", ts,pts);

   printf("\nREFRIGERANTE ENTRADA0  %8.4e  [oC], \n", te);
//   printf("REFRIGERANTE ENTRADA1 . %8.4e  %8.4e   [oC], \n", tmp1);
//   printf("REFRIGERANTE ENTRADA2 . %8.4e  %8.4e   [oC], \n", tmp1);
                            
/*      WRITE(6,28) TMAX,TF,PTF,TC,PTC,TM,PTM,TS,PTS,TE,TMP1,TMP2    
   28 FORMAT(//,' ACIDENTE OU TRANSIENTE COM DURACAO DE',F10.4,    
     & ' SEGUNDOS'      ///    
     & ' T E M P E R A T U R A S   I N I C I A I S'/               
     &                   T25,'  MEDIAS',T40,' DO PICO'/            
     & ' COMBUSTIVEL           =', T25,F8.2,' C', T39,F8.2,' C'/   
     & ' REVESTIMENTO          =', T25,F8.2,' C', T39,F8.2,' C'/   
     & ' REFRIGERANTE MEDIA    =', T25,F8.2,' C', T39,F8.2,' C'/   
     & ' REFRIGERANTE SAIDA    =', T25,F8.2,' C', T39,F8.2,' C'/   
     & ' REFRIGERANTE ENTRADA0 =',T25, F8.2,' C'/                  
     & ' REFRIGERANTE ENTRADA1 =',T25, 1PE12.4,' C'/               
     & ' REFRIGERANTE ENTRADA2 =',T25, 1PE12.4,' C'/)              
C                              
      WRITE(6,914)             
  914 FORMAT(/,'**** TEMPERATURA NA PASTILHA ****',/)              
C                              
      WRITE(6,924) TCL,(  TFR(I),I=1,NN)                           
  924 FORMAT(1X,'PARA CANAL MEDIO = ',6(F12.2))                    
C                              
      WRITE(6,934) PTCL, (PTFR(I), I=1,NN )                        
  934 FORMAT(1X,'PARA CANAL QUENTE = ',6(F12.2))                   
*/


   printf(" TEMPERATURA NA PASTILHA \n" );
   printf(" PARA CANAL MEDIO \n" );
   for(int j=0;j<dv.nn;++j){
      printf("  %i.........%8.4e \n", j+1, tfr[j]);
   }
   printf(" PARA CANAL QUENTE \n" );
   for(int j=0;j<dv.nn;++j){
      printf("  %i.........%8.4e \n", j+1, ptfr[j]);
   }

                        
// INICIO  DO PROCESSO DE CONTAGEM DE TEMPO DO ACIDENTE     
                        
    t     = 0.0;              
    it    = 0;               
    itmax = 12;          

    if(dv.nn>5) itmax=8;                                   

    dv.dtprt=0.95*dv.dtprt;  
    dv.fqprt=0.98*dv.fqprt;  

    ipaf=1;            

label_1: 

    tprt= t;           
    qprt=dv.qa;           

label_2:
    t = t + dt; 
      
    if(t>=dv.tmax) {
         goto label_3;
    }                                
                       
//  CALCULO DA REATIVIDADE E NOVA POTENCIA EM T+DT         
//                     
// GL  RHOT= RHOI+A1*T+A2*T*T+ALFAF*(TF-TF0)+ALFAM*(TM-TM0)
//                     
    rof  = dv.xif*dv.alfaf*(tf-tf0) + (1.-dv.xif)*dv.alfaf*(ptf-ptf0);   
    rom  = dv.xim*dv.alfam*(tm-tm0) + (1.-dv.xim)*dv.alfam*(ptm-ptm0);   
    rhot = dv.rhoi+dv.a1*t+dv.a2*t*t + rof + rom;                    

    cinet( rhot,  dt, c, &(dv.qa), dv.al, dv.beta, beta_t, dv.alam, dv.ngr);

//  CALCULO DE POTENCIA INTEGRADA                                

    qaint=qaint+ dv.qa*dt;    
    dqa=dv.qa/dv.nn;             
    pqa= dv.pfat*dv.qa;          
    pdqa= dv.pfat*dqa;        
                             
// CALCULO DE  NOVAS  T E M P E R A T U R A S                    
//                           
//   C O M B U S T I V E L   
                             
    qex= af[0]*dv.kf*(tcl-tfr[0] )/dr[0];                         
    tcl= tcl + dt*(0.97*dqa - 2.0*qex)/( dmf*dv.cpf);            

    pqex= af[0]*dv.kf*( ptcl - ptfr[0] )/dr[0];                   
    ptcl= ptcl + dt*(0.97*pdqa - 2.0*pqex)/(dmf*dv.cpf);          
                             
                             
    for(int i=0;i<(dv.nn-1);++i){         
       qin= qex;              
       qex= af[i+1]*dv.kf*(tfr[i]-tfr[i+1] )/dr[i+1];                
       tfr0[i]=tfr[i] + dt*(0.97*dqa + qin - qex)/( dmf*dv.cpf );    
                             
       pqin=pqex;             
       pqex= af[i+1]*dv.kf*(ptfr[i] - ptfr[i+1] )/dr[i+1];           
       ptfr0[i]= ptfr[i] + dt*(0.97*pdqa+pqin-pqex)/(dmf*dv.cpf);    
    } 
                             

    qin = qex;                     
    qex = afs*dv.hg*(tfr[dv.nn-1]-tc);                                   
    tfr0[dv.nn-1]= tfr[dv.nn-1] + dt*(0.97*dqa + 2.0*qin - 2.0*qex)/(dmf*dv.cpf );

    pqin= pqex;                   
    pqex= afs*dv.hg*( ptfr[dv.nn-1] - ptc);                                   
    ptfr0[dv.nn-1]= ptfr[dv.nn-1] + dt*(0.97*pdqa + 2.0*pqin - 2.0*pqex)/(dmf*dv.cpf);                   

    for(int i=0; i<dv.nn; ++i){             
         tfr[i] = tfr0[i];              
        ptfr[i] = ptfr0[i];            
    }                             

//    revestimento       


    double hm_=dv.hm;
    
    if(tc>dv.tcvcc) hm_=dv.hdnb;      
    if(ptc>dv.tcvcc) phm= dv.hdnb;   
    qin= qex;                     
    qex= ac*(hm_)*( tc-tm );         
    tc= tc + dt*( qin-qex )/(dv.mc*dv.cpc);                                 

    pqin= pqex;                   
    pqex= ac*phm*(ptc-ptm);       
    ptc= ptc+dt*(pqin-pqex)/(dv.mc*dv.cpc);                                 

//   MODERADOR        
                            
    if(dv.iflux==0){ 
      tm = tm+dt*(0.03*dv.qa+qex)/(dv.mm*dv.cpm);                           
      ts = tm;                 
      te = tm;                 
    }
    else{                                
      poterm= dv.gm*dv.cpm*(ts-te); 
      tmp   = dv.tmp0+dv.tmp1*poterm;  
      te    = tmp-poterm/(2.0*dv.gm*dv.cpm);                                

      if(dv.itemp==0) te=dv.tmp0;

      axm = (dv.mm*dv.cpm)/dt;       
      tm  = tm+(0.03*dv.qa+qex-dv.gm*dv.cpm*(ts-te))/axm;                    
      ts  = 2.0*tm-te;          
    }
                   
    if(dv.iflux==0){        
       ptm = ptm+dt*(0.03*dv.qa+qex)/(dv.mm*dv.cpm);                         
       pts = ptm;               
       pte = ptm;               
    }
    else{                                
       ppterm = dv.gm*dv.cpm*(pts-pte);                                   
       ptmp   = 	ptmp0+dv.tmp1*ppterm;
       pte    = ptmp-ppterm/(2.0*dv.gm*dv.cpm);                              

       if(dv.itemp==0) pte=ptmp0;                                  

       ptm = ptm+(0.03*pqa+pqex-dv.gm*dv.cpm*(pts-pte))/axm;              
       pts = 2.0*ptm-pte;     
    }  
                             
                          
//  CALCULO DE MAXIMA  TEMPERATURA E  TF(MEDIA NO COMBUSTIVEL)
                          
    tfmax  = tcl;         
    tf     = 0.5*tcl;         
    ptfmax = ptcl;       
    ptf    = 0.5*ptcl;      
                       
    for(int i=0;i<dv.nn;++i){                              
       tf = tf + tfr[i];     
       if(tfmax<tfr[i]) 
          tfmax=tfr[i];                       
       ptf = ptf + ptfr[i];   
       if(ptfmax<ptfr[i] ) 
         ptfmax= ptfr[i];                 
    }
// *** 
    tf  = tf- 0.5*tfr[dv.nn-1]; 
    tf  = tf/dv.nn;           
    ptf = ptf-0.5*ptfr[dv.nn-1];                                  
    ptf = ptf/dv.nn;        
             
//  check para print                             

   if((fabs(t-tprt)<dv.dtprt)||(fabs(dv.qa/qprt)<dv.fqprt)||(fabs(qprt/dv.qa)<dv.fqprt))
       goto label_2;                                   

   it=it+1;

//    saida dos   resultados

// reatividade (absoluta) 
     roex = dv.rhoi + dv.a1*t + dv.a2*t*t + dv.a3*t*t*t;  
     rex  = roex*1.e+05;                       
     rhos = rhot*1.e+05;                      
     rofs = rof*1.e+05;                       
     roms = rom*1.e+05;                        

     if(it>=itmax) it=0;                      




printf("\nTEMPO DE TRANSIENTE de %8.4e [s], \n", t);

   printf("POTENCIA ................  %8.4e  [W]  , \n", dv.qa);
   printf("POTENCIA  TERMICA .......  %8.4e  [W]  , \n", poterm);
   printf("ENERGIA LIBERADA ........  %8.4e  [J]  , \n", qaint);
   printf("REATIVIDADE TOTAL .......  %8.4e  [pcm], \n", rhos);
   printf("REATIVIDADE EXTERNA .....  %8.4e  [pcm], \n", rex);
   printf("REATIVIDADE DOPPLER .....  %8.4e  [pcm], \n", rofs);
   printf("REATIVIDADE MODERADOR ...  %8.4e  [pcm], \n", roms);


/*
C                                   
      IF(IPAG.EQ.2) WRITE(6,9111)   
 9111 FORMAT(//)                    
      IF(IPAG.GT.2) WRITE(6,911)    
  911 FORMAT('1')                   
      IF(IPAG.GT.2) IPAG=1          
      WRITE(6,902) T,QA,POTERM,QAINT,RHOS,REX,ROFS,ROMS                 
C                                   
  902 FORMAT(//,' TEMPO DE TRANSIENTE    =  ',T34,1PE10.3,5X,'SEGUNDOS'/
     & ' POTENCIA               =', T32,1PE14.5,2X,'WATTS'/             
     & ' POTENCIA  TERMICA      =', T32,1PE14.5,2X,'WATTS'/             
     & ' ENERGIA LIBERADA       =', T32,1PE14.5,2X,'JOULES'/            
     & ' REATIVIDADE TOTAL      =', T32,1PE14.5,2X,'PCM '/              
     & ' REATIVIDADE EXTERNA    =', T32,1PE14.5,2X,'PCM '/              
     & ' REATIVIDADE DOPPLER    =', T32,1PE14.5,2X,'PCM '/              
     & ' REATIVIDADE MODERADOR  =', T32,1PE14.5,2X,'PCM '/)             
C                                   
      WRITE(6,912) TF,PTF,TC,PTC,TM,PTM,TS,PTS,TE,PTE                   
 912  FORMAT(/,3X,T40,' CANAL MEDIO ',T70,                              
     & ' CANAL  QUENTE   ',/        
     & ' TEMPERATURA DO COMBUSTIVEL           =  ',F8.2,' C',T70,F8.2,  
     & ' C',/                       
     & ' TEMPERATURA DO REVESTIMENTO          =  ',F8.2,' C',T70,F8.2,  
     & ' C',/                       
     & ' TEMPERATURA MEDIA DO REFRIGERANTE    =  ',F8.2,' C',T70,F8.2,  
     & ' C',/                       
     & ' TEMPERATURA SAIDA DO REFRIGERANTE    =  ',F8.2,' C',T70,F8.2,  
     & ' C',/                       
     & ' TEMPERATURA ENTRADA DO REFRIGERANTE  =  ',F8.2,' C',T70,F8.2,  
     & ' C')                        
C                                   
      WRITE(6,9431)                 
 9431 FORMAT(/,'*** TEMPERATURA NA PASTILHA ****',/)                    
      WRITE(6,943)                  
  943 FORMAT(' ',T10,   ' C A N A L    M E D I O')                      
      WRITE(6,904)     TCL,(  TFR(I),I=1,NN)                            
  904 FORMAT(' ',T10,   ' TF(I)=',(T21,1X,6(F12.2)))                    
C                                   
      WRITE(6,944)                  
  944 FORMAT(' ',T10,   ' C A N A L    Q U E N T E')                    
      WRITE(6,904) PTCL, (PTFR(I), I=1,NN )                             

*/


printf("\nTEMPO DE TRANSIENTE de %8.4e [s], \n", t);

   printf("\nTEMPERATURAS NA PASTILHA \n");
   printf("\nMEDIA  PICO\n");

   printf("COMBUSTIVEL ..........  %8.4e  %8.4e   [oC], \n", tf,ptf);
   printf("REVESTIMENTO .........  %8.4e  %8.4e   [oC], \n", tc,ptc);
   printf("REFRIGERANTE MEDIA ...  %8.4e  %8.4e   [oC], \n", tm,ptm);
   printf("REFRIGERANTE SAIDA ...  %8.4e  %8.4e   [oC], \n", ts,pts);

   printf("\nREFRIGERANTE ENTRADA0  %8.4e  [oC], \n", te);

   printf(" TEMPERATURA NA PASTILHA \n" );
   printf(" PARA CANAL MEDIO \n" );
   for(int j=0;j<dv.nn;++j){
      printf("  %i.........%8.4e \n", j+1, tfr[j]);
   }
   printf(" PARA CANAL QUENTE \n" );
   for(int j=0;j<dv.nn;++j){
      printf("  %i.........%8.4e \n", j+1, ptfr[j]);
   }

   if(iprt!=0){       
       pmw  = dv.qa/1.e+6;                 
       ptmw = poterm/1.e+6;            
       veloc= dz/dt;                  
       fccm = pqex/ac/1.e+6;           
       tsc  = ts-tm+dv.tmp0;               
       tec  = te-tm+dv.tmp0;               
    }
/*
      WRITE(2,992) T,PMW,PTMW,FCCM,RHOS,ROMS,ROFS,REX                   
      WRITE(3,992) T,PTF,PTC,TSC,TEC,TM,PRD,Z,VELOC                     
  992 FORMAT(10E12.3)               
*/
    ipag=ipag+1;                   
    goto label_1;                      

label_3: 
/*                                   
 1990 WRITE(6,1992) (I,C(I),I=1,NGR)
 1992 FORMAT(///,2X,40('=')/'0CONCENTRACAO FINAL DOS',                  
     & ' PRECURSORES(WATT)'/'0',40('=')////' GRUPO  CONCENTRACAO'//     
     & (I4,2X,1PE13.5))             
*/                                  

    return 0;
}


