#ifndef __MODEL_H
#define __MODEL_H

#include <math.h>
#include "helper.h"

/* for include in C++ */
#ifdef __cplusplus
extern "C" {
#endif

void cinet(double  rhot, double    dt, double   *c, double *qa, double *al,
           double *beta, double betat, double alam,     int ngr       ) ;


/*
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
              double  p,     // potencia termica  xxxxxxxxxxxxxxxxxxx qa
              double  pfat,  // fator de pico     xxxxxxxxxxxxxxxxxxx pfat
              double *dr,    // raio eff de transferencia de calor  ---------------------
              double  kf,    // condutubilidade combustivel xxxxxxxxxxxxxxxxxxx
              double  hg,    // coef. convecção revestimento (cladding) xxxxxxxxxxxxxxxxxxx
              double  hm,    // coef. convecção moderador xxxxxxxxxxxxxxxxxxx
              double *af,    // area de transferencia combus. (nn) -------------
              double  ac,    // area de transferencia revestimento 
              double  afs,   // area de transferencia pastilha
              double  gm,    // vasao de refri. xxxxxxxxxxxxxxxxxxxxxxx
              double  cpm,   // calor espc. moderador xxxxxxxxxxxxxxxxxxx
              int     nn,    // pontos discretização moderador xxxxxxxxxxxxxxxxxxx
              double *tf0,   // temperatura média combustivel (canal medio)
              double *ptf0); // temperatura média combustivel (canal quente) */

//void thinit ( &tcl0,  tfr,  &tc0, &tm0, &ts0, te, &ptcl0, ptfr, &ptc0, &ptm0, &pts0, &tf0, &ptf0, data *dv);



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
              data   dv  );


void banner();

#ifdef __cplusplus
}
#endif

#endif /* __MODEL_H */
