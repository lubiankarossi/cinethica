/* Example: parse a simple configuration file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inp.h"

#define MAX_P 12

/*
typedef struct
{
    int version;
    const char* name;
    const char* email;
} data;
*/


int split (const char *str, char c, char ***arr)
{
    int count = 1;
    int token_len = 1;
    int i = 0;
    const char *p;
    char *t;

    p = str;
    while (*p != '\0')
    {
        if (*p == c)
            count++;
        p++;
    }

    *arr = (char**) malloc(sizeof(char*) * count);
    if (*arr == NULL)
        exit(1);

    p = str;
    while (*p != '\0')
    {
        if (*p == c)
        {
            (*arr)[i] = (char*) malloc( sizeof(char) * token_len );
            if ((*arr)[i] == NULL)
                exit(1);

            token_len = 0;
            i++;
        }
        p++;
        token_len++;
    }
    (*arr)[i] = (char*) malloc( sizeof(char) * token_len );
    if ((*arr)[i] == NULL)
        exit(1);

    i = 0;
    p = str;
    t = ((*arr)[i]);
    while (*p != '\0')
    {
        if (*p != c && *p != '\0')
        {
            *t = *p;
            t++;
        }
        else
        {
            *t = '\0';
            i++;
            t = ((*arr)[i]);
        }
        p++;
    }

    return count;
}


typedef struct {

 double al  [MAX_P], /* Lambda dos precursores */
        beta[MAX_P]; /* Betas dos precursor    */

 double alam  , /* Tempo de geração de nêutrons. (s)                                */
        rhoi  , /* Reatividade inicial (absoluto)                                   */
        a1    , /* Coeficiente linear da reatividade (s-1)                          */
        a2    , /* Coeficiente quadrático da reatividade (s-2)                      */
        a3    , /* Coeficiente cúbico da reatividade (s-3)                          */
        alfaf , /* Coeficiente de reatividade de combustível (abs/oC)               */
        alfam , /* Coeficiente de reatividade do moderador (abs/oC)                 */
        xim   , /* Fator de ponderação para o cálculo da reatividade do moderador   */
        xif   , /* Fator de ponderação para o cálculo da reatividade do combustível */
        pfat  , /* Fator de pico.                                                   */
        qa    , /* Potência térmcia. (MW)                                           */
        aqint , /* Energia acumulada (MJ) – Utilizar zero.                          */
        cpf   , /* Calor específico do combustível (J/kg°C)                         */
        cpc   , /* Calor específico do revestimento (J/kg°C)                        */
        cpm   , /* Calor específico do moderador (J/kg°C)                           */
        tcvcc , /* Temperatura com convecção com ebulição  (oC)                     */
        gm    , /* Vazão do refrigerante (kg/s)                                     */
        tmp0  , /* Temperatura de entrada do moderador no instante inicial (°C)     */
        tmp1  , /* Temperatura da entrada 1 do refrigerante (oC)                    */
        tmp2  , /* Temperatura da entrada 2 do refrigerante (oC)                    */
        kf    , /* Condutividade do combustível (W/m oC)                            */
        kc    , /* Condutividade do revestimento (W/m oC)                           */
        hg    , /* Coeficiente de convecção do GAP (W/m2 oC)                        */
        hm    , /* Coeficiente de convecção do moderador (W/m2 oC)                  */
        hdnb  , /* Coeficiente de convecção em ebulição (W/m2 oC)                   */
        mf    , /* Massa do combustível (kg).                                       */
        mc    , /* Massa do revestimento (kg).                                      */
        mm    , /* Massa do moderador (kg).                                         */
        rf    , /* Raio do combustível (mm).                                        */
        rc    , /* Raio do revestimento (mm).                                       */
        dc    , /* Espessura do revestimento (mm).                                  */
        lftt  , /* Comprimento total do combustível (m).                            */
        tmax  , /* Tempo máximo do variante (s)                                     */
        dt0   , /* Intervalo de tempo do cálculo (s)                                */
        dtprt , /* Intervalo de tempo para impressão (s)                            */
        fqprt ; /* Intervalo de tempo para a impressão dos dados do canal quente (s)*/
  
 int nn     , /* Número de pontos para o cálculo da distribuição */
     flux   , /* Flag do fluxo de refrigerante.                  */
     itemp  , /* Flag de temperatura de entrada do refrigerante. */
     ngr    ; /* Número de grupos de nêutrons atrasados          */

} data;

static int handler(void* _pdata, const char* section, const char* name,
                   const char* value)
{
    data* pdata = (data*) _pdata;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("neutronic"         , "ngr"       )) {
//       printf("%s - %s - %s \n",section,name,value);
       pdata->ngr = atoi(value);
//       printf("%s - %s - %i \n\n",section,name,pdata->ngr);
    }
    else if (MATCH("neutronic"         , "beta" )) {
//       printf("%s - %s - %s \n",section,name,value);

       char **arr = NULL;

       int c = split(value, ',', &arr);

       for (int i = 0; i < c; i++){
            pdata->beta[i]=atof(arr[i]);
//            printf("string #%d: %lf \n", i, pdata->beta[i]);
       }
       free(arr);

//       printf("\n");

    }
    else if (MATCH("neutronic"         , "alam" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->alam = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->alam);
    }
    else if (MATCH("neutronic"         , "rhoi" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->rhoi = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->rhoi);
    }
    else if (MATCH("neutronic"         , "a1"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->a1 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->a1);
    }
    else if (MATCH("neutronic"         , "a2"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->a2 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->a2);
    }
    else if (MATCH("neutronic"         , "a3"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->a3 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->a3);
    }
    else if (MATCH("neutronic"         , "alfaf")) { 
//       printf("%s - %s - %s \n",section,name,value);
       pdata->alfaf = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->alfaf);
    }
    else if (MATCH("neutronic"         , "alfam")) { 
//       printf("%s - %s - %s \n",section,name,value);
       pdata->alfam = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->alfam);
    }
    else if (MATCH("neutronic"         , "xim"  )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->xim = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->xim);
    }
    else if (MATCH("neutronic"         , "xif"  )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->xif = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->xif);
    }
    else if (MATCH("termo"             , "pfat" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->pfat = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->pfat);
    }
    else if (MATCH("termo"             , "qa"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->qa = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->qa);
    }
    else if (MATCH("termo"             , "aqint")) { 
//       printf("%s - %s - %s \n",section,name,value);
       pdata->aqint = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->aqint);
    }
    else if (MATCH("termo"             , "cpf"  )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->cpf = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->cpf);
    }
    else if (MATCH("termo"             , "cpc"  )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->cpc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->cpc);
    }
    else if (MATCH("termo"             , "cpm"  )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->cpc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->cpc);
    }
    else if (MATCH("termo"             , "tcvcc")) { 
//       printf("%s - %s - %s \n",section,name,value);
       pdata->tcvcc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->tcvcc);
    }
    else if (MATCH("termo"             , "gm"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->gm = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->gm);
    }
    else if (MATCH("termo"             , "tmp0" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->tmp0 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->tmp0);
    }
    else if (MATCH("termo"             , "tmp1" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->tmp1 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->tmp1);
    }
    else if (MATCH("termo"             , "tmp2" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->tmp2 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->tmp2);
    }
    else if (MATCH("termo"             , "kf"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->kf = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->kf);
    }
    else if (MATCH("termo"             , "kc"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->kc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->kc);
    }
    else if (MATCH("termo"             , "hg"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->hg = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->hg);
    }
    else if (MATCH("termo"             , "hm"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->hm = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->hm);
    }
    else if (MATCH("termo"             , "hdnb" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->hdnb = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->hdnb);
    }
    else if (MATCH("geometric-material", "mf"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->mf = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->mf);
    }
    else if (MATCH("geometric-material", "mc"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->mc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->mc);
    }
    else if (MATCH("geometric-material", "mm"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->mm = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->mm);
    }
    else if (MATCH("geometric-material", "rf"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->rf = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->rf);
    }
    else if (MATCH("geometric-material", "rc"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->rc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->rc);
    }
    else if (MATCH("geometric-material", "dc"   )) {    
//       printf("%s - %s - %s \n",section,name,value);
       pdata->dc = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->dc);
    }
    else if (MATCH("geometric-material", "lftt" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->lftt = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->lftt);
    }
    else if (MATCH("parameter"         , "tmax" )) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->tmax = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->tmax);
    }
    else if (MATCH("parameter"         , "dt0"  )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->dt0 = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->dt0);
    }
    else if (MATCH("parameter"         , "dtprt")) { 
//       printf("%s - %s - %s \n",section,name,value);
       pdata->dtprt = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->dtprt);
    }
    else if (MATCH("parameter"         , "fqprt")) { 
//       printf("%s - %s - %s \n",section,name,value);
       pdata->fqprt = atof(value);
//       printf("%s - %s - %lf \n\n",section,name,pdata->fqprt);
    }
    else if (MATCH("geometric-material", "nn"   )) {     
//       printf("%s - %s - %s \n",section,name,value);
       pdata->nn = atoi(value);
//       printf("%s - %s - %i \n\n",section,name,pdata->nn);
    }
    else if (MATCH("parameter"         , "flux" )) {   
//       printf("%s - %s - %s \n",section,name,value);
       pdata->flux = atoi(value);
//       printf("%s - %s - %i \n\n",section,name,pdata->flux);
    }
    else if (MATCH("parameter"         , "itemp")) {  
//       printf("%s - %s - %s \n",section,name,value);
       pdata->itemp = atoi(value);
//       printf("%s - %s - %i \n\n",section,name,pdata->itemp);
    }
    else if (MATCH("neutronic"         , "al"   )) {
//       printf("%s - %s - %s \n",section,name,value);

       char **arr = NULL;

       int c = split(value, ',', &arr);

       for (int i = 0; i < c; i++){
            pdata->al[i]=atof(arr[i]);
//            printf("string #%d: %lf \n", i, pdata->al[i]);
       }
       free(arr);
//       printf("\n");

    }else {
        return 0;  /* unknown section/name, error */
    }
    return 1;    return 1;
}

int main(int argc, char* argv[])
{
    data data_values;

    if (inp_parse("test.inp", handler, &data_values) < 0) {
        printf("Can't load 'test.inp'\n");
        return 1;
    }
 //   printf("Config loaded from 'test.inp': version=%d, name=%s, email=%s\n",
   //     data_values.version, data_values.name, data_values.email);
    return 0;
}
