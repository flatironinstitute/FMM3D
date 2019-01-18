
#include <string.h>
#include "cprini.h"

/*
*
***************************************************************************
*
*        This is the end of the debugging code, and the beginning
*        of the printing functions proper. The code below contains
*        five user-callable functions: cprini, cprin, cprin2, cprinf,
*        cprin_start_stop. Their functions are as follows.
*
*  cprini - communicates to the function cprin8 (NOT a user-callable
*        function) the names of the two files on which the output 
*        will be written. Please note that the name "stdout" is 
*        treated differently from all other names: specifying it 
*        will cause the output to be written on the standard output
*        - screen, most probably.
*  cprin - prints a real (of the type "float") array
*  cprin2 - prints a double precision (of the type "double") array
*  cprinf - prints an integer array
*  cprina - prints a character array, string
*  cprin_start_stop - starts and stops printing to either or both
*        output files
*
***************************************************************************
*
*/



void cprin_init(char *str17,  char *str27)
{
  char *mes;
  float *ap;
  int *afp;
  double *adp;
  char *acp;
  int n, m;
  int itype;
  int i1, i2;
  //int *adp,*ap,itype,i1,i2,mes,afp,n,acp;

  /* This function initializes which files to print to.  */
  /* To print to ONLY one file, set one of str17, str18 to " ".  */
  /* i.e., make the call cprini("stdout"," ") to print ONLY to  */
  /* the screen.  Anything other than " ", such as "    ", will cause   */
  /* a file to be created, in this case "    ".  */

  itype=11;
  cprin_master(mes, ap, afp, adp, acp, m, n, itype, str17, str27, i1, i2);
  return;
}





void cprin_start_stop(int i1, int i2)
{
  char *mes;
  float *ap;
  int *afp;
  double *adp;
  char *acp;
  int n, m;
  int itype;
  char *str17,*str27;

  // . . . stop/resume.  i1 and i2 control printing to file
  // str17 and str27, respectively.  0 means stop printing,
  // and 1 means print.

  itype=21;
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}







void cprin_message(char *mes)
{
  float *ap;
  int *afp;
  char *acp;
  int itype;
  char *str17, *str27;
  int i1, i2;
  double *adp;
  int m, n;
  
  itype = 99;
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}





void cprin_skipline(int n)
{
  float *ap;
  int *afp;
  char *acp;
  int itype;
  char *str17, *str27;
  int i1, i2;
  int status;
  char *mes = NULL;
  double *adp;
  int m;
  
  //  printf("in cprin_skipline, n  = %d\n", n);
  //exit(status);
  
  //int *ap,*afp,itype,str17,str27,i1,i2,acp;
  /*
   *   Print double precision data
   */
  itype = 77;
  
  //printf("in cprin_skipline, itype  = %d\n", itype);
  
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}





void cprinz(char *mes, double _Complex *adp, int n)
{
  float *ap;
  int *afp;
  char *acp;
  int itype;
  char *str17, *str27;
  int i1, i2;
  int status;
  int n2, m;

  //printf("n  = %ld\n", n);
  //exit(status);
  
  //int *ap,*afp,itype,str17,str27,i1,i2,acp;
  /*
   *   Print double precision data
   */
  itype=7;
  n2 = 2*n;
  cprin_master(mes,ap,afp, (double *)adp,
               acp,m,n2,itype,str17,str27,i1,i2);
  return;
}




void cprinf(char *mes, int *afp, int n)
{
  float *ap;
  int itype, i1, i2, m;
  double *adp;
  char *acp, *str17, *str27;
  
  //
  // print integer data
  //
  itype = 2;
  cprin_master(mes, ap, afp, adp, acp, m, n, itype, str17, str27, i1, i2);
  return;
}





void cprind(char *mes, double *adp, int n) {

  float *ap;
  int itype, i1, i2, m, *afp;
  char *acp, *str17, *str27;
  
  //
  // print double precision data
  //
  itype = 3;
  cprin_master(mes, ap, afp, adp, acp, m, n, itype, str17, str27, i1, i2);
  return;
}





void cprind_matrix(char *mes, double *adp, int m, int n)
{
  float *ap;
  int *afp;
  char *acp;
  int itype;
  char *str17, *str27;
  int i1, i2;
  int status;
  int n2;

  //printf("n  = %ld\n", n);
  //exit(status);
  
  //int *ap,*afp,itype,str17,str27,i1,i2,acp;
  /*
   *   Print double precision data
   */
  itype=33;
  cprin_master(mes,ap,afp,adp,acp,m,n,itype,str17,str27,i1,i2);
  return;
}





void cprin_master(char *mes, float *ap, int *afp, double *adp, char *acp, int m,
                  int n, int itype, char *str17, char *str27, int i1, int i2)
{
  static int ifprint1,ifprint2,ifstr1,ifstr2;
  static FILE *st1,*st2;
  int iii;
  
  // If this is the initialization call - open the 
  // files str17, str27 after checking if either is null, i.e. " "

  //printf("in cprin_master, itype = %d\n", itype);
  //return;

  
  if(itype==11)
    {
      ifprint1=0;
      ifprint2=0;
      
      ifstr1=0;
      ifstr2=0;
      
      iii=strcmp(str17," ");
      if(iii != 0) ifprint1=1;
      if(iii != 0) ifstr1=1;
      
      iii=strcmp(str27," ");
      if(iii != 0) ifprint2=1;
      if(iii != 0) ifstr2=1;
      
      if(ifprint1 == 1)
        {
          iii=strcmp(str17,"stdout");
          if(iii != 0) st1=fopen(str17,"w");
          if(iii == 0) st1=stdout;
        }
      
      if(ifprint2 == 1) {
        iii=strcmp(str27,"stdout");
        if(iii != 0) st2=fopen(str27,"w");
        if(iii == 0) st2=stdout;
      }
      
      return;
    }
  
  //
  //   If this is the "stop/resume" call - stop/resume printing 
  //
  if(itype==21){
    if(i1==0) ifprint1=0;
    if((i1 != 0) && (ifstr1 != 0)) ifprint1=1;
    if(i2==0) ifprint2=0;
    if((i2 != 0) && (ifstr2 != 0)) ifprint2=1;
    return;
  }

  //printf("ifprint1 = %d\n", ifprint1);
  //printf("ifprint2 = %d\n", ifprint2);

  
  if(ifprint1 !=0) cprin_all(mes,ap,afp,adp,acp,m,n,itype,st1);
  if(ifprint2 !=0) cprin_all(mes,ap,afp,adp,acp,m,n,itype,st2);
  
  return;
}





void cprin_all(char *mes, float *ap, int *afp, double *adp, char *acp,
               int m, int n, int itype, FILE *str)
{
  int i;
  /*
   *   Process the message
   */

  if(mes) fprintf(str,"%s\n",mes);
  if (itype == 99) return;


  
  
  // Process the double precision data to be printed
  if(itype == 3) {
    for(i=0; i<n; i=i+1) {
      fprintf(str,"  %11.4le", adp[i]);
      if(i%6==5 || i==n-1) fprintf(str,"\n");
    }
    return;
  }


  //
  // print out the double precision matrix
  //
  if(itype ==33) {
    int ijk = 0;
    int j;
    
    //printf("m = %d\n", m);
    //printf("n = %d\n", n);      
    //return;

    
    for(i=0; i<m; i=i+1) {
      for (j=0; j<n; j++) {
        fprintf(str,"  %11.4le", adp[ijk]);
        ijk++;
      }
      fprintf(str,"\n");
    }
    return;
  }

  


  /*
   *   Process the complex double precision data to be printed
   */
  if(itype ==7)
    {
      for(i=0; i<n; i=i+2)
        {
          fprintf(str," (%11.4le %+11.4le)", adp[i], adp[i+1]);
          if(i%6==4 || i==n-2) fprintf(str,"\n");
        }
      return;
    }

/*
*   Process the integer data to be printed
*/
    if(itype ==2)
    {
        for(i=0; i<n; i=i+1)
        {
            fprintf(str," %7d", afp[i]);
            if(i%10==9 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

/*
*   Process the single precision data to be printed
*/
    if(itype ==1)
    {
        for(i=0; i<n; i=i+1)
        {
            fprintf(str,"  %11.4le", ap[i]);
            if(i%6==5 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

/*
*   Process the character data to be printed
*/
    if(itype==4) {
      for(i=0; i<n; i++) {
        fprintf(str,"%c", acp[i]);
        if(i%60==59 || i==n-1) fprintf(str,"\n");
      }
      return;
    }

    
    //
    // insert a line skip
    //
    if(itype == 77) {
      for (i=0; i<n; i++) {
        fprintf(str,"\n");
      }
    }
    return;



}
