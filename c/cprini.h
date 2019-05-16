
#include <stdio.h>


void cprin_init(char *str1, char *str2);

void cprin_master(char *mes, float *ap, int *afp, double *adp, char *acp, int m,
                  int n, int itype, char *str17, char *str27, int i1, int i2);

void cprin_all(char *mes, float *ap, int *afp, double *adp, char *acp,
               int m, int n, int itype, FILE *str);

void cprinf(char *mes, int *ip, int n);

void cprind(char *mes, double *adp, int n);

void cprind_matrix(char *mes, double *adp, int m, int n);

void cprinz(char *mes, double _Complex *adp, int n);

void cprin_message(char *mes);

void cprin_skipline(int n);

void cprin_start_stop(int i1, int i2);

