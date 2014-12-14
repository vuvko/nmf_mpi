#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>
    
enum 
{
    MIN_STR_SIZE = 64
};

char *getline2(FILE *fin);
int inarrayd(int value, int *array, int size);
int arrayidx_str(char *value, char **array, int size);
char *make_param_name(
    char *buf, 
    int size, 
    const char *prefix, 
    const char *name);
void error_open(const char *file);
void error_undefined(const char *param);
void error_invalid(const char *param);

double sqr(double x);
#endif
