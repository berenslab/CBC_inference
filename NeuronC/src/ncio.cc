/* module ncio in program nc */

/* defines io compatible with C++ */

extern "C" {
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
}

int ncfprintf(FILE *stream, const char *format, ...);

#ifndef MAX_PRINTF_STR
#define MAX_PRINTF_STR 16384
#endif

#include <iostream>

std::ostream *nc_stdout = &(std::cout);
std::ostream *nc_stderr = &(std::cerr);

int ncfprintf(FILE *stream, const char *format, ...)
{
  va_list argp;
  va_start(argp, format);

  if ((stream != stdout) && (stream != stderr)) {
    return vfprintf(stream, format, argp);
  }
  
  char buffer[MAX_PRINTF_STR];  
  int nchars = vsprintf(buffer, format, argp);

  if (stream == stdout) {
    *nc_stdout << buffer;
    nc_stdout->flush();
  } else {
    *nc_stderr << buffer;
    nc_stderr->flush();
  }

  return nchars;
}


/* Dummy call to avoid C++ features */

/* to use, uncomment this and comment out C++ version above */

/*
int ncfprintf(FILE *stream, const char *format, ...)
{
  va_list argp;
  va_start(argp, format);

  return vfprintf(stream, format, argp);
}
*/
