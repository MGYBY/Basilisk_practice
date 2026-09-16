/** Simple reader/interpolator for OpenFOAM-style two-column .xy files. */
#ifndef BASILISK_PROFILE1D_H
#define BASILISK_PROFILE1D_H

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int n;
  int capacity;
  double * x;
  double * value;
} Profile1D;

static void profile1d_init (Profile1D * p)
{
  p->n = 0;
  p->capacity = 0;
  p->x = NULL;
  p->value = NULL;
}

static void profile1d_free (Profile1D * p)
{
  free (p->x);
  free (p->value);
  profile1d_init (p);
}

static void profile1d_append (Profile1D * p, double x, double value)
{
  if (p->n == p->capacity) {
    p->capacity = p->capacity ? 2*p->capacity : 128;
    p->x = (double *) realloc (p->x, p->capacity*sizeof(double));
    p->value = (double *) realloc (p->value,
                                    p->capacity*sizeof(double));
    if (!p->x || !p->value) {
      fprintf (stderr, "profile1d: memory allocation failed.\n");
      exit (1);
    }
  }
  if (p->n && x <= p->x[p->n - 1]) {
    fprintf (stderr, "profile1d: abscissae must be strictly increasing.\n");
    exit (1);
  }
  p->x[p->n] = x;
  p->value[p->n] = value;
  p->n++;
}

static void profile1d_read (Profile1D * p, const char * filename)
{
  profile1d_init (p);
  FILE * fp = fopen (filename, "r");
  if (!fp) {
    fprintf (stderr, "Cannot open profile '%s': %s\n", filename,
             strerror(errno));
    exit (1);
  }

  char line[1024];
  while (fgets(line, sizeof(line), fp)) {
    double x, value;
    if (sscanf(line, " ( %lf %lf )", &x, &value) == 2 ||
        sscanf(line, " %lf %lf", &x, &value) == 2)
      profile1d_append (p, x, value);
  }
  fclose (fp);

  if (p->n < 2) {
    fprintf (stderr, "Profile '%s' contains fewer than two points.\n",
             filename);
    exit (1);
  }
}

static double profile1d_eval (const Profile1D * p, double x)
{
  if (x <= p->x[0])
    return p->value[0];
  if (x >= p->x[p->n - 1])
    return p->value[p->n - 1];

  int lo = 0, hi = p->n - 1;
  while (hi - lo > 1) {
    const int mid = (lo + hi)/2;
    if (p->x[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }
  const double r = (x - p->x[lo])/(p->x[hi] - p->x[lo]);
  return (1. - r)*p->value[lo] + r*p->value[hi];
}


#endif
