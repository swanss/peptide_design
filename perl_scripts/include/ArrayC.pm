use Inline C => <<'END_C';

#define _DOUBLE 1
#define _INT    2

long allocateArrayC (int N, int type) {
  if (N <= 0) {
    printf("Error in allocateArrayC: array size must be positive, but is %d\n", N);
    return 0;
  }

  void *arr;
  switch (type) {
    case _DOUBLE:
      arr = (double*) calloc(N, sizeof(double));
      break;
    case _INT:
      arr = (int*) calloc(N, sizeof(int));
      break;
    default:
      printf("Error in allocateArrayC: unrecognized type %d\n", type);
      return 0;
  }
  if (arr == NULL) {
    printf("Error in allocateArrayC: allocation of array of size %d failed...\n", N);
    return 0;
  }

  // no matter what type, cast to long
  return((long) arr);
}

int freeArrayC (long arr_p) {
  free((void*) arr_p);
  return 1;
}

double accessDoubleArrayElementC (long arr_p, int i) {
  double* arr = (double*) arr_p;
  return arr[i];
}

int accessIntArrayElementC (long arr_p, int i) {
  int* arr = (int*) arr_p;
  return arr[i];
}

int setDoubleArrayElementC (long arr_p, int i, double val) {
  double *arr = (double*) arr_p;
  arr[i] = val;
  return 1;
}

int setIntArrayElementC (long arr_p, int i, double val) {
  double *arr = (double*) arr_p;
  arr[i] = val;
  return 1;
}

long resizeArrayC (long arr_p, int newN, int type) {
  int sz;
  void* arr;

  switch (type) {
    case _DOUBLE:
      arr = (double*) realloc((double*) arr_p, newN*sizeof(double));
      break;
    case _INT:
      arr = (int*) realloc((int*) arr_p, newN*sizeof(int));
      break;
    default:
      printf("Error in resizeArrayC: unrecognized type %d\n", type);
      return 0;
  }
  if (arr == NULL) {
    printf("Error in resizeArrayC: reallocation to size %d of type %d failed...\n", newN, type);
    return 0;
  }
  return ((long) arr);
}
END_C

1;
