#include "hermes_common.h"
#include "../../../../testing-core/testing-core.h"
#include <iostream>

using namespace Hermes::Algebra::DenseMatrixOperations;
using namespace Hermes::Solvers;

// Test of linear solvers.
// Read matrix and RHS from a file.

// Max row length in input file.
#define MAX_ROW_LEN  1024

class MatrixEntry
{
public:
  MatrixEntry() { }
  MatrixEntry(int m, int n, double value) {
    this->m = m;
    this->n = n;
    this->value = value;
  }

  int m, n;
  double  value;
};

void show_mat(const char *msg, std::map<unsigned int, MatrixEntry> mp)
{
  std::map<unsigned int, MatrixEntry>::iterator itr;

  std::cout << msg << std::endl;

  for(itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " <<
    (int) itr->second.m << " " <<
    (int) itr->second.n << " " <<
    (double) itr->second.value <<
    std::endl;

  std::cout << std::endl;
}

void show_rhs(const char *msg, std::map<unsigned int, double > mp) {
  std::map<unsigned int, double >::iterator itr;

  std::cout << msg << std::endl;

  for(itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " << (double) itr->second << std::endl;

  std::cout << std::endl;
}

bool testPrint(bool value, const char *msg, bool correct) {
  if(value == correct) {
    return true;
  }
  else {
    return false;
  }
}

bool read_n_nums(char *row, int n, double values[]) {
  int i = 0;
  char delims[] = " \t\n\r";
  char *token = strtok(row, delims);
  while (token != NULL && i < n) {
    double entry_buffer;
    sscanf(token, "%lf", &entry_buffer);
    values[i++] = entry_buffer;

    token = strtok(NULL, delims);
  }

  return (i == n);
}

int read_matrix_and_rhs(char *file_name, int &n, int &nnz, std::map<unsigned int, MatrixEntry>& mat, std::map<unsigned int, double>& rhs)
{
  FILE *file = fopen(file_name, "r");
  if(file == NULL)
    return -1;

  enum EState {
    STATE_N,
    STATE_MATRIX,
    STATE_RHS,
    STATE_NNZ
  }
  state = STATE_N;

  double buffer[4];
  char row[MAX_ROW_LEN];
  while (fgets(row, MAX_ROW_LEN, file) != NULL) {
    switch (state) {
    case STATE_N:
      if(read_n_nums(row, 1, buffer))
      {
        n = (int) buffer[0];
        state = STATE_NNZ;
      }
      break;

    case STATE_NNZ:
      if(read_n_nums(row, 1, buffer))
        nnz = (int) buffer[0];

      state = STATE_MATRIX;
      break;

    case STATE_MATRIX:
      if(read_n_nums(row, 3, buffer))
        mat[mat.size()] = (MatrixEntry ((int) buffer[0], (int) buffer[1], buffer[2]));
      else
        state = STATE_RHS;
      break; //case STATE_MATRIX break.

    case STATE_RHS:
      { // if cplx_2_real is false.
        if(read_n_nums(row, 2, buffer))
          rhs[(int) buffer[0]] = (double) buffer[1];
      }
      break;
    }
  }

  fclose(file);

  return 0;
}

void build_matrix(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, double > &ar_rhs,
                  SparseMatrix<double> *mat, Vector<double> *rhs)
{
  mat->prealloc(n);
  for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
    MatrixEntry &me = it->second;
    mat->pre_add_ij(me.m, me.n);
  }

  mat->alloc();
  for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
    MatrixEntry &me = it->second;
    mat->add(me.m, me.n, me.value);
  }
  mat->finish();

  rhs->alloc(n);
  for (std::map<unsigned int, double >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
    rhs->add(it->first, it->second);
  }
  rhs->finish();
}

void build_matrix_block(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, double > &ar_rhs,
                        SparseMatrix<double> *matrix, Vector<double> *rhs) {
                          matrix->prealloc(n);
                          for (int i = 0; i < n; i++)
                            for (int j = 0; j < n; j++)
                              matrix->pre_add_ij(i, j);

                          matrix->alloc();
                          double  **mat = new_matrix<double>(n, n);
                          int *cols = new int[n];
                          int *rows = new int[n];
                          for (int i = 0; i < n; i++) {
                            cols[i] = i;
                            rows[i] = i;
                          }
                          for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
                            MatrixEntry &me = it->second;
                            mat[me.m][me.n] = me.value;
                          }
                          matrix->add(n, n, mat, rows, cols);
                          matrix->finish();

                          rhs->alloc(n);
                          double  *rs = new double[n];
                          for (std::map<unsigned int, double >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
                            rs[it->first] = it->second;
                          }
                          unsigned int *u_rows = new unsigned int[n];
                          for (int i = 0; i < n; i++)
                            u_rows[i] = rows[i] >= 0 ? rows[i] : 0;
                          rhs->add(n, u_rows, rs);
                          rhs->finish();
}

bool export_and_test_matrix(Matrix<double>* mat, int argc, char *argv[])
{
  bool success = true;
  for(int i = 0; i < 3; i++)
  {
    char* s = new char[1000];

    if(i == 0)
    {
      sprintf(s, "Real-Matrix-%s-%s-plain.mat", argv[1], argv[2]);
      success = mat->export_to_file(s, "A", EXPORT_FORMAT_PLAIN_ASCII) && success;
    }
    if(i == 1)
    {
      sprintf(s, "Real-Matrix-%s-%s-market.mat", argv[1], argv[2]);
      success = mat->export_to_file(s, "A", EXPORT_FORMAT_MATRIX_MARKET) && success;
    }

    if(i == 2)
    {
#ifdef WITH_MATIO
      sprintf(s, "Real-Matrix-%s-%s-matio.mat", argv[1], argv[2]);
      success = mat->export_to_file(s, "A", EXPORT_FORMAT_MATLAB_MATIO) && success;
#else
      break;
#endif
    }

    if(success)
    {
      char* test_s = new char[1000];
#ifdef WIN32
      sprintf(test_s, "win\\%s-stored", s);
#else
      sprintf(test_s, "linux/%s-stored", s);
#endif
      success = Testing::compare_files(s, test_s) && success;
    } 
  }

  return success;
}

bool export_and_test_vector(Vector<double>* vec, int argc, char *argv[])
{
  bool success = true;
  for(int i = 0; i < 3; i++)
  {
    char* s = new char[1000];

    if(i == 0)
    {
      sprintf(s, "Real-Vector-%s-%s-plain.vec", argv[1], argv[2]);
      success = vec->export_to_file(s, "b", EXPORT_FORMAT_PLAIN_ASCII) && success;
    }

    if(i == 1)
    {
      sprintf(s, "Real-Vector-%s-%s-market.vec", argv[1], argv[2]);
      success = vec->export_to_file(s, "b", EXPORT_FORMAT_MATRIX_MARKET) && success;
    }

    if(i == 2)
    {
#ifdef WITH_MATIO
      sprintf(s, "Real-Vector-%s-%s-matio.vec", argv[1], argv[2]);
      success = vec->export_to_file(s, "b", EXPORT_FORMAT_MATLAB_MATIO) && success;
#else
      break;
#endif
    }

    if(success)
    {
      char* test_s = new char[1000];
#ifdef WIN32
      sprintf(test_s, "win\\%s-stored", s);
#else
      sprintf(test_s, "linux/%s-stored", s);
#endif
      success = Testing::compare_files(s, test_s) && success;
    } 
  }

  return success;
}

int main(int argc, char *argv[])
{
  int n;
  int nnz;

  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, double > ar_rhs;

  double* sln;
  switch(atoi(argv[2]))
  {
  case 1:
    if(read_matrix_and_rhs((char*)"in/linsys-1", n, nnz, ar_mat, ar_rhs) != 0)
      throw Hermes::Exceptions::Exception("Failed to read the matrix and rhs.");
    break;
  case 2:
    if(read_matrix_and_rhs((char*)"in/linsys-2", n, nnz, ar_mat, ar_rhs) != 0)
      throw Hermes::Exceptions::Exception("Failed to read the matrix and rhs.");
    break;
  case 3:
    if(read_matrix_and_rhs((char*)"in/linsys-3", n, nnz, ar_mat, ar_rhs) != 0)
      throw Hermes::Exceptions::Exception("Failed to read the matrix and rhs.");
    break;
  }

  SparseMatrix<double> *mat = NULL;
  Vector<double> *rhs = NULL;
  if(strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    mat = new PetscMatrix<double>;
    rhs = new PetscVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
    mat = new PetscMatrix<double>;
    rhs = new PetscVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    mat = new UMFPackMatrix<double>;
    rhs = new UMFPackVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
    mat = new UMFPackMatrix<double>;
    rhs = new UMFPackVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix<double>;
    rhs = new EpetraVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "aztecoo-block") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix<double>;
    rhs = new EpetraVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix<double>;
    rhs = new EpetraVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
  }
#endif
}
  else if(strcasecmp(argv[1], "amesos-block") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix<double>;
    rhs = new EpetraVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
  }
#endif
  }
  else if(strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    mat = new MumpsMatrix<double>;
    rhs = new MumpsVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
    mat = new MumpsMatrix<double>;
    rhs = new MumpsVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }

  if(mat || rhs)
  {
    bool success = mat ? export_and_test_matrix(mat, argc, argv) : true;
    success = rhs ? export_and_test_vector(rhs, argc, argv) && success : success;

    if(success)
    {

      printf("Success!\n");
      return 0;
    }
    else
    {
      printf("Failure!\n");
      return -1;
    } 
  }
}
