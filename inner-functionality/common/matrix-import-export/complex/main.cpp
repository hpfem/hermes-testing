#include "../../../../testing-core/testing-core.h"

using namespace Hermes::Algebra::DenseMatrixOperations;
using namespace Hermes::Solvers;

typedef std::complex<double> complex;

// Test of linear solvers.
// Read matrix and RHS from a file.

// Max row length in input file.
#define MAX_ROW_LEN  1024

class MatrixEntry
{
public:
  MatrixEntry() { }
  MatrixEntry(int m, int n, ::complex  value) {
    this->m = m;
    this->n = n;
    this->value = value;
  }

  int m, n;
  ::complex   value;
};

bool export_and_test_matrix(Matrix<::complex>* mat, int argc, char *argv[])
{
  bool success = true;
  for (int i = 0; i < 3; i++)
  {
    char* s = new char[1000];
    char* s2 = new char[1000];

    if (i == 0)
    {
      sprintf(s, "Complex_Matrix_%s_plain.mat", argv[1]);
      try
      {
        mat->export_to_file(s, "A", EXPORT_FORMAT_PLAIN_ASCII);
      }
      catch (Hermes::Exceptions::MethodNotImplementedException& e)
      {
        e.print_msg();
      }
    }

    if (i == 1)
    {
      sprintf(s, "Complex_Matrix_%s_market.mat", argv[1]);
      try
      {
        mat->export_to_file(s, "A", EXPORT_FORMAT_MATRIX_MARKET);
      }
      catch (Hermes::Exceptions::MethodNotImplementedException& e)
      {
        e.print_msg();
      }
    }

    if (i == 2)
    {
#ifdef WITH_MATIO
      sprintf(s, "Complex_Matrix_%s_matio.mat", argv[1]);
      try
      {
        mat->export_to_file(s, "A", EXPORT_FORMAT_MATLAB_MATIO);
        mat->free();

        mat->import_from_file(s, "A", EXPORT_FORMAT_MATLAB_MATIO);
        sprintf(s2, "Real_%s_matio2.dat", argv[1]);
        mat->export_to_file(s2, "A", EXPORT_FORMAT_MATLAB_MATIO);
        success = Testing::compare_files(s, s2) && success;
      }
      catch (Hermes::Exceptions::MethodNotImplementedException& e)
      {
        e.print_msg();
      }
#else
      break;
#endif
    }
  }

  return success;
}

bool export_and_test_vector(Vector<::complex>* vec, int argc, char *argv[])
{
  bool success = true;
  for (int i = 0; i < 3; i++)
  {
    char* s = new char[1000];
    char* s2 = new char[1000];

    if (i == 0)
    {
      sprintf(s, "Complex_Vector_%s_plain.vec", argv[1]);
      try
      {
        vec->export_to_file(s, "rhs", EXPORT_FORMAT_PLAIN_ASCII);
      }
      catch (Hermes::Exceptions::MethodNotImplementedException& e)
      {
        e.print_msg();
      }
      catch (Hermes::Exceptions::MethodNotOverridenException& e)
      {
        e.print_msg();
      }
    }

    if (i == 1)
    {
      sprintf(s, "Complex_Vector_%s_market.vec", argv[1]);
      try
      {
        vec->export_to_file(s, "rhs", EXPORT_FORMAT_MATRIX_MARKET);
      }
      catch (Hermes::Exceptions::MethodNotImplementedException& e)
      {
        e.print_msg();
      }
      catch (Hermes::Exceptions::MethodNotOverridenException& e)
      {
        e.print_msg();
      }
    }

    if (i == 2)
    {
#ifdef WITH_MATIO
      try
      {
        sprintf(s, "Complex_Vector_%s_matio.vec", argv[1]);
        vec->export_to_file(s, "rhs", EXPORT_FORMAT_MATLAB_MATIO);
        vec->free();

        vec->import_from_file(s, "rhs", EXPORT_FORMAT_MATLAB_MATIO);
        sprintf(s2, "Real_%s_%s_matio2.dat", argv[1], argv[2]);
        vec->export_to_file(s2, "rhs", EXPORT_FORMAT_MATLAB_MATIO);
        success = Testing::compare_files(s, s2) && success;
      }
      catch (Hermes::Exceptions::MethodNotImplementedException& e)
      {
        e.print_msg();
      }
      catch (Hermes::Exceptions::MethodNotOverridenException& e)
      {
        e.print_msg();
      }
#else
      break;
#endif
    }
  }

  return success;
}

void show_mat(const char *msg, std::map<unsigned int, MatrixEntry> mp)
{
  std::map<unsigned int, MatrixEntry>::iterator itr;

  std::cout << msg << std::endl;

  for (itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int)itr->first << ": " <<
    (int)itr->second.m << " " <<
    (int)itr->second.n << " " <<
    (::complex) itr->second.value <<
    std::endl;

  std::cout << std::endl;
}

void show_rhs(const char *msg, std::map<unsigned int, ::complex > mp) {
  std::map<unsigned int, ::complex >::iterator itr;

  std::cout << msg << std::endl;

  for (itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int)itr->first << ": " << (::complex) itr->second << std::endl;

  std::cout << std::endl;
}

bool testPrint(bool value, const char *msg, bool correct) {
  if (value == correct) {
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
  while (token != nullptr && i < n) {
    double entry_buffer;
    sscanf(token, "%lf", &entry_buffer);
    values[i++] = entry_buffer;

    token = strtok(nullptr, delims);
  }

  return (i == n);
}

int read_matrix_and_rhs(char *file_name, int &n, int &nnz,
  std::map<unsigned int, MatrixEntry>& mat, std::map<unsigned int, ::complex >& rhs)
{
  FILE *file = fopen(file_name, "r");
  if (file == nullptr) return -1;

  enum EState {
    STATE_N,
    STATE_MATRIX,
    STATE_RHS,
    STATE_NNZ
  }
  state = STATE_N;

  double buffer[4];
  char row[MAX_ROW_LEN];
  while (fgets(row, MAX_ROW_LEN, file) != nullptr) {
    switch (state) {
    case STATE_N:
      if (read_n_nums(row, 1, buffer)) {
        n = (int)buffer[0];

        state = STATE_NNZ;
      }
      break;

    case STATE_NNZ:
      if (read_n_nums(row, 1, buffer))
        nnz = (int)buffer[0];

      state = STATE_MATRIX;
      break;

    case STATE_MATRIX:
      if (read_n_nums(row, 4, buffer)) {
        ::complex  cmplx_buffer(buffer[2], buffer[3]);
        mat[mat.size()] = (MatrixEntry((int)buffer[0], (int)buffer[1], (::complex) cmplx_buffer));
      }
      else
        state = STATE_RHS;
      break;

    case STATE_RHS:
      if (read_n_nums(row, 3, buffer)) {
        ::complex  cmplx_buffer(buffer[1], buffer[2]);
        rhs[(int)buffer[0]] = (::complex) cmplx_buffer;
      }
      break;
    }
  }

  fclose(file);

  return 0;
}

void build_matrix(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, ::complex > &ar_rhs,
  SparseMatrix<::complex> *mat, Vector<::complex> *rhs)
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
  for (std::map<unsigned int, ::complex >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
    rhs->add(it->first, it->second);
  }
  rhs->finish();
}

void build_matrix_block(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, ::complex > &ar_rhs,
  SparseMatrix<::complex> *matrix, Vector<::complex> *rhs) {
  matrix->prealloc(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      matrix->pre_add_ij(i, j);

  matrix->alloc();
  ::complex * mat = new ::complex[n * n];
  int *cols = new int[n];
  int *rows = new int[n];
  for (int i = 0; i < n; i++) {
    cols[i] = i;
    rows[i] = i;
  }
  for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
    MatrixEntry &me = it->second;
    mat[me.m * n + me.n] = me.value;
  }
  matrix->add(n, n, mat, rows, cols, n);
  matrix->finish();

  rhs->alloc(n);
  ::complex   *rs = new ::complex[n];
  for (std::map<unsigned int, ::complex >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
    rs[it->first] = it->second;
  }
  unsigned int *u_rows = new unsigned int[n];
  for (int i = 0; i < n; i++)
    u_rows[i] = rows[i] >= 0 ? rows[i] : 0;
  rhs->add(n, u_rows, rs);
  rhs->finish();
}

int main(int argc, char *argv[])
{
  int n;
  int nnz;

  char* argv_local[2];
  if (argc < 2)
  {
    argv_local[1] = new char[20];
    sprintf(argv_local[1], "mumps");
  }
  else
  {
    argv_local[1] = argv[1];
  }

  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, ::complex > ar_rhs;
  ::complex * sln;

  if (read_matrix_and_rhs((char*)"in/linsys-cplx-4", n, nnz, ar_mat, ar_rhs) != 0)
    throw Hermes::Exceptions::Exception("Failed to read the matrix and rhs.");

  SparseMatrix<::complex> *mat = nullptr;
  Vector<::complex> *rhs = nullptr;

  if (strcasecmp(argv_local[1], "petsc") == 0) {
#ifdef WITH_PETSC
    mat = new PetscMatrix < ::complex > ;
    rhs = new PetscVector < ::complex > ;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "petsc_block") == 0) {
#ifdef WITH_PETSC
    mat = new PetscMatrix < ::complex > ;
    rhs = new PetscVector < ::complex > ;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    mat = new  CSCMatrix < ::complex > ;
    rhs = new SimpleVector < ::complex > ;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "umfpack_block") == 0) {
#ifdef WITH_UMFPACK
    mat = new CSCMatrix < ::complex > ;
    rhs = new SimpleVector < ::complex > ;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix < ::complex > ;
    rhs = new EpetraVector < ::complex > ;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "aztecoo_block") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix < ::complex > ;
    rhs = new EpetraVector < ::complex > ;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix < ::complex > ;
    rhs = new EpetraVector < ::complex > ;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "amesos_block") == 0) {
#ifdef WITH_TRILINOS
    mat = new EpetraMatrix < ::complex > ;
    rhs = new EpetraVector < ::complex > ;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    mat = new MumpsMatrix < ::complex > ;
    rhs = new SimpleVector < ::complex > ;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if (strcasecmp(argv_local[1], "mumps_block") == 0) {
#ifdef WITH_MUMPS
    mat = new MumpsMatrix < ::complex > ;
    rhs = new SimpleVector < ::complex > ;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }

  if (mat && rhs)
  {
    bool success = export_and_test_matrix(mat, argc, argv_local) && export_and_test_vector(rhs, argc, argv_local);

    if (success)
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
  return 0;
}
