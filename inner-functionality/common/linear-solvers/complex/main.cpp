#include "../../../../testing-core/testing-core.h"

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
  MatrixEntry(int m, int n, ::complex  value) {
    this->m = m;
    this->n = n;
    this->value = value;
  }

  int m, n;
  ::complex   value;
};

void show_mat(const char *msg, std::map<unsigned int, MatrixEntry> mp)
{
  std::map<unsigned int, MatrixEntry>::iterator itr;

  std::cout << msg << std::endl;

  for(itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " <<
    (int) itr->second.m << " " <<
    (int) itr->second.n << " " <<
    (::complex) itr->second.value <<
    std::endl;

  std::cout << std::endl;
}

void show_rhs(const char *msg, std::map<unsigned int, ::complex > mp) {
  std::map<unsigned int, ::complex >::iterator itr;

  std::cout << msg << std::endl;

  for(itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " << (::complex) itr->second << std::endl;

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
  if(file == nullptr) return -1;

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
      if(read_n_nums(row, 1, buffer)) {
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
      if(read_n_nums(row, 4, buffer)) {
        ::complex  cmplx_buffer(buffer[2], buffer[3]);
        mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], (::complex) cmplx_buffer));
      }
      else
        state = STATE_RHS;
      break;

    case STATE_RHS:
      if(read_n_nums(row, 3, buffer)) {
        ::complex  cmplx_buffer(buffer[1], buffer[2]);
        rhs[(int) buffer[0]] = (::complex) cmplx_buffer;
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
                          ::complex  *mat = new ::complex [n * n];
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
                          ::complex   *rs = new ::complex [n];
                          for (std::map<unsigned int, ::complex >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
                            rs[it->first] = it->second;
                          }
                          unsigned int *u_rows = new unsigned int[n];
                          for (int i = 0; i < n; i++)
                            u_rows[i] = rows[i] >= 0 ? rows[i] : 0;
                          rhs->add(n, u_rows, rs);
                          rhs->finish();
}

// Test code.
void solve(LinearMatrixSolver<::complex> &solver, int n)
{
  solver.solve();
  ::complex  *sln = solver.get_sln_vector();
  for (int i = 0; i < n; i++)
    if(sln[i].imag() < 0.0)
      std::cout << std::endl << sln[i].real() << sln[i].imag();
    else
      std::cout << std::endl << sln[i].real() << ' + ' << sln[i].imag();
}

int main(int argc, char *argv[]) {
  int ret = 0;

  int n;
  int nnz;

  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, ::complex > ar_rhs;

  ::complex * sln = nullptr;

  if(read_matrix_and_rhs((char*)"in/linsys-cplx-4", n, nnz, ar_mat, ar_rhs) != 0)
    throw Hermes::Exceptions::Exception("Failed to read the matrix and rhs.");

  if(strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    PetscMatrix<::complex> mat;
    PetscVector<::complex> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearMatrixSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
    PetscMatrix<::complex> mat;
    PetscVector<::complex> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearMatrixSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    CSCMatrix<::complex> mat;
    SimpleVector<::complex> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearMatrixSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
    CSCMatrix<::complex> mat;
    SimpleVector<::complex> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearMatrixSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<::complex> mat;
    EpetraVector<::complex> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "aztecoo-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<::complex> mat;
    EpetraVector<::complex> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<::complex> mat;
    EpetraVector<::complex> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    if(AmesosSolver<::complex>::is_available("Klu")) {
      AmesosSolver<::complex> solver("Klu", &mat, &rhs);
      solve(solver, n);
      sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
    }
#endif
  }
  else if(strcasecmp(argv[1], "amesos-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<::complex> mat;
    EpetraVector<::complex> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    if(AmesosSolver<::complex>::is_available("Klu")) {
      AmesosSolver<::complex> solver("Klu", &mat, &rhs);
      solve(solver, n);
      sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
    }
#endif
  }
  else if(strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix<::complex> mat;
    SimpleVector<::complex> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else if(strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix<::complex> mat;
    SimpleVector<::complex> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver<::complex> solver(&mat, &rhs);
    solve(solver, n);
    sln = new ::complex [mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(::complex));
#endif
  }
  else
    ret = -1;

  if(sln)
  {
    std::cout << sln[0] << sln[1] << sln[2];

    if(std::abs(sln[0] - ::complex (0.800000, -0.600000)) > 1E-6 || std::abs(sln[1] - ::complex (0.470588, -0.882353)) > 1E-6 || std::abs(sln[2] - ::complex (0.486486, -0.918919)) > 1E-6)
      ret = -1;
    else
      ret = 0;

    // Test
    if(ret == -1)
      printf("Failure!\n");
    else
      printf("Success!\n");
  }
  else
    return 0;
}