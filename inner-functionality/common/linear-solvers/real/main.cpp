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
  while (token != nullptr && i < n) {
    double entry_buffer;
    sscanf(token, "%lf", &entry_buffer);
    values[i++] = entry_buffer;

    token = strtok(nullptr, delims);
  }

  return (i == n);
}

int read_matrix_and_rhs(char *file_name, int &n, int &nnz, std::map<unsigned int, MatrixEntry>& mat, std::map<unsigned int, double>& rhs)
{
  FILE *file = fopen(file_name, "r");
  if(file == nullptr)
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
  while (fgets(row, MAX_ROW_LEN, file) != nullptr) {
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

// Test code.
void solve(LinearMatrixSolver<double> &solver, int n)
{
  solver.solve();
}

int main(int argc, char *argv[]) {
  int ret = 0;

  int n;
  int nnz;

  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, double > ar_rhs;

  double* sln = nullptr;
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

  if(strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    PetscMatrix<double> mat;
    PetscVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearMatrixSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
    PetscMatrix<double> mat;
    PetscVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearMatrixSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    CSCMatrix<double> mat;
    SimpleVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearMatrixSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
    CSCMatrix<double> mat;
    SimpleVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearMatrixSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "paralution") == 0) {
#ifdef WITH_PARALUTION
    ParalutionMatrix<double> mat;
    ParalutionVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    Hermes::Solvers::IterativeParalutionLinearMatrixSolver<double> solver(&mat, &rhs);
    if(atoi(argv[2]) != 1)
      solver.set_solver_type(BiCGStab);
    solver.set_precond(new ParalutionPrecond<double>(ILU));
    // Tested as of 13th August 2013.
    solver.set_max_iters(atoi(argv[2]) == 1 ? 1 : atoi(argv[2]) == 2 ? 2 : 5);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "paralution-block") == 0) {
#ifdef WITH_PARALUTION
    ParalutionMatrix<double> mat;
    ParalutionVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    Hermes::Solvers::IterativeParalutionLinearMatrixSolver<double> solver(&mat, &rhs);
    solver.set_precond(new ParalutionPrecond<double>(ILU));
    solver.set_max_iters(10000);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "superlu") == 0) {
#ifdef WITH_SUPERLU
    CSCMatrix<double> mat;
    SimpleVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    Hermes::Solvers::SuperLUSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "superlu-block") == 0) {
#ifdef WITH_SUPERLU
    CSCMatrix<double> mat;
    SimpleVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    Hermes::Solvers::SuperLUSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<double> mat;
    EpetraVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "aztecoo-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<double> mat;
    EpetraVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<double> mat;
    EpetraVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    if(AmesosSolver<double>::is_available("Klu")) {
      AmesosSolver<double> solver("Klu", &mat, &rhs);
      solve(solver, n);
      sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
    }
#endif
  }
  else if(strcasecmp(argv[1], "amesos-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<double> mat;
    EpetraVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    if(AmesosSolver<double>::is_available("Klu")) {
      AmesosSolver<double> solver("Klu", &mat, &rhs);
      solve(solver, n);
      sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
    }
#endif
  }
  else if(strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix<double> mat;
    SimpleVector<double> rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }
  else if(strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix<double> mat;
    SimpleVector<double> rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver<double> solver(&mat, &rhs);
    solve(solver, n);
    sln = new double[mat.get_size()]; memcpy(sln, solver.get_sln_vector(), mat.get_size() * sizeof(double));
#endif
  }

  bool success = true;
  if(sln)
  {
    switch(atoi(argv[2]))
    {
    case 1:
      success = Testing::test_value(sln[0], 4, "sln[0]", 1E-6) && success;
      success = Testing::test_value(sln[1], 2, "sln[1]", 1E-6) && success;
      success = Testing::test_value(sln[2], 3, "sln[2]", 1E-6) && success;
      break;
    case 2:
      success = Testing::test_value(sln[0], 2, "sln[0]", 1E-6) && success;
      success = Testing::test_value(sln[1], 3, "sln[1]", 1E-6) && success;
      success = Testing::test_value(sln[2], 1, "sln[2]", 1E-6) && success;
      success = Testing::test_value(sln[3], -3, "sln[3]", 1E-6) && success;
      success = Testing::test_value(sln[4], -1, "sln[4]", 1E-6) && success;
      break;
    case 3:
      success = Testing::test_value(sln[0], 1, "sln[0]", 1E-6) && success;
      success = Testing::test_value(sln[1], 2, "sln[1]", 1E-6) && success;
      success = Testing::test_value(sln[2], 3, "sln[2]", 1E-6) && success;
      success = Testing::test_value(sln[3], 4, "sln[3]", 1E-6) && success;
      success = Testing::test_value(sln[4], 5, "sln[4]", 1E-6) && success;
      break;
    }

    delete[] sln;

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
  else
    return 0;
}