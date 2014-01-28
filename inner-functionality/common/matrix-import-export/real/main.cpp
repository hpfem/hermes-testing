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
                          double* mat = (double*)calloc(n * n, sizeof(double));
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

template<typename AlgebraicEntity>
bool export_and_test(AlgebraicEntity* ent, int argc, char *argv[])
{
  bool success = true;
  char* s = new char[1000];
  char* s2 = new char[1000];

  for(int i = 0; i < 3; i++)
  {
    if(i == 0)
    {
      sprintf(s, "Real_%s_%s_plain.dat", argv[1], argv[2]);
      ent->export_to_file(s, "A", EXPORT_FORMAT_PLAIN_ASCII);
    }

    if(i == 1)
    {
      sprintf(s, "Real_%s_%s_market.dat", argv[1], argv[2]);
      ent->export_to_file(s, "A", EXPORT_FORMAT_MATRIX_MARKET);
    }

    if(i == 2)
    {
#ifdef WITH_MATIO
      sprintf(s, "Real_%s_%s_matio.dat", argv[1], argv[2]);
      ent->export_to_file(s, "A", EXPORT_FORMAT_MATLAB_MATIO);
      ent->free();

      ent->import_from_file(s, "A", EXPORT_FORMAT_MATLAB_MATIO);
      sprintf(s2, "Real_%s_%s_matio2.dat", argv[1], argv[2]);
      ent->export_to_file("fff", "A", EXPORT_FORMAT_PLAIN_ASCII);
      success = Testing::compare_files(s, s2) && success;
#else
      break;
#endif
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

  SparseMatrix<double> *mat = nullptr;
  Vector<double> *rhs = nullptr;
  if(strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    mat = new PetscMatrix<double>;
    rhs = new PetscVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "petsc_block") == 0) {
#ifdef WITH_PETSC
    mat = new PetscMatrix<double>;
    rhs = new PetscVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    mat = new CSCMatrix<double>;
    rhs = new SimpleVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "umfpack_block") == 0) {
#ifdef WITH_UMFPACK
    mat = new CSCMatrix<double>;
    rhs = new SimpleVector<double>;
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
  else if(strcasecmp(argv[1], "aztecoo_block") == 0) {
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
  else if(strcasecmp(argv[1], "amesos_block") == 0) {
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
    rhs = new SimpleVector<double>;
    build_matrix(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }
  else if(strcasecmp(argv[1], "mumps_block") == 0) {
#ifdef WITH_MUMPS
    mat = new MumpsMatrix<double>;
    rhs = new SimpleVector<double>;
    build_matrix_block(n, ar_mat, ar_rhs, mat, rhs);
#endif
  }

  if(mat || rhs)
  {
    bool success = mat ? export_and_test(mat, argc, argv) : true;
    success = rhs ? export_and_test(rhs, argc, argv) && success : success;

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
