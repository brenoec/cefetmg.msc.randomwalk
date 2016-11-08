#include <iostream>
#include <fstream>
#include <limits>
#include <string>

#include "ranq.cpp"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648

using namespace std;

int main(int argc, char *argv[]) {

  // pseudorandom number generator
  Ranq2 *random = new Ranq2();

  int steps       = 10000;
  int randomwalks = 1000;

  int qcount = 1;

  // q values for q-gaussian distribution
  double q[qcount] = { 1.8 };

  double**** acc    = new double***[2]();
  double**** qgauss = new double***[2]();
  double*** uniform = new double**[2]();

  // matrix initialization
  for (int i = 0; i < 2; i++) {
        acc[i] = new double**[qcount + 1]();
     qgauss[i] = new double**[qcount]();
    uniform[i] = new double*[randomwalks]();

    for (int j = 0; j < qcount + 1; j++) {
      acc[i][j] = new double*[randomwalks]();
      if (j < qcount) {
        qgauss[i][j] = new double*[randomwalks]();
      }

      for (int k = 0; k < randomwalks; k++) {
        acc[i][j][k] = new double[steps]();
          if (j < qcount) {
            qgauss[i][j][k] = new double[steps]();
          }
        uniform[i][k] = new double[steps]();
      }
    }
  }

  // open output files
  ofstream* outputFiles = new ofstream[(qcount+1)*4];
  for (int i = 0; i < (qcount+1)*4; i++) {
    int j;
    int k = floor(i/4);

    string type;
    string path;

    switch (i%2) {
      case 0:
      type = "values";
      break;

      case 1:
      type = "accumulated";
      break;
    }

    switch (i%4) {
      case 0:
      case 1:
      j = 0;
      break;

      case 2:
      case 3:
      j = 1;
      break;
    }

    if (k < qcount) {
      path = "./output/data/" + to_string(q[k]) + "_" + type + "_" + to_string(j) + ".txt";
    }

    else {
      path = "./output/data/uniform__" + type + "_" + to_string(j) + ".txt";
    }

    outputFiles[i].open(path);
  }

  // cout max precision for double
  typedef std::numeric_limits<double> dbl;
  cout.precision(dbl::max_digits10);
  cout << fixed;

  // uniform pseudorandom values alias
  double n1 = 0.0;
  double n2 = 0.0;

  // generate data
  for (int i = 0; i < steps; i++) {

    // step zero, set all to zero and continue to next step
    if (i == 0) {
      for (int j = 0; j < randomwalks; j++) {

        uniform[0][j][i] = 0;
        uniform[1][j][i] = 0;

        for (int k = 0; k < qcount; k++) {

          // gaussian: logarithm
          if (q[k] == 1) {
            qgauss[0][1][j][i] = 0;
            qgauss[1][1][j][i] = 0;
          }
          // q-gaussian: q-logarithm
          else {
            qgauss[0][k][j][i] = 0;
            qgauss[1][k][j][i] = 0;
          }

          // accumulate pseudorandom values for gaussian and q-gaussian
          acc[0][k][j][i] = 0;
          acc[1][k][j][i] = 0;
        }

        // accumulate pseudorandom values for uniform
        acc[0][qcount][j][i] = 0;
        acc[1][qcount][j][i] = 0;

        // write to file
        for (int k = 0; k < qcount; k++) {
          outputFiles[0 + k*4] << qgauss[0][k][j][i] << "\t";
          outputFiles[1 + k*4] << acc[0][k][j][i] << "\t";
          outputFiles[2 + k*4] << qgauss[1][k][j][i] << "\t";
          outputFiles[3 + k*4] << acc[1][k][j][i] << "\t";
        }

        outputFiles[0 + qcount*4] << uniform[0][j][i] << "\t";
        outputFiles[1 + qcount*4] << acc[0][qcount][j][i] << "\t";
        outputFiles[2 + qcount*4] << uniform[1][j][i] << "\t";
        outputFiles[3 + qcount*4] << acc[1][qcount][j][i] << "\t";
      }

      // file breakline
      for (int k = 0; k < (qcount+1)*4; k++) {
        outputFiles[k] << endl;
      }

      continue;
    }

    for (int j = 0; j < randomwalks; j++) {

      // generate pseudorandom values for uniform distribution
      n1 = random->nextd();
      n2 = random->nextd();

      uniform[0][j][i] = (n1 - 0.5) * 2;
      uniform[1][j][i] = (n2 - 0.5) * 2;

      // generate pseudorandom values for q-gaussian distribution
      for (int k = 0; k < qcount; k++) {

        double qline = (1 + q[k]) / (3 - q[k]);

        // gaussian: logarithm
        if (q[k] == 1.0) {
          n1 != 0.0 ? qgauss[0][k][j][i] = sqrt(-2.0*log(n1))*cos(2.0*PI*n2) : qgauss[0][k][j][i] = 0.0;
          n2 != 0.0 ? qgauss[1][k][j][i] = sqrt(-2.0*log(n1))*sin(2.0*PI*n2) : qgauss[1][k][j][i] = 0.0;
        }
        // q-gaussian: q-logarithm
        else {
          n1 != 0.0 ? qgauss[0][k][j][i] = sqrt(-2.0*((pow(n1, 1.0-qline)-1.0)/(1-qline)))*cos(2.0*PI*n2) : qgauss[0][k][j][i] = 0.0;
          n2 != 0.0 ? qgauss[1][k][j][i] = sqrt(-2.0*((pow(n1, 1.0-qline)-1.0)/(1-qline)))*sin(2.0*PI*n2) : qgauss[1][k][j][i] = 0.0;
        }

        // accumulate pseudorandom values for gaussian and q-gaussian
        acc[0][k][j][i] = acc[0][k][j][i-1] + qgauss[0][k][j][i];
        acc[1][k][j][i] = acc[1][k][j][i-1] + qgauss[1][k][j][i];
      }

      // accumulate pseudorandom values for uniform
      acc[0][qcount][j][i] = acc[0][qcount][j][i-1] + uniform[0][j][i];
      acc[1][qcount][j][i] = acc[1][qcount][j][i-1] + uniform[1][j][i];

      // write to file
      for (int k = 0; k < qcount; k++) {
        outputFiles[0 + k*4] << qgauss[0][k][j][i] << "\t";
        outputFiles[1 + k*4] << acc[0][k][j][i] << "\t";
        outputFiles[2 + k*4] << qgauss[1][k][j][i] << "\t";
        outputFiles[3 + k*4] << acc[1][k][j][i] << "\t";
      }

      outputFiles[0 + qcount*4] << uniform[0][j][i] << "\t";
      outputFiles[1 + qcount*4] << acc[0][qcount][j][i] << "\t";
      outputFiles[2 + qcount*4] << uniform[1][j][i] << "\t";
      outputFiles[3 + qcount*4] << acc[1][qcount][j][i] << "\t";
    }

    // file breakline
    for (int k = 0; k < (qcount+1)*4; k++) {
      outputFiles[k] << endl;
    }
  }

  // free memory
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < qcount + 1; j++) {
      for (int k = 0; k < randomwalks; k++) {
        delete [] acc[i][j][k];
        if (j < qcount) {
          delete [] qgauss[i][j][k];
        }
      }
      delete [] acc[i][j];
      if (j < qcount) {
        delete [] qgauss[i][j];
      }
    }
    delete [] acc[i];
    delete [] qgauss[i];
  }

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < randomwalks; j++) {
      delete [] uniform[i][j];
    }
    delete [] uniform[i];
  }

  delete [] acc;
  delete [] qgauss;
  delete [] uniform;

  // close files
  for (int i = 0; i < (qcount+1)*4; i++) {
    outputFiles[i].close();
  }

  delete [] outputFiles;

  cout << ".\n";

  return 0;
}
