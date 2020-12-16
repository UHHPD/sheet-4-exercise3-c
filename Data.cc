#include "Data.hh"

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <cmath>

using namespace std;

Data::Data(const std::string& filename) {
  ifstream file(filename);

  // check if file opened correctly
  if (!file.is_open()) {
    throw std::runtime_error("Error could not open file with name: " +
                             filename);
  }

  int size = -1;
  file >> size;

  // read in data from file: m_siz+1 bin edges
  for (int i = 0; i < size + 1; ++i) {
    double edge;
    file >> edge;
    m_bins.push_back(edge);
  }
  // read in data from file: m_siz bin contents
  for (int i = 0; i < size; ++i) {
    double entries;
    file >> entries;
    m_data.push_back(entries);
  }
  // read in uncertainty from file: m_siz bin errors
  for (int i = 0; i < size; ++i) {
    double uncertainties;
    file >> uncertainties;
    m_errors.push_back(uncertainties);
  }

  // done! close the file
  file.close();

  assertSizes();
};


Data Data::averageDatasets(const Data& in) {
  if(checkCompatibility(in,3)>5) throw runtime_error("Datasets not compatible!!!"); // at what point do we define two datasets to be compatible?
  // I'm gonna assume that all datasets have the same size and bin structure
  Data out;
  out.m_bins = in.m_bins;

  for(int i =0; i<size(); i++) {
    double w1 = 1/pow(error(i),2);
    double w2 = 1/pow(in.error(i),2);
    double y1 = measurement(i);
    double y2 = in.measurement(i);

    double y = (w1*y1+w2*y2)/(w1+w2);
    double dy = 1/sqrt(w1+w2);
    //cout << "y: " << y << endl;
    //cout << "dy: " << dy << endl;
    out.m_data.push_back(y);
    out.m_errors.push_back(dy);
  }
  
  return out;
}


int Data::checkCompatibility(const Data& in, int n) {
  int count = 0;
  for(int i = 0; i < size(); i++) {
    double uncert_this = error(i);
    double uncert_that = in.error(i);
    double uncert = sqrt(pow(error(i),2)+pow(in.error(i),2));
    double deltaY = abs(measurement(i)-in.measurement(i));
    //cout << "uncert: " << uncert << endl;
    //cout << "deltaY: " << deltaY << endl;
    bool compatible = (deltaY < n*uncert);
    //cout << "deltaY>" << n << "*uncert? " << compatible << endl;
    if(!compatible) count++;
  }

  return count; // #count == #of data points incompatible within n std's
}

double Data::fitFunction(double x) {
  double alpha = 0.005;
  double beta = -0.00001;
  double gamma = 0.08;
  double delta = 0.015;
  double f = alpha + beta*x + gamma*exp(-delta*x);
  return f;
}

double Data::chi2perNdf() {
  double chi2 = 0;
  for(int i=0; i<size(); i++) {
    double y_meas = measurement(i);
    double y_pred = fitFunction(binCenter(i));
    double error_exp = error(i);
    double summand = pow(y_meas-y_pred, 2) / pow(error_exp, 2);
    chi2 += summand;
  }
  double ndf = 52;
  return chi2/ndf;
}

void Data::assertSizes() { assert(m_data.size() + 1 == m_bins.size()); }
