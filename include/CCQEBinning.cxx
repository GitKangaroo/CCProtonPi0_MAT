#include "include/CCQEBinning.h"
#include <cmath>
#include <iostream>

// set up bins

std::vector<double> GetQ2BinsGeV( ) {
  std::vector<double> bins;
  double binw = 0.00625;

  bins.push_back(0.0);
  bins.push_back(0.00625);//added
  bins.push_back(0.0125);
  bins.push_back(0.025);
  //bins.push_back(0.0375);
  bins.push_back(0.05);
  bins.push_back(0.1);
  //bins.push_back(0.15);//added
  bins.push_back(0.2);
  //bins.push_back(0.3);//added
  bins.push_back(0.4);
  //bins.push_back(0.6);//added
  bins.push_back(0.8);
  //bins.push_back(1.0);//added
  bins.push_back(1.2);
  bins.push_back(2.0);
  //bins.push_back(4.0);
  //bins.push_back(6.0);
  //bins.push_back(8.0);
  //bins.push_back(10.0);



  return bins;
}

std::vector<double> GetLogQ2BinsGeV( ) {



  std::vector<double> bins = GetQ2BinsGeV();
//  double binw = 0.00625;
//
//  bins.push_back(0.0001);
//  // bins.push_back(0.00625);//added
//  bins.push_back(0.0125);
//  bins.push_back(0.025);
//  //bins.push_back(0.0375);
//  bins.push_back(0.05);
//  bins.push_back(0.1);
//  // bins.push_back(0.15);//added
//  bins.push_back(0.2);
//  // bins.push_back(0.3);//added
//  bins.push_back(0.4);
//  //bins.push_back(0.6);//added
//  bins.push_back(0.8);
//  // bins.push_back(1.0);//added
//  bins.push_back(1.2);
//  bins.push_back(2.0);
//  bins.push_back(4.0);
//  //bins.push_back(6.0);
//  // bins.push_back(8.0);
//  //bins.push_back(10.0);

  int n = bins.size();
  std::vector<double> logbins = bins;
  for (int i = 0 ; i < n; i++){
    logbins[i] = std::log10(bins[i]);
    std::cout << " q2 " << logbins[i] << std::endl;
  }
  logbins[0]=-3.;

  return logbins;
}
