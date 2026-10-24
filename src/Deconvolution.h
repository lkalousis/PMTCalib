
#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "iostream"
#include "iomanip"
#include <math.h>

#include "TObject.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TRandom.h"
#include "TCanvas.h"

#include "SPEResponse.h"
#include "SPEFitter.h"

#include "fftw3.h"

class Deconvolution : public TObject
{
 private:
  
  Double_t Q0;
  Double_t dQ0;
  
  Double_t s0;
  Double_t ds0;
  
  Int_t N;
  Double_t x1;
  Double_t x2;
  Double_t wbin;
  
 public:
  
  Deconvolution();
  
  virtual ~Deconvolution();
  
  Deconvolution( Double_t _Q0, Double_t _dQ0, Double_t _s0, Double_t _ds0 );
    
  void CleanUps( TH1D *h );
  Double_t fftPhase( Double_t vy, Double_t vz );
  void Deconvolute( Double_t _Q0, Double_t _s0 );
  TH1D* Run( TH1D* h, Int_t ntoys );
  
  TGraph* grRe;
  TGraph* grIm;
  
  Double_t mu_est;

  TH1D *h1;
  TH1D *h2;

  TH1D *h3;
  TH1D *h4;
  
  TH1D *h5;
    
  ClassDef( Deconvolution, 1 )
        
};

#endif
