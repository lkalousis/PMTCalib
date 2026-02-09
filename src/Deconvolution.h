
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

  Double_t x1;

  Float_t mu_bf;
  
 public:
  
  Deconvolution();
  
  virtual ~Deconvolution();
  
  Deconvolution( Double_t _Q0, Double_t _dQ0, Double_t _s0, Double_t _ds0 );

  TH1D* CleanUps( TH1D *h );
  Double_t fftPhase( Double_t vy, Double_t vz );

  TH1D* Deconvolute( TH1D* h, Double_t _Q0, Double_t _s0, Double_t _mu );
  Float_t FindMu( TH1D* h, Double_t _Q0, Double_t _s0 );
  TH1D* RunSingle( TH1D* h, Double_t _Q0, Double_t _s0 );
  TH1D* Run( TH1D* h, Int_t ntoys );    
    
  ClassDef( Deconvolution, 1 )
        
};

#endif
