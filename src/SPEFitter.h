
#ifndef SPEFITTER_H
#define SPEFITTER_H

#include "TObject.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"

#include "DFTmethod.h"
#include "PMTModel.h"

#include "TMinuit.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/FunctionMinimum.h" 
#include "Minuit2/MnMigrad.h" 
#include "Minuit2/MnUserParameters.h" 
#include "Minuit2/MnPrint.h" 
#include "Minuit2/FCNBase.h" 
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


class SPEFitter : public TObject
{
 private:

  DFTmethod dft;
  //ROOT::Minuit2::Minuit2Minimizer *mFFT;
    
  PMTModel mod;
  ROOT::Minuit2::Minuit2Minimizer *mMOD;
  
  TF1 *ped_func;
  
 public:
  
  SPEFitter();
  
  virtual ~SPEFitter();

  ROOT::Minuit2::Minuit2Minimizer *mFFT;

  
  Int_t fit_status;
  
  Double_t vals[20];
  Double_t errs[20];

  Double_t ndof;
  Double_t chi2r;

  Double_t FindMu( TH1 *hspec, Double_t _Q0, Double_t _s0 );
  Double_t FindG( TH1 *hspec, Double_t _Q0, Double_t _mu );
  
  void SetDFTmethod( DFTmethod _dft );
  void FitwDFTmethod( TH1 *hspec );
  
  void SetPMTModel( PMTModel _mod );
  void FitwPMTModel( TH1 *hspec );
  
  ClassDef( SPEFitter, 1 )
    
};

#endif
