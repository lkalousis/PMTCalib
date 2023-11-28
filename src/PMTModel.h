
#ifndef PMTMODEL_H
#define PMTMODEL_H

#include "TObject.h"
#include "TGraph.h"
#include "TF1.h"

#include "PMType.h"

class PMTModel : public TObject
{
 private:

  Int_t nbins;
  
  Double_t xmin;
  Double_t xmax;

  Double_t step;
          
 public:
  
  PMTModel();
  
  virtual ~PMTModel();

  PMTModel( Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Model _modtype );

  PMType::Model modtype;
  
  Int_t nparams;
  Double_t params[20]={-1.0};
  
  Double_t wbin;
  
  void SetParams( Double_t _params[] );
  Double_t GetValue( Double_t xx );
  
  Double_t F1( Double_t xx ); // SIMPLE GAUSS 1
  Double_t F2( Double_t xx ); // SIMPLE GAUSS 2
  Double_t F3( Double_t xx ); // ANALYTICAL GAUSS 2
  Double_t F4( Double_t xx ); // EXPLICIT GAUSS 2
  
  TGraph* GetGraph();
  TGraph* GetGraphN( Int_t n );
  TGraph* GetGraphN2( Int_t n );
  
  ClassDef( PMTModel, 1 )
    
};

#endif
