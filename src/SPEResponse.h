
#ifndef SPERESPONSE_H
#define SPERESPONSE_H

#include "TObject.h"
#include "TMath.h"
#include "TF1.h"

#include "PMType.h"

class SPEResponse : public TObject
{
 private:

  PMType::Response spetype;
    
      
 public:
  
  SPEResponse();
  
  virtual ~SPEResponse();

  SPEResponse( PMType::Response _spetype, Double_t _params[] );
  
  Int_t nparams;

  void SetParams( Double_t _params[] );
  
  Double_t GetValue( Double_t xx );
    
  Double_t GenQ();

  TF1 *spefunc;

  Double_t params[10]={-1.0};
  
  ClassDef( SPEResponse, 1 )
    
};

#endif
