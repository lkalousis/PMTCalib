
#include "Occupancy.h"

using namespace std;


ClassImp( Occupancy )

Occupancy::Occupancy()
{}

Occupancy::~Occupancy()
{}

Occupancy::Occupancy( Double_t _Q0, Double_t _s0 )
{    
  Q0 = _Q0;
  s0 = _s0;
        
}

Double_t Occupancy::Gauss1( Float_t x )
{
  Double_t result = 0.0;
    
  Double_t arg = 0.0; 
  if ( s0!=0.0 ) arg = ( x - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero." << endl;
  result = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg*arg );

  return result;

}

Float_t Occupancy::FindG( TH1D* _h, Float_t f )
{
  Int_t n1 = 0;
  Float_t x1[6000];
  Float_t y1[6000];
  
  x1[n1] = _h->GetXaxis()->GetBinLowEdge(1);
  y1[n1] = 0.0;
  n1++;
  
  Double_t sum1 = 0.0;
  Int_t nbins = _h->GetXaxis()->GetNbins();
  Float_t wbin = _h->GetXaxis()->GetBinWidth(1);

  for ( Int_t i=1; i<=nbins; i++ )
    {
      sum1 += Gauss1( _h->GetXaxis()->GetBinCenter(i) )*wbin;
      x1[n1] = _h->GetXaxis()->GetBinUpEdge(i);
      y1[n1] = sum1;
      n1++;
      
    }
  
  TGraph *gr1 = new TGraph( n1, y1, x1 );
  Float_t cut = gr1->Eval( f );

  Int_t n2 = 0;
  Float_t x2[6000];
  Float_t y2[6000];

  x2[n2] = _h->GetXaxis()->GetBinLowEdge(1);
  y2[n2] = 0.0;
  n2++;

  Double_t sum2 = 0.0;
  
  for ( Int_t i=1; i<=nbins; i++ )
    {
      sum2 += _h->GetBinContent(i);
      x2[n2] = _h->GetXaxis()->GetBinUpEdge(i);
      y2[n2] = sum2;
      n2++;
      
    }

  TGraph *gr2 = new TGraph( n2, x2, y2 );
  Float_t Ncut = gr2->Eval( cut );
  Float_t N0 = Ncut/f;
  Double_t Ntot = _h->Integral();  
  Float_t mu = -TMath::Log( N0/Ntot );;

  delete gr1;
  delete gr2;

  Double_t G = ( _h->GetMean()-Q0 )/mu;

  return G;
  
}
