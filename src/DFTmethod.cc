
#include "DFTmethod.h"

using namespace std;


ClassImp( DFTmethod )

DFTmethod::DFTmethod()
{}

DFTmethod::~DFTmethod()
{}

DFTmethod::DFTmethod( Int_t _nbins, Double_t _xmin, Double_t _xmax, SPEResponse _spef )
{
  nbins = _nbins;

  xmin = _xmin;
  xmax = _xmax;

  step = ( xmax-xmin )/( 1.0*nbins*1.0 );
  
  spef = _spef;
    
  N = 2*nbins+60; 
  M = N/2+1;

  
  xvalues.clear();
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xmin + 1.0*i*step;
      
      xvalues.push_back( xx );
      
    }

  edge = xmin;
  
    
}

void DFTmethod::CalculateValues()
{
  fftw_plan FWfftBG;
  fftw_plan FWfftSG;
  
  Double_t wfinBG[N]; fftw_complex wfoutBG[M];
  Double_t wfinSG[N]; fftw_complex wfoutSG[M];
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xvalues.at( i ) - edge;
      
      Double_t arg = 0.0; 
      if ( s0!=0.0 ) arg = ( xx - Q0 + edge )/s0;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t yy = 1.0/( sqrt( 2.0 * TMath::Pi() ) * s0 ) * TMath::Exp( -0.5*arg*arg );
      wfinBG[i] = yy;
      
      wfinSG[i] = spef.GetValue( xx );

    }

  FWfftBG = fftw_plan_dft_r2c_1d( N, wfinBG, wfoutBG, FFTW_ESTIMATE );
  fftw_execute( FWfftBG );
  fftw_destroy_plan( FWfftBG );

  FWfftSG = fftw_plan_dft_r2c_1d( N, wfinSG, wfoutSG, FFTW_ESTIMATE );
  fftw_execute( FWfftSG );
  fftw_destroy_plan( FWfftSG );


  fftw_complex wfout[M];
  Double_t fftout[N];
  
  for ( UInt_t i=0; i<M; i++ )
    {
      Double_t amp_BG = sqrt( pow( wfoutBG[i][0], 2.0 )+pow( wfoutBG[i][1], 2.0 ) );
      Double_t ph_BG = fftPhase( wfoutBG[i][1], wfoutBG[i][0] );
      
      Double_t ReS = wfoutSG[i][0];
      Double_t ImS = wfoutSG[i][1];
      
      double ph = ( ph_BG + mu*ImS*step );
      
      wfout[i][0] = amp_BG*TMath::Exp( mu*ReS*step ) * TMath::Cos( ph );
      wfout[i][1] = amp_BG*TMath::Exp( mu*ReS*step ) * TMath::Sin( ph );
      
    }
  
  fftw_plan BWfft;
  BWfft = fftw_plan_dft_c2r_1d( N, wfout, fftout, FFTW_ESTIMATE );
  fftw_execute( BWfft );
  fftw_destroy_plan( BWfft );
  
  yvalues.clear();
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t yy = Norm * wbin * TMath::Exp( -1.0*mu ) * fftout[i]/( 1.0*N*1.0 );
      yvalues.push_back( yy );
      
    }
  
  
  Double_t x[nbins];
  Double_t y[nbins];

  for ( Int_t i=0; i<nbins; i++ )
    {
      x[i] = xvalues.at( i );
      y[i] = yvalues.at( i );
      
    }

  gr = new TGraph( nbins, x, y );
  
  return;
  
}

Double_t DFTmethod::GetValue( Double_t xx )
{
  Double_t y_ = gr->Eval( xx ); 

  return y_;

}

TGraph* DFTmethod::GetGraph()
{
  CalculateValues();
  
  Double_t x[nbins];
  Double_t y[nbins];
  
  for ( Int_t i=0; i<nbins; i++ )
    {
      x[i] = xvalues.at( i );

      Double_t y_ = GetValue( x[i] );
      
      if ( y_<1.0e-10 ) y[i] = 1.e-4;
      else y[i] = y_;
      
    }
  
  TGraph *_gr = new TGraph( nbins, x, y );
  
  _gr->SetLineWidth( 2 );
  int cc = kAzure+1;
  _gr->SetLineColor( cc );
  _gr->SetMarkerColor( cc );
  _gr->SetMarkerSize( 0.1 );
    
  return _gr;
  
}

TGraph* DFTmethod::GetGraphN( Int_t n )
{
  CalculateValues();
  
  fftw_plan FWfftBG;
  fftw_plan FWfftSG;
  
  Double_t wfinBG[N]; fftw_complex wfoutBG[M];
  Double_t wfinSG[N]; fftw_complex wfoutSG[M];
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = xvalues.at( i ) - edge;
      
      Double_t arg = 0.0; 
      if ( s0!=0.0 ) arg = ( xx - Q0 + edge )/s0;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t yy = 1.0/( sqrt( 2.0 * TMath::Pi() ) * s0 ) * TMath::Exp( -0.5*arg*arg );
      wfinBG[i] = yy;
      
      wfinSG[i] = spef.GetValue( xx );

    }

  FWfftBG = fftw_plan_dft_r2c_1d( N, wfinBG, wfoutBG, FFTW_ESTIMATE );
  fftw_execute( FWfftBG );
  fftw_destroy_plan( FWfftBG );

  FWfftSG = fftw_plan_dft_r2c_1d( N, wfinSG, wfoutSG, FFTW_ESTIMATE );
  fftw_execute( FWfftSG );
  fftw_destroy_plan( FWfftSG );


  fftw_complex wfout[M];
  Double_t fftout[N];
  
  for ( UInt_t i=0; i<M; i++ )
    {
      Double_t amp_BG = sqrt( pow( wfoutBG[i][0], 2.0 )+pow( wfoutBG[i][1], 2.0 ) );
      Double_t ph_BG = fftPhase( wfoutBG[i][1], wfoutBG[i][0] );

      Double_t amp_SG = sqrt( pow( wfoutSG[i][0], 2.0 )+pow( wfoutSG[i][1], 2.0 ) );
      Double_t ph_SG = fftPhase( wfoutSG[i][1], wfoutSG[i][0] );
            
      double ph = ( ph_BG + 1.0*n*ph_SG*step );
      
      wfout[i][0] = amp_BG*pow( amp_SG, 1.0*n*1.0 )*TMath::Cos( ph );
      wfout[i][1] = amp_BG*pow( amp_SG, 1.0*n*1.0 )*TMath::Sin( ph );
      
    }
  
  fftw_plan BWfft;
  BWfft = fftw_plan_dft_c2r_1d( N, wfout, fftout, FFTW_ESTIMATE );
  fftw_execute( BWfft );
  fftw_destroy_plan( BWfft );
  
  
  Double_t x[nbins];
  Double_t y[nbins];
  
  for ( Int_t i=0; i<nbins; i++ )
    {
      x[i] = xvalues.at( i );
      Double_t y_ = Norm * wbin * TMath::Exp( -1.0*mu )/TMath::Factorial( n ) * pow( mu, 1.0*n*1.0 ) * fftout[i]/( 1.0*N*1.0 );

      if ( y_<1.0e-10 ) y[i] = 1.e-4;
      else y[i] = y_;
            
    }

  TGraph *_gr = new TGraph( nbins, x, y );

  _gr->SetLineWidth( 2 );
  _gr->SetLineStyle( 3 );
  _gr->SetLineColor( kBlack );
  _gr->SetMarkerColor( kBlack );
  _gr->SetMarkerSize( 0.1 );
  
  return _gr;
  
}

Double_t DFTmethod::fftPhase( Double_t vy, Double_t vz )
{
  Double_t thetayz = -999.0;
  
  Double_t pi = TMath::Pi();
  
  
  if ( vz>0 && vy>0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); }

  else if ( vz<0 && vy>0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); thetayz=pi-thetayz; }

  else if ( vz<0 && vy<0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); thetayz=thetayz+pi; }

  else if ( vz>0 && vy<0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); thetayz=2.0*pi-thetayz; }

  else if ( vz==0 && vy>0 ) { thetayz=pi/2.0; }

  else if ( vz==0 && vy<0 ) { thetayz=3.0*pi/2.0; }

  else if ( vz>0 && vy==0 ) { thetayz=0.0; }

  else if ( vz<0 && vy==0 ) { thetayz=pi; }
  
  thetayz = fmod( thetayz, 2.0*pi );
  
  return thetayz;

}
