
#include "Deconvolution.h"

using namespace std;


ClassImp( Deconvolution )

Deconvolution::Deconvolution()
{
  h1 = NULL;
  h2= NULL;
  
  h3 = NULL;
  h4 = NULL;

  h5 = NULL;
  
}

Deconvolution::~Deconvolution()
{
  delete h1;
  delete h2;

  delete h3;
  delete h4;
  
  delete h5;
  
}

Deconvolution::Deconvolution( Double_t _Q0, Double_t _dQ0, Double_t _s0, Double_t _ds0 )
{    
  Q0 = _Q0;
  dQ0 = _dQ0;
  
  s0 = _s0;
  ds0 = _ds0;

  h1 = NULL;
  h2= NULL;

  h3 = NULL;
  h4 = NULL;

  h5 = NULL;
  
}

void Deconvolution::CleanUps( TH1D *h )
{
  h1 = (TH1D*)h->Clone( "h1" );
    
  Int_t nbins = h1->GetXaxis()->GetNbins();
  
  for ( Int_t i=1; i<=nbins; i++ )
    {
      Double_t xx = h1->GetXaxis()->GetBinCenter(i);
      if ( xx<Q0-5.25*s0 ) h1->SetBinContent(i,0.0);
      else break;
      
    }
  
  Int_t n1 = -1;
  for ( Int_t i=1; i<=nbins; i++ )
    {
      Double_t yy = h1->GetBinContent(i);
      if ( yy!=0.0 ) { n1 = i; break; }
      
    }
 
  Int_t n2 = -1;
  for ( Int_t i=nbins; i>=1; i-- )
    {
      Double_t yy = h1->GetBinContent(i);
      if ( yy!=0.0 ) { n2 = i; break; }
      
    }
    
  N = n2-n1+1; 
  x1 = h1->GetXaxis()->GetBinCenter(n1);
  wbin = h1->GetBinWidth(1);
  x2 = -wbin/2.0+N*wbin;  
  
  h2 = new TH1D( "h2", "SPE Project; Charge in A.U.; Entries", N, -wbin/2.0, x2 );
  h3 = new TH1D( "h3", "SPE project; Charge in A.U.; Entries", N, -wbin/2.0, x2 );
  h4 = new TH1D( "h4", "SPE project; Charge in A.U.; Entries", N, -wbin/2.0, x2 );
  h5 = new TH1D( "h5", "SPE project; Charge in A.U.; Entries", N, -wbin/2.0, x2 );
  
  for ( Int_t i=1; i<=nbins; i++ )
    {
      Double_t xx = h1->GetXaxis()->GetBinCenter(i)-x1;
      Double_t yy = h1->GetBinContent(i);
      h2->Fill( xx, yy );
      
    }
     
  return;

}

Double_t Deconvolution::fftPhase( Double_t vy, Double_t vz )
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
  if ( thetayz>pi ) thetayz -= 2.0*pi;
  
  return thetayz;

}

void Deconvolution::Deconvolute( Double_t _Q0, Double_t _s0 )
{
  _Q0 = _Q0 - x1;
  
  Int_t M = N/2+1;
  Double_t a = N*wbin;

  fftw_plan FWfftR;
  Double_t wfinR[N];
  fftw_complex wfoutR[M];
  
  for ( Int_t i=0; i<N; i++ ) wfinR[i] = h3->GetBinContent(i+1);
  
  FWfftR = fftw_plan_dft_r2c_1d( N, wfinR, wfoutR, FFTW_ESTIMATE );
  fftw_execute( FWfftR );
  fftw_destroy_plan( FWfftR );
  
  Double_t cut1 = 1.0/_s0;
  Double_t cut2 = 1.0/_s0;
  
  Int_t n_t1 = 0;
  Double_t xx1[M];
  Double_t ReSID[M];

  Int_t n_t2 = 0;
  Double_t xx2[M];
  Double_t ImSID[M];

  Double_t C1;
  Double_t D;
  
  for ( Int_t i=0; i<M; i++ )
    {
      Double_t Re = wfoutR[i][0];
      Double_t Im = wfoutR[i][1];
      
      Double_t amp = sqrt( pow( Re, 2.0 )+pow( Im, 2.0 ) );
      Double_t phiR = fftPhase( Im, Re );
      
      Double_t k = i*2.0*TMath::Pi()/a;
            
      Double_t phi = phiR + _Q0*k;
      phi = fftPhase( TMath::Sin(phi), TMath::Cos(phi) );
      if ( phi>TMath::Pi() ) phi -= 2.0*TMath::Pi();

      if ( k<=cut1 )
	{
	  xx1[n_t1] = k;
	  ReSID[n_t1] = TMath::Exp( 0.5*pow( _s0*k, 2.0 ) )*amp*TMath::Cos( phi );
	  C1 = ReSID[n_t1];
	  n_t1++;
	  
	}
      else
	{
	  xx1[n_t1] = k;
	  ReSID[n_t1] = C1; //TMath::Exp( -_mu );
	  n_t1++;

	}
      
      if ( k<=cut2 )
	{
	  xx2[n_t2] = k;
	  ImSID[n_t2] = TMath::Exp( 0.5*pow( _s0*k, 2.0 ) )*amp*TMath::Sin( phi );
	  D = ImSID[n_t2] * k;
	  n_t2++;
	  
	}
      else
	{
	  xx2[n_t2] = k;
	  ImSID[n_t2] = D/k;
	  n_t2++;
	  
	}
            
    }

  grRe = new TGraph( n_t1, xx1, ReSID );
  grIm = new TGraph( n_t2, xx2, ImSID );

  mu_est = -TMath::Log( C1 );


  fftw_complex wfout1[M];

  for ( Int_t i=0; i<M; i++ )
    {
      Double_t k = i*2.0*TMath::Pi()/a;
      
      wfout1[i][0] = grRe->Eval(k);
      wfout1[i][1] = grIm->Eval(k);

    }

  fftw_complex wfout[M];
  Double_t fft[N];

  for ( Int_t i=0; i<M; i++ )
    {
      Double_t Re = wfout1[i][0];
      Double_t Im = wfout1[i][1];

      Double_t amp = sqrt( pow( Re, 2.0 )+pow( Im, 2.0 ) );
      Double_t phi = fftPhase( Im, Re );
      if ( phi>TMath::Pi() ) phi -= 2.0*TMath::Pi();
      
      wfout[i][0] = mu_est + TMath::Log( amp );
      wfout[i][1] = phi;
                        
    }

  fftw_plan BW;
  BW = fftw_plan_dft_c2r_1d( N, wfout, fft, FFTW_ESTIMATE );
  fftw_execute( BW );
  fftw_destroy_plan( BW );

  h4->Reset();
   
  for ( Int_t i=0; i<N; i++ )
    {
      Double_t xx = h4->GetBinCenter( i+1 );
      Double_t yy = fft[i]/(1.0*N*1.0);
      
      if ( yy>1.0e-10 ) h4->Fill( xx, yy );
      else h4->Fill( xx, 0.0 );
      
    }
  
  h4->SetBinContent( 1, 2.0*h4->GetBinContent(1) );
  h4->SetBinError( 1, 2.0*h4->GetBinError(1) );
  
  return;
    
}

TH1D* Deconvolution::Run( TH1D* h, Int_t ntoys )
{
  cout << " Deconvolution starts ... " << endl;
  cout << " " << endl;
  
  CleanUps( h );
  gRandom->SetSeed( 0 );
  
  Double_t xx[6000] = { 0.0 };
  Double_t yy[6000] = { 0.0 };
  Double_t eyy[6000] = { 0.0 };
    
  for ( Int_t i=0; i<N; i++ ) xx[i] = h2->GetXaxis()->GetBinCenter(i+1);
      
  for ( Int_t i=0; i<ntoys; i++ )
    {
      Double_t Qprime = gRandom->Gaus( Q0, dQ0 );
      Double_t sprime = gRandom->Gaus( s0, ds0 );
      
      //h3->Reset();
      for ( Int_t j=0; j<N; j++ )
	{
	  Double_t _n = h2->GetBinContent( j+1 );
	  Double_t _nprime = gRandom->Poisson( _n );
	  	  
	  h3->SetBinContent( j+1, _nprime );
	  	  
	}

      Double_t sum = h3->Integral();
      h3->Scale( 1.0/sum );
      
      Deconvolute( Qprime, sprime );
      for ( Int_t j=0; j<N; j++ )  yy[j] += h4->GetBinContent(j+1);
      for ( Int_t j=0; j<N; j++ ) eyy[j] += h4->GetBinContent(j+1)*h4->GetBinContent(j+1);
      
      cout << "\r Progress : " << Form( "%.1f", i*1.0/ntoys*100.0 ) << " %";
      std::cout.flush();
      
    }

  cout << "\r Progress 2 : " << 100 << " %";
  std::cout.flush();
  cout << "" << endl;
  
  for ( Int_t j=0; j<N; j++ ) yy[j]  = yy[j]/(1.0*ntoys);
  for ( Int_t j=0; j<N; j++ ) eyy[j] = sqrt( eyy[j]/(1.0*ntoys) - yy[j]*yy[j] );
  
  for ( Int_t j=0; j<N; j++ )
    {
      h5->SetBinContent( j+1, yy[j] );
      h5->SetBinError( j+1, eyy[j] );
      
    }
  
  h5->SetMinimum( 0.0 );
  h5->SetMarkerStyle( 20 );
  h5->SetMarkerSize( 0.75 );
  h5->SetLineColor( kBlack );
  h5->SetMarkerColor( kBlack );
  
  cout << " " << endl;
  cout << " ... and now it ends !" << endl;
  cout << " " << endl;
  
  return h5;
    
}

