
#include "Deconvolution.h"

using namespace std;


ClassImp( Deconvolution )

Deconvolution::Deconvolution()
{}

Deconvolution::~Deconvolution()
{}

Deconvolution::Deconvolution( Double_t _Q0, Double_t _dQ0, Double_t _s0, Double_t _ds0, Double_t _mu )
{    
  Q0 = _Q0;
  dQ0 = _dQ0;
  
  s0 = _s0;
  ds0 = _ds0;

  mu = _mu;
    
}

TH1D* Deconvolution::CleanUps( TH1D *h )
{
  TH1D *h1 = (TH1D*)h->Clone( "h1" );
    
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
    
  UInt_t N = n2-n1+1; //cout << N << endl;
  x1 = h1->GetXaxis()->GetBinCenter(n1);
  Double_t wbin = h1->GetBinWidth(1);
  Double_t x2 = -wbin/2.0+N*wbin;  
  TH1D *h2 = new TH1D( "h2", "Project", N, -wbin/2.0, x2 );

  for ( Int_t i=1; i<=nbins; i++ )
    {
      Double_t xx = h1->GetXaxis()->GetBinCenter(i)-x1;
      Double_t yy = h1->GetBinContent(i);

      h2->Fill( xx, yy );
      
    }
  
  Double_t sum = h2->Integral();
  h2->Scale( 1.0/sum );
  //Q0 = Q0 - x1;
  //cout << N << endl;
  
  delete h1;
  return h2;

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

TH1D* Deconvolution::Deconvolute( TH1D* h, Double_t _Q0, Double_t _s0, Double_t _mu )
{
  TH1D *h2 = CleanUps( h );
  UInt_t N = h2->GetXaxis()->GetNbins();
  Double_t wbin = h2->GetBinWidth(1);
  _Q0 = _Q0 - x1;
  
  UInt_t M = N/2+1;
  Double_t a = N*wbin;

  fftw_plan FWfftR;
  Double_t wfinR[N];
  fftw_complex wfoutR[M];
  
  for ( UInt_t i=0; i<N; i++ ) wfinR[i] = h2->GetBinContent(i+1);

  FWfftR = fftw_plan_dft_r2c_1d( N, wfinR, wfoutR, FFTW_ESTIMATE );
  fftw_execute( FWfftR );
  fftw_destroy_plan( FWfftR );

  Double_t cut1 = 1.0/_s0;
  //Double_t cut2 = 1.5/_s0;

  Int_t n_t1 = 0;
  Double_t xx1[M];
  Double_t ReSID[M];

  Int_t n_t2 = 0;
  Double_t xx2[M];
  Double_t ImSID[M];
  //fftw_complex wfout2[M];
  
  for ( UInt_t i=0; i<M; i++ )
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
	  ReSID[n_t1] = TMath::Exp( _mu + 0.5*pow( _s0*k, 2.0 ) )*amp*TMath::Cos( phi );
	  n_t1++;

	}

      if ( k<=cut1 )
	{
	  xx2[n_t2] = k;
	  ImSID[n_t2] = TMath::Exp( _mu + 0.5*pow( _s0*k, 2.0 ) )*amp*TMath::Sin( phi );
	  n_t2++;

	}

      //wfout2[i][0] = TMath::Exp( _mu )*amp*TMath::Cos( phi );
      //wfout2[i][1] = TMath::Exp( _mu )*amp*TMath::Sin( phi );
      
    }

  Double_t k0 = (M-1)*2.0*TMath::Pi()/a;

  xx1[n_t1] = k0;
  ReSID[n_t1] = 1.0;
  n_t1++;

  xx2[n_t2] = k0;
  ImSID[n_t2] = 0.0;
  n_t2++;
    
  TGraph *grRe = new TGraph( n_t1, xx1, ReSID );
  TGraph *grIm = new TGraph( n_t2, xx2, ImSID );
  
  fftw_complex wfout1[M];

  for ( UInt_t i=0; i<M; i++ )
    {
      Double_t k = i*2.0*TMath::Pi()/a;
      
      if (k<=cut1)
      wfout1[i][0] = grRe->Eval(k);
      else wfout1[i][0] = 1.0;

      //if (k<2.0/_s0)
      wfout1[i][1] = grIm->Eval(k);
      //else wfout1[i][1] = 0.0;
      
    }
  /*
  Double_t fftout1[N];
  Double_t fftout2[N];
  
  fftw_plan BWfft1;
  BWfft1 = fftw_plan_dft_c2r_1d( N, wfout1, fftout1, FFTW_ESTIMATE );
  fftw_execute( BWfft1 );
  fftw_destroy_plan( BWfft1 );

  fftw_plan BWfft2;
  BWfft2 = fftw_plan_dft_c2r_1d( N, wfout2, fftout2, FFTW_ESTIMATE );
  fftw_execute( BWfft2 );
  fftw_destroy_plan( BWfft2 );

  Double_t x2 = -wbin/2.0+N*wbin;  
  TH1D *h3 = new TH1D( "h3", "SPE project; Charge; Entries", N, -wbin/2.0, x2 );
  TH1D *h4 = new TH1D( "h4", "SPE project; Charge; Entries", N, -wbin/2.0, x2 );
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = h2->GetBinCenter( i+1 );
      Double_t yy = fftout1[i]/(1.0*N*1.0);
      h3->Fill( xx, yy );

      Double_t yy1 = fftout2[i]/(1.0*N*1.0);
      h4->Fill( xx, yy1 );
                  
    }
  */
  /*
  fftw_plan FWfft6;
  Double_t wfin6[N];
  fftw_complex wfout6[M];
  
  TH1D *h4 = new TH1D( "h4", "SPE project; Charge; Entries", n, x1, x2 );

  for ( Int_t i=0; i<n; i++ )
    {
      Double_t xx = h2->GetBinCenter( i+1 );
      Double_t yy1 = h2->GetBinContent( i+1 );
      Double_t yy2 = h3->GetBinContent( i+1 );

      Double_t yy = 0.0;
      
      if ( i<=n/2 ) //xx<30.0*s0 )
	yy=yy1;
      else if ( x2-5.0*s0<xx ) yy=0.0;
      else yy = yy2;

      h4->Fill( xx, yy );
      wfin6[i] = yy;

    }

  c1->SetLogy(1);
  
  hID->Draw("HIST");
  
  h4->SetLineColor(kBlue);
  h4->Draw("HIST,SAME");
  
  //c1->Update();
  //c1->WaitPrimitive();

  FWfft6 = fftw_plan_dft_r2c_1d( N, wfin6, wfout6, FFTW_ESTIMATE );
  fftw_execute( FWfft6 );
  fftw_destroy_plan( FWfft6 );
  */
  
  fftw_complex wfout[M];
  Double_t fft[N];

  for ( UInt_t i=0; i<M; i++ )
    {
      Double_t Re = wfout1[i][0];
      Double_t Im = wfout1[i][1];

      Double_t amp = sqrt( pow( Re, 2.0 )+pow( Im, 2.0 ) );
      Double_t phi = fftPhase( Im, Re );
      if ( phi>TMath::Pi() ) phi -= 2.0*TMath::Pi();

      wfout[i][0] = TMath::Log( amp );
      wfout[i][1] = phi;
                  
    }

  fftw_plan BW;
  BW = fftw_plan_dft_c2r_1d( N, wfout, fft, FFTW_ESTIMATE );
  fftw_execute( BW );
  fftw_destroy_plan( BW );

  Double_t x2 = -wbin/2.0+N*wbin;  
  TH1D *h5 = new TH1D( "h5", "SPE project; Charge; Entries", N, -wbin/2.0, x2 );
  
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = h5->GetBinCenter( i+1 );
      Double_t yy = fft[i]/(1.0*N*1.0);
      
      if ( yy>1.0e-10 )
      h5->Fill( xx, yy );
      else h5->Fill( xx, 0.0 );
      
    }
  
  delete grRe;
  delete grIm;
    
  delete h2;

  return h5;
  
}

Double_t Deconvolution::GridMu( TH1D *h, Double_t _Q0, Double_t _s0 )
{
  Double_t tsum = 666.0;
      
  Double_t lim1 = 0.7*mu;
  Double_t lim2 = 1.3*mu;

  Float_t step = 0.01;
  Int_t n1 = (lim2-lim1)/step;
    
  for ( Int_t k=0; k<=n1; k++ )
    {
      Double_t mu_k = lim1 + 1.0*k*step;
      //cout << mu_k << endl;
      
      TH1D* h3 = Deconvolute( h, _Q0, _s0, mu_k );
      	  
      Double_t sum1 = h3->Integral();
      Double_t sum2 = mu_k;
      Double_t sum = TMath::Abs( sum1-sum2 )/sum2;
                        
      if ( sum<tsum )
	{
	  mu_bf = sum2;
	  tsum = sum;

	}
            
      delete h3;

    }
  
  return mu_bf;

}

TH1D* Deconvolution::RunSingle( TH1D* h, Double_t _Q0, Double_t _s0 )
{
  Float_t rnd = GridMu( h, _Q0, _s0 );
  TH1D* h3 = Deconvolute( h, _Q0, _s0, rnd );
  //cout << "-> mu " << rnd << endl;
  //cout << " " << endl;
  return h3;

}

TH1D* Deconvolution::Run( TH1D* h, Int_t ntoys )
{
  cout << " Deconvolution starts ... " << endl;
  cout << " " << endl;
  //cout << "[" << flush;
  
  Int_t nbins = h->GetXaxis()->GetNbins();
  TH1D *hnew = (TH1D*)h->Clone( "hnew" );
  
  Double_t xx[6000] = { 0.0 };
  Double_t yy[6000] = { 0.0 };
  Double_t eyy[6000] = { 0.0 };

  TH1D *h4 = Deconvolute( h, Q0, s0, 0.0 );
  Int_t n = h4->GetXaxis()->GetNbins();
  Double_t wbin = h4->GetBinWidth(1);
  
  for ( Int_t i=0; i<n; i++ ) xx[i] = h4->GetXaxis()->GetBinCenter(i+1);
  delete h4;
    
  for ( Int_t i=0; i<ntoys; i++ )
    {
      //cout << " Run : " << i << endl;
      //cout << " " << endl;

      Double_t Qprime = gRandom->Gaus( Q0, dQ0 );
      Double_t sprime = gRandom->Gaus( s0, ds0 );

      hnew->Reset();
      
      for ( Int_t j=0; j<nbins; j++ )
	{
	  Double_t _n = h->GetBinContent( j+1 );
	  Double_t _nprime = gRandom->Poisson( _n );
	  if ( _n!=0 && _nprime==0 ) _nprime = 1.0e-19;
	  hnew->SetBinContent( j+1, _nprime );
	  
	}
      
      TH1D* h4 = RunSingle( hnew, Qprime, sprime );
      for ( Int_t j=0; j<n; j++ )  yy[j] += h4->GetBinContent(j+1);
      for ( Int_t j=0; j<n; j++ ) eyy[j] += h4->GetBinContent(j+1)*h4->GetBinContent(j+1);

      delete h4;

      cout << "\r Progress: " << Form( "%.1f", i*1.0/ntoys*100.0 ) << " %";
      std::cout.flush();
      
    }

  for ( Int_t j=0; j<n; j++ ) yy[j]  = yy[j]/(1.0*ntoys);
  for ( Int_t j=0; j<n; j++ ) eyy[j] = sqrt( eyy[j]/(1.0*ntoys) - yy[j]*yy[j] );

  TH1D *h5 = new TH1D( "h5", "Project", n, -wbin/2.0, -wbin/2+1.0*n*wbin );
  
  for ( Int_t j=0; j<n; j++ )
    {
      h5->SetBinContent( j+1, yy[j] );
      h5->SetBinError( j+1, eyy[j] );
      
    }
  
  delete hnew;
  
  h5->SetMinimum( 0.0 );
  h5->SetMarkerStyle( 20 );
  h5->SetMarkerSize( 0.75 );
  h5->SetLineColor( kBlack );
  h5->SetMarkerColor( kBlack );

  cout << " " << endl;
  cout << " " << endl;
  cout << " ... and now it ends !" << endl;
 
  
  return h5;
    
}

