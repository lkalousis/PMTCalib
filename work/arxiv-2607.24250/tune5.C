
void tune5()
{
  gROOT->Reset();
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetCanvasColor( 0 );
  gStyle->SetPadColor( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetFrameBorderMode( 0 );
  gStyle->SetTitleColor( 0 );   
  gStyle->SetTitleFillColor( 0 );  
  gStyle->SetTitleBorderSize( 0 );
  gStyle->SetTitleX( 0.10 );
  gStyle->SetTitleY( 0.98 );
  gStyle->SetTitleFont( 22, "" );
  gStyle->SetTitleSize( 0.055, ""  );
  gStyle->SetStatColor( 0 );
  gStyle->SetStatFont( 22 );
  gStyle->SetStatBorderSize( 1 );
  gStyle->SetStatX( 0.90 );
  gStyle->SetStatY( 0.90 );
  gStyle->SetStatFontSize( 0.04 );
  gStyle->SetOptStat( 1110 );
  gStyle->SetTitleFont( 22, "XYZ"  );
  gStyle->SetTitleSize( 0.05, "XYZ"  );
  gStyle->SetTitleColor( kBlack, "XYZ"  );
  gStyle->SetTitleAlign(13);
  gStyle->SetLabelFont( 22, "XYZ"  );
  gStyle->SetLabelSize( 0.04, "XYZ"  );
  gStyle->SetOptStat( 1110 ); 
  gStyle->SetOptFit( 0 ); 
  
  gStyle->SetPalette( 1 ); // Standard
  //gStyle->SetPalette( 51, 0 ); // DeepSea
  //gStyle->SetPalette( 53, 0 ); // Black body
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
  gStyle->SetNumberContours( NCont );

  gROOT->ForceStyle();

  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetGridx(1);
  c1->SetGridy(1);

  
  Double_t x[6] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };
  
  Double_t y1[6];
  Double_t ey1[6];

  Double_t y2[6];
  Double_t ey2[6];
    
  for ( Int_t i=0; i<6; i++ )
    {
      TFile *f1 = new TFile( Form( "out/test-occ2-%0.2f.root", x[i] ), "read" );
      TH1D *hG1 = (TH1D*)f1->Get( "hG" );
      
      y1[i] = hG1->GetRMS();
      ey1[i] = hG1->GetRMSError();

      TFile *f2 = new TFile( Form( "out/test-fit2-%0.2f.root", x[i] ), "read" );
      TH1D *hG2 = (TH1D*)f2->Get( "hG" );

      y2[i] = hG2->GetRMS();
      ey2[i] = hG2->GetRMSError();
                  
    }


  TGraphErrors *gr1 = new TGraphErrors( 6, x, y1, 0, ey1 );
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.75);
  gr1->SetLineWidth(2.0);
  gr1->SetTitle( "; #mu; Resolution of #Delta G  in % " );
  gr1->GetXaxis()->SetRangeUser( 0.0, 3.5 );
  gr1->SetMinimum(0.0);
  gr1->SetMaximum(8.0);
  gr1->Draw( "APEZL" );

  TGraphErrors *gr2 = new TGraphErrors( 6, x, y2, 0, ey2 );
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.75);
  gr2->SetLineWidth(2.0);
  gr2->SetMarkerColor(50);
  gr2->SetLineColor(50);
  gr2->Draw( "PEZL,SAME" );
    
  TLegend *leg = new TLegend( 0.55, 0.67, 0.85, 0.85 );
  leg->SetTextFont( 22 );
  //leg->SetBorderSize( 0 );
  leg->AddEntry( gr1, "Occupancy", "PL" );
  leg->AddEntry( gr2, "Fit", "PL" );
  leg->SetFillColor( 0 ); 
  leg->Draw();

  
}
