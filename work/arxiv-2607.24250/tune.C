
void tune()
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

  
  Double_t x[18] = { 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85 };
  
  Double_t y1[18];
  Double_t ey1[18];

  Double_t y2[18];
  Double_t ey2[18];

  Double_t y3[18];
  Double_t ey3[18];
  
  for ( Int_t i=0; i<18; i++ )
    {
      TFile *f1 = new TFile( Form( "out/test-occ1-%0.2f-%0.2f.root", 1.00, x[i] ), "read" );
      TH1D *hG1 = (TH1D*)f1->Get( "hG" );

      y1[i] = hG1->GetMean();
      ey1[i] = hG1->GetMeanError();

      TFile *f2 = new TFile( Form( "out/test-occ1-%0.2f-%0.2f.root", 2.00, x[i] ), "read" );
      TH1D *hG2 = (TH1D*)f2->Get( "hG" );

      y2[i] = hG2->GetMean();
      ey2[i] = hG2->GetMeanError();

      TFile *f3 = new TFile( Form( "out/test-occ1-%0.2f-%0.2f.root", 3.00, x[i] ), "read" );
      TH1D *hG3 = (TH1D*)f3->Get( "hG" );

      y3[i] = hG3->GetMean();
      ey3[i] = hG3->GetMeanError();
            
    }


  TGraphErrors *gr1 = new TGraphErrors( 18, x, y1, 0, ey1 );
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.75);
  gr1->SetLineWidth(2.0);
  gr1->SetTitle( "; Threshold fraction (f); Mean of #Delta G  in % " );
  gr1->GetXaxis()->SetRangeUser( 0.0, 0.86 );
  gr1->SetMinimum(0.0);
  gr1->SetMaximum(10.0);
  gr1->Draw( "APEZL" );

  TGraphErrors *gr2 = new TGraphErrors( 18, x, y2, 0, ey2 );
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.75);
  gr2->SetLineWidth(2.0);
  gr2->SetMarkerColor(50);
  gr2->SetLineColor(50);
  gr2->Draw( "PEZL,SAME" );

  TGraphErrors *gr3 = new TGraphErrors( 18, x, y3, 0, ey3 );
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.75);
  gr3->SetLineWidth(2.0);
  gr3->SetMarkerColor(29);
  gr3->SetLineColor(29);
  gr3->Draw( "PEZL,SAME" );
  
  TLegend *leg = new TLegend( 0.55, 0.67, 0.85, 0.85 );
  leg->SetTextFont( 22 );
  //leg->SetBorderSize( 0 );
  leg->AddEntry( gr1, "#mu = 1", "PL" );
  leg->AddEntry( gr2, "#mu = 2", "PL" );
  leg->AddEntry( gr3, "#mu = 3", "PL" );
  leg->SetFillColor( 0 ); 
  leg->Draw();

  
}
