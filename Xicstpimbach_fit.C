#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include <iostream>
#include <TCanvas.h>
#include <TLine.h>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooVoigtian.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit ;

using std::cout;
using std::endl;
using std::ofstream;

void Xicstpimbach_fit_matt()
{

  TString afs_home = "/Users/needham/andrew/";

  TChain treeXipi("DecayTree");
  TChain tree_lumi("GetIntegratedLuminosity/LumiTuple");

  //TString Xic_cuts("p_ProbNNp>0.8 & K1_ProbNNk>0.8 & K2_ProbNNk>0.8 & pip_ProbNNpi>0.8 & Xic_IPCHI2_OWNPV<16 & Xic_ENDVERTEX_CHI2/Xic_ENDVERTEX_NDOF<3 & Xic_M>2460 & Xic_M<2480 & Omegacst_PT>5000");

  //TString pibach_cuts = "pipbach_ProbNNpi>0.7 & pipbach_IPCHI2_OWNPV < 9 & Omegacst_CHI2NDOF_DTF_Xic_PV<3 & Omegacst_CHI2NDOF_DTFXic_PV > 0 & pipbach_PT>500";

  // TString pimbach_cuts = "pimbach_ProbNNpi>0.5 ";
   TString pimbach_cuts = "pimbach_ProbNNpi>0.7 & pimbach_IPCHI2_OWNPV < 9 & pimbach_PT>500";
  TString Xic0pip_cuts = "Omegacst_M12_DTF_Xic_PV>2642. & Omegacst_M12_DTF_Xic_PV<2648.7";
  //  TString Xic0pip_cuts = "Omegacst_M12_DTF_Xic_PV>2640. & Omegacst_M12_DTF_Xic_PV<2650";

  treeXipi.Add(afs_home + "Xicst_Xic0pipi_ALLCUT.root");
  tree_lumi.Add(afs_home + "Xicst_Xic0pipi_ALL.root");

  // -----------------------Fitting ----------------------------------//

  RooRealVar M("M", "M", 2785, 3085); // Mass range for Xicstpim plot&fit



  // Building threshold function for bkg
  RooRealVar m0("m0", "m(Xicst)+m(pim)", 2785);
  RooRealVar alpha("alpha", "alpha",0.6,  0.1, 1.);
  RooRealVar gamma("gamma", "gamma",-0.003, -0.01, 0.);
  RooGenericPdf threshold("threshold", "threshold", "pow((M-m0),alpha)*exp(gamma*(M-m0))", RooArgList(M, m0, alpha, gamma));
  m0.setConstant(kTRUE);

  RooRealVar a("a","a", 0.1, 0.0001, 1);


  RooRealVar Xicstpim1_M("Xicstpim1_M", "Xicstpim1_M", 2820,2815,2825);
  RooRealVar Xicstpim1_width("Xicstpim1_sigma","Xicstpim1_sigma",2.54);
  //  RooRealVar Xicstpim1_sigma("Xicstpim1_width","Xicstpim1_width",3.,1.,5.);
  RooFormulaVar Xicstpim1_sigma("Xicstpim1_sigma", "Xicstpim1_sigma", "a*sqrt(Xicstpim1_M-m0)", RooArgList(a,Xicstpim1_M, m0));


  RooVoigtian voigt1("voigt1", "voigt1", M, Xicstpim1_M, Xicstpim1_width, Xicstpim1_sigma);
  //Xicstpim1_width.setConstant(kTRUE);

  RooRealVar Xicstpim2_M("Xicstpim2_M", "Xicstpim2_M", 2969.4, 2955, 2985);
  RooRealVar Xicstpim2_width("Xicstpim2_width", "Xicstpim2_width", 28.1);
  RooFormulaVar Xicstpim2_sigma("Xicstpim2_sigma", "Xicstpim1_sigma", "a*sqrt(Xicstpim2_M-m0)", RooArgList(a,Xicstpim2_M, m0));
  RooVoigtian voigt2("voigt2", "voigt2", M, Xicstpim2_M, Xicstpim2_width, Xicstpim2_sigma);

  RooRealVar Xicstpim3_M("Xicstpim3_M", "Xicstpim3_M", 2923.04, 2915, 2929); Xicstpim3_M.setConstant(true);
  RooRealVar Xicstpim3_width("Xicstpim3_width", "Xicstpim3_width", 7.1, 0., 20.); Xicstpim3_width.setConstant(true);
  RooFormulaVar Xicstpim3_sigma("Xicstpim3_sigma", "Xicstpim3_sigma", "a*sqrt(Xicstpim3_M-m0)", RooArgList(a,Xicstpim3_M, m0));
  RooVoigtian voigt3("voigt3", "voigt3", M, Xicstpim3_M, Xicstpim3_width, Xicstpim1_sigma);


  // Add the components
  RooRealVar voigt1frac("voigt1frac", "voigt1frac",0.5, 0., 1.);
  RooRealVar voigt2frac("voigt2frac", "voigt2frac",0.2, 0., 1.);
  RooRealVar voigt3frac("voigt3frac", "voigt3frac",0.0, 0., 1.); // voigt3frac.setConstant(true);
  RooAddPdf sigXicstpim("sigXicstpim", "voigt1+voigt2+voigt3+threshold", RooArgList(voigt1, voigt2, voigt3, threshold), RooArgList(voigt1frac, voigt2frac, voigt3frac));

  //---------------------------Sideband background ---------------------------------------//

  TString lowerXicSideband_cuts = "Omegacst_M12_DTF_Xic_PV<2630 & Omegacst_M12_DTF_Xic_PV>2610";
  TString upperXicSideband_cuts = "Omegacst_M12_DTF_Xic_PV>2700 & Omegacst_M12_DTF_Xic_PV<2740";
  TString middleXicSideband_cuts = "Omegacst_M12_DTF_Xic_PV>2655 & Omegacst_M12_DTF_Xic_PV<2665";

  RooRealVar alpha2("alpha2", "alpha2", 1.,0., 1.5);
  RooRealVar beta2("beta2", "beta2", -7, -10, 10);
  //RooRealVar gamma2("gamma2", "gamma2",-4e-03, -10, 10);
  RooRealVar gamma2("gamma2", "gamma2",-4e-03,-1, 1);
  //RooRealVar alpha2("alpha2", "alpha2",0.001 , -0.01, 0.01);
  //RooRealVar gamma2("gamma2", "gamma2",0.2, 0, 1);

  RooGenericPdf sideband_threshold("sideband threshold", "sideband threshold", "pow((M-m0),alpha2)*exp(-beta2+(gamma2*(M-m0)))", RooArgList(M, m0, alpha2, beta2, gamma2));
 //---------------------------------------Plotting on TCanvas --------------------------------------//

  TLine *line = new TLine();
  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(0);
  // -----------------------------------------Plot Xicstpimbach_M ----------------------------------//

  TH1F* h_Xicstpim = new TH1F("h_Xicstpim", "m(#Xi^{*+}_{c}#pi^{-})", 50, 2785, 3085);
  TH1F* h_XicstpimBKG = new TH1F("h_XicstpimBKG","m(#Xi^{*+}_{c}#pi^{-})", 100, 2785, 3085);

  treeXipi.Draw("Omegacst_M_DTF_Xic_PV >> h_Xicstpim", pimbach_cuts +"&"+ Xic0pip_cuts);
  // treeXipi.Draw("Omegacst_M_DTF_Xic_PV >> h_Xicstpim", pimbach_cuts +"&"+ Xic0pip_cuts);
  //treeXipi.Draw("Omegacst_M_DTF_Xic_PV >> h_XicstpimBKG",pimbach_cuts +"&"+ (upperXicSideband_cuts +"|"+ lowerXicSideband_cuts + "|" + middleXicSideband_cuts));
  treeXipi.Draw("Omegacst_M_DTF_Xic_PV >> h_XicstpimBKG",pimbach_cuts +"&"+ upperXicSideband_cuts +"|"+ lowerXicSideband_cuts);

  RooDataHist Xicstpim_data("Xicstpim_data", "Xicstpim_data", M, h_Xicstpim);
  RooFitResult* res = sigXicstpim.fitTo(Xicstpim_data, Save());

  std::cout << "likelihood " << res->minNll() << std::endl;

  // RooDataHist sideband_Xicstpim_data("sideband_Xicstpim_data", "sideband_Xicstpim_data", M, h_XicstpimBKG);
  //sideband_threshold.chi2FitTo(sideband_Xicstpim_data);

  RooPlot* xframeXicstpim = M.frame(Title("Fitted #Xi^{*+}_{c}#pi^{-} TEST; m(#Xi_{c}^{*+}#pi^{-}) [MeV]"));

  Xicstpim_data.plotOn(xframeXicstpim);
  sigXicstpim.plotOn(xframeXicstpim);
  // sideband_threshold.plotOn(xframeXicstpim, LineStyle(kDashed));

  c1->SetWindowSize(1400., 800.);
  xframeXicstpim->Draw();

  c1->SaveAs("Xicstpimbach_TEST.pdf");

  // ---------------------------Get integrated luminosity -------------------------//

  Double_t  IntegratedLuminosity;
  tree_lumi.SetBranchAddress("IntegratedLuminosity",&IntegratedLuminosity);
  Double_t TotintegLumi = 0.;
  Int_t nevent = tree_lumi.GetEntries();
  for (Int_t i=0;i<nevent;i++)
  {
    tree_lumi.GetEvent(i);
    //read complete accepted event in memory
    TotintegLumi = TotintegLumi+IntegratedLuminosity;
  }
  cout<<"Integrated Luminosity = "<<TotintegLumi<<endl;

}
