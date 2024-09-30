#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
// for fake rate
#include "QCDFakeRate.h"
// for prescales
//#include "Run2PhotonTriggerPrescales.h"
#include "PrescaleProvider.h"
#include <set>
#include <limits>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

  analysisClass::~analysisClass(){}

bool CheckTriggerThreshold(const float thresh1, const float thresh2, const float hltPhotonPt, const float offlineElectronPt)
{
      if(hltPhotonPt >= thresh1  && hltPhotonPt < thresh2)
        return true;
      else if(offlineElectronPt >= thresh1 && offlineElectronPt < thresh2)
        return true;
      return false;
}

bool CheckTriggerThresholds(const float thresh1, const float thresh2, const float hltPhotonPt1, const float offlineElectronPt1,
    const float hltPhotonPt2, const float offlineElectronPt2)
{
  return CheckTriggerThreshold(thresh1, thresh2, hltPhotonPt1, offlineElectronPt1) ||
    CheckTriggerThreshold(thresh1, thresh2, hltPhotonPt2, offlineElectronPt2);
}

bool CheckTriggerThreshold(const float thresh1, const float thresh2, const float matchingObjectPt)
{
      if(matchingObjectPt >= thresh1  && matchingObjectPt < thresh2)
        return true;
      return false;
}

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop() begins" <<std::endl;   

  //--------------------------------------------------------------------------
  // Final selection mass points
  //--------------------------------------------------------------------------
  //const int n_lq_mass = 11;
  //int LQ_MASS[n_lq_mass] = {
  //  1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000
  //};

  const int n_lq_mass = 28;
  int LQ_MASS[n_lq_mass] = { 
    300,  400,  500,  600,
    700,  800,  900,  1000,
    1100, 1200, 1300, 1400,
    1500, 1600, 1700, 1800,
    1900, 2000, 2100, 2200,
    2300, 2400, 2500, 2600,
    2700, 2800, 2900, 3000,
    //3500, 4000
  };

  // LQ650 only 2012
  //const int n_lq_mass = 1;
  //int LQ_MASS[n_lq_mass] = { 650 };

  std::vector<bool> passed_vector;

  char cut_name[100];
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------

  fillAllPreviousCuts              ( true  ) ;
  fillAllOtherCuts                 ( true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;

  bool do_roi_plots = false;

  //--------------------------------------------------------------------------
  // Get pre-cut values
  //--------------------------------------------------------------------------
  // eta boundaries

  double eleEta_bar            = getPreCutValue1("eleEta_bar");
  double eleEta_end1_min       = getPreCutValue1("eleEta_end1");
  double eleEta_end1_max       = getPreCutValue2("eleEta_end1");
  double eleEta_end2_min       = getPreCutValue1("eleEta_end2");
  double eleEta_end2_max       = getPreCutValue2("eleEta_end2");

  // prescales
	PrescaleProvider psProv(getPreCutString1("PrescaleProviderInfo"));

  //--------------------------------------------------------------------------
  // Analysis year
  //--------------------------------------------------------------------------
  string analysisYear = getPreCutString1("AnalysisYear");
  int analysisYearInt = -1;
  if(analysisYear.find("2016") != string::npos)
      analysisYearInt = 2016;
  else {
    analysisYearInt = getPreCutValue1("AnalysisYear");
    analysisYear = to_string(analysisYearInt);
  }

  //--------------------------------------------------------------------------
  // electron ID
  //--------------------------------------------------------------------------
  std::string electronIDType     = getPreCutString1("electronIDType");
  if(electronIDType != "HEEP" && electronIDType != "EGMLoose") {
    STDOUT("electronIDType=" << electronIDType << " is unknown! Please implement it in the analysisClass code. Exiting.");
    exit(-5);
  }

  //--------------------------------------------------------------------------
  // QCD Fake Rate loading part
  //--------------------------------------------------------------------------
  std::string qcdFileName = getPreCutString1("QCDFakeRateFilename");
  std::vector<std::string> regionVec;
  if(analysisYearInt < 2018)
    regionVec = {"2Jet_TrkIsoHEEP7vsHLTPt_PAS"};
  else
    regionVec = {
      "2Jet_TrkIsoHEEP7vsHLTPt_pre319077",
      "2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077",
      "2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"
    };
  QCDFakeRate qcdFR(qcdFileName, "fr2D_", regionVec, true); // look up from hist

  std::string qcdCalcType  = getPreCutString1("QCDCalcType");
  bool doSingleFR = false;
  if(qcdCalcType == "single")
    doSingleFR = true;
  else if(qcdCalcType != "double")
    throw runtime_error("Could not understand value in cut file for precut QCDCalcType: " + qcdCalcType);
  //--------------------------------------------------------------------------
  // B-tag stuff
  //--------------------------------------------------------------------------
  std::string btagAlgo = getPreCutString1("BTagAlgo");
  std::string btagWP = getPreCutString1("BTagWP");
  double btagCut = getPreCutValue1("BTagCutValue");

  //--------------------------------------------------------------------------
  // BDT weight file
  //--------------------------------------------------------------------------
  std::string bdtWeightFileName = "";
  bool evaluateBDT = false;
  if(hasPreCut("BDTWeightFileName")) {
    bdtWeightFileName = getPreCutString1("BDTWeightFileName");
    evaluateBDT = true;
  }
  else if(hasPreCutMatch("BDTWeightFile")) {
    evaluateBDT = true;
  }
  else if(hasPreCut("EvaluateBDT")) {
    std::string evalBDT = getPreCutString1("EvaluateBDT");
    if(evalBDT == "true" || evalBDT == "True")
      evaluateBDT = true;
  }

  //--------------------------------------------------------------------------
  // Create hists
  //--------------------------------------------------------------------------

  CreateUserHist( "EventCount"            ,    1    , 0      , 1	 );    

  CreateUserHist( "FakeRateEffective"     ,    50   , 0      , 1	 );    
  CreateUserHist( "FakeRateEffective_PassNEle"     ,    50   , 0      , 1	 );    
  CreateUserHist( "MinPrescale"           ,    1000 , 0      , 1000	 );    

  CreateUserHist( "Trigger0OrOffline1Match" ,    2, 0, 2);
  CreateUserHist( "Trigger0OrOffline1Match_PassingTrigger" ,    2, 0, 2);
  CreateUserHist( "Trigger0OrOffline1Match_PAS" ,    2, 0, 2);
  CreateUserHist( "Trigger_OfflineMatch_PassingTrigger_Ele1Pt" ,    100, 0, 1000);
  GetUserHist<TH1>("Trigger0OrOffline1Match").GetXaxis()->SetBinLabel(1, "trigger");
  GetUserHist<TH1>("Trigger0OrOffline1Match").GetXaxis()->SetBinLabel(2, "offline");
  GetUserHist<TH1>("Trigger0OrOffline1Match_PassingTrigger").GetXaxis()->SetBinLabel(1, "trigger");
  GetUserHist<TH1>("Trigger0OrOffline1Match_PassingTrigger").GetXaxis()->SetBinLabel(2, "offline");
  GetUserHist<TH1>("Trigger0OrOffline1Match_PAS").GetXaxis()->SetBinLabel(1, "trigger");
  GetUserHist<TH1>("Trigger0OrOffline1Match_PAS").GetXaxis()->SetBinLabel(2, "offline");
  GetUserHist<TH1>("Trigger_OfflineMatch_PassingTrigger_Ele1Pt").GetXaxis()->SetBinLabel(1, "trigger");
  GetUserHist<TH1>("Trigger_OfflineMatch_PassingTrigger_Ele1Pt").GetXaxis()->SetBinLabel(2, "offline");

  CreateUserHist( "M_j1j3_PAS"            ,    200 , 0       , 2000	 );    
  CreateUserHist( "M_j2j3_PAS"            ,    200 , 0       , 2000	 ); 
  CreateUserHist( "M_e1j3_PAS"            ,    200 , 0       , 2000	 );    
  CreateUserHist( "M_e2j3_PAS"            ,    200 , 0       , 2000	 ); 
  CreateUserHist( "M_eejjj_PAS"           ,    500 , 0       , 5000	 ); 

  //CreateUserHist( "M_j1j3_PASandMee100"   ,    200 , 0       , 2000	 );    
  //CreateUserHist( "M_j2j3_PASandMee100"   ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "M_e1j3_PASandMee100"   ,    200 , 0       , 2000	 );    
  //CreateUserHist( "M_e2j3_PASandMee100"   ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "M_eejjj_PASandMee100"  ,    500 , 0       , 5000	 ); 


  CreateUserHist( "sTfrac_Jet1_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserHist( "sTfrac_Jet2_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserHist( "sTfrac_Ele1_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserHist( "sTfrac_Ele2_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserHist( "sTfrac_Jet_PAS"        ,   100  ,  0.0    , 1.0      );
  CreateUserHist( "sTfrac_Ele_PAS"        ,   100  ,  0.0    , 1.0      );
  //CreateUserHist( "sTfrac_Jet1_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  //CreateUserHist( "sTfrac_Jet2_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  //CreateUserHist( "sTfrac_Ele1_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  //CreateUserHist( "sTfrac_Ele2_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  //CreateUserHist( "sTfrac_Jet_PASandMee100"        ,   100  ,  0.0    , 1.0      );
  //CreateUserHist( "sTfrac_Ele_PASandMee100"        ,   100  ,  0.0    , 1.0      );
  CreateUserHist( "nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nMuon_PAS"             ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nJet_PAS"              ,    10  , -0.5    , 9.5      );
  //CreateUserHist( "nJet_PASandMee100"        ,    10  , -0.5    , 9.5      );
  CreateUserHist( "Pt1stEle_PAS"	   , 	100 , 0       , 1000     ); 
  //CreateUserHist( "Pt1stEle_PASandMee100" , 	100 , 0       , 1000     ); 
  CreateUserHist( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "SCEta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "DeltaEtaEleTrk1stEle_Presel", 400, -0.5,   0.5 );
  CreateUserHist( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt2ndEle_PAS"	   , 	300 , 0       , 3000     ); 
  //CreateUserHist( "Pt2ndEle_PASandMee100" , 	300 , 0       , 3000     ); 
  CreateUserHist( "Eta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "SCEta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "DeltaEtaEleTrk2ndEle_Presel", 400, -0.5,   0.5 );
  CreateUserHist( "Phi2ndEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "Charge2ndEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "EleChargeSum_PAS"         ,    3   , -2.5    , 2.5  );
  //CreateUserHist( "EleChargeSum_PASandMee100",    3   , -2.5    , 2.5  );
  CreateUserHist( "MET_PAS"               ,    200 , 0       , 1000	 ); 
  CreateUserHist( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt1stJet_PAS"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Pt2ndJet_PAS"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Pt3rdJet_PAS"          ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "Pt1stJet_PASandMee100" ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "Pt2ndJet_PASandMee100" ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Eta3rdJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Phi3rdJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sTjet_PAS"             ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sTlep_PASandMee100"    ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sTjet_PASandMee100"    ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sT_PAS"                ,    300   , 0       , 3000	 ); 
  CreateUserHist( "sT_zjj_PAS"            ,    300   , 0       , 3000	  ); 
  //CreateUserHist( "sT_zjj_PASandMee100"   ,    300   , 0       , 3000	  ); 
  //CreateUserHist( "sT_PASandMee100"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee110"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee120"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee130"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee140"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee150"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee160"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee170"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee180"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee190"       ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "sT_PASandMee200"       ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "Mjj_PASandMee100"	   ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mee_PAS"		   ,    1000 , 0       , 2000	 ); 
  //CreateUserHist( "Mee_PASandST445"       ,    2000 , 0       , 2000	 ); 
  //CreateUserHist( "MTenu_PAS"             ,    200 , 0       , 1000	 ); 
  CreateUserHist( "Me1j1_PAS"             ,    200 , 0       , 2000	 ); 
  //CreateUserHist( "Me1j1_PASandMee100"    ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me1j2_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j1_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j2_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me1j_selected_PAS"     ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j_selected_PAS"     ,    200 , 0       , 2000     );
  CreateUserHist( "Mej_selected_min_PAS"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_max_PAS"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_asym_PAS"          ,    50  , 0       , 1   ); 
  CreateUserHist( "Mej_minmax_PAS"        ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_avg_PAS"  ,    200 , 0       , 2000     );
  //CreateUserHist( "Mej_selected_avg_PASandMee100"  ,    200 , 0       , 2000     );
  CreateUserHist( "Mejj_PAS"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meej_PAS"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meejj_PAS"             ,    400 , 0       , 4000     );
  // muon kinematics
  CreateUserHist( "Pt1stMuon_PAS"	   , 	100 , 0       , 1000     ); 
  CreateUserHist( "Eta1stMuon_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stMuon_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt2ndMuon_PAS"	   , 	300 , 0       , 3000     ); 
  CreateUserHist( "Eta2ndMuon_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi2ndMuon_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 


  CreateUserHist( "Ptj1j2j3_PAS"                    ,    400 , 0       , 4000     );
  CreateUserHist( "Ptj1j2_PAS"                      ,    400 , 0       , 4000     );
  CreateUserHist( "Ptj2j3_PAS"                      ,    400 , 0       , 4000     );
  CreateUserHist( "Ptj1j3_PAS"                      ,    400 , 0       , 4000     );

  CreateUserHist( "Ptee_Minus_Ptj1j2_PAS"           ,    200 , -500    , 500      );
  CreateUserHist( "Ptee_Minus_Ptj1j2j3_PAS"         ,    200 , -500    , 500      );


  //CreateUserHist( "Ptj1j2j3_PASandMee100"           ,    200 , 0       , 2000     );
  //CreateUserHist( "Ptj1j2_PASandMee100"             ,    200 , 0       , 2000     );
  //CreateUserHist( "Ptj2j3_PASandMee100"             ,    200 , 0       , 2000     );
  //CreateUserHist( "Ptj1j3_PASandMee100"             ,    200 , 0       , 2000     );

  //CreateUserHist( "Ptee_Minus_Ptj1j2_PASandMee100"  ,    200 , -500    , 500      );
  //CreateUserHist( "Ptee_Minus_Ptj1j2j3_PASandMee100",    200 , -500    , 500      );

  CreateUserHist( "Ptee_PAS"              ,    200 , 0       , 2000     );
  //CreateUserHist( "Ptee_PASandMee100"     ,    200 , 0       , 2000     );

  CreateUserHist( "nVertex_PAS"                     ,    101   , -0.5   , 100.5	 ) ; 
  //CreateUserHist( "nVertex_PASandMee100"            ,    101   , -0.5   , 100.5	 ) ; 

  CreateUserHist( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
  CreateUserHist( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
  CreateUserHist( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
  CreateUserHist( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
  CreateUserHist( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserHist( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserHist( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserHist( "minDR_ZJet_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

  CreateUserHist( "DR_ZJet1_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserHist( "DR_ZJet2_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 


  CreateUserHist2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
  CreateUserHist2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;


  CreateUserHist2D( "MeeVsST_PAS"                 ,     200, 0, 2000, 200, 0, 2000) ;
  //CreateUserHist2D( "MeeVsST_PASandMee100"        ,     200, 0, 2000, 200, 0, 2000) ;


  //CreateUserHist( "Mee_80_100_Preselection", 200, 60, 120 );
  //CreateUserHist( "Mee_70_110_Preselection", 200, 60, 120 );
  //CreateUserHist( "Mee_70_110_ST600_Preselection", 200, 60, 120 );

  CreateUserHist( "Mee_EBEB_PAS"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEE_PAS"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EEEE_PAS"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EB_PAS" 		   ,    60 , 60       , 120	 ); 
  CreateUserHist( "Mee_End2End2_PAS" 		   ,    60 , 60       , 120	 ); 

  //CreateUserHist( "Mee_EBEB_80_100_PAS"	   ,    60 , 60       , 120	 ); 
  //CreateUserHist( "Mee_EBEE_80_100_PAS"	   ,    60 , 60       , 120	 ); 
  //CreateUserHist( "Mee_EEEE_80_100_PAS"	   ,    60 , 60       , 120	 ); 
  //CreateUserHist( "Mee_EB_80_100_PAS" 	     	   ,    60 , 60       , 120	 ); 

  //CreateUserHist( "Mee_EBEB_70_110_PAS"	   ,    60 , 60       , 120	 ); 
  //CreateUserHist( "Mee_EBEE_70_110_PAS"	   ,    60 , 60       , 120	 ); 
  //CreateUserHist( "Mee_EEEE_70_110_PAS"	   ,    60 , 60       , 120	 ); 
  //CreateUserHist( "Mee_EB_70_110_PAS" 	     	   ,    60 , 60       , 120	 ); 

  CreateUserHist( "PileupWeight"   , 100, -10, 10 );
  CreateUserHist( "GeneratorWeight", 100, -2.0 , 2.0 );

  CreateUserHist("CorrIsolation_1stEle_PAS"                 , 200,-25.0 ,  25.0  ); CreateUserHist("CorrIsolation_2ndEle_PAS"                 , 200,-25.0 ,  25.0  );
  CreateUserHist("DeltaEtaTrkSC_1stEle_PAS"                 , 200, -0.01,   0.01 ); CreateUserHist("DeltaEtaTrkSC_2ndEle_PAS"                 , 200, -0.01,   0.01 );
  CreateUserHist("EcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserHist("EcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
  CreateUserHist("HcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserHist("HcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
  CreateUserHist("TrkIsolation_1stEle_PAS"                  , 200,  0.0,    5.0  ); CreateUserHist("TrkIsolation_2ndEle_PAS"                  , 200,  0.0,    5.0  );
  CreateUserHist("HasMatchedPhot_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserHist("HasMatchedPhot_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
  CreateUserHist("HoE_1stEle_PAS"                           , 200,  0.0 ,   0.05 ); CreateUserHist("HoE_2ndEle_PAS"                           , 200,  0.0 ,   0.05 );
  CreateUserHist("LeadVtxDistXY_1stEle_PAS"                 , 200, -0.05,   0.05 ); CreateUserHist("LeadVtxDistXY_2ndEle_PAS"                 , 200, -0.05,   0.05 );
  CreateUserHist("LeadVtxDistZ_1stEle_PAS"                  , 200, -0.2 ,   0.2  ); CreateUserHist("LeadVtxDistZ_2ndEle_PAS"                  , 200, -0.2 ,   0.2  );
  CreateUserHist("MissingHits_1stEle_PAS"                   , 2  , -0.5,    1.5  ); CreateUserHist("MissingHits_2ndEle_PAS"                   , 2  , -0.5,    1.5  );
  CreateUserHist("Full5x5SigmaIEtaIEta_Barrel_1stEle_PAS"   , 200,  0.0,    0.04 ); CreateUserHist("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PAS"   , 200,  0.0,    0.04 );
  CreateUserHist("Full5x5SigmaIEtaIEta_Endcap_1stEle_PAS"   , 200,  0.0,    0.1  ); CreateUserHist("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PAS"   , 200,  0.0,    0.1  );

  //CreateUserHist("CorrIsolation_1stEle_PASandMee100"        , 200,-25.0 ,  25.0  ); CreateUserHist("CorrIsolation_2ndEle_PASandMee100"        , 200,-25.0 ,  25.0  );
  //CreateUserHist("DeltaEtaTrkSC_1stEle_PASandMee100"        , 200, -0.01,   0.01 ); CreateUserHist("DeltaEtaTrkSC_2ndEle_PASandMee100"        , 200, -0.01,   0.01 );
  //CreateUserHist("EcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserHist("EcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
  //CreateUserHist("HcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserHist("HcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
  //CreateUserHist("TrkIsolation_1stEle_PASandMee100"         , 200,  0.0,    5.0  ); CreateUserHist("TrkIsolation_2ndEle_PASandMee100"         , 200,  0.0,    5.0  );
  //CreateUserHist("HasMatchedPhot_1stEle_PASandMee100"       , 2,   -0.5 ,   1.5  ); CreateUserHist("HasMatchedPhot_2ndEle_PASandMee100"       , 2,   -0.5 ,   1.5  );
  //CreateUserHist("HoE_1stEle_PASandMee100"                  , 200,  0.0 ,   0.05 ); CreateUserHist("HoE_2ndEle_PASandMee100"                  , 200,  0.0 ,   0.05 );
  //CreateUserHist("LeadVtxDistXY_1stEle_PASandMee100"        , 200, -0.05,   0.05 ); CreateUserHist("LeadVtxDistXY_2ndEle_PASandMee100"        , 200, -0.05,   0.05 );
  //CreateUserHist("LeadVtxDistZ_1stEle_PASandMee100"         , 200, -0.2 ,   0.2  ); CreateUserHist("LeadVtxDistZ_2ndEle_PASandMee100"         , 200, -0.2 ,   0.2  );
  //CreateUserHist("MissingHits_1stEle_PASandMee100"          , 2  , -0.5,    1.5  ); CreateUserHist("MissingHits_2ndEle_PASandMee100"          , 2  , -0.5,    1.5  );
  //CreateUserHist("Full5x5SigmaIEtaIEta_Barrel_1stEle_PASandMee100" , 200,  0.0,    0.02 ); CreateUserHist("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PASandMee100" , 200,  0.0,    0.02 );
  //CreateUserHist("Full5x5SigmaIEtaIEta_Endcap_1stEle_PASandMee100" , 200,  0.0,    0.1  ); CreateUserHist("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PASandMee100" , 200,  0.0,    0.1  );

  if(do_roi_plots) {
    CreateUserHist( "M_j1j3_ROI"            ,    200 , 0       , 2000	 );    
    CreateUserHist( "M_j2j3_ROI"            ,    200 , 0       , 2000	 ); 
    CreateUserHist( "M_e1j3_ROI"            ,    200 , 0       , 2000	 );    
    CreateUserHist( "M_e2j3_ROI"            ,    200 , 0       , 2000	 ); 
    CreateUserHist( "M_eejjj_ROI"           ,    500 , 0       , 5000	 ); 
    CreateUserHist( "sTfrac_Jet1_ROI"       ,   100  ,  0.0    , 1.0      );
    CreateUserHist( "sTfrac_Jet2_ROI"       ,   100  ,  0.0    , 1.0      );
    CreateUserHist( "sTfrac_Ele1_ROI"       ,   100  ,  0.0    , 1.0      );
    CreateUserHist( "sTfrac_Ele2_ROI"       ,   100  ,  0.0    , 1.0      );
    CreateUserHist( "sTfrac_Jet_ROI"        ,   100  ,  0.0    , 1.0      );
    CreateUserHist( "sTfrac_Ele_ROI"        ,   100  ,  0.0    , 1.0      );
    CreateUserHist( "nJet_ROI"              ,    10  , -0.5    , 9.5      );
    CreateUserHist( "Pt1stEle_ROI"	   , 	100 , 0       , 1000     ); 
    CreateUserHist( "Pt2ndEle_ROI"	   , 	300 , 0       , 3000     ); 
    CreateUserHist( "EleChargeSum_ROI"         ,    3   , -2.5    , 2.5  );
    CreateUserHist( "MET_ROI"               ,    200 , 0       , 1000	 ); 
    CreateUserHist( "Pt1stJet_ROI"          ,    200 , 0       , 2000	 ); 
    CreateUserHist( "Pt2ndJet_ROI"          ,    200 , 0       , 2000	 ); 
    CreateUserHist( "sTlep_ROI"             ,    200 , 0       , 2000	 ); 
    CreateUserHist( "sTjet_ROI"             ,    200 , 0       , 2000	 ); 
    CreateUserHist( "sT_zjj_ROI"                      ,    200   , 0       , 2000	  ); 
    CreateUserHist( "sT_ROI"                ,    200 , 0       , 2000	 ); 
    CreateUserHist( "Mjj_ROI"		   ,    200 , 0       , 2000	 ); 
    CreateUserHist( "Mee_ROI"		   ,    200 , 0       , 2000	 ); 
    CreateUserHist( "Me1j1_ROI"             ,    200 , 0       , 2000	 ); 
    CreateUserHist( "Mej_selected_avg_ROI"  ,    200 , 0       , 2000     );
    CreateUserHist( "Meejj_ROI"             ,    400 , 0       , 4000     );
    CreateUserHist( "Mejj_ROI"              ,    400 , 0       , 4000     );
    CreateUserHist( "Meej_ROI"              ,    400 , 0       , 4000     );
    CreateUserHist( "Eta1stJet_ROI"                   ,    100   , -5      , 5	  ); 
    CreateUserHist( "Eta2ndJet_ROI"                   ,    100   , -5      , 5	  ); 
    CreateUserHist( "Eta1stEle_ROI"	             , 	100    , -5      , 5	  ); 
    CreateUserHist( "Eta2ndEle_ROI"	             , 	100    , -5      , 5	  ); 

    CreateUserHist( "Phi1stJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
    CreateUserHist( "Phi2ndJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
    CreateUserHist( "Phi1stEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
    CreateUserHist( "Phi2ndEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
    CreateUserHist( "Ptj1j2j3_ROI"                    ,    200 , 0       , 2000     );
    CreateUserHist( "Ptj1j2_ROI"                      ,    200 , 0       , 2000     );
    CreateUserHist( "Ptj2j3_ROI"                      ,    200 , 0       , 2000     );
    CreateUserHist( "Ptj1j3_ROI"                      ,    200 , 0       , 2000     );
    CreateUserHist( "Ptee_Minus_Ptj1j2_ROI"           ,    200 , -500    , 500      );
    CreateUserHist( "Ptee_Minus_Ptj1j2j3_ROI"         ,    200 , -500    , 500      );
    CreateUserHist( "Ptee_ROI"              ,    200 , 0       , 2000     );
    CreateUserHist( "nVertex_ROI"                     ,    101   , -0.5   , 100.5	 ) ; 
    CreateUserHist( "minDR_ZJet_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    CreateUserHist( "DR_ZJet1_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    CreateUserHist2D( "MeeVsST_ROI"                 ,     200, 0, 2000, 200, 0, 2000) ;
    CreateUserHist( "DR_ZJet2_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    CreateUserHist("CorrIsolation_1stEle_ROI"                 , 200,-25.0 ,  25.0  ); CreateUserHist("CorrIsolation_2ndEle_ROI"                 , 200,-25.0 ,  25.0  );
    CreateUserHist("DeltaEtaTrkSC_1stEle_ROI"                 , 200, -0.01,   0.01 ); CreateUserHist("DeltaEtaTrkSC_2ndEle_ROI"                 , 200, -0.01,   0.01 );
    CreateUserHist("E1x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserHist("E1x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
    CreateUserHist("E2x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserHist("E2x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
    CreateUserHist("EcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserHist("EcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
    CreateUserHist("HcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserHist("HcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
    CreateUserHist("TrkIsolation_1stEle_ROI"                  , 200,  0.0,    5.0  ); CreateUserHist("TrkIsolation_2ndEle_ROI"                  , 200,  0.0,    5.0  );
    CreateUserHist("HasMatchedPhot_1stEle_ROI"                , 2,   -0.5 ,   1.5  ); CreateUserHist("HasMatchedPhot_2ndEle_ROI"                , 2,   -0.5 ,   1.5  );
    CreateUserHist("HoE_1stEle_ROI"                           , 200,  0.0 ,   0.05 ); CreateUserHist("HoE_2ndEle_ROI"                           , 200,  0.0 ,   0.05 );
    CreateUserHist("LeadVtxDistXY_1stEle_ROI"                 , 200, -0.05,   0.05 ); CreateUserHist("LeadVtxDistXY_2ndEle_ROI"                 , 200, -0.05,   0.05 );
    CreateUserHist("LeadVtxDistZ_1stEle_ROI"                  , 200, -0.2 ,   0.2  ); CreateUserHist("LeadVtxDistZ_2ndEle_ROI"                  , 200, -0.2 ,   0.2  );
    CreateUserHist("MissingHits_1stEle_ROI"                   , 2  , -0.5,    1.5  ); CreateUserHist("MissingHits_2ndEle_ROI"                   , 2  , -0.5,    1.5  );
    CreateUserHist("SigmaIEtaIEta_Barrel_1stEle_ROI"          , 200,  0.0,    0.02 ); CreateUserHist("SigmaIEtaIEta_Barrel_2ndEle_ROI"          , 200,  0.0,    0.02 );
    CreateUserHist("SigmaIEtaIEta_Endcap_1stEle_ROI"          , 200,  0.0,    0.1  ); CreateUserHist("SigmaIEtaIEta_Endcap_2ndEle_ROI"          , 200,  0.0,    0.1  );
  }

  // for scale factor dependence studies
  //CreateUserHist( "Mee_NJetEq2_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetEq3_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetEq4_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetEq5_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetEq6_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetEq7_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetGeq3_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_NJetGeq4_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT300To500_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT500To750_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT750To1250_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT1250ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin100To200_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin200To300_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin300To400_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin400To500_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin500To650_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin650ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
  // with zero B-tags
  CreateUserHist( "Mee_PAS_noBtaggedJets"		       ,      1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEB_PAS_noBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEE_PAS_noBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EEEE_PAS_noBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  //CreateUserHist( "Mee_sT300To500_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT500To750_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT750To1250_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT1250ToInf_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin100To200_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin200To300_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin300To400_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin400To500_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin500To650_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin650ToInf_PAS_noBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  CreateUserHist( "nElectron_noBtaggedJets"         ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nMuon_noBtaggedJets"             ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nJet_noBtaggedJets"              ,    10  , -0.5    , 9.5      );
  CreateUserHist( "Pt1stEle_noBtaggedJets"	   , 	100 , 0       , 1000     ); 
  CreateUserHist( "Eta1stEle_noBtaggedJets"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stEle_noBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt2ndEle_noBtaggedJets"	   , 	300 , 0       , 3000     ); 
  CreateUserHist( "Eta2ndEle_noBtaggedJets"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi2ndEle_noBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Charge1stEle_noBtaggedJets"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "Charge2ndEle_noBtaggedJets"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "MET_noBtaggedJets"               ,    200 , 0       , 1000	 ); 
  CreateUserHist( "METPhi_noBtaggedJets"		   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt1stJet_noBtaggedJets"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Pt2ndJet_noBtaggedJets"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Eta1stJet_noBtaggedJets"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Eta2ndJet_noBtaggedJets"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stJet_noBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Phi2ndJet_noBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "sTlep_noBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sTjet_noBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sT_noBtaggedJets"                ,    300   , 0       , 3000	 ); 
  CreateUserHist( "sT_zjj_noBtaggedJets"            ,    300   , 0       , 3000	  ); 
  CreateUserHist( "Mjj_noBtaggedJets"		   ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mee_noBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "MTenu_noBtaggedJets"             ,    200 , 0       , 1000	 ); 
  CreateUserHist( "Me1j1_noBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me1j2_noBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j1_noBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j2_noBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mej_selected_min_noBtaggedJets"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_max_noBtaggedJets"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_minmax_noBtaggedJets"        ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_avg_noBtaggedJets"  ,    200 , 0       , 2000     );
  CreateUserHist( "Mejj_noBtaggedJets"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meej_noBtaggedJets"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meejj_noBtaggedJets"             ,    400 , 0       , 4000     );
  // with >= 1 B-tags
  CreateUserHist( "Mee_PAS_gteOneBtaggedJet"		       ,      1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEB_PAS_gteOneBtaggedJet"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEE_PAS_gteOneBtaggedJet"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EEEE_PAS_gteOneBtaggedJet"		   ,    1000 , 0       , 2000	 ); 
  //CreateUserHist( "Mee_sT300To500_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT500To750_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT750To1250_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT1250ToInf_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin100To200_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin200To300_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin300To400_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin400To500_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin500To650_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin650ToInf_PAS_gteOneBtaggedJet"		             ,    200   , 0       , 2000	  ); 
  CreateUserHist( "nElectron_gteOneBtaggedJet"         ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nMuon_gteOneBtaggedJet"             ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nJet_gteOneBtaggedJet"              ,    10  , -0.5    , 9.5      );
  CreateUserHist( "Pt1stEle_gteOneBtaggedJet"	   , 	100 , 0       , 1000     ); 
  CreateUserHist( "Eta1stEle_gteOneBtaggedJet"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stEle_gteOneBtaggedJet"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt2ndEle_gteOneBtaggedJet"	   , 	300 , 0       , 3000     ); 
  CreateUserHist( "Eta2ndEle_gteOneBtaggedJet"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi2ndEle_gteOneBtaggedJet"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Charge1stEle_gteOneBtaggedJet"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "Charge2ndEle_gteOneBtaggedJet"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "MET_gteOneBtaggedJet"               ,    200 , 0       , 1000	 ); 
  CreateUserHist( "METPhi_gteOneBtaggedJet"		   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt1stJet_gteOneBtaggedJet"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Pt2ndJet_gteOneBtaggedJet"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Eta1stJet_gteOneBtaggedJet"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Eta2ndJet_gteOneBtaggedJet"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stJet_gteOneBtaggedJet"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Phi2ndJet_gteOneBtaggedJet"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "sTlep_gteOneBtaggedJet"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sTjet_gteOneBtaggedJet"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sT_gteOneBtaggedJet"                ,    300   , 0       , 3000	 ); 
  CreateUserHist( "sT_zjj_gteOneBtaggedJet"            ,    300   , 0       , 3000	  ); 
  CreateUserHist( "Mjj_gteOneBtaggedJet"		   ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mee_gteOneBtaggedJet"		   ,     1000 , 0       , 2000	 ); 
  CreateUserHist( "MTenu_gteOneBtaggedJet"             ,    200 , 0       , 1000	 ); 
  CreateUserHist( "Me1j1_gteOneBtaggedJet"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me1j2_gteOneBtaggedJet"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j1_gteOneBtaggedJet"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j2_gteOneBtaggedJet"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mej_selected_min_gteOneBtaggedJet"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_max_gteOneBtaggedJet"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_minmax_gteOneBtaggedJet"        ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_avg_gteOneBtaggedJet"  ,    200 , 0       , 2000     );
  CreateUserHist( "Mejj_gteOneBtaggedJet"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meej_gteOneBtaggedJet"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meejj_gteOneBtaggedJet"             ,    400 , 0       , 4000     );
  // with >= 2 B-tags
  CreateUserHist( "Mee_PAS_gteTwoBtaggedJets"		       ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEB_PAS_gteTwoBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EBEE_PAS_gteTwoBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  CreateUserHist( "Mee_EEEE_PAS_gteTwoBtaggedJets"		   ,    1000 , 0       , 2000	 ); 
  //CreateUserHist( "Mee_sT300To500_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT500To750_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT750To1250_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_sT1250ToInf_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin100To200_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin200To300_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin300To400_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin400To500_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin500To650_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "Mee_MejMin650ToInf_PAS_gteTwoBtaggedJets"		             ,    200   , 0       , 2000	  ); 
  CreateUserHist( "nElectron_gteTwoBtaggedJets"         ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nMuon_gteTwoBtaggedJets"             ,    5   , -0.5    , 4.5      );
  CreateUserHist( "nJet_gteTwoBtaggedJets"              ,    10  , -0.5    , 9.5      );
  CreateUserHist( "Pt1stEle_gteTwoBtaggedJets"	   , 	100 , 0       , 1000     ); 
  CreateUserHist( "Eta1stEle_gteTwoBtaggedJets"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stEle_gteTwoBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt2ndEle_gteTwoBtaggedJets"	   , 	300 , 0       , 3000     ); 
  CreateUserHist( "Eta2ndEle_gteTwoBtaggedJets"	   , 	100 , -5      , 5	 ); 
  CreateUserHist( "Phi2ndEle_gteTwoBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Charge1stEle_gteTwoBtaggedJets"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "Charge2ndEle_gteTwoBtaggedJets"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserHist( "MET_gteTwoBtaggedJets"               ,    200 , 0       , 1000	 ); 
  CreateUserHist( "METPhi_gteTwoBtaggedJets"		   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Pt1stJet_gteTwoBtaggedJets"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Pt2ndJet_gteTwoBtaggedJets"          ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Eta1stJet_gteTwoBtaggedJets"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Eta2ndJet_gteTwoBtaggedJets"         ,    100 , -5      , 5	 ); 
  CreateUserHist( "Phi1stJet_gteTwoBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "Phi2ndJet_gteTwoBtaggedJets"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserHist( "sTlep_gteTwoBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sTjet_gteTwoBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "sT_gteTwoBtaggedJets"                ,    300   , 0       , 3000	 ); 
  CreateUserHist( "sT_zjj_gteTwoBtaggedJets"            ,    300   , 0       , 3000	  ); 
  CreateUserHist( "Mjj_gteTwoBtaggedJets"		   ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mee_gteTwoBtaggedJets"		   ,     1000 , 0       , 2000	 ); 
  CreateUserHist( "MTenu_gteTwoBtaggedJets"             ,    200 , 0       , 1000	 ); 
  CreateUserHist( "Me1j1_gteTwoBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me1j2_gteTwoBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j1_gteTwoBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Me2j2_gteTwoBtaggedJets"             ,    200 , 0       , 2000	 ); 
  CreateUserHist( "Mej_selected_min_gteTwoBtaggedJets"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_max_gteTwoBtaggedJets"  ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_minmax_gteTwoBtaggedJets"        ,    200 , 0       , 2000     ); 
  CreateUserHist( "Mej_selected_avg_gteTwoBtaggedJets"  ,    200 , 0       , 2000     );
  CreateUserHist( "Mejj_gteTwoBtaggedJets"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meej_gteTwoBtaggedJets"              ,    400 , 0       , 4000     );
  CreateUserHist( "Meejj_gteTwoBtaggedJets"             ,    400 , 0       , 4000     );

  // bkg control region plots
  CreateUserHistWithSysts( "Mee_BkgControlRegion"		                  ,    1000   , 0       , 2000	  ); 
  CreateUserHistWithSysts( "Mee_BkgControlRegion_gteOneBtaggedJet"		,    1000 , 0       , 2000	 ); 
  CreateUserHistWithSysts( "Mee_BkgControlRegion_gteTwoBtaggedJets"		,    1000 , 0       , 2000	 ); 
  CreateUserHistWithSysts( "Mee_EB_BkgControlRegion"		         ,    1000 , 0       , 2000	 ); 
  CreateUserHistWithSysts( "Mee_EBEB_BkgControlRegion"		       ,    1000   , 0       , 2000	 ); 
  CreateUserHistWithSysts( "Mee_EBEE_BkgControlRegion"		       ,    1000   , 0       , 2000	 ); 
  CreateUserHistWithSysts( "Mee_EEEE_BkgControlRegion"		       ,    1000   , 0       , 2000	 ); 
  CreateUserHistWithSysts( "Mee_End2End2_BkgControlRegion"		   ,    1000   , 0       , 2000	 ); 
  CreateUserHist2D("MeeVsNJet_BkgControlRegion", 10, -0.5, 9.5, 1000, 0, 2000);
  CreateUserHist2D("MeeVsST_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsSTjet_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsSTlep_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMejMin_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMejMax_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMeejj_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe1j1_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe1j2_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe2j1_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe2j2_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsEle1Pt_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsEle2Pt_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet1Pt_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet2Pt_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet3Pt_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsPFMETType1Pt_BkgControlRegion", 1000, 0, 2000, 1000, 0, 2000);
  // >= 1 b-tag
  CreateUserHist2D("MeeVsNJet_BkgControlRegion_gteOneBtaggedJet", 10, -0.5, 9.5, 1000, 0, 2000);
  CreateUserHist2D("MeeVsST_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsSTjet_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsSTlep_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMejMin_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMejMax_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMeejj_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe1j1_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe1j2_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe2j1_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe2j2_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsEle1Pt_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsEle2Pt_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet1Pt_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet2Pt_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet3Pt_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsPFMETType1Pt_BkgControlRegion_gteOneBtaggedJet", 1000, 0, 2000, 1000, 0, 2000);
  // >= 2 b-tags
  CreateUserHist2D("MeeVsNJet_BkgControlRegion_gteTwoBtaggedJets", 10, -0.5, 9.5, 1000, 0, 2000);
  CreateUserHist2D("MeeVsST_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsSTjet_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsSTlep_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMejMin_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMejMax_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMeejj_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe1j1_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe1j2_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe2j1_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsMe2j2_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsEle1Pt_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsEle2Pt_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet1Pt_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet2Pt_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsJet3Pt_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  CreateUserHist2D("MeeVsPFMETType1Pt_BkgControlRegion_gteTwoBtaggedJets", 1000, 0, 2000, 1000, 0, 2000);
  //// test opt
  //CreateUserHist( "Mee_sT2000_PAS"		             ,    200   , 0       , 2000	  ); 
  //CreateUserHist( "OptBinLQ600", 200, 0, 2000);
  //CreateUserHist( "PtZforOptBin600",1000,0,2000);
  //CreateUserHist( "PtEEforOptBin600",1000,0,2000);
  //CreateUserHist2D( "PtEEVsZPtForOptBin600",1000,0,2000,1000,0,2000);
  //CreateUserHist( "OptBinLQ650", 200, 0, 2000);
  //CreateUserHist( "OptBinLQ700", 200, 0, 2000);
  //CreateUserHist( "OptBinLQ600_noWeight", 200, 0, 2000);
  //CreateUserHist( "OptBinLQ650_noWeight", 200, 0, 2000);
  //CreateUserHist( "OptBinLQ700_noWeight", 200, 0, 2000);
  // 3D opt cut space
  //CreateUserTH3D( "OptimizationCutSpace", 200, 0, 2000, 200, 0, 2000, 200, 0, 2000);


  //--------------------------------------------------------------------------
  // Final selection plots
  //--------------------------------------------------------------------------
  bool doFinalSelections = false;
  // check if there is a final Mej specific in cutfile for any LQ mass
  for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
    int lq_mass = LQ_MASS[i_lq_mass];
    //TODO FIXME; hack for now
    //sprintf(cut_name, "min_M_ej_LQ%d"   , lq_mass );
    sprintf(cut_name, "BDTOutput_TrainRegion_LQ%d"   , lq_mass );
    CreateUserHist(cut_name,2000,-1,1.001);
    sprintf(cut_name, "BDTOutput_noWeight_TrainRegion_LQ%d"   , lq_mass );
    CreateUserHist(cut_name,2000,-1,1.001);
    sprintf(cut_name, "BDTOutput_LQ%d"   , lq_mass );
    if(hasCut(cut_name))
      doFinalSelections = true;
  }
  if(doFinalSelections && !evaluateBDT) {
    STDOUT("ERROR: No BDTWeightFileName specified, yet asked to do final selections. Quitting here.");
    exit(-5);
  }
  // now, we must have an Mej cut and optimization must be off to have final selections enabled
  doFinalSelections = doFinalSelections && !isOptimizationEnabled();

  char plot_name[100];

  for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
    int lq_mass = LQ_MASS[i_lq_mass];
    sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); CreateUserHist ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); CreateUserHist ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); CreateUserHist ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); CreateUserHist ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); CreateUserHist ( plot_name, 30  , 0 , 3000 );
    sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); CreateUserHist ( plot_name, 40  , 0 , 2000 );
    sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); CreateUserHist2D ( plot_name, 150  , 0 , 3000, 150  , 0 , 3000 );
    sprintf(plot_name, "DR_Ele1Jet1_LQ%d"            , lq_mass ); CreateUserHist ( plot_name, 
        getHistoNBins("DR_Ele1Jet1"), 
        getHistoMin  ("DR_Ele1Jet1"), 
        getHistoMax  ("DR_Ele1Jet1"));

    sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200,-25.0 ,  25.0  );
    sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200, -0.01,   0.01 );
    sprintf(plot_name, "EcalIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "HcalIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "TrkIsolation_1stEle_LQ%d"         , lq_mass ); CreateUserHist( plot_name , 200,  0.0,    5.0  );
    sprintf(plot_name, "HasMatchedPhot_1stEle_LQ%d"       , lq_mass ); CreateUserHist( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "HoE_1stEle_LQ%d"                  , lq_mass ); CreateUserHist( plot_name , 200,  0.0 ,   0.05 );
    sprintf(plot_name, "LeadVtxDistXY_1stEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200, -0.05,   0.05 );
    sprintf(plot_name, "LeadVtxDistZ_1stEle_LQ%d"         , lq_mass ); CreateUserHist( plot_name , 200, -0.2 ,   0.2  );
    sprintf(plot_name, "MissingHits_1stEle_LQ%d"          , lq_mass ); CreateUserHist( plot_name , 2  , -0.5,    1.5  );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d" , lq_mass ); CreateUserHist( plot_name , 200,  0.0,    0.02 );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d" , lq_mass ); CreateUserHist( plot_name , 200,  0.0,    0.1  );

    sprintf(plot_name, "CorrIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200,-25.0 ,  25.0  );
    sprintf(plot_name, "DeltaEtaTrkSC_2ndEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200, -0.01,   0.01 );
    sprintf(plot_name, "EcalIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "HcalIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "TrkIsolation_2ndEle_LQ%d"         , lq_mass ); CreateUserHist( plot_name , 200,  0.0,    5.0  );
    sprintf(plot_name, "HasMatchedPhot_2ndEle_LQ%d"       , lq_mass ); CreateUserHist( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "HoE_2ndEle_LQ%d"                  , lq_mass ); CreateUserHist( plot_name , 200,  0.0 ,   0.05 );
    sprintf(plot_name, "LeadVtxDistXY_2ndEle_LQ%d"        , lq_mass ); CreateUserHist( plot_name , 200, -0.05,   0.05 );
    sprintf(plot_name, "LeadVtxDistZ_2ndEle_LQ%d"         , lq_mass ); CreateUserHist( plot_name , 200, -0.2 ,   0.2  );
    sprintf(plot_name, "MissingHits_2ndEle_LQ%d"          , lq_mass ); CreateUserHist( plot_name , 2  , -0.5,    1.5  );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ%d" , lq_mass ); CreateUserHist( plot_name , 200,  0.0,    0.02 );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ%d" , lq_mass ); CreateUserHist( plot_name , 200,  0.0,    0.1  );

    sprintf(plot_name, "EleChargeSum_LQ%d"        , lq_mass ); CreateUserHist( plot_name ,      3 , -2.5    , 2.5      );
    sprintf(plot_name, "sTfrac_Jet1_LQ%d"         , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Jet2_LQ%d"         , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Ele1_LQ%d"         , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Ele2_LQ%d"         , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Jet_LQ%d"          , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Ele_LQ%d"          , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "nJet_LQ%d"                , lq_mass ); CreateUserHist( plot_name ,     10 , -0.5    , 9.5      );
    sprintf(plot_name, "Pt1stEle_LQ%d"	           , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi1stEle_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,     60 , -3.1416 , +3.1416  ); 
    sprintf(plot_name, "Pt2ndEle_LQ%d"	           , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta2ndEle_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi2ndEle_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,     60 , -3.1416 , 3.1416   ); 
    sprintf(plot_name, "MET_LQ%d"                 , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Pt1stJet_LQ%d"            , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Pt2ndJet_LQ%d"            , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta1stJet_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Eta2ndJet_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi1stJet_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,     60 , -3.1416 , 3.1416   ); 
    sprintf(plot_name, "Phi2ndJet_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,     60 , -3.1416 , 3.1416   ); 
    sprintf(plot_name, "sTlep_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "sTjet_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "sT_zjj_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "Mjj_LQ%d"		   , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "Meejj_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    400 ,  0.0    , 4000     );
    sprintf(plot_name, "Mejj_LQ%d"                , lq_mass ); CreateUserHist( plot_name ,    400 ,  0.0    , 4000     );
    sprintf(plot_name, "Meej_LQ%d"                , lq_mass ); CreateUserHist( plot_name ,    400 ,  0.0    , 4000     );
    sprintf(plot_name, "Ptj1j2j3_LQ%d"            , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptj1j2_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptj2j3_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptj1j3_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptee_Minus_Ptj1j2_LQ%d"   , lq_mass ); CreateUserHist( plot_name ,    200 , -500    ,  500     );
    sprintf(plot_name, "Ptee_Minus_Ptj1j2j3_LQ%d" , lq_mass ); CreateUserHist( plot_name ,    200 , -500    ,  500     );
    sprintf(plot_name, "Ptee_LQ%d"                , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me1j1_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me1j2_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me2j1_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me2j2_LQ%d"               , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "M_j1j3_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );    
    sprintf(plot_name, "M_j2j3_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "M_e1j3_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     );    
    sprintf(plot_name, "M_e2j3_LQ%d"              , lq_mass ); CreateUserHist( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "M_eejjj_LQ%d"             , lq_mass ); CreateUserHist( plot_name ,    500 ,  0.0    , 5000     ); 
    sprintf(plot_name, "nVertex_LQ%d"             , lq_mass ); CreateUserHist( plot_name ,    101 , -0.5    ,  100.5   ); 
    sprintf(plot_name, "MeeVsST_LQ%d"             , lq_mass ); CreateUserHist2D( plot_name ,    200 ,  0.0, 2000, 200, 0, 2000) ;
    sprintf(plot_name, "minDR_ZJet_LQ%d"          , lq_mass ); CreateUserHist( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    sprintf(plot_name, "DR_ZJet1_LQ%d"            , lq_mass ); CreateUserHist( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    sprintf(plot_name, "DR_ZJet2_LQ%d"            , lq_mass ); CreateUserHist( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    // muon kinematics
    sprintf(plot_name, "Pt1stMuon_LQ%d"	           , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta1stMuon_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi1stMuon_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,     60 , -3.1416 , +3.1416  ); 
    sprintf(plot_name, "Pt2ndMuon_LQ%d"	           , lq_mass ); CreateUserHist( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta2ndMuon_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi2ndMuon_LQ%d"           , lq_mass ); CreateUserHist( plot_name ,     60 , -3.1416 , 3.1416   ); 

  }
  // for SF at final selections
  CreateUserHist( "Mee_70_110_LQ300", 200, 60, 120 );
  CreateUserHist( "Mee_70_110_LQ600", 200, 60, 120 );
  CreateUserHist( "Mee_70_110_LQ800", 200, 60, 120 );
  CreateUserHist( "Mee_70_110_LQ900", 200, 60, 120 );
  CreateUserHist( "Mee_70_110_LQ1000", 200, 60, 120 );

  //--------------------------------------------------------------------------
  // Add new skim tree branches if needed
  //--------------------------------------------------------------------------
  double fakeRateEffective = 0;
  if(!hasBranch("FakeRateEffective"))
    addSkimTreeBranch("FakeRateEffective", &fakeRateEffective, "FakeRateEffective/D");
  double min_prescale = 1;
  if(!hasBranch("MinPrescale"))
    addSkimTreeBranch("MinPrescale", &min_prescale, "MinPrescale/D");
  double eventWeight = 1;
  if(!hasBranch("EventWeight"))
    addSkimTreeBranch("EventWeight", &eventWeight, "EventWeight/D");

  map<string, int> hltPathToL1PrescaleValueMap;
  bool checkForDifferingL1Prescales = false; // only works when running on all data in a single year in a single job
  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;

  //--------------------------------------------------------------------------
  // Loop over the chain
  //--------------------------------------------------------------------------
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);

    //--------------------------------------------------------------------------
    // Tricky part: refine Weight branch to sign of gen weight for powhegMiNNLO
    //  and EventTriggerScaleFactor: not used for QCD
    //--------------------------------------------------------------------------
    float weight = 1.0;
    resetSkimTreeBranchAddress("Weight", &weight);
    float eventTriggerScaleFactor = 1.0;
    resetSkimTreeBranchAddress("EventTriggerScaleFactor", &eventTriggerScaleFactor);
    
    //------------------------------------------------------------------
    // Print progress
    //------------------------------------------------------------------
    if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

    //--------------------------------------------------------------------------
    // Reset the cuts
    //--------------------------------------------------------------------------

    resetCuts();

    //--------------------------------------------------------------------------
    // Get information about gen-level reweighting (should be for Sherpa only)
    //--------------------------------------------------------------------------
    std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());

    float gen_weight = readerTools_->ReadValueBranch<Float_t>("Weight");
    if ( isData() ) gen_weight = 1.0;
    // special handling for powhegMiNNLO samples
    if(current_file_name.find("powhegMiNNLO") != std::string::npos) gen_weight = TMath::Sign(1, gen_weight);
    weight = gen_weight;
    //if ( isData && Ele2_ValidFrac > 998. ){
    //  gen_weight = 0.0;
    //  if      (  60.0 < M_e1e2 < 120. ) gen_weight = 0.61;
    //  else if ( 120.0 < M_e1e2 < 200. ) gen_weight = 0.42;
    //  else if ( 200.0 < M_e1e2        ) gen_weight = 0.42;
    //}

    //std::cout << "Gen weight = " << gen_weight << "; pileup_weight = " << pileup_weight << std::endl;

    //--------------------------------------------------------------------------
    // Do pileup re-weighting
    //--------------------------------------------------------------------------
    float pileup_weight = readerTools_->ReadValueBranch<Float_t>("puWeight");
    if ( isData() ) pileup_weight = 1.0;
    gen_weight*=pileup_weight;

    //--------------------------------------------------------------------------
    // Get information about prefire reweighting
    //--------------------------------------------------------------------------
    float prefire_weight = 1.0;
    if(analysisYearInt < 2018 && hasBranch("PrefireWeight") && !isData())
      prefire_weight = readerTools_->ReadValueBranch<Float_t>("PrefireWeight");
    //prefire_weight = readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Nom");
    gen_weight*=prefire_weight;

    // SIC remove March 2018
    //// TopPt reweight
    //// only valid for powheg
    //if(current_file_name.find("TT_") != std::string::npos) {
    //  gen_weight*=TopPtWeight;
    //}

    //// apply PDF rescale for LQToDEle/LQToBEle 2016 only
    //if(analysisYearInt==2016) {
    //  if(current_file_name.find("LQToBEle") != std::string::npos || current_file_name.find("LQToDEle") != std::string::npos ) {
    //    gen_weight*=readerTools_->ReadArrayBranch<Float_t>("LHEPdfWeight", 0);
    //    //std::cout << "INFO: Applying LHEPdfWeight=" << readerTools_->ReadArrayBranch<Float_t>("LHEPdfWeight", 0) << " for run: " << run << " ls: " << ls << " event: " << event << std::endl;
    //  }
    //}

    //std::cout << "prefire weight = " << prefire_weight << std::endl;

    // Electron scale factors for MC only
    if(!isData()) {
      float recoHeepSF = 1.0;
      bool verbose = false;
      float recoSFEle1 = readerTools_->ReadValueBranch<Float_t>("Ele1_RecoSF");
      float recoSFEle2 = readerTools_->ReadValueBranch<Float_t>("Ele2_RecoSF");
      // HEEP ID
      float heepSFEle1 = readerTools_->ReadValueBranch<Float_t>("Ele1_HEEPSF");
      float heepSFEle2 = readerTools_->ReadValueBranch<Float_t>("Ele2_HEEPSF");
      // EGM loose ID
      float looseSFEle1 = readerTools_->ReadValueBranch<Float_t>("Ele1_EGMLooseIDSF");
      float looseSFEle2 = readerTools_->ReadValueBranch<Float_t>("Ele2_EGMLooseIDSF");
      // z vtx (2017)
      float zVtxSF =  readerTools_->ReadValueBranch<Float_t>("ZVtxSF");
      gen_weight*=zVtxSF;
      gen_weight*=recoSFEle1*recoSFEle2*heepSFEle1*heepSFEle2;
      // EWK NLO
      if(current_file_name.find("DYJetsTo") != std::string::npos) {
        gen_weight*=readerTools_->ReadValueBranch<Float_t>("EWKNLOCorrection");
      }
    }

    // run ls event
    unsigned int run = readerTools_->ReadValueBranch<UInt_t>("run");
    unsigned int ls = readerTools_->ReadValueBranch<UInt_t>("ls");
    unsigned long long int event = readerTools_->ReadValueBranch<ULong64_t>("event");

    //--------------------------------------------------------------------------
    // Find the right prescale for this event
    //--------------------------------------------------------------------------
    min_prescale = 1;
    bool passTrigger = false;
    std::string triggerName;
    std::string triggerMatch = "none";

    float hltPhoton1Pt = readerTools_->ReadValueBranch<Float_t>("Ele1_MatchedHLTriggerObjectPt");
    //float hltPhoton2Pt = readerTools_->ReadValueBranch<Float_t>("Ele2_MatchedHLTriggerObjectPt");
    float Ele1_Pt = readerTools_->ReadValueBranch<Float_t>("Ele1_Pt");
    float Ele2_Pt = readerTools_->ReadValueBranch<Float_t>("Ele2_Pt");
    float Ele3_Pt = readerTools_->ReadValueBranch<Float_t>("Ele3_Pt");
    float leadElePt = Ele1_Pt;
    
    vector<string> triggersToConsider;
    vector<int> photonThresholds;
    if(analysisYearInt==2016)
      photonThresholds = {22, 30, 36, 50, 75, 90, 120, 175};
    else if(analysisYearInt==2017)
      photonThresholds = {25, 33, 50, 75, 90, 120, 150, 175, 200};
    else if(analysisYearInt==2018)
      photonThresholds = {33, 50, 75, 90, 120, 150, 175, 200};
    //for(const auto& thresh : photonThresholds)
    //  triggersToConsider.push_back("Photon"+to_string(thresh));

    for(unsigned int idx = 0; idx < photonThresholds.size(); ++idx) {
      int lowerThresh = photonThresholds[idx];
      int upperThresh = (idx+1) >= photonThresholds.size() ? std::numeric_limits<int>::max() : photonThresholds[idx+1];
      std::string trigName = "Photon"+to_string(lowerThresh);
      if(readerTools_->ReadValueBranch<Float_t>("H_"+trigName) > 0.1) {
       if(CheckTriggerThreshold(lowerThresh, upperThresh,  hltPhoton1Pt)) {
         passTrigger = true;
         triggerName = trigName;
         triggerMatch = "online";
       }
       else if(CheckTriggerThreshold(lowerThresh, upperThresh,  Ele1_Pt)) {
         // we no longer select these events, but we keep track of whether they _would_ have been kept, at least as far as the trigger matching is concerned
         //passTrigger = true;
         //triggerName = trigName;
         triggerMatch = "offline";
       }
      }
    }

    if(isData() && passTrigger) {
      //min_prescale = run2PhotonTriggerPrescales.LookupPrescale(analysisYear,triggerName);
      int hltPrescale = psProv.hltPrescale("HLT_"+triggerName+"_v", run, ls);
      //std::cout << "INFO: Found HLT prescale = " << hltPrescale << " for trigger " << triggerName << endl;
      int l1Prescale = 1;
      std::string l1Seed = "";
      if(triggerName == "Photon22") {
        l1Seed = "L1_SingleEG18";
        l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
      }
      else if(triggerName == "Photon25" || triggerName == "Photon30" || triggerName == "Photon36") {
        l1Seed = "L1_SingleEG26";
        l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
      }
      else if(triggerName == "Photon33") {
        if(analysisYearInt == 2017) {
          l1Seed = "L1_SingleEG26";
          l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
        }
        if(analysisYearInt == 2018) {
          l1Seed = "L1_SingleEG26er2p5";
          l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
          if(checkForDifferingL1Prescales && hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end()) {
            if(hltPathToL1PrescaleValueMap[triggerName] != l1Prescale) {
              std::cout << "WARN: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << " while previously we stored a prescale = " << hltPathToL1PrescaleValueMap[triggerName] << "; hlt prescale = " << hltPrescale << 
                "; run = " << run << " ls = " << ls << "; psColumn = " << psProv.getRunInfo(run)->psColumn(ls) << "; l1 menu=" << psProv.getRunInfo(run)->l1Menu() <<
                "; hlt menu=" << psProv.getRunInfo(run)->hltMenu() << "; trig mode=" << psProv.getRunInfo(run)->triggerMode() << std::endl;
            }
          }
          if(l1Prescale <= 0) {
            l1Seed = "L1_SingleEG26";
            l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
            if(checkForDifferingL1Prescales && hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end()) {
              if(hltPathToL1PrescaleValueMap[triggerName] != l1Prescale) {
                std::cout << "WARN: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << " while previously we stored a prescale = " << hltPathToL1PrescaleValueMap[triggerName] << "; hlt prescale = " << hltPrescale << 
                  "; run = " << run << " ls = " << ls << "; psColumn = " << psProv.getRunInfo(run)->psColumn(ls) << "; l1 menu=" << psProv.getRunInfo(run)->l1Menu() <<
                  "; hlt menu=" << psProv.getRunInfo(run)->hltMenu() << "; trig mode=" << psProv.getRunInfo(run)->triggerMode() << std::endl;
              }
            }
          }
          if(hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end())
            hltPathToL1PrescaleValueMap[triggerName] = l1Prescale;
        }
      }
      else if(triggerName == "Photon50" || triggerName == "Photon75" || triggerName == "Photon90" || triggerName == "Photon120" ||
          triggerName == "Photon150" || (analysisYearInt > 2016 && triggerName == "Photon175") ) {
        l1Seed = "L1_SingleEG40";
        l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
        if(analysisYearInt == 2018) {
          l1Seed = "L1_SingleEG42er2p5";
          l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
          if(checkForDifferingL1Prescales && hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end()) {
            if(hltPathToL1PrescaleValueMap[triggerName] != l1Prescale) {
              std::cout << "WARN: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << " while previously we stored a prescale = " << hltPathToL1PrescaleValueMap[triggerName] << "; hlt prescale = " << hltPrescale << 
                "; run = " << run << " ls = " << ls << "; psColumn = " << psProv.getRunInfo(run)->psColumn(ls) << "; l1 menu=" << psProv.getRunInfo(run)->l1Menu() <<
                "; hlt menu=" << psProv.getRunInfo(run)->hltMenu() << "; trig mode=" << psProv.getRunInfo(run)->triggerMode() << std::endl;
            }
          }
          if(l1Prescale <= 0) {
            l1Seed = "L1_SingleEG42";
            l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
            if(checkForDifferingL1Prescales && hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end()) {
              if(hltPathToL1PrescaleValueMap[triggerName] != l1Prescale) {
                std::cout << "WARN: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << " while previously we stored a prescale = " << hltPathToL1PrescaleValueMap[triggerName] << "; hlt prescale = " << hltPrescale << 
                  "; run = " << run << " ls = " << ls << "; psColumn = " << psProv.getRunInfo(run)->psColumn(ls) << "; l1 menu=" << psProv.getRunInfo(run)->l1Menu() <<
                  "; hlt menu=" << psProv.getRunInfo(run)->hltMenu() << "; trig mode=" << psProv.getRunInfo(run)->triggerMode() << std::endl;
              }
            }
            if(l1Prescale <= 0) {
              l1Seed = "L1_SingleEG45er2p5";
              l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
              if(checkForDifferingL1Prescales && hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end()) {
                if(hltPathToL1PrescaleValueMap[triggerName] != l1Prescale) {
                  std::cout << "WARN: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << " while previously we stored a prescale = " << hltPathToL1PrescaleValueMap[triggerName] << "; hlt prescale = " << hltPrescale << 
                    "; run = " << run << " ls = " << ls << "; psColumn = " << psProv.getRunInfo(run)->psColumn(ls) << "; l1 menu=" << psProv.getRunInfo(run)->l1Menu() <<
                    "; hlt menu=" << psProv.getRunInfo(run)->hltMenu() << "; trig mode=" << psProv.getRunInfo(run)->triggerMode() << std::endl;
                }
              }
            }
          }
          if(checkForDifferingL1Prescales && hltPathToL1PrescaleValueMap.find(triggerName) != hltPathToL1PrescaleValueMap.end())
            hltPathToL1PrescaleValueMap[triggerName] = l1Prescale;
        }
      }
      else if(analysisYearInt == 2016 && triggerName == "Photon175") {
        l1Prescale = 1;
        hltPrescale = 1;
      }
      else if((analysisYearInt == 2017 || analysisYearInt == 2018) && triggerName == "Photon200") {
        l1Prescale = 1;
        hltPrescale = 1;
      }
      if(l1Prescale <= 0 || hltPrescale <= 0) {
        const RunInfo* runInfo = psProv.getRunInfo(run);
        auto psCol = runInfo ? to_string(runInfo->psColumn(ls)) : "unknown";
        auto l1Menu = runInfo ? runInfo->l1Menu() : "unknown";
        auto hltMenu = runInfo ? runInfo->hltMenu() : "unknown";
        auto trigMode = runInfo ? runInfo->triggerMode() : "unknown";
        std::cout << "WARN: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << "; hlt prescale = " << hltPrescale << 
          "; run = " << run << " ls = " << ls << "; psColumn = " << psCol << "; l1 menu=" << l1Menu <<
          "; hlt menu=" << hltMenu << "; trig mode=" << trigMode << std::endl;
      }
      assert(l1Prescale > 0);
      assert(hltPrescale > 0);
      min_prescale = l1Prescale * hltPrescale;
      //std::cout << "INFO: Found L1 prescale = " << l1Prescale << " for L1 seed " << l1Seed << endl;
      //if(hltPhotonPt < 100)
      //  std::cout << "INFO: found prescale " << min_prescale << " for trigger name " << triggerName << " with hltPhotonPt = " << hltPhotonPt << ", for year: " << analysisYear << std::endl;
      //if(min_prescale <= 0)
      //  passTrigger = false;
    }
    //else {
    //  std::vector<std::string> triggerNames {"H_Photon22", "H_Photon30", "H_Photon36", "H_Photon50", "H_Photon75", "H_Photon90", "H_Photon120", "H_Photon175"};
    //  std::cout << "Photon with highest hltPt = " << hltPhotonPt << " did not pass any of the single photon triggers. run: " << run << " ls: " << ls << " event: " << event << std::endl;
    //  for(const auto& trigName : triggerNames)
    //    std::cout << "\t" << trigName << " --> " << readerTools_->ReadValueBranch<Float_t>(trigName) << std::endl;
    //  
    //}

    //// test output
    //std::cout << "Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;
    //std::cout << "Ele1_hltPhotonPt: " << hltPhoton1Pt << std::endl;
    ////std::cout << "Ele2_hltPhotonPt: " << hltPhoton2Pt << std::endl;
    ////std::cout << "Ele_hltPhotonPt: " << hltPhotonPt << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; trigger=" << triggerName << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "Passed Photon22 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon22")  > 0.1) << std::endl;
    //std::cout << "Passed Photon22 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon25")  > 0.1) << std::endl;
    //std::cout << "Passed Photon30 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon30")  > 0.1) << std::endl;
    //std::cout << "Passed Photon33 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon33")  > 0.1) << std::endl;
    //std::cout << "Passed Photon36 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon36")  > 0.1) << std::endl;
    //std::cout << "Passed Photon50 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon50")  > 0.1) << std::endl;
    //std::cout << "Passed Photon75 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon75")  > 0.1) << std::endl;
    //std::cout << "Passed Photon90 ? "  << (readerTools_->ReadValueBranch<Float_t>("H_Photon90")  > 0.1) << std::endl;
    //std::cout << "Passed Photon120 ? " << (readerTools_->ReadValueBranch<Float_t>("H_Photon120") > 0.1) << std::endl;
    //std::cout << "Passed Photon175 ? " << (readerTools_->ReadValueBranch<Float_t>("H_Photon175") > 0.1) << std::endl;
    //std::cout << "Passed Photon200 ? " << (readerTools_->ReadValueBranch<Float_t>("H_Photon200") > 0.1) << std::endl;
    //// test output
    //--------------------------------------------------------------------------
    // Fill HLT -- analysis trigger
    //--------------------------------------------------------------------------
    //bool passTrigger = readerTools_->ReadValueBranch<bool>("PassTrigger");
    //if ( isData() ) {
    //  passTrigger = false;
    //  std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
    //  if(current_file_name.find("SinglePhoton") != std::string::npos) {
    //    if(analysisYearInt==2016) {
    //      if (readerTools_->ReadValueBranch<float>("H_Photon175") == 1) // take events triggered by Photon175 only plus those triggered by Photon175 AND Ele27/Ele115
    //        passTrigger = true;
    //    }
    //    else {
    //      if (readerTools_->ReadValueBranch<float>("H_Photon200") == 1) // take events triggered by Photon200 only plus those triggered by Photon200 AND Ele35
    //        passTrigger = true;
    //    }
    //  }
    //  else if(current_file_name.find("SingleElectron") != std::string::npos) {
    //    if(analysisYearInt==2016) {
    //      if (readerTools_->ReadValueBranch<float>("H_Photon175") != 1 && 
    //          //(readerTools_->ReadValueBranch<float>("H_Ele27_WPTight") == 1 || readerTools_->ReadValueBranch<float>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele27 OR Ele115
    //          //(readerTools_->ReadValueBranch<float>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele115
    //        (readerTools_->ReadValueBranch<float>("H_Ele27_WPTight") == 1) )
    //          passTrigger = true;
    //    }
    //    else {
    //      if (readerTools_->ReadValueBranch<float>("H_Photon200") != 1 && 
    //          readerTools_->ReadValueBranch<float>("H_Ele35_WPTight") == 1 ) // take events triggered only by Ele35
    //        passTrigger = true;
    //    }
    //  }
    //  else if(analysisYearInt==2018) {
    //    if (readerTools_->ReadValueBranch<float>("H_Photon200") == 1 ||
    //        //readerTools_->ReadValueBranch<float>("H_Ele32_WPTight") == 1 ||
    //        //readerTools_->ReadValueBranch<float>("H_Ele115_CIdVT_GsfIdT") == 1) // take events triggered by Photon200 OR Ele32 OR Ele115
    //      readerTools_->ReadValueBranch<float>("H_Ele32_WPTight") == 1) // take events triggered by Photon200 OR Ele32
    //        passTrigger = true;
    //  }
    //}

    //--------------------------------------------------------------------------
    // Make this a QCD fake rate calculation
    //--------------------------------------------------------------------------
    float Ele1_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele1_SCEta");
    float Ele2_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele2_SCEta");
    float Ele3_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele3_SCEta");
    //float Ele1_PtHeep = readerTools_->ReadValueBranch<Float_t>("Ele1_PtHeep");
    //float Ele2_PtHeep = readerTools_->ReadValueBranch<Float_t>("Ele2_PtHeep");
    float Ele1_Phi = readerTools_->ReadValueBranch<Float_t>("Ele1_Phi");
    float Ele2_Phi = readerTools_->ReadValueBranch<Float_t>("Ele2_Phi");
    float Ele3_Phi = readerTools_->ReadValueBranch<Float_t>("Ele3_Phi");

    double fakeRate1 = -1;
    double fakeRate2 = -1;
    double fakeRate3 = -1;
    if(analysisYearInt < 2018) {
      fakeRate1 = qcdFR.GetFakeRate(Ele1_Pt, "2Jet", Ele1_SCEta);
      fakeRate2 = qcdFR.GetFakeRate(Ele2_Pt, "2Jet", Ele2_SCEta);
      fakeRate3 = qcdFR.GetFakeRate(Ele3_Pt, "2Jet", Ele3_SCEta);
    }
    else {
      fakeRate1 = qcdFR.GetFakeRate(Ele1_Pt, Ele1_SCEta, Ele1_Phi, run);
      fakeRate2 = qcdFR.GetFakeRate(Ele2_Pt, Ele2_SCEta, Ele2_Phi, run);
      fakeRate3 = qcdFR.GetFakeRate(Ele3_Pt, Ele3_SCEta, Ele3_Phi, run);
    }

    //--------------------------------------------------------------------------
    // Finally have the effective fake rate
    //--------------------------------------------------------------------------
    std::string idSuffix = "PassHEEPID";
    if(electronIDType == "EGMLoose")
      idSuffix = "PassEGMLooseID";
    int nEle_store = readerTools_->ReadValueBranch<Int_t>("nEle_store");
    int nEle_ptCut = readerTools_->ReadValueBranch<Int_t>("nEle_ptCut");
    bool Ele1_PassID = nEle_store >= 1 ? readerTools_->ReadValueBranch<Bool_t>("Ele1_"+idSuffix) : false;
    bool Ele2_PassID = nEle_store >= 2 ? readerTools_->ReadValueBranch<Bool_t>("Ele2_"+idSuffix) : false;
    bool Ele3_PassID = nEle_store >= 3 ? readerTools_->ReadValueBranch<Bool_t>("Ele3_"+idSuffix) : false;
    bool passIDRequirements = false;
    std::set<int> electronIndicesUsed;
    fakeRateEffective = 0;

    if(doSingleFR) {
      if(Ele1_PassID) {
        if(nEle_store >= 2 && !Ele2_PassID) {
          passIDRequirements = true;
          fakeRateEffective += fakeRate2/(1-fakeRate2);
          electronIndicesUsed.insert(1);
          electronIndicesUsed.insert(2);
        }
        if(nEle_store >= 3 && !Ele3_PassID) {
          passIDRequirements = true;
          fakeRateEffective += fakeRate3/(1-fakeRate3);
          electronIndicesUsed.insert(1);
          electronIndicesUsed.insert(3);
        }
      }
      if(Ele2_PassID) {
        if(!Ele1_PassID) {
          passIDRequirements = true;
          fakeRateEffective += fakeRate1/(1-fakeRate1);
          electronIndicesUsed.insert(2);
          electronIndicesUsed.insert(1);
        }
        if(nEle_store >= 3 && !Ele3_PassID) {
          passIDRequirements = true;
          fakeRateEffective += fakeRate3/(1-fakeRate3);
          electronIndicesUsed.insert(2);
          electronIndicesUsed.insert(3);
        }
      }
      if(Ele3_PassID) {
        if(!Ele1_PassID) {
          passIDRequirements = true;
          fakeRateEffective += fakeRate1/(1-fakeRate1);
          electronIndicesUsed.insert(3);
          electronIndicesUsed.insert(1);
        }
        if(!Ele2_PassID) {
          passIDRequirements = true;
          fakeRateEffective += fakeRate2/(1-fakeRate2);
          electronIndicesUsed.insert(3);
          electronIndicesUsed.insert(2);
        }
      }
    }
    else {
      fakeRateEffective = fakeRate1/(1-fakeRate1);
      passIDRequirements = !Ele1_PassID && !Ele2_PassID;
      if(nEle_store >= 2) {
        fakeRateEffective *= fakeRate2/(1-fakeRate2);
      }
      electronIndicesUsed.insert(1);
      electronIndicesUsed.insert(2);
      //XXX do we consider electron 3 here?
    }
    int indexOfHighestPtEleUsed = 1;
    int indexOfSecondHighestPtEleUsed = 2;
    if(passIDRequirements) {
      indexOfHighestPtEleUsed = *electronIndicesUsed.begin();
      indexOfSecondHighestPtEleUsed = *std::next(electronIndicesUsed.begin(), 1);
    }
    std::string ele1KeyName = "Ele"+to_string(indexOfHighestPtEleUsed);
    std::string ele2KeyName = "Ele"+to_string(indexOfSecondHighestPtEleUsed);

    //if(1-fakeRate1 <= 0)
    //{
    //  cout << "ERROR: Found fakeRate1: " << fakeRate1 << " for SCEta=" << Ele1_SCEta << " SCEt="
    //    << Ele1_SCEnergy/cosh(Ele1_SCEta) << "=" << Ele1_SCEnergy << "/" << 
    //    cosh(Ele1_SCEta) << endl;
    //}
    //FIXME: add error on fake rate as well
    // double eFakeRateEffective = fakeRateEffective * sqrt (  ( eFakeRate1 / fakeRate1 ) * ( eFakeRate1 / fakeRate1 ) +
    //					     ( eFakeRate2 / fakeRate2 ) * ( eFakeRate2 / fakeRate2 ) );
    double eFakeRateEffective = 0.0;

    eventWeight = fakeRateEffective * min_prescale * gen_weight;
    //cout << "event=" << event << " ls=" << ls << ", original weight = " << readerTools_->ReadValueBranch<Float_t>("Weight") << ", &fakeRateEffective=" << &fakeRateEffective << ", fakeRateEffective=" << fakeRateEffective << ", min_prescale=" << min_prescale << ", product = " << fakeRateEffective*min_prescale << endl;

    //--------------------------------------------------------------------------
    // What kind of event is this?
    //   - Barrel
    //   - Endcap 1 (eta < 2.0)
    //   - Endcap 2 (eta > 2.0) 
    //--------------------------------------------------------------------------

    bool ele1_isBarrel  = false;
    bool ele1_isEndcap1 = false;
    bool ele1_isEndcap2 = false;
    bool ele2_isBarrel  = false;
    bool ele2_isEndcap1 = false;
    bool ele2_isEndcap2 = false;

    // redefine eta/phi/pt
    Ele1_SCEta = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_SCEta");
    Ele2_SCEta = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_SCEta");
    //Ele1_PtHeep = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_PtHeep");
    //Ele2_PtHeep = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_PtHeep");
    Ele1_Pt = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_Pt");
    Ele2_Pt = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_Pt");
    Ele1_Phi = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_Phi");
    Ele2_Phi = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_Phi");
    float Ele1_Eta = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_Eta");
    float Ele2_Eta = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_Eta");

    if( fabs( Ele1_Eta  ) < eleEta_bar )        ele1_isBarrel  = true;
    if( fabs( Ele1_Eta  ) > eleEta_end1_min &&
        fabs( Ele1_Eta  ) < eleEta_end1_max )   ele1_isEndcap1 = true;
    if( fabs( Ele1_Eta  ) > eleEta_end2_min &&
        fabs( Ele1_Eta  ) < eleEta_end2_max )   ele1_isEndcap2 = true;

    if( fabs( Ele2_Eta  ) < eleEta_bar )        ele2_isBarrel  = true;
    if( fabs( Ele2_Eta  ) > eleEta_end1_min &&
        fabs( Ele2_Eta  ) < eleEta_end1_max )   ele2_isEndcap1 = true;
    if( fabs( Ele2_Eta  ) > eleEta_end2_min &&
        fabs( Ele2_Eta  ) < eleEta_end2_max )   ele2_isEndcap2 = true;

    bool ele1_isEndcap = ( ele1_isEndcap1 || ele1_isEndcap2 ) ;
    bool isEnd2_1 = ( fabs(Ele1_Eta) > eleEta_end2_min &&
        fabs(Ele1_Eta) < eleEta_end2_max ) ;
    bool ele2_isEndcap = ( ele2_isEndcap1 || ele2_isEndcap2 ) ;
    bool isEnd2_2 = ( fabs(Ele2_Eta) > eleEta_end2_min &&
        fabs(Ele2_Eta) < eleEta_end2_max ) ;

    bool isEBEB = ( ele1_isBarrel && ele2_isBarrel ) ;
    bool isEBEE = ( ( ele1_isBarrel && ele2_isEndcap ) ||
        ( ele2_isBarrel && ele1_isEndcap ) );
    bool isEEEE = ( ele1_isEndcap && ele2_isEndcap ) ;
    bool isEB   = ( isEBEB || isEBEE ) ;
    bool isEnd2End2 = isEnd2_1 && isEnd2_2;

    //--------------------------------------------------------------------------
    // Fill variables
    //--------------------------------------------------------------------------
    FillUserHist("EventCount"           , 1                   , 1   ); 
    FillUserHist("FakeRateEffective"    , fakeRateEffective); 
    FillUserHist("MinPrescale"          , min_prescale     ); 
    int trigMatchValToFill = -1;
    if(triggerMatch != "none") {
      if(triggerMatch == "online")
        trigMatchValToFill = 0;
      if(triggerMatch == "offline")
        trigMatchValToFill = 1;
    }
    FillUserHist("Trigger0OrOffline1Match", trigMatchValToFill);
    if(passTrigger) {
      FillUserHist("Trigger0OrOffline1Match_PassingTrigger", trigMatchValToFill);
      if(trigMatchValToFill == 1)
        FillUserHist("Trigger_OfflineMatch_PassingTrigger_Ele1Pt", leadElePt);
    }
    //if(fakeRateEffective * min_prescale != 1.0)
    //  std::cout << "!!!!!THIS EVENT HAD fakeRateEffective * min_prescale != 1.0: " << fakeRateEffective * min_prescale << std::endl;
    //if(fakeRateEffective >= 1) {
    //  std::cout << "!!!!!EVENT " << jentry << " HAD fakeRateEffective = " << fakeRateEffective << " and min_prescale = " << min_prescale << std::endl;
    //  std::cout.precision(0);
    //  std::cout << fixed <<  "Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;
    //  std::cout.precision(2);
    //  cout << "\tFound fakeRate1: " << fakeRate1 << " for SCEta=" << Ele1_SCEta << " SCEt="
    //    << Ele1_SCEnergy/cosh(Ele1_SCEta) << "=" << Ele1_SCEnergy << "/" << 
    //    cosh(Ele1_SCEta) << endl;
    //  cout << "\tFound fakeRate2: " << fakeRate2 << " for SCEta=" << Ele2_SCEta << " SCEt="
    //    << Ele2_SCEnergy/cosh(Ele2_SCEta) << "=" << Ele2_SCEnergy << "/" << 
    //    cosh(Ele2_SCEta) << endl;
    //}
    //if(min_prescale != 1.0) {
    //if(passIDRequirements) {
    //  std::cout << "!!!!!EVENT " << jentry << " HAD fakeRateEffective=" << fakeRateEffective << " and min_prescale= " << min_prescale << "; Ele1_MatchedHLTriggerObjectPt=" << readerTools_->ReadValueBranch<Float_t>("Ele1_MatchedHLTriggerObjectPt") <<  std::endl;
    //  std::cout.precision(0);
    //  std::cout << fixed <<  "Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;
    //  std::cout.precision(2);
    //  std::cout << "\tEle1_MatchedHLTriggerObjectPt: " << readerTools_->ReadValueBranch<Float_t>("Ele1_MatchedHLTriggerObjectPt") << std::endl;
    //  std::cout << "\tEle2_MatchedHLTriggerObjectPt: " << readerTools_->ReadValueBranch<Float_t>("Ele2_MatchedHLTriggerObjectPt") << std::endl;
    //  std::cout << "\tEle3_MatchedHLTriggerObjectPt: " << readerTools_->ReadValueBranch<Float_t>("Ele3_MatchedHLTriggerObjectPt") << std::endl;
    //  std::cout << "\tEle1_Pt: " << readerTools_->ReadValueBranch<Float_t>("Ele1_Pt") <<"; Ele1_Eta: " << readerTools_->ReadValueBranch<Float_t>("Ele1_Eta") <<  "; Ele1_Phi: " << readerTools_->ReadValueBranch<Float_t>("Ele1_Phi") << "; passID=" << Ele1_PassID << std::endl;
    //  std::cout << "\tEle2_Pt: " << readerTools_->ReadValueBranch<Float_t>("Ele2_Pt") <<"; Ele2_Eta: " << readerTools_->ReadValueBranch<Float_t>("Ele2_Eta") <<  "; Ele2_Phi: " << readerTools_->ReadValueBranch<Float_t>("Ele2_Phi") << "; passID=" << Ele2_PassID << std::endl;
    //  std::cout << "\tEle3_Pt: " << readerTools_->ReadValueBranch<Float_t>("Ele3_Pt") <<"; Ele3_Eta: " << readerTools_->ReadValueBranch<Float_t>("Ele3_Eta") <<  "; Ele3_Phi: " << readerTools_->ReadValueBranch<Float_t>("Ele3_Phi") << "; passID=" << Ele3_PassID << std::endl;
    //  if(analysisYearInt==2016) {
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon22")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon22"  ? "[selected] " : "") <<  "prescale of H_Photon22 = " << psProv.hltPrescale("HLT_Photon22_v", run, ls)  << std::endl;
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon30")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon30"  ? "[selected] " : "") <<  "prescale of H_Photon30 = " << psProv.hltPrescale("HLT_Photon30_v", run, ls)  << std::endl;
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon36")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon36"  ? "[selected] " : "") <<  "prescale of H_Photon36 = " << psProv.hltPrescale("HLT_Photon36_v", run, ls)  << std::endl;
    //  }
    //  if(analysisYearInt==2017) {
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon25")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon25"  ? "[selected] " : "") <<  "prescale of H_Photon25 = " << psProv.hltPrescale("HLT_Photon25_v", run, ls)  << std::endl;
    //  }
    //  if(analysisYearInt > 2016)
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon33")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon33"  ? "[selected] " : "") <<  "prescale of H_Photon33 = " << psProv.hltPrescale("HLT_Photon33_v", run, ls)  << std::endl;
    //  std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon50")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon50"  ? "[selected] " : "") <<  "prescale of H_Photon50 = " << psProv.hltPrescale("HLT_Photon50_v", run, ls)  << std::endl;
    //  std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon75")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon75"  ? "[selected] " : "") <<  "prescale of H_Photon75 = " << psProv.hltPrescale("HLT_Photon75_v", run, ls)  << std::endl;
    //  std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon90")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon90"  ? "[selected] " : "") <<  "prescale of H_Photon90 = " << psProv.hltPrescale("HLT_Photon90_v", run, ls)  << std::endl;
    //  std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon120")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon120" ? "[selected] " : "") <<  "prescale of H_Photon120= " << psProv.hltPrescale("HLT_Photon120_v", run, ls) << std::endl;
    //  if(analysisYearInt > 2016)
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon150")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon150" ? "[selected] " : "") <<  "prescale of H_Photon150= " << psProv.hltPrescale("HLT_Photon150_v", run, ls) << std::endl;
    //  std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon175")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon175" ? "[selected] " : "") <<  "prescale of H_Photon175= " << psProv.hltPrescale("HLT_Photon175_v", run, ls) << std::endl;
    //  if(analysisYearInt > 2016)
    //    std::cout << "\t" << (readerTools_->ReadValueBranch<Float_t>("H_Photon200")   > 0.1 ? "[passed] " : "") << (triggerName == "Photon200" ? "[selected] " : "") <<  "prescale of H_Photon200= " << psProv.hltPrescale("HLT_Photon200_v", run, ls) << std::endl;
    //}

    // reweighting
    fillVariableWithValue ( "Reweighting", 1, fakeRateEffective * min_prescale * gen_weight ) ; 

    // LHE cuts
    bool passLHECuts = true;
    // for jet-pt bin stitching
    //if(current_file_name.find("DYJetsToLL_M-50_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_newPMX_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_ext_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_ext_newPMX_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_ext1_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_ext1_newPMX_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_ext2_amcatnloFXFX") != std::string::npos ||
    //    current_file_name.find("DYJetsToLL_M-50_ext2_newPMX_amcatnloFXFX") != std::string::npos) {
    //  passLHECuts = false;
    //  if(readerTools_->ReadValueBranch<Float_t>("LHE_NpNLO") == 0) passLHECuts = true; // if zero jets, take it
    //  else if(readerTools_->ReadValueBranch<Float_t>("LHE_Vpt") < 50) passLHECuts = true; // if Z/gamma Pt < 50 GeV, take it
    //}
    // pt-binned sample stitching, Jul. 2022
    // see: https://cms-talk.web.cern.ch/t/bug-in-ul-pt-binned-dy-samples/11639
    if(current_file_name.find("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX") != std::string::npos) {
      passLHECuts = false;
      if(readerTools_->ReadValueBranch<Float_t>("LHE_Vpt") == 0)
        passLHECuts = true; // if Z/gamma Pt == 0 GeV, take it
    }
    // for stitching with MG LO HT-binned samples
    if(current_file_name.find("DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM") != std::string::npos) {
      passLHECuts = false;
      if(readerTools_->ReadValueBranch<Float_t>("LHE_HTIncoming") < 70)
        passLHECuts = true;
    }
    // for stitching with powhegMiNNLO mass-binned samples, Aug. 2023
    if(current_file_name.find("DYJetsToEE_M-50_massWgtFix_TuneCP5") != std::string::npos) {
      passLHECuts = false;
      // construct mass of dielectron system
      TLorentzVector e1, e2;
      float LHEEle1_Pt = readerTools_->ReadValueBranch<Float_t>("LHEElectron1_Pt");
      float LHEEle2_Pt = readerTools_->ReadValueBranch<Float_t>("LHEElectron2_Pt");
      float LHEEle1_Eta = readerTools_->ReadValueBranch<Float_t>("LHEElectron1_Eta");
      float LHEEle2_Eta = readerTools_->ReadValueBranch<Float_t>("LHEElectron2_Eta");
      float LHEEle1_Phi = readerTools_->ReadValueBranch<Float_t>("LHEElectron1_Phi");
      float LHEEle2_Phi = readerTools_->ReadValueBranch<Float_t>("LHEElectron2_Phi");
      e1.SetPtEtaPhiM ( LHEEle1_Pt, LHEEle1_Eta, LHEEle1_Phi, 0.0 );
      e2.SetPtEtaPhiM ( LHEEle2_Pt, LHEEle2_Eta, LHEEle2_Phi, 0.0 );
      float dielectron_mass = (e1 + e2).M();
      if(dielectron_mass > 50 && dielectron_mass < 100) // next mass bin starts at 100 GeV and previous ends at 50 GeV
        passLHECuts = true;
    }
    fillVariableWithValue("PassLHECuts",passLHECuts, fakeRateEffective * min_prescale * gen_weight);

    // JSON variable
    fillVariableWithValue ("PassJSON", readerTools_->ReadValueBranch<Bool_t>("PassJSON"), fakeRateEffective * min_prescale * gen_weight); 

    //--------------------------------------------------------------------------
    // Fill noise filters
    //--------------------------------------------------------------------------
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , readerTools_->ReadValueBranch<Bool_t>("PassGlobalSuperTightHalo2016Filter")     , fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassGoodVertices"                   , readerTools_->ReadValueBranch<Bool_t>("PassGoodVertices")                       , fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassHBHENoiseFilter"                , readerTools_->ReadValueBranch<Bool_t>("PassHBHENoiseFilter")                    , fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassHBHENoiseIsoFilter"             , readerTools_->ReadValueBranch<Bool_t>("PassHBHENoiseIsoFilter")                 , fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassBadEESupercrystalFilter"        , readerTools_->ReadValueBranch<Bool_t>("PassBadEESupercrystalFilter")            , fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassEcalDeadCellTrigPrim"           , readerTools_->ReadValueBranch<Bool_t>("PassEcalDeadCellTrigPrim")               , fakeRateEffective * min_prescale * gen_weight);
    // not recommended
    //fillVariableWithValue("PassChargedCandidateFilter"         , int(readerTools_->ReadValueB<Float_t>("PassChargedCandidateFilter")            == 1), fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassBadPFMuonFilter"                , readerTools_->ReadValueBranch<Bool_t>("PassBadPFMuonFilter")                    , fakeRateEffective * min_prescale * gen_weight);
    fillVariableWithValue("PassBadPFMuonDzFilter"              , readerTools_->ReadValueBranch<Bool_t>("PassBadPFMuonDzFilter")                  , fakeRateEffective * min_prescale * gen_weight);
    // EcalBadCalibV2 for 2017, 2018
    if(analysisYearInt > 2016)
      fillVariableWithValue("PassEcalBadCalibFilter"         , readerTools_->ReadValueBranch<Bool_t>("PassEcalBadCalibFilter")               , fakeRateEffective * min_prescale * gen_weight);
    else
      fillVariableWithValue("PassEcalBadCalibFilter"         , 1                                                                                , fakeRateEffective * min_prescale * gen_weight);

    // Electrons
    int PassNEle = 0;
    //if ( nEle_ptCut == 2 ) PassNEle = 1;
    if ( nEle_ptCut >= 2 ) PassNEle = 1;
    // we only look at events that have at least two loose electrons with pT > 35 GeV

    //FIXME fix the counting to avoid the below
    //// nEle_store is the number of loose eles with pt>10
    //if(nEle_store==3)
    //{
    //  // if 3 stored, 2 leads must pass Pt cut and third must fail
    //  if(Ele3_Pt<50 && Ele1_Pt>=50 && Ele2_Pt>=50)
    //    PassNEle=1;
    //}
    //else if(nEle_store==2)
    //{
    //  // if 2 stored, OK
    //  if(Ele1_Pt>=50 && Ele2_Pt>=50)
    //      PassNEle=1;
    //}

    double M_ej_avg = 0;
    double M_ej_min = 0;
    double M_ej_max = 0;
    double M_ej_asym = 0;

    // Muons
    int nMuon_ptCut = readerTools_->ReadValueBranch<Int_t>("nMuon_ptCut");
    double Muon1_Pt = readerTools_->ReadValueBranch<Float_t>("Muon1_Pt");
    double Muon1_Eta = readerTools_->ReadValueBranch<Float_t>("Muon1_Eta");
    double Muon1_Phi = readerTools_->ReadValueBranch<Float_t>("Muon1_Phi");
    int PassNMuon = 0;
    if ( nMuon_ptCut == 0 ) PassNMuon = 1;

    //std::cout << "INFO: Call PassHLT with value = " << passTrigger << " for trigger = " << triggerName << " (triggerMatch=" << triggerMatch << "), hltPhoton1Pt=" << hltPhoton1Pt << " and weight=" << fakeRateEffective * min_prescale * gen_weight << std::endl;
    fillVariableWithValue ( "PassHLT"                        , passTrigger             , fakeRateEffective * min_prescale * gen_weight ) ;
    fillVariableWithValue("nEle_hltMatched",-1, fakeRateEffective * min_prescale * gen_weight ) ;
    fillVariableWithValue("nJet_hltMatched",-1, fakeRateEffective * min_prescale * gen_weight ) ;


    // Jets
    int nJet_ptCut = readerTools_->ReadValueBranch<Int_t>("nJet_ptCut");
    int nJet_store = readerTools_->ReadValueBranch<Int_t>("nJet_store");
    float Jet1_Pt = readerTools_->ReadValueBranch<Float_t>("Jet1_Pt");
    float Jet1_Eta = readerTools_->ReadValueBranch<Float_t>("Jet1_Eta");
    float Jet1_Phi = readerTools_->ReadValueBranch<Float_t>("Jet1_Phi");
    float Jet2_Pt = readerTools_->ReadValueBranch<Float_t>("Jet2_Pt");
    float Jet2_Eta = readerTools_->ReadValueBranch<Float_t>("Jet2_Eta");
    float Jet2_Phi = readerTools_->ReadValueBranch<Float_t>("Jet2_Phi");
    float DR_Jet1Jet2 = readerTools_->ReadValueBranch<Float_t>("DR_Jet1Jet2");
    float Jet3_Pt = readerTools_->ReadValueBranch<Float_t>("Jet3_Pt");
    float Jet3_Eta = readerTools_->ReadValueBranch<Float_t>("Jet3_Eta");
    float Jet3_Phi = readerTools_->ReadValueBranch<Float_t>("Jet3_Phi");
    fillVariableWithValue(   "nJet"                          , nJet_ptCut      , fakeRateEffective * min_prescale * gen_weight ) ;
    if ( nJet_store >= 1 ) { 						                
      fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt         , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Jet1_Pt_skim"                  , Jet1_Pt         , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta        , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Jet1_Phi"                      , Jet1_Phi        , fakeRateEffective * min_prescale * gen_weight ) ;
    }
    else
      fillVariableWithValue( "Jet1_Pt"                       , -1.0         , fakeRateEffective * min_prescale * gen_weight ) ;

    //--------------------------------------------------------------------------
    // Recalculate some variables
    //--------------------------------------------------------------------------

    TLorentzVector e1, j1, e2, j2,j3, mu, met;
    TLorentzVector eejj, e1e2mu;
    TLorentzVector eej, ejj, ee;
    TLorentzVector e1j3, e2j3, j1j3, j2j3, j1j2, j1j2j3, eejjj;
    TLorentzVector e1j1, e1j2, e2j1, e2j2;

    e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
    e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
    j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
    j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
    mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
    float PFMET_Type1_Pt = readerTools_->ReadValueBranch<Float_t>("PFMET_Type1_Pt");
    float PFMET_Type1_Phi = readerTools_->ReadValueBranch<Float_t>("PFMET_Type1_Phi");
    met.SetPtEtaPhiM ( PFMET_Type1_Pt, 0.0, PFMET_Type1_Phi, 0.0 );

    eejj = e1 + e2 + j1 + j2 ; 
    eej  = e1 + e2 + j1;
    ejj  = e1 + j1 + j2;
    ee   = e1 + e2;
    j1j2 = j1 + j2;
    e1j1  = e1 + j1;
    e1j2  = e1 + j2;
    e2j1  = e2 + j1;
    e2j2  = e2 + j2;

    double min_DR_EleJet = 1000;
    double DR_Ele1Jet3 = 0;
    double DR_Ele2Jet3 = 0;
    double DR_Ele1Ele2 = e1.DeltaR( e2 );

    double M_eejj = eejj.M();
    double M_eej  = eej.M();
    double M_ejj  = ejj.M();
    double M_e1e2 = ee.M();
    double M_e1j1 = e1j1.M();
    double M_e1j2 = e1j2.M();
    double M_e2j1 = e2j1.M();
    double M_e2j2 = e2j2.M();
    double MT_Ele1MET = sqrt ( 2.0 * e1.Pt() * met.Pt() * ( 1.0 - cos ( met.DeltaPhi(e1))));

    double Pt_j1j2 = j1j2.Pt();
    double Pt_j1j3;
    double Pt_j2j3;
    double Pt_j1j2j3;
    double Pt_e1e2 = ee.Pt();

    double M_j1j3, M_j2j3;
    double M_e1j3, M_e2j3;
    double M_eejjj;

    double min_DeltaR_Zj = 999.;
    if ( ee.DeltaR ( j1 ) < min_DeltaR_Zj ) min_DeltaR_Zj = ee.DeltaR ( j1 );
    if ( ee.DeltaR ( j2 ) < min_DeltaR_Zj ) min_DeltaR_Zj = ee.DeltaR ( j2 );
    double DR_ZJ1 = ee.DeltaR ( j1 );
    double DR_ZJ2 = ee.DeltaR ( j2 );

    if ( nJet_ptCut > 2 ) { 
      j3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );

      e1j3 = e1 + j3;
      e2j3 = e2 + j3;
      j1j3 = j1 + j3;
      j2j3 = j2 + j3;
      j1j2j3 = j1 + j2 + j3;
      eejjj = e1 + e2 + j1 + j2 + j3;

      DR_Ele1Jet3 = e1.DeltaR( j3 );
      DR_Ele2Jet3 = e2.DeltaR( j3 );
      M_e1j3 = e1j3.M();
      M_e2j3 = e2j3.M();
      M_j1j3 = j1j3.M();
      M_j2j3 = j2j3.M();
      Pt_j1j3 = j1j3.Pt();
      Pt_j2j3 = j2j3.Pt();
      Pt_j1j2j3 = j1j2j3.Pt();
      M_eejjj = eejjj.M();
    }

    double sT_zjj = Pt_e1e2 + Jet1_Pt + Jet2_Pt;
    double sT_eejj = e1.Pt() + e2.Pt() + j1.Pt() + j2.Pt();

    // Muons
    fillVariableWithValue(   "PassNMuon"                     , PassNMuon               , fakeRateEffective * min_prescale * gen_weight ) ;

    // Electrons
    fillVariableWithValue(   "PassIDRequirements"            , passIDRequirements      , fakeRateEffective * min_prescale * gen_weight ) ;
    fillVariableWithValue(   "PassNEle"                      , PassNEle                , fakeRateEffective * min_prescale * gen_weight ) ;
    bool Ele1_PassHEEPEta = readerTools_->ReadValueBranch<Bool_t>(ele1KeyName+"_PassHEEPGsfEleSCEtaMultiRangeCut");
    bool Ele2_PassHEEPEta = readerTools_->ReadValueBranch<Bool_t>(ele2KeyName+"_PassHEEPGsfEleSCEtaMultiRangeCut");
    if ( nEle_store >= 1 ) {
      fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt            , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Ele1_Pt_skim"                  , Ele1_Pt            , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Ele1_Eta",            Ele1_Eta, fakeRateEffective * min_prescale * gen_weight  ) ;
      fillVariableWithValue( "Ele1_Phi",            Ele1_Phi, fakeRateEffective * min_prescale * gen_weight  ) ;
      //fillVariableWithValue( "Ele1_PassHEEPSCEtaCut"                      , Ele1_PassHEEPEta            , fakeRateEffective * min_prescale * gen_weight ) ;
    }										        
    //if ( nEle_store >= 2 ) { 							        
      fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt            , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Ele2_Pt_skim"                  , Ele2_Pt            , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Ele2_Eta",            Ele2_Eta, fakeRateEffective * min_prescale * gen_weight  ) ;
      fillVariableWithValue( "Ele2_Phi",            Ele2_Phi, fakeRateEffective * min_prescale * gen_weight  ) ;
      //fillVariableWithValue( "Ele2_PassHEEPSCEtaCut"                      , Ele2_PassHEEPEta            , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Pt_e1e2"                       , Pt_e1e2             , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Pt_e1e2_skim"                       , Pt_e1e2             , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "M_e1e2"                        , M_e1e2                  , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "M_e1e2_skim"                        , M_e1e2                  , fakeRateEffective * min_prescale * gen_weight ) ;
      //fillVariableWithValue( "M_e1e2_opt"                    , M_e1e2                  , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "M_e1e2_bkgCR"     , M_e1e2 , fakeRateEffective * min_prescale * gen_weight ) ;
    //}

    // DeltaR
    float DR_Ele1Jet1 = readerTools_->ReadValueBranch<Float_t>("DR_"+ele1KeyName+"Jet1");
    float DR_Ele2Jet1 = readerTools_->ReadValueBranch<Float_t>("DR_"+ele2KeyName+"Jet1");
    float DR_Ele1Jet2 = readerTools_->ReadValueBranch<Float_t>("DR_"+ele1KeyName+"Jet2");
    float DR_Ele2Jet2 = readerTools_->ReadValueBranch<Float_t>("DR_"+ele2KeyName+"Jet2");
    if ( nEle_store >= 2 && nJet_store >= 1) {
      fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1             , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1             , fakeRateEffective * min_prescale * gen_weight ) ;
      if(nJet_store >= 2) {
        fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2             , fakeRateEffective * min_prescale * gen_weight ) ;
        fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2             , fakeRateEffective * min_prescale * gen_weight ) ;
      }
    }

    if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
    if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
    if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
    if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
    if (nJet_ptCut > 2 ) {
      if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
      if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
    }

    //--------------------------------------------------------------------------
    // Calculate electron-jet pair mass values
    //--------------------------------------------------------------------------
    if ( nJet_store >= 2 ) { 
      fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt         , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Jet2_Pt_skim"                  , Jet2_Pt         , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta        , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "Jet2_Phi"                      , Jet2_Phi        , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2             , fakeRateEffective * min_prescale * gen_weight ) ;
      if ( nEle_store >= 2 && nJet_store >= 2) {
        if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
          M_ej_avg = (M_e1j1 + M_e2j2) / 2.0;
          if    ( M_e1j1 < M_e2j2 ) { M_ej_min = M_e1j1; M_ej_max = M_e2j2; }
          else                      { M_ej_min = M_e2j2; M_ej_max = M_e1j1; }
        }
        else { 
          M_ej_avg = (M_e1j2 + M_e2j1) / 2.0;
          if    ( M_e1j2 < M_e2j1 ) { M_ej_min = M_e1j2; M_ej_max = M_e2j1; }
          else                      { M_ej_min = M_e2j1; M_ej_max = M_e1j2; }
        }
        // compute Mej_asym
        M_ej_asym = (M_ej_max - M_ej_min)/(M_ej_max + M_ej_min);
      }
    }
    else
      fillVariableWithValue( "Jet2_Pt"                       , -1.0         , fakeRateEffective * min_prescale * gen_weight ) ;

    fillVariableWithValue( "Jet3_Pt"    , Jet3_Pt     , fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "Jet3_Eta"   , Jet3_Eta    , fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "Jet3_Phi"   , Jet3_Phi    , fakeRateEffective * min_prescale * gen_weight  ) ;
    //if ( nEle_store >= 2 && nJet_store >= 2) {
      // SIC recompute sT using PtHeep. FIXME: this is now being done in skims
      //sT_eejj = Ele1_PtHeep+Ele2_PtHeep+Jet1_Pt+Jet2_Pt;
      fillVariableWithValue( "sT_eejj"                       , sT_eejj                 , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "sT_eejj_bkgCR"    , sT_eejj , fakeRateEffective * min_prescale * gen_weight ) ;
      fillVariableWithValue( "sT_eejj_skim"    , sT_eejj , fakeRateEffective * min_prescale * gen_weight ) ;
    //}
    fillVariableWithValue( "M_e1j1", M_e1j1, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "M_e1j2", M_e1j2, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "M_e2j1", M_e2j1, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "M_e2j2", M_e2j2, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "Masym", M_ej_asym, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "MejMin", M_ej_min, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "MejMax", M_ej_max, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "Meejj", M_eejj, fakeRateEffective * min_prescale * gen_weight  ) ;
    //fillVariableWithValue( "PFMET_opt", PFMET_Type1_Pt, fakeRateEffective * min_prescale * gen_weight  ) ;
    fillVariableWithValue( "PFMET_Type1_Pt", PFMET_Type1_Pt,  fakeRateEffective * min_prescale * gen_weight ) ;
    fillVariableWithValue( "PFMET_Type1_Phi", PFMET_Type1_Phi, fakeRateEffective * min_prescale * gen_weight ) ;

    // Dummy variables
    fillVariableWithValue ("skim_selection", 1, fakeRateEffective * min_prescale * gen_weight ); 
    fillVariableWithValue ("preselection", 1, fakeRateEffective * min_prescale * gen_weight ); 
    fillVariableWithValue ("trainingSelection", 1, fakeRateEffective * min_prescale * gen_weight ); 


    //--------------------------------------------------------------------------
    // Fill bjet variables
    //--------------------------------------------------------------------------
    // require at least 1 b-tagged jet for TTBar control region
    // require zero b-tagged jets for DYJets control region
    // and then apply the b-tag scale factors (2-SF in the veto case)
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
    // for using event weights, we follow: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1c_Event_reweighting_using_scale

    double Jet1_btagDisc = readerTools_->ReadValueBranch<Float_t>("Jet1_btag"+btagAlgo);
    double Jet2_btagDisc = readerTools_->ReadValueBranch<Float_t>("Jet2_btag"+btagAlgo);
    double Jet3_btagDisc = readerTools_->ReadValueBranch<Float_t>("Jet3_btag"+btagAlgo);
    double Jet4_btagDisc = readerTools_->ReadValueBranch<Float_t>("Jet4_btag"+btagAlgo);
    double Jet5_btagDisc = readerTools_->ReadValueBranch<Float_t>("Jet5_btag"+btagAlgo);
    int Jet1_flavor = readerTools_->ReadValueBranch<Int_t>("Jet1_HadronFlavor");
    int Jet2_flavor = readerTools_->ReadValueBranch<Int_t>("Jet2_HadronFlavor");
    int Jet3_flavor = readerTools_->ReadValueBranch<Int_t>("Jet3_HadronFlavor");
    int Jet4_flavor = readerTools_->ReadValueBranch<Int_t>("Jet4_HadronFlavor");
    int Jet5_flavor = readerTools_->ReadValueBranch<Int_t>("Jet5_HadronFlavor");
    std::string sfSuffix = btagWP+btagAlgo;
    double Jet1_btagSF = Jet1_flavor==0 ? readerTools_->ReadValueBranch<Float_t>("Jet1_btagSF"+sfSuffix+"Incl") : readerTools_->ReadValueBranch<Float_t>("Jet1_btagSF"+sfSuffix+"Comb");
    double Jet2_btagSF = Jet2_flavor==0 ? readerTools_->ReadValueBranch<Float_t>("Jet2_btagSF"+sfSuffix+"Incl") : readerTools_->ReadValueBranch<Float_t>("Jet2_btagSF"+sfSuffix+"Comb");
    double Jet3_btagSF = Jet3_flavor==0 ? readerTools_->ReadValueBranch<Float_t>("Jet3_btagSF"+sfSuffix+"Incl") : readerTools_->ReadValueBranch<Float_t>("Jet3_btagSF"+sfSuffix+"Comb");
    double Jet4_btagSF = Jet4_flavor==0 ? readerTools_->ReadValueBranch<Float_t>("Jet4_btagSF"+sfSuffix+"Incl") : readerTools_->ReadValueBranch<Float_t>("Jet4_btagSF"+sfSuffix+"Comb");
    double Jet5_btagSF = Jet5_flavor==0 ? readerTools_->ReadValueBranch<Float_t>("Jet5_btagSF"+sfSuffix+"Incl") : readerTools_->ReadValueBranch<Float_t>("Jet5_btagSF"+sfSuffix+"Comb");

    float weightZeroBJets = 1.0;
    float weightAtLeastOneBJet = 1.0;
    float weightAtLeastTwoBJets = 1.0;
    float weightAtLeastTwoBJetsOneBtagBin = 1.0;
    float weightZeroBJetsBeyondLeadingTwo = 1.0;
    float weightAtLeastOneBJetBeyondLeadingTwo = 1.0;
    float weightZeroBJetsUpShift = 1.0;
    float weightAtLeastOneBJetUpShift = 1.0;
    float weightZeroBJetsDownShift = 1.0;
    float weightAtLeastOneBJetDownShift = 1.0;

    double discArray[5] = {Jet1_btagDisc, Jet2_btagDisc, Jet3_btagDisc, Jet4_btagDisc, Jet5_btagDisc};
    if(!isData())
    {
      float weightAtLeastTwoBJetsOneBtagBin = 0.0;
      // calculate and apply scale factors to MC only
      double sfArray[5] = {Jet1_btagSF, Jet2_btagSF, Jet3_btagSF, Jet4_btagSF, Jet5_btagSF };
      for(unsigned int iJet = 0; iJet < 5; ++iJet) {
        if (discArray[iJet] > btagCut) {
          weightZeroBJets*=(1-sfArray[iJet]);
          float tmpWeight = 1.0;
          for(unsigned int jJet = 0; jJet < 5; ++jJet) {
            if (discArray[jJet] > btagCut && jJet != iJet)
              tmpWeight*=(1-sfArray[jJet]);
          }
          weightAtLeastTwoBJetsOneBtagBin+=tmpWeight*sfArray[iJet];
        }
      }
      weightAtLeastOneBJet = 1 - weightZeroBJets;
      weightAtLeastTwoBJets = 1 - weightZeroBJets - weightAtLeastTwoBJetsOneBtagBin;
      //
      //if ( Jet3_btagDisc > btagCut ) weightZeroBJetsBeyondLeadingTwo*=(1-Jet3_btagSF);
      //if ( Jet4_btagDisc > btagCut ) weightZeroBJetsBeyondLeadingTwo*=(1-Jet4_btagSF);
      //if ( Jet5_btagDisc > btagCut ) weightZeroBJetsBeyondLeadingTwo*=(1-Jet5_btagSF);
      //weightAtLeastOneBJetBeyondLeadingTwo = 1 - weightZeroBJetsBeyondLeadingTwo;
      //
      //if ( Jet1_btagDisc > btagCut ) weightZeroBJetsUpShift*=(1-Jet1_btagSF);
      //if ( Jet2_btagDisc > btagCut ) weightZeroBJetsUpShift*=(1-Jet2_btagSF);
      //if ( Jet3_btagDisc > btagCut ) weightZeroBJetsUpShift*=(1-Jet3_btagSF);
      //if ( Jet4_btagDisc > btagCut ) weightZeroBJetsUpShift*=(1-Jet4_btagSF);
      //if ( Jet5_btagDisc > btagCut ) weightZeroBJetsUpShift*=(1-Jet5_btagSF);
      //weightAtLeastOneBJetUpShift = 1 - weightZeroBJetsUpShift;
      ////
      //if ( Jet1_btagDisc > btagCut ) weightZeroBJetsDownShift*=(1-Jet1_btagSF);
      //if ( Jet2_btagDisc > btagCut ) weightZeroBJetsDownShift*=(1-Jet2_btagSF);
      //if ( Jet3_btagDisc > btagCut ) weightZeroBJetsDownShift*=(1-Jet3_btagSF);
      //if ( Jet4_btagDisc > btagCut ) weightZeroBJetsDownShift*=(1-Jet4_btagSF);
      //if ( Jet5_btagDisc > btagCut ) weightZeroBJetsDownShift*=(1-Jet5_btagSF);
      //weightAtLeastOneBJetDownShift = 1 - weightZeroBJetsDownShift;
    }
    int nBJet_ptCut  = 0;
    int nBJet_ptCut_beyondLeadingTwo = 0;

    for(unsigned int iJet = 0; iJet < 5; ++iJet) {
      if (discArray[iJet] > btagCut)
        nBJet_ptCut++;
    }

    if ( Jet3_btagDisc > btagCut ) nBJet_ptCut_beyondLeadingTwo++;
    if ( Jet4_btagDisc > btagCut ) nBJet_ptCut_beyondLeadingTwo++;
    if ( Jet5_btagDisc > btagCut ) nBJet_ptCut_beyondLeadingTwo++;

    //std::cout << "INFO: weightAtLeastOneBJet=" << weightAtLeastOneBJet << "; while weightZeroBJets=" << weightZeroBJets << "; this event has " << nJet_store << " stored jets, with " << 
    // nBJet_medium_ptCut << " passing the medium Btag cut." << std::endl;

    //--------------------------------------------------------------------------
    // Fill mass-dependent meejj cuts
    //--------------------------------------------------------------------------
    if(doFinalSelections) {
      for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ) { 
        int lq_mass = LQ_MASS[i_lq_mass];
        sprintf(cut_name, "MeejjLQ%d", lq_mass );
        fillVariableWithValue( cut_name, M_eejj, gen_weight * pileup_weight  ) ;
      }
    }

    //--------------------------------------------------------------------------
    // Evaluate the cuts
    //--------------------------------------------------------------------------

    evaluateCuts();

    //--------------------------------------------------------------------------
    // Did we at least pass the noise filtes?
    //--------------------------------------------------------------------------

    bool passed_minimum = ( passedAllPreviousCuts("PassEcalBadCalibFilter") && passedCut ("PassEcalBadCalibFilter"));

    //--------------------------------------------------------------------------
    // Did we make it to the background control region?
    //--------------------------------------------------------------------------
    bool bkgControlRegion = passedAllPreviousCuts("M_e1e2_bkgCR") && passedCut("M_e1e2_bkgCR");

    //--------------------------------------------------------------------------
    // Did we pass preselection?
    //--------------------------------------------------------------------------

    bool passed_preselection = ( passedAllPreviousCuts("preselection") && passedCut ("preselection") );
    std::string preselectionCut = "preselection";

    //--------------------------------------------------------------------------
    // Are we in the region of interest?
    //--------------------------------------------------------------------------

    bool passed_region_of_interest = bool (
        passed_preselection && M_e1e2 > 170. && sT_eejj > 900.0
        );

    //--------------------------------------------------------------------------
    // Did we pass any final selection cuts?
    //--------------------------------------------------------------------------

    passed_vector.clear();
    if(doFinalSelections)
    {
      for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
        int lq_mass = LQ_MASS[i_lq_mass];
        //sprintf(cut_name, "M_e1e2_LQ%d", lq_mass );
        //sprintf(cut_name, "min_M_ej_LQ%d", lq_mass ); // this is actually the last cut in the cut file...!
        // TODO FIXME the right way; hack for now
        sprintf(cut_name, "BDTOutput_LQ%d", lq_mass ); // this is actually the last cut in the cut file...!
        bool decision = bool ( passedAllPreviousCuts(cut_name) && passedCut (cut_name));
        passed_vector.push_back (decision);
      }
    }

    //--------------------------------------------------------------------------
    // Fill plots with no selection applied
    //--------------------------------------------------------------------------

    FillUserHist( "PileupWeight"   , -1.0 );
    FillUserHist( "GeneratorWeight", -1.0 );

    bool passed_nele = ( passedAllPreviousCuts("PassNEle") && passedCut ("PassNEle") );
    if(passed_nele)
    {
      FillUserHist("FakeRateEffective_PassNEle"    , fakeRateEffective); 
      //if(fakeRateEffective < 4e-4)
      //{
      //  std::cout << "!!!!!THIS EVENT PASSING NELE_PTCUT HAD fakeRateEffective = " << fakeRateEffective << std::endl;
      //  std::cout << "\tFound fakeRate1: " << fakeRate1 << " for SCEta=" << Ele1_SCEta << " Pt=" << Ele1_Pt << " --> term1=" << fakeRate1/(1-fakeRate1) << std::endl;
      //  std::cout << "\tFound fakeRate2: " << fakeRate2 << " for SCEta=" << Ele2_SCEta << " Pt=" << Ele2_Pt << " --> term2=" << fakeRate2/(1-fakeRate2) << std::endl;
      //}
    }

    //--------------------------------------------------------------------------
    // Fill background control region plots
    //--------------------------------------------------------------------------
    if(bkgControlRegion) {
      FillUserHist("Mee_BkgControlRegion"	                ,    M_e1e2,    fakeRateEffective * min_prescale * gen_weight , "M_e1e2_bkgCR");
      if(nBJet_ptCut>=1) {
        FillUserHist( "Mee_BkgControlRegion_gteOneBtaggedJet"      , M_e1e2,  pileup_weight * gen_weight * weightAtLeastOneBJet, "preselection" ) ;
        FillUserHist2D("MeeVsNJet_BkgControlRegion_gteOneBtaggedJet", nJet_ptCut, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsST_BkgControlRegion_gteOneBtaggedJet", sT_eejj, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsSTjet_BkgControlRegion_gteOneBtaggedJet", Jet1_Pt + Jet2_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsSTlep_BkgControlRegion_gteOneBtaggedJet", Ele1_Pt + Ele2_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMejMin_BkgControlRegion_gteOneBtaggedJet", M_ej_min, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMejMax_BkgControlRegion_gteOneBtaggedJet", M_ej_max, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMeejj_BkgControlRegion_gteOneBtaggedJet", M_eejj, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMe1j1_BkgControlRegion_gteOneBtaggedJet", M_e1j1, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMe1j2_BkgControlRegion_gteOneBtaggedJet", M_e1j2, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMe2j1_BkgControlRegion_gteOneBtaggedJet", M_e2j1, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsMe2j2_BkgControlRegion_gteOneBtaggedJet", M_e2j2, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsEle1Pt_BkgControlRegion_gteOneBtaggedJet", Ele1_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsEle2Pt_BkgControlRegion_gteOneBtaggedJet", Ele2_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsJet1Pt_BkgControlRegion_gteOneBtaggedJet", Jet1_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsJet2Pt_BkgControlRegion_gteOneBtaggedJet", Jet2_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsJet3Pt_BkgControlRegion_gteOneBtaggedJet", Jet3_Pt, M_e1e2, pileup_weight * gen_weight);
        FillUserHist2D("MeeVsPFMETType1Pt_BkgControlRegion_gteOneBtaggedJet", PFMET_Type1_Pt, M_e1e2, pileup_weight * gen_weight);
        if(nBJet_ptCut>=2) {
          FillUserHist( "Mee_BkgControlRegion_gteTwoBtaggedJets"      , M_e1e2,  pileup_weight * gen_weight * weightAtLeastTwoBJets, "preselection" ) ;
          FillUserHist2D("MeeVsNJet_BkgControlRegion_gteTwoBtaggedJets", nJet_ptCut, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsST_BkgControlRegion_gteTwoBtaggedJets", sT_eejj, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsSTjet_BkgControlRegion_gteTwoBtaggedJets", Jet1_Pt + Jet2_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsSTlep_BkgControlRegion_gteTwoBtaggedJets", Ele1_Pt + Ele2_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMejMin_BkgControlRegion_gteTwoBtaggedJets", M_ej_min, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMejMax_BkgControlRegion_gteTwoBtaggedJets", M_ej_max, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMeejj_BkgControlRegion_gteTwoBtaggedJets", M_eejj, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMe1j1_BkgControlRegion_gteTwoBtaggedJets", M_e1j1, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMe1j2_BkgControlRegion_gteTwoBtaggedJets", M_e1j2, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMe2j1_BkgControlRegion_gteTwoBtaggedJets", M_e2j1, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsMe2j2_BkgControlRegion_gteTwoBtaggedJets", M_e2j2, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsEle1Pt_BkgControlRegion_gteTwoBtaggedJets", Ele1_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsEle2Pt_BkgControlRegion_gteTwoBtaggedJets", Ele2_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsJet1Pt_BkgControlRegion_gteTwoBtaggedJets", Jet1_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsJet2Pt_BkgControlRegion_gteTwoBtaggedJets", Jet2_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsJet3Pt_BkgControlRegion_gteTwoBtaggedJets", Jet3_Pt, M_e1e2, pileup_weight * gen_weight);
          FillUserHist2D("MeeVsPFMETType1Pt_BkgControlRegion_gteTwoBtaggedJets", PFMET_Type1_Pt, M_e1e2, pileup_weight * gen_weight);
        }
      }
      if      ( isEB   ) FillUserHist( "Mee_EB_BkgControlRegion"  , M_e1e2, fakeRateEffective * min_prescale * gen_weight, "M_e1e2_bkgCR" ); 
      if      ( isEBEB ) FillUserHist( "Mee_EBEB_BkgControlRegion", M_e1e2, fakeRateEffective * min_prescale * gen_weight, "M_e1e2_bkgCR" ); 
      else if ( isEBEE ) FillUserHist( "Mee_EBEE_BkgControlRegion", M_e1e2, fakeRateEffective * min_prescale * gen_weight, "M_e1e2_bkgCR" ); 
      else if ( isEEEE ) FillUserHist( "Mee_EEEE_BkgControlRegion", M_e1e2, fakeRateEffective * min_prescale * gen_weight, "M_e1e2_bkgCR" ); 
      if      ( isEnd2End2 ) FillUserHist( "Mee_End2End2_BkgControlRegion", M_e1e2, fakeRateEffective * min_prescale * gen_weight, "M_e1e2_bkgCR" ); 
      // scale factor dependence histos
      FillUserHist2D("MeeVsNJet_BkgControlRegion", nJet_ptCut, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsST_BkgControlRegion", sT_eejj, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsSTjet_BkgControlRegion", Jet1_Pt + Jet2_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsSTlep_BkgControlRegion", Ele1_Pt + Ele2_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMejMin_BkgControlRegion", M_ej_min, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMejMax_BkgControlRegion", M_ej_max, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMeejj_BkgControlRegion", M_eejj, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMe1j1_BkgControlRegion", M_e1j1, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMe1j2_BkgControlRegion", M_e1j2, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMe2j1_BkgControlRegion", M_e2j1, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsMe2j2_BkgControlRegion", M_e2j2, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsEle1Pt_BkgControlRegion", Ele1_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsEle2Pt_BkgControlRegion", Ele2_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsJet1Pt_BkgControlRegion", Jet1_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsJet2Pt_BkgControlRegion", Jet2_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsJet3Pt_BkgControlRegion", Jet3_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
      FillUserHist2D("MeeVsPFMETType1Pt_BkgControlRegion", PFMET_Type1_Pt, M_e1e2, fakeRateEffective * min_prescale * gen_weight);
    }

    //--------------------------------------------------------------------------
    // Fill preselection plots
    //--------------------------------------------------------------------------

    if ( passed_preselection ) {
      //XXX SIC DEBUG
      //if(hltPhotonPt < 100) {
      //  std::cout << "INFO: found prescale " << min_prescale << " for trigger name " << triggerName << " with hltPhotonPt = " << hltPhotonPt << ", for year: " << analysisYear << std::endl;
      //  std::cout << "The weight of this event: fakeRateEffective=" << fakeRateEffective << " x min_prescale=" << min_prescale << " = " << fakeRateEffective*min_prescale << std::endl;
      //  std::cout << "run: " << run << " ls: " << ls << " event: " << event << ", file = " << current_file_name << std::endl;
      //}
      //
      //std::cout << "We have electrons which look like (ptHEEP,eta,phi): (" << Ele1_PtHeep << "," << Ele1_Eta << "," << Ele1_Phi << "); FR = " << qcdFakeRate.GetFakeRate(Ele1_Eta,Ele2_PtHeep) << std::endl;
      //std::cout << "We have electrons which look like (pt,eta,phi): (" << Ele1_Pt << "," << Ele1_Eta << "," << Ele1_Phi << "); FR = " << qcdFakeRate.GetFakeRate(Ele1_Eta,Ele2_Pt) << std::endl;
      //std::cout << "Used fake rate=" << fakeRate1 << std::endl;
      ////float fakeRate2 = qcdFakeRate.GetFakeRate(Ele2_Eta,Ele2_PtHeep);
      //std::cout << "We have electrons which look like (pt,eta,phi): (" << Ele2_Pt << "," << Ele2_Eta << "," << Ele2_Phi << ")" << std::endl;
      //if(nEle_store > 2)
      //  std::cout << "We have electrons which look like (pt,eta,phi): (" << Ele3_Pt << "," << Ele3_Eta << "," << Ele3_Phi << ")" << std::endl;
      FillUserHist("Trigger0OrOffline1Match_PAS", trigMatchValToFill, min_prescale * gen_weight * fakeRateEffective);

      //--------------------------------------------------------------------------
      // Electron quality histograms (preselection)
      //--------------------------------------------------------------------------
      double Ele1_CorrIsolation        = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_CorrIsolation"); 
      double Ele1_DeltaEtaTrkSC        = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_DeltaEtaTrkSC"); 
      double Ele1_EcalIsolation        = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_EcalIsolation"); 
      double Ele1_HcalIsolation        = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_HcalIsolation"); 
      double Ele1_TrkIsolation         = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_TrkIsolation"); 
      bool Ele1_HasMatchedPhot         = readerTools_->ReadValueBranch<Bool_t>(ele1KeyName+"_HasMatchedPhot"); 
      double Ele1_HoE                  = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_HoE"); 
      double Ele1_LeadVtxDistXY        = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_LeadVtxDistXY"); 
      double Ele1_LeadVtxDistZ         = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_LeadVtxDistZ"); 
      int Ele1_MissingHits          = readerTools_->ReadValueBranch<Int_t>(ele1KeyName+"_MissingHits")        ; 
      double Ele1_SCEta                = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_SCEta");
      double Ele1_Full5x5SigmaIEtaIEta = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_Full5x5SigmaIEtaIEta");
      int Ele1_Charge               = readerTools_->ReadValueBranch<Int_t>(ele1KeyName+"_Charge");
      //double Ele1_PtHeep               = readerTools_->ReadValueBranch<Float_t>(ele1KeyName+"_PtHeep");

      FillUserHist("CorrIsolation_1stEle_PAS"         , Ele1_CorrIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("DeltaEtaTrkSC_1stEle_PAS"         , Ele1_DeltaEtaTrkSC                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("EcalIsolation_1stEle_PAS"         , Ele1_EcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("HcalIsolation_1stEle_PAS"         , Ele1_HcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("TrkIsolation_1stEle_PAS"          , Ele1_TrkIsolation                   , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("HasMatchedPhot_1stEle_PAS"        , Ele1_HasMatchedPhot                 , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("HoE_1stEle_PAS"                   , Ele1_HoE                            , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("LeadVtxDistXY_1stEle_PAS"         , Ele1_LeadVtxDistXY                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("LeadVtxDistZ_1stEle_PAS"          , Ele1_LeadVtxDistZ                   , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("MissingHits_1stEle_PAS"           , Ele1_MissingHits                    , min_prescale * gen_weight * fakeRateEffective   ); 
      if ( fabs(Ele1_SCEta) < eleEta_bar ) { 
        FillUserHist("Full5x5SigmaIEtaIEta_Barrel_1stEle_PAS", Ele1_Full5x5SigmaIEtaIEta                  , min_prescale * gen_weight * fakeRateEffective   ); 
      }
      else if ( fabs(Ele1_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
        FillUserHist("Full5x5SigmaIEtaIEta_Endcap_1stEle_PAS", Ele1_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      }

      double Ele2_CorrIsolation        = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_CorrIsolation"); 
      double Ele2_DeltaEtaTrkSC        = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_DeltaEtaTrkSC"); 
      double Ele2_EcalIsolation        = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_EcalIsolation"); 
      double Ele2_HcalIsolation        = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_HcalIsolation"); 
      double Ele2_TrkIsolation         = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_TrkIsolation"); 
      bool Ele2_HasMatchedPhot         = readerTools_->ReadValueBranch<Bool_t>(ele2KeyName+"_HasMatchedPhot")     ; 
      double Ele2_HoE                  = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_HoE"); 
      double Ele2_LeadVtxDistXY        = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_LeadVtxDistXY"); 
      double Ele2_LeadVtxDistZ         = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_LeadVtxDistZ"); 
      int Ele2_MissingHits          = readerTools_->ReadValueBranch<Int_t>(ele2KeyName+"_MissingHits")        ; 
      double Ele2_Full5x5SigmaIEtaIEta = readerTools_->ReadValueBranch<Float_t>(ele2KeyName+"_Full5x5SigmaIEtaIEta"); 
      int Ele2_Charge               = readerTools_->ReadValueBranch<Int_t>(ele2KeyName+"_Charge");

      FillUserHist("CorrIsolation_2ndEle_PAS"         , Ele2_CorrIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("DeltaEtaTrkSC_2ndEle_PAS"         , Ele2_DeltaEtaTrkSC                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("EcalIsolation_2ndEle_PAS"         , Ele2_EcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("HcalIsolation_2ndEle_PAS"         , Ele2_HcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("TrkIsolation_2ndEle_PAS"          , Ele2_TrkIsolation                   , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("HasMatchedPhot_2ndEle_PAS"        , Ele2_HasMatchedPhot                 , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("HoE_2ndEle_PAS"                   , Ele2_HoE                            , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("LeadVtxDistXY_2ndEle_PAS"         , Ele2_LeadVtxDistXY                  , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("LeadVtxDistZ_2ndEle_PAS"          , Ele2_LeadVtxDistZ                   , min_prescale * gen_weight * fakeRateEffective   ); 
      FillUserHist("MissingHits_2ndEle_PAS"           , Ele2_MissingHits                    , min_prescale * gen_weight * fakeRateEffective   ); 
      if ( fabs(Ele2_SCEta) < eleEta_bar ) { 
        FillUserHist("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PAS", Ele2_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      }
      else if ( fabs(Ele2_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
        FillUserHist("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PAS", Ele2_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      }

      //--------------------------------------------------------------------------
      // Preselection histograms
      //--------------------------------------------------------------------------
      double M_j1j2 = readerTools_->ReadValueBranch<Float_t>("M_j1j2");
      double Muon2_Pt = readerTools_->ReadValueBranch<Float_t>("Muon2_Pt");
      double Muon2_Eta = readerTools_->ReadValueBranch<Float_t>("Muon2_Eta");
      double Muon2_Phi = readerTools_->ReadValueBranch<Float_t>("Muon2_Phi");
      int nVertex = readerTools_->ReadValueBranch<Int_t>("nVertex");

      FillUserHist( "Ptj1j2_PAS"           , Pt_j1j2                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist( "Ptee_Minus_Ptj1j2_PAS", Pt_e1e2 - Pt_j1j2              , min_prescale * gen_weight * fakeRateEffective ) ;

      FillUserHist("minDR_EleJet_PAS"     , min_DR_EleJet                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("EleChargeSum_PAS"     , Ele1_Charge + Ele2_Charge, min_prescale * gen_weight * fakeRateEffective ) ;

      FillUserHist("nElectron_PAS"        , nEle_ptCut                    , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("nMuon_PAS"            , nMuon_ptCut                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("nJet_PAS"             , nJet_ptCut                 , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Pt1stEle_PAS"	   , Ele1_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Eta1stEle_PAS"	   , Ele1_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("SCEta1stEle_PAS"	   , Ele1_SCEta                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Phi1stEle_PAS"	   , Ele1_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Pt2ndEle_PAS"	   , Ele2_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Eta2ndEle_PAS"	   , Ele2_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("SCEta2ndEle_PAS"	   , Ele2_SCEta                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Phi2ndEle_PAS"	   , Ele2_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Charge1stEle_PAS"	   , Ele1_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Charge2ndEle_PAS"	   , Ele2_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("MET_PAS"              , PFMET_Type1_Pt                  , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("METPhi_PAS"	   , PFMET_Type1_Phi                 , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Pt1stJet_PAS"         , Jet1_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Pt2ndJet_PAS"         , Jet2_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Eta1stJet_PAS"        , Jet1_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Eta2ndJet_PAS"        , Jet2_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Phi1stJet_PAS"	   , Jet1_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Phi2ndJet_PAS"	   , Jet2_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
      if(nJet_store > 2) {
        FillUserHist("Pt3rdJet_PAS"         , Jet3_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta3rdJet_PAS"        , Jet3_Eta                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi3rdJet_PAS"     , Jet3_Phi                    , min_prescale * gen_weight * fakeRateEffective ) ;
      }
      FillUserHist("sTlep_PAS"            , Ele1_Pt + Ele2_Pt        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("sTjet_PAS"            , Jet1_Pt + Jet2_Pt  , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("sT_PAS"               , sT_eejj                            , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("sT_zjj_PAS"           , sT_zjj                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Mjj_PAS"		   , M_j1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Mee_PAS"		   , M_e1e2                             , min_prescale * gen_weight * fakeRateEffective ) ;
      //FillUserHist("MTenu_PAS"            , MT_Ele1MET                         , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Me1j1_PAS"            , M_e1j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Me1j2_PAS"            , M_e1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Me2j1_PAS"            , M_e2j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Me2j2_PAS"            , M_e2j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Ptee_PAS"             , Pt_e1e2                            , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("nVertex_PAS"          , nVertex                            , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2                        , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Meejj_PAS"            , M_eejj                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Meej_PAS"             , M_eej                              , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Mejj_PAS"             , M_ejj                              , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("minDR_ZJet_PAS"       , min_DeltaR_Zj                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_ZJet1_PAS"         , DR_ZJ1                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("DR_ZJet2_PAS"         , DR_ZJ2                             , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Mej_selected_avg_PAS" , M_ej_avg                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
      FillUserHist("Mej_selected_min_PAS" , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
      FillUserHist("Mej_selected_max_PAS" , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
      FillUserHist("Mej_minmax_PAS"       , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
      FillUserHist("Mej_minmax_PAS"       , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
      FillUserHist("Mej_asym_PAS"         , M_ej_asym                        , min_prescale * gen_weight * fakeRateEffective );	   
      // muon kinematics
      FillUserHist("Pt1stMuon_PAS"	   , Muon1_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Eta1stMuon_PAS"	   , Muon1_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Phi1stMuon_PAS"	   , Muon1_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Pt2ndMuon_PAS"	   , Muon2_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Eta2ndMuon_PAS"	   , Muon2_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
      FillUserHist("Phi2ndMuon_PAS"	   , Muon2_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;

      FillUserHist2D("MeeVsST_PAS" , M_e1e2, sT_eejj, min_prescale * gen_weight * fakeRateEffective ) ;	   
      // scale factor dependence histos
      //if ( nJet_ptCut == 2 )
      //  FillUserHist("Mee_NJetEq2_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if( nJet_ptCut == 3 )
      //  FillUserHist("Mee_NJetEq3_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if( nJet_ptCut == 4 )
      //  FillUserHist("Mee_NJetEq4_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if( nJet_ptCut == 5 )
      //  FillUserHist("Mee_NJetEq5_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if( nJet_ptCut == 6 )
      //  FillUserHist("Mee_NJetEq6_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if( nJet_ptCut == 7 )
      //  FillUserHist("Mee_NJetEq7_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      ////
      //if ( nJet_ptCut >= 3 )
      //  FillUserHist("Mee_NJetGeq3_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //if ( nJet_ptCut >= 4 )
      //  FillUserHist("Mee_NJetGeq4_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      ////
      //if (sT_eejj >= 300 && sT_eejj < 500)
      //  FillUserHist("Mee_sT300To500_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (sT_eejj >= 500 && sT_eejj < 750)
      //  FillUserHist("Mee_sT500To750_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (sT_eejj >= 750 && sT_eejj < 1250)
      //  FillUserHist("Mee_sT750To1250_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (sT_eejj >= 1250)
      //  FillUserHist("Mee_sT1250ToInf_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      ////
      //if (M_ej_min >= 100 && M_ej_min < 200)
      //  FillUserHist("Mee_MejMin100To200_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (M_ej_min >= 200 && M_ej_min < 300)
      //  FillUserHist("Mee_MejMin200To300_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (M_ej_min >= 300 && M_ej_min < 400)
      //  FillUserHist("Mee_MejMin300To400_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (M_ej_min >= 400 && M_ej_min < 500)
      //  FillUserHist("Mee_MejMin400To500_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (M_ej_min >= 500 && M_ej_min < 650)
      //  FillUserHist("Mee_MejMin500To650_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //else if (M_ej_min >= 650)
      //  FillUserHist("Mee_MejMin650ToInf_PAS", M_e1e2                         , min_prescale * gen_weight * fakeRateEffective );
      //-------------------------------------------------------------------------- 
      // no b tags
      //-------------------------------------------------------------------------- 
      if((isData() && nBJet_ptCut==0) || !isData()) {
        FillUserHist("nElectron_noBtaggedJets"        , nEle_ptCut                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("nMuon_noBtaggedJets"            , nMuon_ptCut                        , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("nJet_noBtaggedJets"             , nJet_ptCut                 , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt1stEle_noBtaggedJets"	   , Ele1_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta1stEle_noBtaggedJets"	   , Ele1_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi1stEle_noBtaggedJets"	   , Ele1_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt2ndEle_noBtaggedJets"	   , Ele2_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta2ndEle_noBtaggedJets"	   , Ele2_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi2ndEle_noBtaggedJets"	   , Ele2_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Charge1stEle_noBtaggedJets"	   , Ele1_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Charge2ndEle_noBtaggedJets"	   , Ele2_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("MET_noBtaggedJets"              , PFMET_Type1_Pt                  , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("METPhi_noBtaggedJets"	   , PFMET_Type1_Phi                 , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt1stJet_noBtaggedJets"         , Jet1_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt2ndJet_noBtaggedJets"         , Jet2_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta1stJet_noBtaggedJets"        , Jet1_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta2ndJet_noBtaggedJets"        , Jet2_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi1stJet_noBtaggedJets"	   , Jet1_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi2ndJet_noBtaggedJets"	   , Jet2_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sTlep_noBtaggedJets"            , Ele1_Pt + Ele2_Pt        , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sTjet_noBtaggedJets"            , Jet1_Pt + Jet2_Pt  , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sT_noBtaggedJets"               , sT_eejj                            , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sT_zjj_noBtaggedJets"           , sT_zjj                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mjj_noBtaggedJets"		   , M_j1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mee_noBtaggedJets"		   , M_e1e2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("MTenu_noBtaggedJets"            , MT_Ele1MET                         , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me1j1_noBtaggedJets"            , M_e1j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me1j2_noBtaggedJets"            , M_e1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me2j1_noBtaggedJets"            , M_e2j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me2j2_noBtaggedJets"            , M_e2j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mej_selected_min_noBtaggedJets" , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_selected_max_noBtaggedJets" , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_minmax_noBtaggedJets"       , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_minmax_noBtaggedJets"       , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_selected_avg_noBtaggedJets" , M_ej_avg                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mejj_noBtaggedJets"             , M_ejj                              , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Meej_noBtaggedJets"             , M_eej                              , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Meejj_noBtaggedJets"            , M_eejj                             , min_prescale * gen_weight * fakeRateEffective ) ;

        FillUserHist( "Mee_PAS_noBtaggedJets"      , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightZeroBJets ) ;
        if      ( isEBEB ) FillUserHist( "Mee_EBEB_PAS_noBtaggedJets"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightZeroBJets ); 
        else if ( isEBEE ) FillUserHist( "Mee_EBEE_PAS_noBtaggedJets"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightZeroBJets ); 
        else if ( isEEEE ) FillUserHist( "Mee_EEEE_PAS_noBtaggedJets"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightZeroBJets ); 

        //if (sT_eejj >= 300 && sT_eejj < 500)
        //  FillUserHist("Mee_sT300To500_PAS_noBtaggedJets", M_e1e2      , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (sT_eejj >= 500 && sT_eejj < 750)
        //  FillUserHist("Mee_sT500To750_PAS_noBtaggedJets", M_e1e2      , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (sT_eejj >= 750 && sT_eejj < 1250)
        //  FillUserHist("Mee_sT750To1250_PAS_noBtaggedJets", M_e1e2     , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (sT_eejj >= 1250)
        //  FillUserHist("Mee_sT1250ToInf_PAS_noBtaggedJets", M_e1e2     , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );

        //if (M_ej_min >= 100 && M_ej_min < 200)
        //  FillUserHist("Mee_MejMin100To200_PAS_noBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (M_ej_min >= 200 && M_ej_min < 300)
        //  FillUserHist("Mee_MejMin200To300_PAS_noBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (M_ej_min >= 300 && M_ej_min < 400)
        //  FillUserHist("Mee_MejMin300To400_PAS_noBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (M_ej_min >= 400 && M_ej_min < 500)
        //  FillUserHist("Mee_MejMin400To500_PAS_noBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (M_ej_min >= 500 && M_ej_min < 650)
        //  FillUserHist("Mee_MejMin500To650_PAS_noBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
        //else if (M_ej_min >= 650)
        //  FillUserHist("Mee_MejMin650ToInf_PAS_noBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightZeroBJets );
      }
      if(nBJet_ptCut>=1) {
        FillUserHist("nElectron_gteOneBtaggedJet"        , nEle_ptCut                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("nMuon_gteOneBtaggedJet"            , nMuon_ptCut                        , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("nJet_gteOneBtaggedJet"             , nJet_ptCut                 , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt1stEle_gteOneBtaggedJet"	   , Ele1_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta1stEle_gteOneBtaggedJet"	   , Ele1_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi1stEle_gteOneBtaggedJet"	   , Ele1_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt2ndEle_gteOneBtaggedJet"	   , Ele2_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta2ndEle_gteOneBtaggedJet"	   , Ele2_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi2ndEle_gteOneBtaggedJet"	   , Ele2_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Charge1stEle_gteOneBtaggedJet"	   , Ele1_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Charge2ndEle_gteOneBtaggedJet"	   , Ele2_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("MET_gteOneBtaggedJet"              , PFMET_Type1_Pt                  , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("METPhi_gteOneBtaggedJet"	   , PFMET_Type1_Phi                 , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt1stJet_gteOneBtaggedJet"         , Jet1_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt2ndJet_gteOneBtaggedJet"         , Jet2_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta1stJet_gteOneBtaggedJet"        , Jet1_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta2ndJet_gteOneBtaggedJet"        , Jet2_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi1stJet_gteOneBtaggedJet"	   , Jet1_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi2ndJet_gteOneBtaggedJet"	   , Jet2_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sTlep_gteOneBtaggedJet"            , Ele1_Pt + Ele2_Pt        , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sTjet_gteOneBtaggedJet"            , Jet1_Pt + Jet2_Pt  , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sT_gteOneBtaggedJet"               , sT_eejj                            , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sT_zjj_gteOneBtaggedJet"           , sT_zjj                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mjj_gteOneBtaggedJet"		   , M_j1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mee_gteOneBtaggedJet"		   , M_e1e2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("MTenu_gteOneBtaggedJet"            , MT_Ele1MET                         , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me1j1_gteOneBtaggedJet"            , M_e1j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me1j2_gteOneBtaggedJet"            , M_e1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me2j1_gteOneBtaggedJet"            , M_e2j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me2j2_gteOneBtaggedJet"            , M_e2j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mej_selected_min_gteOneBtaggedJet" , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_selected_max_gteOneBtaggedJet" , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_minmax_gteOneBtaggedJet"       , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_minmax_gteOneBtaggedJet"       , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_selected_avg_gteOneBtaggedJet" , M_ej_avg                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mejj_gteOneBtaggedJet"             , M_ejj                              , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Meej_gteOneBtaggedJet"             , M_eej                              , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Meejj_gteOneBtaggedJet"            , M_eejj                             , min_prescale * gen_weight * fakeRateEffective ) ;

        FillUserHist( "Mee_PAS_gteOneBtaggedJet"      , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet ) ;
        if      ( isEBEB ) FillUserHist( "Mee_EBEB_PAS_gteOneBtaggedJet"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet ); 
        else if ( isEBEE ) FillUserHist( "Mee_EBEE_PAS_gteOneBtaggedJet"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet ); 
        else if ( isEEEE ) FillUserHist( "Mee_EEEE_PAS_gteOneBtaggedJet"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet ); 

        //if (sT_eejj >= 300 && sT_eejj < 500)
        //  FillUserHist("Mee_sT300To500_PAS_gteOneBtaggedJet", M_e1e2      , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (sT_eejj >= 500 && sT_eejj < 750)
        //  FillUserHist("Mee_sT500To750_PAS_gteOneBtaggedJet", M_e1e2      , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (sT_eejj >= 750 && sT_eejj < 1250)
        //  FillUserHist("Mee_sT750To1250_PAS_gteOneBtaggedJet", M_e1e2     , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (sT_eejj >= 1250)
        //  FillUserHist("Mee_sT1250ToInf_PAS_gteOneBtaggedJet", M_e1e2     , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );

        //if (M_ej_min >= 100 && M_ej_min < 200)
        //  FillUserHist("Mee_MejMin100To200_PAS_gteOneBtaggedJet", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (M_ej_min >= 200 && M_ej_min < 300)
        //  FillUserHist("Mee_MejMin200To300_PAS_gteOneBtaggedJet", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (M_ej_min >= 300 && M_ej_min < 400)
        //  FillUserHist("Mee_MejMin300To400_PAS_gteOneBtaggedJet", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (M_ej_min >= 400 && M_ej_min < 500)
        //  FillUserHist("Mee_MejMin400To500_PAS_gteOneBtaggedJet", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (M_ej_min >= 500 && M_ej_min < 650)
        //  FillUserHist("Mee_MejMin500To650_PAS_gteOneBtaggedJet", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
        //else if (M_ej_min >= 650)
        //  FillUserHist("Mee_MejMin650ToInf_PAS_gteOneBtaggedJet", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastOneBJet );
      }
      if(nBJet_ptCut>=2) {
        FillUserHist("nElectron_gteTwoBtaggedJets"        , nEle_ptCut                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("nMuon_gteTwoBtaggedJets"            , nMuon_ptCut                        , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("nJet_gteTwoBtaggedJets"             , nJet_ptCut                 , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt1stEle_gteTwoBtaggedJets"	   , Ele1_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta1stEle_gteTwoBtaggedJets"	   , Ele1_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi1stEle_gteTwoBtaggedJets"	   , Ele1_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt2ndEle_gteTwoBtaggedJets"	   , Ele2_Pt                       , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta2ndEle_gteTwoBtaggedJets"	   , Ele2_Eta                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi2ndEle_gteTwoBtaggedJets"	   , Ele2_Phi                      , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Charge1stEle_gteTwoBtaggedJets"	   , Ele1_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Charge2ndEle_gteTwoBtaggedJets"	   , Ele2_Charge                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("MET_gteTwoBtaggedJets"              , PFMET_Type1_Pt                  , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("METPhi_gteTwoBtaggedJets"	   , PFMET_Type1_Phi                 , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt1stJet_gteTwoBtaggedJets"         , Jet1_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Pt2ndJet_gteTwoBtaggedJets"         , Jet2_Pt                    , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta1stJet_gteTwoBtaggedJets"        , Jet1_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Eta2ndJet_gteTwoBtaggedJets"        , Jet2_Eta                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi1stJet_gteTwoBtaggedJets"	   , Jet1_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Phi2ndJet_gteTwoBtaggedJets"	   , Jet2_Phi                   , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sTlep_gteTwoBtaggedJets"            , Ele1_Pt + Ele2_Pt        , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sTjet_gteTwoBtaggedJets"            , Jet1_Pt + Jet2_Pt  , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sT_gteTwoBtaggedJets"               , sT_eejj                            , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("sT_zjj_gteTwoBtaggedJets"           , sT_zjj                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mjj_gteTwoBtaggedJets"		   , M_j1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mee_gteTwoBtaggedJets"		   , M_e1e2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("MTenu_gteTwoBtaggedJets"            , MT_Ele1MET                         , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me1j1_gteTwoBtaggedJets"            , M_e1j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me1j2_gteTwoBtaggedJets"            , M_e1j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me2j1_gteTwoBtaggedJets"            , M_e2j1                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Me2j2_gteTwoBtaggedJets"            , M_e2j2                             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Mej_selected_min_gteTwoBtaggedJets" , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_selected_max_gteTwoBtaggedJets" , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_minmax_gteTwoBtaggedJets"       , M_ej_min                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_minmax_gteTwoBtaggedJets"       , M_ej_max                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mej_selected_avg_gteTwoBtaggedJets" , M_ej_avg                           , min_prescale * gen_weight * fakeRateEffective ) ;	   
        FillUserHist("Mejj_gteTwoBtaggedJets"             , M_ejj                              , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Meej_gteTwoBtaggedJets"             , M_eej                              , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist("Meejj_gteTwoBtaggedJets"            , M_eejj                             , min_prescale * gen_weight * fakeRateEffective ) ;

        FillUserHist( "Mee_PAS_gteTwoBtaggedJets"      , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets ) ;
        if      ( isEBEB ) FillUserHist( "Mee_EBEB_PAS_gteTwoBtaggedJets"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets ); 
        else if ( isEBEE ) FillUserHist( "Mee_EBEE_PAS_gteTwoBtaggedJets"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets ); 
        else if ( isEEEE ) FillUserHist( "Mee_EEEE_PAS_gteTwoBtaggedJets"		   , M_e1e2,  min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets ); 

        //if (sT_eejj >= 300 && sT_eejj < 500)
        //  FillUserHist("Mee_sT300To500_PAS_gteTwoBtaggedJets", M_e1e2      , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (sT_eejj >= 500 && sT_eejj < 750)
        //  FillUserHist("Mee_sT500To750_PAS_gteTwoBtaggedJets", M_e1e2      , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (sT_eejj >= 750 && sT_eejj < 1250)
        //  FillUserHist("Mee_sT750To1250_PAS_gteTwoBtaggedJets", M_e1e2     , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (sT_eejj >= 1250)
        //  FillUserHist("Mee_sT1250ToInf_PAS_gteTwoBtaggedJets", M_e1e2     , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );

        //if (M_ej_min >= 100 && M_ej_min < 200)
        //  FillUserHist("Mee_MejMin100To200_PAS_gteTwoBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (M_ej_min >= 200 && M_ej_min < 300)
        //  FillUserHist("Mee_MejMin200To300_PAS_gteTwoBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (M_ej_min >= 300 && M_ej_min < 400)
        //  FillUserHist("Mee_MejMin300To400_PAS_gteTwoBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (M_ej_min >= 400 && M_ej_min < 500)
        //  FillUserHist("Mee_MejMin400To500_PAS_gteTwoBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (M_ej_min >= 500 && M_ej_min < 650)
        //  FillUserHist("Mee_MejMin500To650_PAS_gteTwoBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
        //else if (M_ej_min >= 650)
        //  FillUserHist("Mee_MejMin650ToInf_PAS_gteTwoBtaggedJets", M_e1e2  , min_prescale * gen_weight * fakeRateEffective * weightAtLeastTwoBJets );
      }
      //FillUserTH3D("OptimizationCutSpace",sT_eejj,M_ej_min,M_e1e2, min_prescale * gen_weight * fakeRateEffective );

      //--------------------------------------------------------------------------
      // Mass-pairing histograms at preselection
      //--------------------------------------------------------------------------

      if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
        FillUserHist("Me1j_selected_PAS"   , M_e1j1,         min_prescale * gen_weight * fakeRateEffective );	   
        FillUserHist("Me2j_selected_PAS"   , M_e2j2,         min_prescale * gen_weight * fakeRateEffective );	   
        FillUserHist2D("Me1jVsMe2j_selected" , M_e1j1, M_e2j2, min_prescale * gen_weight * fakeRateEffective );
        FillUserHist2D("Me1jVsMe2j_rejected" , M_e1j2, M_e2j1, min_prescale * gen_weight * fakeRateEffective );
      }
      else {
        FillUserHist("Me1j_selected_PAS"   , M_e1j2,         min_prescale * gen_weight * fakeRateEffective );	   
        FillUserHist("Me2j_selected_PAS"   , M_e2j1,         min_prescale * gen_weight * fakeRateEffective );	   
        FillUserHist2D("Me1jVsMe2j_selected" , M_e1j2, M_e2j1, min_prescale * gen_weight * fakeRateEffective );
        FillUserHist2D("Me1jVsMe2j_rejected" , M_e1j1, M_e2j2, min_prescale * gen_weight * fakeRateEffective );
      }

      //--------------------------------------------------------------------------
      // Preselection + N(Jet) > 2 
      //--------------------------------------------------------------------------

      if ( nJet_ptCut > 2 ){ 
        FillUserHist( "M_e1j3_PAS"  , M_e1j3, min_prescale * gen_weight * fakeRateEffective  ); 
        FillUserHist( "M_e2j3_PAS"  , M_e2j3, min_prescale * gen_weight * fakeRateEffective  ); 
        FillUserHist( "M_j1j3_PAS"  , M_j1j3, min_prescale * gen_weight * fakeRateEffective  ); 
        FillUserHist( "M_j2j3_PAS"  , M_j2j3, min_prescale * gen_weight * fakeRateEffective  ); 
        FillUserHist( "M_eejjj_PAS" , M_eejjj,min_prescale * gen_weight * fakeRateEffective  ); 

        FillUserHist( "Ptj1j2j3_PAS"            , Pt_j1j2j3           , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist( "Ptj2j3_PAS"              , Pt_j2j3             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist( "Ptj1j3_PAS"              , Pt_j1j3             , min_prescale * gen_weight * fakeRateEffective ) ;
        FillUserHist( "Ptee_Minus_Ptj1j2j3_PAS" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * gen_weight * fakeRateEffective ) ;
      }

      //--------------------------------------------------------------------------
      // Preselection + event type (EBEB, EEEB, EEEE, etc)
      //--------------------------------------------------------------------------

      if      ( isEB   ) FillUserHist( "Mee_EB_PAS"  , M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      if      ( isEBEB ) FillUserHist( "Mee_EBEB_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      else if ( isEBEE ) FillUserHist( "Mee_EBEE_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      else if ( isEEEE ) FillUserHist( "Mee_EEEE_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      if      ( isEnd2End2 ) FillUserHist( "Mee_End2End2_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective ); 

      //--------------------------------------------------------------------------
      // Preselection + high ST plot
      //--------------------------------------------------------------------------

      //if ( sT_eejj > 445. ) FillUserHist( "Mee_PASandST445", M_e1e2, min_prescale * gen_weight * fakeRateEffective ) ;

      //--------------------------------------------------------------------------
      // High M(ee) plots
      //--------------------------------------------------------------------------

      ////if ( M_e1e2 > 100. ) FillUserHist("sT_PASandMee100"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 110. ) FillUserHist("sT_PASandMee110"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 120. ) FillUserHist("sT_PASandMee120"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 130. ) FillUserHist("sT_PASandMee130"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 140. ) FillUserHist("sT_PASandMee140"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 150. ) FillUserHist("sT_PASandMee150"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 160. ) FillUserHist("sT_PASandMee160"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 170. ) FillUserHist("sT_PASandMee170"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 180. ) FillUserHist("sT_PASandMee180"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 190. ) FillUserHist("sT_PASandMee190"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 
      //if ( M_e1e2 > 200. ) FillUserHist("sT_PASandMee200"   , sT_eejj , min_prescale * gen_weight * fakeRateEffective  ); 


      //if ( M_e1e2 > 100. ) { 

      ////  FillUserHist("CorrIsolation_1stEle_PASandMee100"         , Ele1_CorrIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("DeltaEtaTrkSC_1stEle_PASandMee100"         , Ele1_DeltaEtaTrkSC                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("EcalIsolation_1stEle_PASandMee100"         , Ele1_EcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("HcalIsolation_1stEle_PASandMee100"         , Ele1_HcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("TrkIsolation_1stEle_PASandMee100"          , Ele1_TrkIsolation                   , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("HasMatchedPhot_1stEle_PASandMee100"        , Ele1_HasMatchedPhot                 , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("HoE_1stEle_PASandMee100"                   , Ele1_HoE                            , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("LeadVtxDistXY_1stEle_PASandMee100"         , Ele1_LeadVtxDistXY                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("LeadVtxDistZ_1stEle_PASandMee100"          , Ele1_LeadVtxDistZ                   , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("MissingHits_1stEle_PASandMee100"           , Ele1_MissingHits                    , min_prescale * gen_weight * fakeRateEffective   ); 
      //  if ( fabs(Ele1_SCEta) < eleEta_bar ) { 
      ////    FillUserHist("Full5x5SigmaIEtaIEta_Barrel_1stEle_PASandMee100", Ele1_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      //  }
      //  else if ( fabs(Ele1_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
      ////    FillUserHist("Full5x5SigmaIEtaIEta_Endcap_1stEle_PASandMee100", Ele1_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      //  }

      ////  FillUserHist("CorrIsolation_2ndEle_PASandMee100"         , Ele2_CorrIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("DeltaEtaTrkSC_2ndEle_PASandMee100"         , Ele2_DeltaEtaTrkSC                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("EcalIsolation_2ndEle_PASandMee100"         , Ele2_EcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("HcalIsolation_2ndEle_PASandMee100"         , Ele2_HcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("TrkIsolation_2ndEle_PASandMee100"          , Ele2_TrkIsolation                   , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("HasMatchedPhot_2ndEle_PASandMee100"        , Ele2_HasMatchedPhot                 , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("HoE_2ndEle_PASandMee100"                   , Ele2_HoE                            , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("LeadVtxDistXY_2ndEle_PASandMee100"         , Ele2_LeadVtxDistXY                  , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("LeadVtxDistZ_2ndEle_PASandMee100"          , Ele2_LeadVtxDistZ                   , min_prescale * gen_weight * fakeRateEffective   ); 
      ////  FillUserHist("MissingHits_2ndEle_PASandMee100"           , Ele2_MissingHits                    , min_prescale * gen_weight * fakeRateEffective   ); 
      //  if ( fabs(Ele2_SCEta) < eleEta_bar ) { 
      ////    FillUserHist("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PASandMee100", Ele2_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      //  }
      //  else if ( fabs(Ele2_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
      ////    FillUserHist("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PASandMee100", Ele2_Full5x5SigmaIEtaIEta    , min_prescale * gen_weight * fakeRateEffective   ); 
      //  }

      ////  FillUserHist("Me1j1_PASandMee100"           , M_e1j1                              , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Ptee_PASandMee100"            , Pt_e1e2                             , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist2D("MeeVsST_PASandMee100" , M_e1e2, sT_eejj, min_prescale * gen_weight * fakeRateEffective ) ;	   
      ////  FillUserHist("sT_zjj_PASandMee100"          , sT_zjj                              , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("nVertex_PASandMee100"         , nVertex                             , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sT_PASandMee100"              , sT_eejj                             , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("EleChargeSum_PASandMee100"    , Ele1_Charge + Ele2_Charge , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("nJet_PASandMee100"            , nJet_ptCut                  , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTlep_PASandMee100"           , Ele1_Pt    + Ele2_Pt      , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTjet_PASandMee100"           , Jet1_Pt + Jet2_Pt   , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Mjj_PASandMee100"             , M_j1j2                              , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Pt1stEle_PASandMee100"        , Ele1_Pt                        , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Pt2ndEle_PASandMee100"        , Ele2_Pt                        , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Pt1stJet_PASandMee100"        , Jet1_Pt                     , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Pt2ndJet_PASandMee100"        , Jet2_Pt                     , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Mej_selected_avg_PASandMee100", M_ej_avg                            , min_prescale * gen_weight * fakeRateEffective ) ;

      ////  FillUserHist("sTfrac_Jet1_PASandMee100"     , Jet1_Pt / sT_eejj                       , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTfrac_Jet2_PASandMee100"     , Jet2_Pt / sT_eejj                       , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTfrac_Ele1_PASandMee100"     , Ele1_Pt / sT_eejj                          , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTfrac_Ele2_PASandMee100"     , Ele2_Pt / sT_eejj                          , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTfrac_Jet_PASandMee100"      , ( Jet1_Pt + Jet2_Pt ) / sT_eejj , min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("sTfrac_Ele_PASandMee100"      , ( Ele1_Pt + Ele2_Pt ) / sT_eejj       , min_prescale * gen_weight * fakeRateEffective ) ;

      ////  FillUserHist("Ptj1j2_PASandMee100"            , Pt_j1j2                        ,  min_prescale * gen_weight * fakeRateEffective ) ;
      ////  FillUserHist("Ptee_Minus_Ptj1j2_PASandMee100" , Pt_e1e2 - Pt_j1j2              ,  min_prescale * gen_weight * fakeRateEffective ) ;

      //  if ( nJet_ptCut > 2 ) { 	 
      ////    FillUserHist( "M_j1j3_PASandMee100" , M_j1j3, min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "M_j2j3_PASandMee100" , M_j2j3, min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "M_e1j3_PASandMee100" , M_e1j3, min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "M_e2j3_PASandMee100" , M_e2j3, min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "M_eejjj_PASandMee100", M_eejjj,min_prescale * gen_weight * fakeRateEffective ) ;

      ////    FillUserHist( "Ptj1j2j3_PASandMee100"            , Pt_j1j2j3           , min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "Ptj2j3_PASandMee100"              , Pt_j2j3             , min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "Ptj1j3_PASandMee100"              , Pt_j1j3             , min_prescale * gen_weight * fakeRateEffective ) ;
      ////    FillUserHist( "Ptee_Minus_Ptj1j2j3_PASandMee100" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * gen_weight * fakeRateEffective ) ;
      //  }
      //}

      //--------------------------------------------------------------------------
      // Preselection + M(ee) normalization region plots
      //--------------------------------------------------------------------------

      //if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
      //  FillUserHist("Mee_80_100_Preselection", M_e1e2, min_prescale * gen_weight * fakeRateEffective ) ;
      //  if      ( isEBEB ) FillUserHist( "Mee_EBEB_80_100_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //  else if ( isEBEE ) FillUserHist( "Mee_EBEE_80_100_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //  else if ( isEEEE ) FillUserHist( "Mee_EEEE_80_100_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //  if      ( isEB   ) FillUserHist( "Mee_EB_80_100_PAS"  , M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //}

      //if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
      //  FillUserHist("Mee_70_110_Preselection", M_e1e2, min_prescale * gen_weight * fakeRateEffective ) ;
      //  if ( sT_eejj > 600 ) 	 FillUserHist("Mee_70_110_ST600_Preselection", M_e1e2, min_prescale * gen_weight * fakeRateEffective ) ;
      //  if      ( isEBEB ) FillUserHist( "Mee_EBEB_70_110_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //  else if ( isEBEE ) FillUserHist( "Mee_EBEE_70_110_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //  else if ( isEEEE ) FillUserHist( "Mee_EEEE_70_110_PAS", M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //  if      ( isEB   ) FillUserHist( "Mee_EB_70_110_PAS"  , M_e1e2, min_prescale * gen_weight * fakeRateEffective  ); 
      //}

      //--------------------------------------------------------------------------
      // BDT preselection plots
      //--------------------------------------------------------------------------
      if(evaluateBDT) {
        for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
          int lq_mass = LQ_MASS[i_lq_mass];
          sprintf(cut_name, "BDTOutput_LQ%d", lq_mass );
          float bdtOutput = getVariableValue(cut_name);
          sprintf(cut_name, "BDTOutput_TrainRegion_LQ%d", lq_mass );
          FillUserHist(cut_name, bdtOutput, min_prescale * gen_weight * fakeRateEffective );
          sprintf(cut_name, "BDTOutput_noWeight_TrainRegion_LQ%d", lq_mass );
          FillUserHist(cut_name, bdtOutput );
        }
      }

      //--------------------------------------------------------------------------
      // Region of interest plots
      //-------------------------------------------------------------------------- 

      if ( do_roi_plots && passed_region_of_interest ) { 

        FillUserHist("CorrIsolation_1stEle_ROI"         , Ele1_CorrIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("DeltaEtaTrkSC_1stEle_ROI"         , Ele1_DeltaEtaTrkSC                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("EcalIsolation_1stEle_ROI"         , Ele1_EcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("HcalIsolation_1stEle_ROI"         , Ele1_HcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("TrkIsolation_1stEle_ROI"          , Ele1_TrkIsolation                   , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("HasMatchedPhot_1stEle_ROI"        , Ele1_HasMatchedPhot                 , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("HoE_1stEle_ROI"                   , Ele1_HoE                            , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("LeadVtxDistXY_1stEle_ROI"         , Ele1_LeadVtxDistXY                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("LeadVtxDistZ_1stEle_ROI"          , Ele1_LeadVtxDistZ                   , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("MissingHits_1stEle_ROI"           , Ele1_MissingHits                    , min_prescale * gen_weight * fakeRateEffective   ); 
        if ( fabs(Ele1_SCEta) < eleEta_bar ) { 
          FillUserHist("SigmaIEtaIEta_Barrel_1stEle_ROI", Ele1_Full5x5SigmaIEtaIEta                  , min_prescale * gen_weight * fakeRateEffective   ); 
        }
        else if ( fabs(Ele1_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) > eleEta_end2_max ){
          FillUserHist("SigmaIEtaIEta_Endcap_1stEle_ROI", Ele1_Full5x5SigmaIEtaIEta                  , min_prescale * gen_weight * fakeRateEffective   ); 
        }

        FillUserHist("CorrIsolation_2ndEle_ROI"         , Ele2_CorrIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("DeltaEtaTrkSC_2ndEle_ROI"         , Ele2_DeltaEtaTrkSC                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("EcalIsolation_2ndEle_ROI"         , Ele2_EcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("HcalIsolation_2ndEle_ROI"         , Ele2_HcalIsolation                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("TrkIsolation_2ndEle_ROI"          , Ele2_TrkIsolation                   , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("HasMatchedPhot_2ndEle_ROI"        , Ele2_HasMatchedPhot                 , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("HoE_2ndEle_ROI"                   , Ele2_HoE                            , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("LeadVtxDistXY_2ndEle_ROI"         , Ele2_LeadVtxDistXY                  , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("LeadVtxDistZ_2ndEle_ROI"          , Ele2_LeadVtxDistZ                   , min_prescale * gen_weight * fakeRateEffective   ); 
        FillUserHist("MissingHits_2ndEle_ROI"           , Ele2_MissingHits                    , min_prescale * gen_weight * fakeRateEffective   ); 
        if ( fabs(Ele2_Eta) < eleEta_bar ) { 
          FillUserHist("SigmaIEtaIEta_Barrel_2ndEle_ROI", Ele2_Full5x5SigmaIEtaIEta                  , min_prescale * gen_weight * fakeRateEffective   ); 
        }
        else if ( fabs(Ele2_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
          FillUserHist("SigmaIEtaIEta_Endcap_2ndEle_ROI", Ele2_Full5x5SigmaIEtaIEta                  , min_prescale * gen_weight * fakeRateEffective   ); 
        }

        FillUserHist("Me1j1_ROI"           , M_e1j1                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Ptee_ROI"            , Pt_e1e2                                        , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Eta1stJet_ROI"       , Jet1_Eta                               , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Eta2ndJet_ROI"       , Jet2_Eta                               , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Eta1stEle_ROI"	    , Ele1_SCEta                                  , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Eta2ndEle_ROI"	    , Ele2_SCEta                                  , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Phi1stJet_ROI"       , Jet1_Phi                               , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Phi2ndJet_ROI"       , Jet2_Phi                               , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Phi1stEle_ROI"	    , Ele1_Phi                                  , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Phi2ndEle_ROI"	    , Ele2_Phi                                  , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist2D("MeeVsST_ROI"         , M_e1e2                                , sT_eejj, min_prescale * gen_weight * fakeRateEffective );	   
        FillUserHist("Mee_ROI"		    , M_e1e2                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sT_zjj_ROI"          , sT_zjj                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("nVertex_ROI"         , nVertex                                        , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("EleChargeSum_ROI"    , Ele1_Charge + Ele2_Charge            , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("nJet_ROI"            , nJet_ptCut                             , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Mej_selected_avg_ROI", M_ej_avg                                       , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Meejj_ROI"           , M_eejj                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Meej_ROI"            , M_eej                                          , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Mejj_ROI"            , M_ejj                                          , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("minDR_ZJet_ROI"      , min_DeltaR_Zj                                  , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("DR_ZJet1_ROI"        , DR_ZJ1                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("DR_ZJet2_ROI"        , DR_ZJ2                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("MET_ROI"             , PFMET_Type1_Pt                              , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Mjj_ROI"             , M_j1j2                                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sT_ROI"              , sT_eejj                                        , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTlep_ROI"           , Ele1_Pt    + Ele2_Pt                 , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTjet_ROI"           , Jet1_Pt + Jet2_Pt              , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Pt1stEle_ROI"        , Ele1_Pt                                   , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Pt2ndEle_ROI"        , Ele2_Pt                                   , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Pt1stJet_ROI"        , Jet1_Pt                                , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Pt2ndJet_ROI"        , Jet2_Pt                                , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTfrac_Jet1_ROI"     , Jet1_Pt / sT_eejj                      , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTfrac_Jet2_ROI"     , Jet2_Pt / sT_eejj                      , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTfrac_Ele1_ROI"     , Ele1_Pt / sT_eejj                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTfrac_Ele2_ROI"     , Ele2_Pt / sT_eejj                         , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTfrac_Jet_ROI"      , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("sTfrac_Ele_ROI"      , ( Ele1_Pt + Ele2_Pt )       / sT_eejj, min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Ptj1j2_ROI"            , Pt_j1j2                                      , min_prescale * gen_weight * fakeRateEffective );
        FillUserHist("Ptee_Minus_Ptj1j2_ROI" , Pt_e1e2 - Pt_j1j2                            , min_prescale * gen_weight * fakeRateEffective );

        if ( nJet_ptCut > 2 ) { 	 
          FillUserHist( "M_e1j3_ROI" , M_e1j3, min_prescale * gen_weight * fakeRateEffective ) ;
          FillUserHist( "M_e2j3_ROI" , M_e2j3, min_prescale * gen_weight * fakeRateEffective ) ;
          FillUserHist( "M_j1j3_ROI" , M_j1j3, min_prescale * gen_weight * fakeRateEffective ) ;
          FillUserHist( "M_j2j3_ROI" , M_j2j3, min_prescale * gen_weight * fakeRateEffective ) ;
          FillUserHist( "M_eejjj_ROI", M_eejjj,min_prescale * gen_weight * fakeRateEffective ) ;
          FillUserHist( "Ptj1j2j3_ROI"            , Pt_j1j2j3           , min_prescale * gen_weight * fakeRateEffective );
          FillUserHist( "Ptj2j3_ROI"              , Pt_j2j3             , min_prescale * gen_weight * fakeRateEffective );
          FillUserHist( "Ptj1j3_ROI"              , Pt_j1j3             , min_prescale * gen_weight * fakeRateEffective );
          FillUserHist( "Ptee_Minus_Ptj1j2j3_ROI" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * gen_weight * fakeRateEffective );
        }
      }

      //-------------------------------------------------------------------------- 
      // Final selection plots
      //-------------------------------------------------------------------------- 

      if(doFinalSelections)
      {
        for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
          int lq_mass = LQ_MASS[i_lq_mass];
          bool pass = passed_vector[i_lq_mass];
          if ( !pass ) continue;
          //std::cout << fixed <<  "Passed Final Selection: Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;

          sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); FillUserHist ( plot_name, M_ej_avg          , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); FillUserHist ( plot_name, M_ej_min          , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); FillUserHist ( plot_name, M_ej_max          , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserHist ( plot_name, M_ej_min          , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserHist ( plot_name, M_ej_max          , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); FillUserHist ( plot_name, sT_eejj           , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); FillUserHist ( plot_name, M_e1e2            , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "DR_Ele1Jet1_LQ%d"            , lq_mass ); FillUserHist ( plot_name, DR_Ele1Jet1       , min_prescale * gen_weight * fakeRateEffective);
          sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); FillUserHist2D ( plot_name, M_ej_min, M_ej_max, min_prescale * gen_weight * fakeRateEffective);

          sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele1_CorrIsolation             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele1_DeltaEtaTrkSC             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "EcalIsolation_1stEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele1_EcalIsolation             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "HcalIsolation_1stEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele1_HcalIsolation             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "TrkIsolation_1stEle_LQ%d"       , lq_mass );   FillUserHist(plot_name,  Ele1_TrkIsolation              , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "HasMatchedPhot_1stEle_LQ%d"     , lq_mass );   FillUserHist(plot_name,  Ele1_HasMatchedPhot            , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "HoE_1stEle_LQ%d"                , lq_mass );   FillUserHist(plot_name,  Ele1_HoE                       , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistXY_1stEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele1_LeadVtxDistXY             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistZ_1stEle_LQ%d"       , lq_mass );   FillUserHist(plot_name,  Ele1_LeadVtxDistZ              , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "MissingHits_1stEle_LQ%d"        , lq_mass );   FillUserHist(plot_name,  Ele1_MissingHits               , min_prescale * gen_weight * fakeRateEffective ); 

          if ( fabs(Ele1_Eta) < eleEta_bar ) { 
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d", lq_mass ); FillUserHist( plot_name , Ele1_Full5x5SigmaIEtaIEta , min_prescale * gen_weight * fakeRateEffective    ); 
          }
          else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d", lq_mass ); FillUserHist( plot_name , Ele1_Full5x5SigmaIEtaIEta , min_prescale * gen_weight * fakeRateEffective    ); 
          }

          sprintf(plot_name, "CorrIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele2_CorrIsolation             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "DeltaEtaTrkSC_2ndEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele2_DeltaEtaTrkSC             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "EcalIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele2_EcalIsolation             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "HcalIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele2_HcalIsolation             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "TrkIsolation_2ndEle_LQ%d"       , lq_mass );   FillUserHist(plot_name,  Ele2_TrkIsolation              , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "HasMatchedPhot_2ndEle_LQ%d"     , lq_mass );   FillUserHist(plot_name,  Ele2_HasMatchedPhot            , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "HoE_2ndEle_LQ%d"                , lq_mass );   FillUserHist(plot_name,  Ele2_HoE                       , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistXY_2ndEle_LQ%d"      , lq_mass );   FillUserHist(plot_name,  Ele2_LeadVtxDistXY             , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistZ_2ndEle_LQ%d"       , lq_mass );   FillUserHist(plot_name,  Ele2_LeadVtxDistZ              , min_prescale * gen_weight * fakeRateEffective ); 
          sprintf(plot_name, "MissingHits_2ndEle_LQ%d"        , lq_mass );   FillUserHist(plot_name,  Ele2_MissingHits               , min_prescale * gen_weight * fakeRateEffective ); 

          if ( fabs(Ele2_Eta) < eleEta_bar ) { 
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ%d", lq_mass ); FillUserHist( plot_name , Ele2_Full5x5SigmaIEtaIEta , min_prescale * gen_weight * fakeRateEffective    ); 
          }
          else if ( fabs(Ele2_Eta) > eleEta_end2_min && fabs(Ele2_Eta) < eleEta_end2_max ){
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ%d", lq_mass ); FillUserHist( plot_name , Ele2_Full5x5SigmaIEtaIEta , min_prescale * gen_weight * fakeRateEffective    ); 
          }

          sprintf(plot_name, "Me1j1_LQ%d"             , lq_mass ); FillUserHist( plot_name , M_e1j1                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Me1j2_LQ%d"             , lq_mass ); FillUserHist( plot_name , M_e1j2                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Me2j1_LQ%d"             , lq_mass ); FillUserHist( plot_name , M_e2j1                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Me2j2_LQ%d"             , lq_mass ); FillUserHist( plot_name , M_e2j2                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Ptee_LQ%d"              , lq_mass ); FillUserHist( plot_name , Pt_e1e2                                        , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Eta1stJet_LQ%d"         , lq_mass ); FillUserHist( plot_name , Jet1_Eta                               , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Eta2ndJet_LQ%d"         , lq_mass ); FillUserHist( plot_name , Jet2_Eta                               , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Eta1stEle_LQ%d"         , lq_mass ); FillUserHist( plot_name , Ele1_Eta                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Eta2ndEle_LQ%d"         , lq_mass ); FillUserHist( plot_name , Ele2_Eta                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Phi1stJet_LQ%d"         , lq_mass ); FillUserHist( plot_name , Jet1_Phi                               , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Phi2ndJet_LQ%d"         , lq_mass ); FillUserHist( plot_name , Jet2_Phi                               , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Phi1stEle_LQ%d"         , lq_mass ); FillUserHist( plot_name , Ele1_Phi                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Phi2ndEle_LQ%d"         , lq_mass ); FillUserHist( plot_name , Ele2_Phi                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "MeeVsST_LQ%d"           , lq_mass ); FillUserHist2D( plot_name , M_e1e2, sT_eejj                                , min_prescale * gen_weight * fakeRateEffective );	   
          sprintf(plot_name, "sT_zjj_LQ%d"            , lq_mass ); FillUserHist( plot_name , sT_zjj                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "nVertex_LQ%d"           , lq_mass ); FillUserHist( plot_name , nVertex                                        , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "nJet_LQ%d"              , lq_mass ); FillUserHist( plot_name , nJet_ptCut                             , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "EleChargeSum_LQ%d"      , lq_mass ); FillUserHist( plot_name , Ele1_Charge + Ele2_Charge            , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Meejj_LQ%d"             , lq_mass ); FillUserHist( plot_name , M_eejj                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Meej_LQ%d"              , lq_mass ); FillUserHist( plot_name , M_eej                                          , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Mejj_LQ%d"              , lq_mass ); FillUserHist( plot_name , M_ejj                                          , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Mjj_LQ%d"               , lq_mass ); FillUserHist( plot_name , M_j1j2                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "minDR_ZJet_LQ%d"        , lq_mass ); FillUserHist( plot_name , min_DeltaR_Zj                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "DR_ZJet1_LQ%d"          , lq_mass ); FillUserHist( plot_name , DR_ZJ1                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "DR_ZJet2_LQ%d"          , lq_mass ); FillUserHist( plot_name , DR_ZJ2                                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "MET_LQ%d"               , lq_mass ); FillUserHist( plot_name , PFMET_Type1_Pt                              , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTlep_LQ%d"             , lq_mass ); FillUserHist( plot_name , Ele1_Pt + Ele2_Pt                    , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTjet_LQ%d"             , lq_mass ); FillUserHist( plot_name , Jet1_Pt + Jet2_Pt              , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Pt1stEle_LQ%d"          , lq_mass ); FillUserHist( plot_name , Ele1_Pt                                   , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Pt2ndEle_LQ%d"          , lq_mass ); FillUserHist( plot_name , Ele2_Pt                                   , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Pt1stJet_LQ%d"          , lq_mass ); FillUserHist( plot_name , Jet1_Pt                                , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Pt2ndJet_LQ%d"          , lq_mass ); FillUserHist( plot_name , Jet2_Pt                                , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Jet1_LQ%d"       , lq_mass ); FillUserHist( plot_name , Jet1_Pt / sT_eejj                      , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Jet2_LQ%d"       , lq_mass ); FillUserHist( plot_name , Jet2_Pt / sT_eejj                      , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Ele1_LQ%d"       , lq_mass ); FillUserHist( plot_name , Ele1_Pt / sT_eejj                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Ele2_LQ%d"       , lq_mass ); FillUserHist( plot_name , Ele2_Pt / sT_eejj                         , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Jet_LQ%d"        , lq_mass ); FillUserHist( plot_name , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Ele_LQ%d"        , lq_mass ); FillUserHist( plot_name , ( Ele1_Pt + Ele2_Pt ) / sT_eejj      , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Ptj1j2_LQ%d"            , lq_mass ); FillUserHist( plot_name , Pt_j1j2                                        , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Ptee_Minus_Ptj1j2_LQ%d" , lq_mass ); FillUserHist( plot_name , Pt_e1e2 - Pt_j1j2                              , min_prescale * gen_weight * fakeRateEffective );
          // muon kinematics
          sprintf(plot_name, "Pt1stMuon_LQ%d"          , lq_mass ); FillUserHist( plot_name , Muon1_Pt                                   , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Pt2ndMuon_LQ%d"          , lq_mass ); FillUserHist( plot_name , Muon2_Pt                                   , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Eta1stMuon_LQ%d"         , lq_mass ); FillUserHist( plot_name , Muon1_Eta                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Eta2ndMuon_LQ%d"         , lq_mass ); FillUserHist( plot_name , Muon2_Eta                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Phi1stMuon_LQ%d"         , lq_mass ); FillUserHist( plot_name , Muon1_Phi                                  , min_prescale * gen_weight * fakeRateEffective );
          sprintf(plot_name, "Phi2ndMuon_LQ%d"         , lq_mass ); FillUserHist( plot_name , Muon2_Phi                                  , min_prescale * gen_weight * fakeRateEffective );

        } // End final selection

        if( hasCut("sT_eejj_LQ300") && passedCut("sT_eejj_LQ300") && passedCut("min_M_ej_LQ300"))
          FillUserHist("Mee_70_110_LQ300", M_e1e2 , min_prescale * gen_weight * fakeRateEffective );
        if( hasCut("sT_eejj_LQ600") && passedCut("sT_eejj_LQ600") && passedCut("min_M_ej_LQ600"))
          FillUserHist("Mee_70_110_LQ600", M_e1e2 , min_prescale * gen_weight * fakeRateEffective );
        if( hasCut("sT_eejj_LQ800") && passedCut("sT_eejj_LQ800") && passedCut("min_M_ej_LQ800"))
          FillUserHist("Mee_70_110_LQ800", M_e1e2 , min_prescale * gen_weight * fakeRateEffective );
        if( hasCut("sT_eejj_LQ900") && passedCut("sT_eejj_LQ900") && passedCut("min_M_ej_LQ900"))
          FillUserHist("Mee_70_110_LQ900", M_e1e2 , min_prescale * gen_weight * fakeRateEffective );
        if( hasCut("sT_eejj_LQ1000") && passedCut("sT_eejj_LQ1000") && passedCut("min_M_ej_LQ1000"))
          FillUserHist("Mee_70_110_LQ1000", M_e1e2 , min_prescale * gen_weight * fakeRateEffective );

      } // End do final selections

    } // End preselection 
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

