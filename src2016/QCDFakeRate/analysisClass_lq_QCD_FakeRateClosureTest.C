#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

//#include <TFile.h>
//#include <TH1F.h>
//#include <TF1.h>
//#include <TFitResultPtr.h>
//#include <TFitResult.h>
//
//#include "qcdFitter.h"
// for scale factors
#include "ElectronScaleFactors.C"
// 2016 trigger efficiency
#include "TriggerEfficiency2016.h"
// for prescales
#include "Run2PhotonTriggerPrescales.h"
#include "PrescaleProvider.h"
#include <set>
//#include "HistoReader.h"
#include "QCDFakeRate.h"

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}
   
void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillAllPreviousCuts              (  true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Get pre-cut values
   //--------------------------------------------------------------------------

   // eta boundaries
   double eleEta_bar_max = getPreCutValue1("eleEta_bar");
   double eleEta_end_min = getPreCutValue1("eleEta_end1");
   double eleEta_end_max = getPreCutValue2("eleEta_end2");


   // eta boundaries
   double eleEta_bar            = getPreCutValue1("eleEta_bar");
   double eleEta_end1_min       = getPreCutValue1("eleEta_end1");
   double eleEta_end1_max       = getPreCutValue2("eleEta_end1");
   double eleEta_end2_min       = getPreCutValue1("eleEta_end2");
   double eleEta_end2_max       = getPreCutValue2("eleEta_end2");


   // override the fake rate?
   double fakeRate_override = getPreCutValue1("fakeRate_override");
   bool override_fakeRate = ( fakeRate_override > 0.0 );
   
   //prescales
   PrescaleProvider psProv("/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ul-analysis-inputs/prescales/2016/triggerData2016");
   Run2PhotonTriggerPrescales run2PhotonTriggerPrescales;
   
   //--------------------------------------------------------------------------
   // Analysis year
   //--------------------------------------------------------------------------
   std::string getAnalysisYear = getPreCutString1("AnalysisYear");
   int analysisYear;
   std::string analysisYearStr;
   if (getAnalysisYear.find("pre") != string::npos ){
      analysisYear = 2016;
      analysisYearStr = "2016preVFP";
   }
   else if (getAnalysisYear.find("post") != string::npos){
      analysisYear = 2016;
      analysisYearStr = "2016postVFP";
   }
   else if (getAnalysisYear.find("17") != string::npos){
      analysisYear = 2017;
      analysisYearStr = "2017";
   }
   else if (getAnalysisYear.find("18") != string::npos){
      analysisYear = 2018;
      analysisYearStr = "2018";
   }
   else{
      std::cout<<"ERROR: cannot determine analysis year from cutfile"<<std::endl;
   }

   //--------------------------------------------------------------------------
   // QCD Fake Rate loading part
   //--------------------------------------------------------------------------
   std::string qcdFileName = getPreCutString1("QCDFakeRateFileName");
   //HistoReader qcdFakeRateReader(qcdFileName,"fr2D_1Jet_TrkIsoHEEP7vsHLTPt_PAS","fr2D_1Jet_TrkIsoHEEP7vsHLTPt_PAS",true,false);
   std::vector<std::string> regionVec;
   if(analysisYear != 2018) regionVec = {"1Jet_TrkIsoHEEP7vsHLTPt_PAS"};
   else{
     regionVec = {
       "1Jet_TrkIsoHEEP7vsHLTPt_pre319077",
       "1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077",
       "1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077",
     };
   } 
   QCDFakeRate qcdFakeRateReader(qcdFileName, "fr2D_", regionVec, true);

   //--------------------------------------------------------------------------
   // reco scale factors
   //--------------------------------------------------------------------------
   std::string recoSFFileName = getPreCutString1("RecoSFFileName");
   std::unique_ptr<HistoReader> recoScaleFactorReader = std::unique_ptr<HistoReader>(new HistoReader(recoSFFileName,"EGamma_SF2D","EGamma_SF2D",true,false));

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------
   //Since we have two electrons, we have three possible combinations of eta regions they could be in:
   //both in the barrel (BB), both in the endcap (EE), or one in each (BE).
   //I also want hists for the total from all regions, which is what the i==3 case is for
   std::string region;
   for (int i=0; i<4; i++){
     if (i==0){region = "_BB";}
     else if (i==1){region = "_BE";}
     else if (i==2){region = "_EE";}
     else {region="";}
    
     CreateUserTH1D( "nElectron_PAS"+region          ,    5   , -0.5    , 4.5      );
     CreateUserTH1D( "nMuon_PAS"+region              ,    5   , -0.5    , 4.5      );
     CreateUserTH1D( "nJet_PAS"+region               ,    10  , -0.5    , 9.5      );
     CreateUserTH1D( "Pt1stEle_PAS"+region           , 	  100 , 0       , 1000     ); 
     CreateUserTH1D( "Eta1stEle_PAS"+region          ,    100 , -5      , 5	   ); 
     CreateUserTH1D( "Phi1stEle_PAS"+region	     ,    60  , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "Pt2ndEle_PAS"+region	     , 	  100 , 0       , 1000     ); 
     CreateUserTH1D( "Eta2ndEle_PAS"+region	     , 	  100 , -5      , 5	   ); 
     CreateUserTH1D( "Phi2ndEle_PAS"+region	     , 	  60  , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "Charge1stEle_PAS"+region	     , 	  2   , -1.0001 , 1.0001   ); 
     CreateUserTH1D( "Charge2ndEle_PAS"+region	     , 	  2   , -1.0001 , 1.0001   ); 
     CreateUserTH1D( "MET_PAS"+region                ,    200 , 0       , 1000	   ); 
     CreateUserTH1D( "METPhi_PAS"+region	     , 	  60  , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "Pt1stJet_PAS"+region           ,    100 , 0       , 1000	   ); 
     CreateUserTH1D( "Eta1stJet_PAS"+region          ,    100 , -5      , 5	   ); 
     CreateUserTH1D( "Phi1stJet_PAS"+region	     , 	  60  , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "sT_PAS"+region                 ,    200 , 0       , 2000	   ); 
     CreateUserTH1D( "Mee_PAS"+region		     ,    200 , 0       , 2000	   ); 
     CreateUserTH1D( "Me1j1_PAS"+region		     ,    200 , 0       , 2000	   ); 
     CreateUserTH1D( "Me2j1_PAS"+region		     ,    200 , 0       , 2000	   ); 
     CreateUserTH1D( "Meejj_PAS"+region              ,    200 , 0       , 2000     );
     CreateUserTH1D( "Ptee_PAS"+region               ,    200 , 0       , 2000     );
     CreateUserTH1D( "Mt_MET_Ele1_PAS"+region        ,    200 , 0       , 2000     );
     CreateUserTH1D( "Mt_MET_Ele2_PAS"+region        ,    200 , 0       , 2000     );
     CreateUserTH1D( "HT"+region                     ,    200 , 0       , 2000     );
     CreateUserTH1D( "Pt1stEle_tight"+region         ,    100 , 0       , 1000     );
     CreateUserTH1D( "Mee_tight"+region              ,    200 , 0       , 2000     );
     CreateUserTH1D( "Me1j1_tight"+region            ,    200 , 0       , 2000     );
     CreateUserTH1D( "sT_tight"+region               ,    200 , 0       , 2000     );
     CreateUserTH1D( "Pt2ndEle_tight"+region         ,    100 , 0       , 1000     );
     CreateUserTH1D( "Me2j1_tight"+region            ,    200 , 0       , 2000     );
     CreateUserTH1D( "MET_tight"+region              ,    200 , 0       , 1000     );
     CreateUserTH1D( "HT_tight"+region               ,    200 , 0       , 2000     );
     CreateUserTH1D( "Mt_MET_Ele1_tight"+region      ,    200 , 0       , 2000     );
     CreateUserTH1D( "Mt_MET_Ele2_tight"+region      ,    200 , 0       , 2000     );
     CreateUserTH1D( "Phi1stEle_tight"+region        ,    60  , -3.1416 , +3.1416  );
     CreateUserTH1D( "Phi2ndEle_tight"+region        ,    60  , -3.1416 , +3.1416  );
     CreateUserTH1D( "METPhi_tight"+region           ,    60  , -3.1416 , +3.1416  );
     //error squared for fake rate:
     CreateUserTH1D("errFRsq_Pt1stEle"+region      ,    100 , 0       , 1000);
     CreateUserTH1D("errFRsq_Mee"+region           ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Me1j1"+region         ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_sT"+region            ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_MET"+region           ,    200 , 0       , 1000);
     CreateUserTH1D("errFRsq_Me2j1"+region         ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Mt_MET_Ele1"+region   ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Mt_MET_Ele2"+region   ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Pt2ndEle"+region      ,    100 , 0       , 1000);
     CreateUserTH1D("errFRsq_Pt1stEle_tight"+region      ,    100 , 0       , 1000);
     CreateUserTH1D("errFRsq_Mee_tight"+region           ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Me1j1_tight"+region         ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_sT_tight"+region            ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_MET_tight"+region           ,    200 , 0       , 1000);
     CreateUserTH1D("errFRsq_Me2j1_tight"+region         ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Pt2ndEle_tight"+region      ,    100 , 0       , 1000);
     CreateUserTH1D("errFRsq_Mt_MET_Ele1_tight"+region   ,    200 , 0       , 2000);
     CreateUserTH1D("errFRsq_Mt_MET_Ele2_tight"+region   ,    200 , 0       , 2000);

     //for the 1P1F regions:
     CreateUserTH1D( "HEEPEle_Pt"+region             ,    100 , 0       , 1000     );
     CreateUserTH1D( "LooseEle_Pt"+region            ,    100 , 0       , 1000     );		                           
		                           
     CreateUserTH1D( "nVertex_PAS"+region           ,    31   , -0.5   , 30.5	 ) ; 
		                           
     CreateUserTH1D( "DR_Ele1Jet1_PAS"+region	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
     CreateUserTH1D( "DR_Ele2Jet1_PAS"+region	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   }

     //control region for scale factors
     CreateUserTH1D("MeeControlReg" , 200 , 0 , 2000);  

     //HT spectrum study
     CreateUserTH1D("HT_fullSpectrum" , 200 , 0 ,2000);
     CreateUserTH1D("fakeRateEffective", 220, 0 ,1.1); 
     CreateUserTH1D("min_prescale", 100 , 0.5, 1.5);
     CreateUserTH1D("pileup_weight", 100, 0.5, 1.5);
     CreateUserTH1D("gen_weight", 200, 0, 12);
     CreateUserTH1D("electronScaleFactors", 100, 0.5, 1.5);
     CreateUserTH1D("totalWeight", 200, 0, 15);
     CreateUserTH1D("Mee_HTStudyRegion", 200, 0, 2000);
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
     //-----------------------------------------------------------------
     // Print progress
     //-----------------------------------------------------------------
     if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass:Loop(): jentry = " << jentry << "/" << nentries << std::endl;
     //if (jentry==150) break; //run over just a few events for troubleshooting
     //// run ls event
     unsigned int run = readerTools_->ReadValueBranch<UInt_t>("run");
     unsigned int ls = readerTools_->ReadValueBranch<UInt_t>("ls");
     unsigned long long int event = readerTools_->ReadValueBranch<ULong64_t>("event");
     //std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;
     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int passedJSON = passJSON ( run,
         readerTools_->ReadValueBranch<UInt_t>("ls"),
         isData() ) ;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     double pileup_weight = readerTools_->ReadValueBranch<Float_t>("puWeight");
     //std::cout<<"pileup_weight = "<<pileup_weight<<std::endl;
     double pileup_weight_only = pileup_weight;
     if ( isData() ) pileup_weight = 1.0;
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = readerTools_->ReadValueBranch<Float_t>("Weight");
     double gen_weight_copy = gen_weight;
     if ( isData() ) gen_weight = 1.0;
     std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
     if(current_file_name.find("powhegMiNNLO") != std::string::npos) gen_weight = TMath::Sign(1, gen_weight);
     //if ( isData && Ele2_ValidFrac > 998. ){
     //  gen_weight = 0.0;
     //  if      (  60.0 < M_e1e2 < 120. ) gen_weight = 0.61;
     //  else if ( 120.0 < M_e1e2 < 200. ) gen_weight = 0.42;
     //  else if ( 200.0 < M_e1e2        ) gen_weight = 0.42;
     //}

     // std::cout << "Gen weight = " << int ( 1.0 / gen_weight ) << std::endl;
     //std::cout << "Gen weight = " << gen_weight << "; isData? " << isData() << std::endl;

     //std::cout<<"current file = "<<current_file_name<<std::endl;
     // SIC remove March 2018
     //// TopPt reweight
     //// only valid for powheg
     //if(current_file_name.find("TT_") != std::string::npos) {
     //  gen_weight*=TopPtWeight;
     //}
    
     //--------------------------------------------------------------------------
     // Electron scale factors for MC only
     //--------------------------------------------------------------------------
     float recoHeepSF = 1.0;
     if(!isData()) {
       float ele1PtUncorr = readerTools_->ReadValueBranch<Float_t>("Ele1_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle1_ECorr");
       float ele2PtUncorr = readerTools_->ReadValueBranch<Float_t>("Ele2_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle2_ECorr");

       float heepSFEle1 = 1.0;
       float heepSFEle2 = 1.0;
       float recoSFEle1 = 1.0;
       float recoSFEle2 = 1.0;
       float zVtxSF = 1.0;
       bool verbose = false;

         if(analysisYear==2017) {
           zVtxSF = ElectronScaleFactors2017::zVtxSF;
         }
	 //std::cout<<"zVtxSF = "<<zVtxSF<<std::endl;
         if(readerTools_->ReadValueBranch<Int_t>("nEle_store")>=1){
           recoSFEle1 = readerTools_->ReadValueBranch<Float_t>("Ele1_RecoSF");
         }
         if (readerTools_->ReadValueBranch<Int_t>("nEle_store") >=2){
           recoSFEle2 = readerTools_->ReadValueBranch<Float_t>("Ele2_RecoSF");
         }
         if(readerTools_->ReadValueBranch<Bool_t>("Ele1_PassHEEPID") == true){
           heepSFEle1 = readerTools_->ReadValueBranch<Float_t>("Ele1_HEEPSF");
         }
         if(readerTools_->ReadValueBranch<Bool_t>("Ele2_PassHEEPID") == true){
           heepSFEle2 = readerTools_->ReadValueBranch<Float_t>("Ele2_HEEPSF");
         }
         //std::cout<<"HEEPSF ele1: "<<heepSFEle1<<"; HEEPSF ele2: "<<heepSFEle2<<std::endl;
         //std::cout<<"recoSF ele1: "<<recoSFEle1<<"; recoSF ele2: "<<recoSFEle2<<std::endl;
         recoHeepSF *= zVtxSF*recoSFEle1*recoSFEle2*heepSFEle1*heepSFEle2;
         //std::cout<<"recoHEEPSF = "<<recoHeepSF<<std::endl;
       gen_weight*=recoHeepSF;
       // stick the gen_weight in with the pileup_weight
       pileup_weight*=gen_weight;
     }

     //--------------------------------------------------------------------------
     // Trigger
     //--------------------------------------------------------------------------
     // Find the right prescale for this event
     double min_prescale = 1;
     int passTrigger = 0;
     std::string triggerName = "";
     double Ele1_hltPhotonPt = readerTools_->ReadValueBranch<Float_t>("Ele1_MatchedHLTriggerObjectPt");
     bool passSinglePhoton = false;
     bool passDoubleEle = false;

     if ( Ele1_hltPhotonPt > 0.0 ) {
       if(analysisYear==2016) {
         if (run >= 276453 && run <= 278822){
           passDoubleEle = readerTools_->ReadValueBranch<Float_t>("H_DoubleEle33_CIdL_GsfIdVL") > 0.1;
      	 }else{
           passDoubleEle = readerTools_->ReadValueBranch<Float_t>("H_DoubleEle33_CIdL_MW") > 0.1;
         }
	 passSinglePhoton = readerTools_->ReadValueBranch<Float_t>("H_Photon175")  > 0.1 && Ele1_hltPhotonPt >= 175.;
       }
       else if(analysisYear==2017) {
         passDoubleEle = readerTools_->ReadValueBranch<Float_t>("H_DoubleEle33_CIdL_MW") > 0.1;
         passSinglePhoton = readerTools_->ReadValueBranch<Float_t>("H_Photon200")  > 0.1 && Ele1_hltPhotonPt >= 200.;
       }
       else if(analysisYear==2018) {
         passSinglePhoton = readerTools_->ReadValueBranch<Float_t>("H_Photon200")  > 0.1 && Ele1_hltPhotonPt >= 200.;
	 passDoubleEle = readerTools_->ReadValueBranch<Float_t>("H_DoubleEle25_CIdL_MW") > 0.1;
       }
       
       if(current_file_name.find("SinglePhoton") != std::string::npos) {
	 if(passSinglePhoton & !passDoubleEle) { passTrigger = 1; triggerName = "Photon"; }
	 if(passSinglePhoton && passDoubleEle) { passTrigger = 1; triggerName = "Photon_and_DoubleEle"; }
       }
       else if(current_file_name.find("SingleElectron") != std::string::npos) {
	 if(passDoubleEle && !passSinglePhoton) { passTrigger = 1; triggerName = "DoubleEle"; }
       }
       else{//EGamma datasets for 2018, or MC 
         if(passSinglePhoton & !passDoubleEle) {passTrigger = 1; triggerName = "Photon";}
	 if(passDoubleEle & !passSinglePhoton) {passTrigger = 1; triggerName = "DoubleEle";}
	 if(passDoubleEle && passSinglePhoton) {passTrigger = 1; triggerName = "Photon_and_DoubleEle";}
       }
     }
     if(isData() && passTrigger) {//if we use photon175 only it's unprescaled so all of this becomes unecessary
       //std::cout << "INFO: lookup trigger name " << triggerName << " for year: " << year << std::endl;
      //if (triggerName.find("Photon") != std::string::npos){
        //min_prescale = run2PhotonTriggerPrescales.LookupPrescale(analysisYearStr,triggerName);
      //}
      //int hltPrescale = psProv.hltPrescale("HLT_"+triggerName+"_v", run, ls);
      //int l1Prescale = 1;
      //std::string l1Seed = "";
      //if(triggerName == "Photon22") {
        //l1Seed = "L1_SingleEG18";
        //l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
      //}
      //else if(triggerName == "Photon25" || triggerName == "Photon30" || triggerName == "Photon36") {
        //l1Seed = "L1_SingleEG26";
        //l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
      //}
      //else if(triggerName == "Photon33") {
        //l1Seed = "L1_SingleEG26er2p5";
        //l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
      //}
      //else if(triggerName == "Photon50" || triggerName == "Photon75" || triggerName == "Photon90" || triggerName == "Photon120" ||
        //  triggerName == "Photon150" || (analysisYear > 2016 && triggerName == "Photon175") ) {
        //int eg34Prescale = psProv.l1Prescale("L1_SingleEG34", run, ls);
      //l1Seed = "L1_SingleEG40";
        //if(analysisYear == 2018)
          //l1Seed = "SingleEG42er2p5";
        //l1Prescale = psProv.l1Prescale(l1Seed, run, ls);
      //}
      //else if(analysisYear == 2016 && triggerName == "Photon175") {
        //l1Prescale = 1;
        //hltPrescale = 1;
      //}
      //if(l1Prescale <= 0 || hltPrescale <= 0)
        //std::cout << "INFO: " << triggerName << ": l1 seed = " << l1Seed << " has prescale = " << l1Prescale << "; hlt prescale = " << hltPrescale << 
          //"; run = " << run << " ls = " << ls << "; psColumn = " << psProv.getRunInfo(run)->psColumn(ls) << "; l1 menu=" << psProv.getRunInfo(run)->l1Menu() <<
          //"; hlt menu=" << psProv.getRunInfo(run)->hltMenu() << "; trig mode=" << psProv.getRunInfo(run)->triggerMode() << std::endl;
      //min_prescale = l1Prescale * hltPrescale;
      if(min_prescale <= 0)
        passTrigger = false;
     }
     //std::cout<<"trigger name: "<<triggerName<<", prescale: "<<min_prescale<<", pass trigger: "<<passTrigger<<", is data?: "<<isData()<<std::endl;

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
     
     double Ele1_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele1_SCEta");
     double Ele2_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele2_SCEta");
     if( fabs( Ele1_SCEta  ) < eleEta_bar )        ele1_isBarrel  = true;
     if( fabs( Ele1_SCEta  ) > eleEta_end1_min &&
         fabs( Ele1_SCEta  ) < eleEta_end1_max )   ele1_isEndcap1 = true;
     if( fabs( Ele1_SCEta  ) > eleEta_end2_min &&
         fabs( Ele1_SCEta  ) < eleEta_end2_max )   ele1_isEndcap2 = true;

     if( fabs( Ele2_SCEta  ) < eleEta_bar )        ele2_isBarrel  = true;
     if( fabs( Ele2_SCEta  ) > eleEta_end1_min &&
         fabs( Ele2_SCEta  ) < eleEta_end1_max )   ele2_isEndcap1 = true;
     if( fabs( Ele2_SCEta  ) > eleEta_end2_min &&
         fabs( Ele2_SCEta  ) < eleEta_end2_max )   ele2_isEndcap2 = true;

     bool ele1_isEndcap = ( ele1_isEndcap1 || ele1_isEndcap2 ) ;
     bool ele2_isEndcap = ( ele2_isEndcap1 || ele2_isEndcap2 ) ;

     bool isBB = ( ele1_isBarrel && ele2_isBarrel ) ;
     bool isEB = ( ( ele1_isBarrel && ele2_isEndcap ) || ( ele2_isBarrel && ele1_isEndcap ) );
     bool isEE = ( ele1_isEndcap && ele2_isEndcap ) ;

     //--------------------------------------------------------------------------
     // Make this a QCD fake rate calculation
     //--------------------------------------------------------------------------
     // LooseEle Pt is the uncorrected SCEt
     double Ele1_Pt = readerTools_->ReadValueBranch<Float_t>("Ele1_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle1_ECorr");
     double Ele2_Pt = readerTools_->ReadValueBranch<Float_t>("Ele2_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle2_ECorr");
     double Ele1_Phi = readerTools_->ReadValueBranch<Float_t>("Ele1_Phi");
     double Ele2_Phi = readerTools_->ReadValueBranch<Float_t>("Ele2_Phi");

     bool verboseFakeRateCalc = false;
     double fakeRate1 = -1;
     double fakeRate2 = -1;
     double eFakeRate1 = 0;
     double eFakeRate2 = 0;
     if (analysisYear != 2018){
       fakeRate1 = qcdFakeRateReader.GetFakeRate(Ele1_Pt,"",Ele1_SCEta);
       fakeRate2 = qcdFakeRateReader.GetFakeRate(Ele2_Pt, "", Ele2_SCEta);
       eFakeRate1 = qcdFakeRateReader.GetFakeRateError(Ele1_Pt,"",Ele1_SCEta);
       eFakeRate2 = qcdFakeRateReader.GetFakeRateError(Ele2_Pt, "", Ele2_SCEta);
     }else{
       fakeRate1 = qcdFakeRateReader.GetFakeRate(Ele1_Pt, Ele1_SCEta, Ele1_Phi, run);
       fakeRate2 = qcdFakeRateReader.GetFakeRate(Ele2_Pt, Ele2_SCEta, Ele2_Phi, run);
       eFakeRate1 = qcdFakeRateReader.GetFakeRateError(Ele1_Pt, Ele1_SCEta, Ele1_Phi, run);
       eFakeRate2 = qcdFakeRateReader.GetFakeRateError(Ele2_Pt, Ele2_SCEta, Ele2_Phi, run);
     }
     //std::cout<<"fake rate 1 = "<<fakeRate1<<" +/- "<<eFakeRate1<<std::endl;

     //--------------------------------------------------------------------------
     // Finally have the effective fake rate
     //--------------------------------------------------------------------------

     // add error on fake rate as well
     double fakeRateEffective  = fakeRate1/(1-fakeRate1); // require loose electron to fail HEEP ID
     double eFakeRateEff1 = eFakeRate1 / ( (1-fakeRate1)*(1-fakeRate1) );
     // (eFakeRate1) * dFakeRateEffective / dFakeRate1

     //if(1-fakeRate1 <= 0)
     //{
     //  cout << "ERROR: Found fakeRate1: " << fakeRate1 << " for SCEta=" << LooseEle1_SCEta << " SCEt="
     //    << LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) << "=" << LooseEle1_SCEnergy << "/" << 
     //    cosh(LooseEle1_SCEta) << endl;
     //}
     double nEle_store = readerTools_->ReadValueBranch<Int_t>("nEle_store");
     double eFakeRateEff2;
     if ( nEle_store < 2 ) {fakeRate2 = 0; eFakeRate2 = 0;}     

     fakeRateEffective += fakeRate2/(1-fakeRate2);
     eFakeRateEff2 = eFakeRate2 / ( (1-fakeRate2)*(1-fakeRate2) );

     double eFakeRateEffectiveSQ = fakeRateEffective * fakeRateEffective + eFakeRateEff1 * eFakeRateEff1 + eFakeRateEff2 * eFakeRateEff2 ;
     //Error on f(x,y) = A(x) + B(y) = sqrt( (err_x * dA/dx)^2 + (err_y * dB/dy)^2 )
     //For N events, the uncertainty will be [ sumOverNEvents(eFakeRateEffectiveSQ_i) ]^(1/2)
     double eFakeRateEffective = sqrt(eFakeRateEffectiveSQ); 
     //double eFakeRateEffective = sqrt(fakeRateEffective);
     //std::cout<<"FR ele1 = "<<fakeRate1<<" +/- "<<eFakeRate1<<"; FR ele2 = "<<fakeRate2<<" +/- "<<eFakeRate2<<std::endl;
     //std::cout<<"  FR effective = FR1/(1 - FR1) + FR2/(1 - FR2) = "<<fakeRateEffective<<" +/- "<<eFakeRateEffective<<std::endl;
     //std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
     //--------------------------------------------------------------------------
     // User has the option to use a flat fake rate (e.g. 1.0 = no fake rate)
     //--------------------------------------------------------------------------
     
     if ( override_fakeRate ) {
       fakeRateEffective = fakeRate_override;
       eFakeRateEffectiveSQ = 0;
     }
     if(fakeRateEffective != 1 && !isData()) {std::cout<<"ERROR: found fake rate != 1 for MC"<<std::endl;}

     //--------------------------------------------------------------------------
     // How many loose electrons have HEEP ID?
     //--------------------------------------------------------------------------

     int nPass = 0;
     if ( readerTools_->ReadValueBranch<Bool_t>("Ele1_PassHEEPID") == true ) nPass ++;
     if ( readerTools_->ReadValueBranch<Bool_t>("Ele2_PassHEEPID") == true ) nPass ++;

     //--------------------------------------------------------------------------
     // Calculate a few missing variables
     //--------------------------------------------------------------------------

     double Ele1_Eta = readerTools_->ReadValueBranch<Float_t>("Ele1_Eta");
     double Ele2_Eta = readerTools_->ReadValueBranch<Float_t>("Ele2_Eta");
     //double Ele1_Phi = readerTools_->ReadValueBranch<Float_t>("Ele1_Phi");
     //double Ele2_Phi = readerTools_->ReadValueBranch<Float_t>("Ele2_Phi");
     double Ele1_Charge = readerTools_->ReadValueBranch<Int_t>("Ele1_Charge");
     double Ele2_Charge = readerTools_->ReadValueBranch<Int_t>("Ele2_Charge");
     double Jet1_Pt = readerTools_->ReadValueBranch<Float_t>("Jet1_Pt");
     double Jet2_Pt = readerTools_->ReadValueBranch<Float_t>("Jet2_Pt");
     double Jet1_Eta = readerTools_->ReadValueBranch<Float_t>("Jet1_Eta");
     double Jet2_Eta = readerTools_->ReadValueBranch<Float_t>("Jet2_Eta");
     double Jet1_Phi = readerTools_->ReadValueBranch<Float_t>("Jet1_Phi");
     double Jet2_Phi = readerTools_->ReadValueBranch<Float_t>("Jet2_Phi");
     double nMuon_ptCut = readerTools_->ReadValueBranch<Int_t>("nMuon_ptCut");
     TLorentzVector loose_ele1, loose_ele2 , jet1, jet2;
     loose_ele1.SetPtEtaPhiM ( Ele1_Pt , Ele1_Eta , Ele1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( Ele2_Pt , Ele2_Eta , Ele2_Phi , 0.0 );
     jet1.SetPtEtaPhiM       ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
     jet2.SetPtEtaPhiM       ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );

     TLorentzVector loose_e1e2 = loose_ele1 + loose_ele2;
     TLorentzVector j1j2 = jet1 + jet2;

     TLorentzVector e1j1 = loose_ele1 + jet1;
     TLorentzVector e2j1 = loose_ele2 + jet1;

     //--------------------------------------------------------------------------
     // Now fill cut values
     // DON'T use the pileup weight ... it's included by default
     // DO    use the prescale.  It's already 1.0 for MC.
     // DO    use the effective fake rate.  It will only mean anything when you run
     //       over data, though, so be sure to override it to 1.0 
     //       if you run over Monte Carlo
     //--------------------------------------------------------------------------


     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------
     double PFMET_Type1_Pt  = readerTools_->ReadValueBranch<Float_t>("PFMET_Type1_Pt");
     double PFMET_Type1_Phi  = readerTools_->ReadValueBranch<Float_t>("PFMET_Type1_Phi");
     double nEle_ptCut = readerTools_->ReadValueBranch<Int_t>("nEle_ptCut");
     double nJet_ptCut = readerTools_->ReadValueBranch<Int_t>("nJet_store"); //ptCut");
     double nJet_store = readerTools_->ReadValueBranch<Int_t>("nJet_store");
     double nVertex = readerTools_->ReadValueBranch<Int_t>("nVertex");
     // reweighting
     fillVariableWithValue ( "Reweighting", 1, pileup_weight * min_prescale * fakeRateEffective) ; 
     //stitching in missing DY events: https://cms-talk.web.cern.ch/t/bug-in-ul-pt-binned-dy-samples/11639
     bool passLHECuts = true;
     double Vpt = -1;
     if(current_file_name.find("DYJetsToLL_M-50_TuneCP5") != std::string::npos) {
       //std::cout<<"found inclusive DYJets file"<<std::endl;
       passLHECuts = false;
       //Vpt = (readerTools_->ReadValueBranch<Float_t>("LHE_Vpt"));
       //std::cout<<"LHE_Vpt = "<<Vpt<<std::endl;
       if(readerTools_->ReadValueBranch<Float_t>("LHE_Vpt") == 0){
         passLHECuts = true;
         //std::cout<<"evaluated statement if LHE_Vpt == 0 as true"<<std::endl;
       }
     }

     //DYJets NNLO sample stitching
     if(current_file_name.find("DYJetsToEE_M-50_massWgtFix_TuneCP5") != std::string::npos) {
       passLHECuts = false;
       TLorentzVector e1, e2;
       float LHEEle1_Pt = readerTools_->ReadValueBranch<Float_t>("LHEElectron1_Pt");
       float LHEEle2_Pt = readerTools_->ReadValueBranch<Float_t>("LHEElectron2_Pt");
       float LHEEle1_Eta = readerTools_->ReadValueBranch<Float_t>("LHEElectron1_Eta");
       float LHEEle2_Eta = readerTools_->ReadValueBranch<Float_t>("LHEElectron2_Eta");
       float LHEEle1_Phi = readerTools_->ReadValueBranch<Float_t>("LHEElectron1_Phi");
       float LHEEle2_Phi = readerTools_->ReadValueBranch<Float_t>("LHEElectron2_Phi");
       if(LHEEle2_Pt > 0) std::cout<<"Ele2 pt, eta, phi = "<<LHEEle2_Pt<<", "<<LHEEle2_Eta<<", "<<LHEEle2_Phi<<std::endl;
       e1.SetPtEtaPhiM ( LHEEle1_Pt, LHEEle1_Eta, LHEEle1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( LHEEle2_Pt, LHEEle2_Eta, LHEEle2_Phi, 0.0 );
       float dielectron_mass = (e1 + e2).M();
       if(dielectron_mass>0) std::cout<<"LHE Mee = "<<dielectron_mass<<std::endl;
       if(dielectron_mass > 50 && dielectron_mass < 100)
         passLHECuts = true;       
     }

     float HT = readerTools_->ReadValueBranch<Float_t>("LHE_HT");
     //WJets HT binned sample stitching
     if (current_file_name.find("WJetsToLNu_TuneCP5_13TeV-madgraphMLM") != std::string::npos) {
       passLHECuts = false;
       if (HT < 80) {
         passLHECuts = true;
         //std::cout<<"found HT = "<<HT<<" in WJets inclusive sample, set passLHECuts to true"<<std::endl;
       }
     }
     if (current_file_name.find("WJetsToLNu_HT-70To100") != std::string::npos) {
       passLHECuts = false;
       if (HT >= 80) {
         passLHECuts = true;
       }
     }

     //std::cout<<"set passLHECuts to "<<passLHECuts<<std::endl;
     //std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
     fillVariableWithValue("PassLHECuts",passLHECuts, pileup_weight* min_prescale * fakeRateEffective);
     // JSON variable
     fillVariableWithValue( "PassJSON" , passedJSON, pileup_weight * min_prescale * fakeRateEffective) ; 

     fillVariableWithValue ( "PassHLT", passTrigger, pileup_weight * min_prescale * fakeRateEffective ) ;

     // Fill noise filters
     // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
     // we filled these at skim time
     fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Bool_t>("PassGlobalSuperTightHalo2016Filter")     == 1), fakeRateEffective * min_prescale * pileup_weight);
     fillVariableWithValue("PassGoodVertices"                   , int(readerTools_->ReadValueBranch<Bool_t>("PassGoodVertices")                       == 1), fakeRateEffective * min_prescale * pileup_weight);
     fillVariableWithValue("PassHBHENoiseFilter"                , int(readerTools_->ReadValueBranch<Bool_t>("PassHBHENoiseFilter")                    == 1), fakeRateEffective * min_prescale * pileup_weight);
     fillVariableWithValue("PassHBHENoiseIsoFilter"             , int(readerTools_->ReadValueBranch<Bool_t>("PassHBHENoiseIsoFilter")                 == 1), fakeRateEffective * min_prescale * pileup_weight);
     // eBadScFilter not suggested for MC
     if(isData())
       fillVariableWithValue("PassBadEESupercrystalFilter"      , int(readerTools_->ReadValueBranch<Bool_t>("PassBadEESupercrystalFilter")            == 1), fakeRateEffective * min_prescale*pileup_weight);
     else
       fillVariableWithValue("PassBadEESupercrystalFilter"      , 1                                                                                          , fakeRateEffective * min_prescale*pileup_weight);
     fillVariableWithValue("PassEcalDeadCellTrigPrim"           , int(readerTools_->ReadValueBranch<Bool_t>("PassEcalDeadCellTrigPrim")               == 1), fakeRateEffective * min_prescale*pileup_weight);
     // not recommended
     //fillVariableWithValue("PassChargedCandidateFilter"         , int(readerTools_->ReadValueBranch<Double_t>("PassChargedCandidateFilter")             == 1), fakeRateEffective * min_prescale);
     fillVariableWithValue("PassBadPFMuonFilter"                , int(readerTools_->ReadValueBranch<Bool_t>("PassBadPFMuonFilter")                    == 1), fakeRateEffective * min_prescale*pileup_weight);
     // EcalBadCalibV2 for 2017, 2018
     if(analysisYear > 2016)
       fillVariableWithValue("PassEcalBadCalibV2Filter"         , int(readerTools_->ReadValueBranch<Bool_t>("PassEcalBadCalibFilter")               == 1), fakeRateEffective * min_prescale*pileup_weight);
     else
       fillVariableWithValue("PassEcalBadCalibV2Filter"         , 1                                                                                          , fakeRateEffective * min_prescale*pileup_weight);

     // MET
     fillVariableWithValue ( "PFMET"  , PFMET_Type1_Pt, pileup_weight * min_prescale * fakeRateEffective ) ;
     
     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon_ptCut           , pileup_weight * min_prescale  * fakeRateEffective ); 
										      			      
     // Electrons
     fillVariableWithValue(   "nEleLoose"                     , nEle_ptCut     , pileup_weight * min_prescale  * fakeRateEffective );
     fillVariableWithValue(   "nEleTight"                     , nPass            , pileup_weight * min_prescale  * fakeRateEffective );
     if ( nEle_ptCut >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt   , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta  , pileup_weight * min_prescale  * fakeRateEffective ) ;
     }
     if ( nEle_ptCut >= 2 ) { 
       fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt   , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Ele2_Eta"                      , Ele2_Eta  , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "M_e1e2"                        , loose_e1e2.M()   , pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "M_e1e2_220"                    , loose_e1e2.M()   , pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "Pt_e1e2"                       , loose_e1e2.Pt()  , pileup_weight * min_prescale  * fakeRateEffective );
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJet_ptCut , pileup_weight * min_prescale  * fakeRateEffective );
     if ( nJet_store >= 1 ) {
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt   , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta , pileup_weight * min_prescale  * fakeRateEffective ) ;
     }

     // DeltaR
     if ( nEle_ptCut >= 2 && nJet_store >= 1) {
       double sT_eej = Ele1_Pt + Ele2_Pt + Jet1_Pt ;
       fillVariableWithValue( "DR_Ele1Jet1"                   , loose_ele1.DeltaR ( jet1 ), pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "DR_Ele2Jet1"                   , loose_ele2.DeltaR ( jet1 ), pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "sT_eej_200"               , sT_eej                   , pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "sT_eej_400"               , sT_eej                   , pileup_weight * min_prescale  * fakeRateEffective );
     //  fillVariableWithValue( "sT_eej_850"               , sT_eej                   , pileup_weight * min_prescale  * fakeRateEffective );
     } 

     // transverse mass Mt^2 = 2 * Ele_Pt * MET * (1 - cos( Ele_phi - MET_phi )) 
     double Mt_MET_Ele1 = sqrt( 2 * Ele1_Pt * PFMET_Type1_Pt * (1 - cos( Ele1_Phi - PFMET_Type1_Phi )));
     double Mt_MET_Ele2 = sqrt( 2 * Ele2_Pt * PFMET_Type1_Pt * (1 - cos( Ele2_Phi - PFMET_Type1_Phi )));

     fillVariableWithValue( "closureTestControlReg"   , 1           , pileup_weight * min_prescale  * fakeRateEffective );
     fillVariableWithValue( "closureTestPreselection" , 1                , pileup_weight * min_prescale  * fakeRateEffective );  
     fillVariableWithValue( "closureTestBDTSelection" , 1           , pileup_weight * min_prescale  * fakeRateEffective );
     fillVariableWithValue( "HTSpectrumStudy"         , 1           , pileup_weight * min_prescale  * fakeRateEffective );
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     // DO    use the pileup weight.  It's equal to 1.0 for data.  
     // DO    use the min_prescale.  It's equal to 1.0 for Monte Carlo
     // DO    use the effective fake rate.  It will only mean anything when you run
     //       over data, though, so be sure to override it to 1.0 
     //       if you run over Monte Carlo
     //--------------------------------------------------------------------------     

     bool passed_controlRegion = ( passedAllPreviousCuts("closureTestControlReg") && passedCut("closureTestControlReg"));
     bool passed_preselection = ( passedAllPreviousCuts("closureTestPreselection") && passedCut("closureTestPreselection") );
     bool passed_tightSelection = ( passedAllPreviousCuts("closureTestBDTSelection") && passedCut("closureTestBDTSelection") );
     bool passed_HTSelection = (passedAllPreviousCuts("HTSpectrumStudy") && passedCut("HTSpectrumStudy"));

     if (passed_HTSelection) {
       FillUserTH1D("HT_fullSpectrum", HT, pileup_weight * min_prescale * fakeRateEffective ); 
       FillUserTH1D("min_prescale", min_prescale, 1);
       FillUserTH1D("fakeRateEffective", fakeRateEffective, 1);
       FillUserTH1D("pileup_weight", pileup_weight_only, 1);
       FillUserTH1D("gen_weight", gen_weight_copy, 1);
       FillUserTH1D("electronScaleFactors", recoHeepSF, 1);
       FillUserTH1D("totalWeight", pileup_weight*min_prescale*fakeRateEffective,1);
       FillUserTH1D("Mee_HTStudyRegion", loose_e1e2.M(), pileup_weight * min_prescale * fakeRateEffective );
     }

     if (passed_controlRegion) { FillUserTH1D("MeeControlReg", loose_e1e2.M() ,pileup_weight * min_prescale * fakeRateEffective ); }  

     if ( passed_preselection ) {
       //std::cout<<"FR = "<<fakeRateEffective<<" +/- "<<sqrt(eFakeRateEffectiveSQ)<<std::endl;
       double sT_eej = Ele1_Pt + Ele2_Pt + Jet1_Pt ;
       //fill these for all events
       FillUserTH1D("nElectron_PAS"        , nEle_ptCut                , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nMuon_PAS"            , nMuon_ptCut               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nJet_PAS"             , nJet_ptCut                , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stEle_PAS"	   , Ele1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stEle_PAS"	   , Ele1_Eta                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stEle_PAS"	   , Ele1_Phi                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt2ndEle_PAS"	   , Ele2_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta2ndEle_PAS"	   , Ele2_Eta                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi2ndEle_PAS"	   , Ele2_Phi                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge1stEle_PAS"	   , Ele1_Charge               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge2ndEle_PAS"	   , Ele2_Charge               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("MET_PAS"              , PFMET_Type1_Pt            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METPhi_PAS"	   , PFMET_Type1_Phi           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stJet_PAS"         , Jet1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stJet_PAS"        , Jet1_Eta                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stJet_PAS"	   , Jet1_Phi                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("sT_PAS"               , sT_eej                    , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mee_PAS"		   , loose_e1e2.M()            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Ptee_PAS"             , loose_e1e2.Pt()           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nVertex_PAS"          , nVertex                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , loose_ele1.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , loose_ele2.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me1j1_PAS"            , e1j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me2j1_PAS"            , e2j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("HT"                   , HT                        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mt_MET_Ele1_PAS"      , Mt_MET_Ele1               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mt_MET_Ele2_PAS"      , Mt_MET_Ele2               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("errFRsq_Pt1stEle"     , Ele1_Pt                   , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Mee"          , loose_e1e2.M()            , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Me1j1"        , e1j1.M()                  , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_sT"           , sT_eej                    , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Pt2ndEle"     , Ele2_Pt                   , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Me2j1"        , e2j1.M()                  , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_MET"          , PFMET_Type1_Pt            , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Mt_MET_Ele1"   , Mt_MET_Ele1               , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Mt_MET_Ele2"   , Mt_MET_Ele2               , eFakeRateEffectiveSQ );
       //For the 1P1F region:
       if(override_fakeRate){
         if(readerTools_->ReadValueBranch<Bool_t>("Ele1_PassHEEPID") == true && readerTools_->ReadValueBranch<Bool_t>("Ele2_PassHEEPID") == false){
            FillUserTH1D("HEEPEle_Pt"      , Ele1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
            FillUserTH1D("LooseEle_Pt"     , Ele2_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
         }
         if(readerTools_->ReadValueBranch<Bool_t>("Ele1_PassHEEPID") == false && readerTools_->ReadValueBranch<Bool_t>("Ele2_PassHEEPID") == true){
            FillUserTH1D("HEEPEle_Pt"      , Ele2_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
            FillUserTH1D("LooseEle_Pt"     , Ele1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
         }
       }
       //fill these for only the region that this event is in
       //it is possible for none of these to be true because there is a gap between the barrel and endcap
       //regions. If at least one ele falls into that gap then we aren't in one of these three regions. 
       std::string region;
       if (isBB){region = "_BB";}
       else if (isEB){region = "_BE";}
       else if (isEE){region = "_EE";}
       else { std::cout<<"skip filling eta reg hists for event no. "<<jentry<<" in file "<<current_file_name<<std::endl; continue; }
       FillUserTH1D("nElectron_PAS"+region        , nEle_ptCut                , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nMuon_PAS"+region            , nMuon_ptCut               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nJet_PAS"+region             , nJet_ptCut                , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stEle_PAS"+region         , Ele1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stEle_PAS"+region        , Ele1_Eta                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stEle_PAS"+region        , Ele1_Phi                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt2ndEle_PAS"+region         , Ele2_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta2ndEle_PAS"+region        , Ele2_Eta                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi2ndEle_PAS"+region        , Ele2_Phi                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge1stEle_PAS"+region     , Ele1_Charge               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge2ndEle_PAS"+region     , Ele2_Charge               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("MET_PAS"+region              , PFMET_Type1_Pt            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METPhi_PAS"+region           , PFMET_Type1_Phi           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stJet_PAS"+region         , Jet1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stJet_PAS"+region        , Jet1_Eta                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stJet_PAS"+region        , Jet1_Phi                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("sT_PAS"+region               , sT_eej                    , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mee_PAS"+region              , loose_e1e2.M()            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Ptee_PAS"+region             , loose_e1e2.Pt()           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nVertex_PAS"+region          , nVertex                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele1Jet1_PAS"+region      , loose_ele1.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele2Jet1_PAS"+region      , loose_ele2.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me1j1_PAS"+region            , e1j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me2j1_PAS"+region            , e2j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("HT"+region                   , HT                        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mt_MET_Ele1_PAS"+region      , Mt_MET_Ele1               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mt_MET_Ele2_PAS"+region      , Mt_MET_Ele2               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("errFRsq_Pt1stEle"+region     , Ele1_Pt                   , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Mee"+region          , loose_e1e2.M()            , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Me1j1"+region        , e1j1.M()                  , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_sT"+region           , sT_eej                    , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Pt2ndEle"+region     , Ele2_Pt                   , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Me2j1"+region        , e2j1.M()                  , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_MET"+region          , PFMET_Type1_Pt            , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Mt_MET_Ele1"+region   , Mt_MET_Ele1               , eFakeRateEffectiveSQ );
       FillUserTH1D("errFRsq_Mt_MET_Ele2"+region   , Mt_MET_Ele2               , eFakeRateEffectiveSQ );
       if(override_fakeRate){
         if(readerTools_->ReadValueBranch<Bool_t>("Ele1_PassHEEPID") == true && readerTools_->ReadValueBranch<Bool_t>("Ele2_PassHEEPID") == false){
            FillUserTH1D("HEEPEle_Pt"+region      , Ele1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
            FillUserTH1D("LooseEle_Pt"+region     , Ele2_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
         }
         if(readerTools_->ReadValueBranch<Bool_t>("Ele1_PassHEEPID") == false && readerTools_->ReadValueBranch<Bool_t>("Ele2_PassHEEPID") == true){
            FillUserTH1D("HEEPEle_Pt"+region      , Ele2_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
            FillUserTH1D("LooseEle_Pt"+region     , Ele1_Pt                   , pileup_weight * min_prescale * fakeRateEffective );
         }
       }
     }
     if(passed_tightSelection){
         double sT_eej = Ele1_Pt + Ele2_Pt + Jet1_Pt ;
         FillUserTH1D("Pt1stEle_tight"          , Ele1_Pt , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Mee_tight"        , loose_e1e2.M() , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Me1j1_tight"             , e1j1.M(), pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("sT_tight"                 , sT_eej , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Pt2ndEle_tight"          , Ele2_Pt , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Me2j1_tight"             , e2j1.M(), pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("MET_tight"        ,  PFMET_Type1_Pt, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("HT_tight"                , HT      , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Mt_MET_Ele1_tight"      , Mt_MET_Ele1, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Mt_MET_Ele2_tight"      , Mt_MET_Ele2, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Phi1stEle_tight"         , Ele1_Phi, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Phi2ndEle_tight"         , Ele2_Phi, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("METPhi_tight"     , PFMET_Type1_Phi, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("errFRsq_Pt1stEle_tight"     , Ele1_Pt                   , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Mee_tight"          , loose_e1e2.M()            , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Me1j1_tight"        , e1j1.M()                  , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_sT_tight"           , sT_eej                    , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Pt2ndEle_tight"     , Ele2_Pt                   , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Me2j1_tight"        , e2j1.M()                  , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_MET_tight"          , PFMET_Type1_Pt            , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Mt_MET_Ele1_tight"  , Mt_MET_Ele1               , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Mt_MET_Ele2_tight"  , Mt_MET_Ele2               , eFakeRateEffectiveSQ );
         std::string region;
         if (isBB){region = "_BB";}
         else if (isEB){region = "_BE";}
         else if (isEE){region = "_EE";}
         else { std::cout<<"skip filling eta reg hists for event no. "<<jentry<<" in file "<<current_file_name<<std::endl; continue; }
         FillUserTH1D("Pt1stEle_tight"+region                , Ele1_Pt , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Mee_tight"+region        , loose_e1e2.M() , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Me1j1_tight"+region             , e1j1.M(), pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("sT_tight"+region                 , sT_eej , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Pt2ndEle_tight"+region          , Ele2_Pt , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Me2j1_tight"+region             , e2j1.M(), pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("MET_tight"+region        ,  PFMET_Type1_Pt, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("HT_tight"+region                , HT      , pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Mt_MET_Ele1_tight"+region      , Mt_MET_Ele1, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Mt_MET_Ele2_tight"+region      , Mt_MET_Ele2, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Phi1stEle_tight"+region         , Ele1_Phi, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("Phi2ndEle_tight"+region         , Ele2_Phi, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("METPhi_tight"+region     , PFMET_Type1_Phi, pileup_weight * min_prescale * fakeRateEffective );
         FillUserTH1D("errFRsq_Pt1stEle_tight"+region     , Ele1_Pt                   , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Mee_tight"+region          , loose_e1e2.M()            , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Me1j1_tight"+region        , e1j1.M()                  , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_sT_tight"+region           , sT_eej                    , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Pt2ndEle_tight"+region     , Ele2_Pt                   , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Me2j1_tight"+region        , e2j1.M()                  , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_MET_tight"+region          , PFMET_Type1_Pt            , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Mt_MET_Ele1_tight"+region  , Mt_MET_Ele1               , eFakeRateEffectiveSQ );
         FillUserTH1D("errFRsq_Mt_MET_Ele2_tight"+region  , Mt_MET_Ele2               , eFakeRateEffectiveSQ );
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

