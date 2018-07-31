// -*- C++ -*-
//
// Package:    MuonObjectInfoExtractorToCsv
// Class:      MuonObjectInfoExtractorToCsv
// 
/**\class MuonObjectInfoExtractorToCsv MuonObjectInfoExtractorToCsv.cc PhysicsObjectsInfo/MuonObjectInfoExtractorToCsv/src/MuonObjectInfoExtractorToCsv.cc

 Description: [Example on how to extract physics information from a CMS EDM Muon Collection]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Edgar F. Carrera Jarrin (ecarrera@cern.ch)
//         Created:  Wed Jul  4 13:38:41 CEST 2018
// $Id$
//
// Notes:
// 
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//classes included to extract muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 

//classes included to extract tracking for the muons
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//additional classes for storage, containers and operations
#include<vector>
#include<string>
#include "TFile.h"
#include "TTree.h"
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<sstream>

// for trigger analysis

//Following the HLTEventAnalyzerAOD.h, 
//include the following headers:
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Also include headers from  HLTEventAnalyzerAOD.cc
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include <cassert>


//
// class declaration
//

class MuonObjectInfoExtractorToCsv : public edm::EDAnalyzer {
   public:
      explicit MuonObjectInfoExtractorToCsv(const edm::ParameterSet&);
      ~MuonObjectInfoExtractorToCsv();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);

      virtual void analyze(const edm::Event&, const edm::EventSetup&);


      //the follwing are not being used here

      virtual void beginJob() ;

      virtual void endJob() ;

      virtual void endRun(edm::Run const&, edm::EventSetup const&);

      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);



   private:

      
      // from HLTEventAnalyzerAOD.h
      /// module config parameters
      std::string   processName_;
      std::string   triggerName_;
      std::string   datasetName_;
      edm::InputTag triggerResultsTag_;
      edm::InputTag triggerEventTag_;


      // additional class data memebers
      // these are actually the containers where we will store
      // the trigger information
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
      HLTConfigProvider hltConfig_;
     
      
 //declare a function to do the muon analysis
      void analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons);
  //function to store info in csv
  void dumpMuonsToCsv();
  //declare the input tag for the muons collection to be used (read from cofiguration)
  edm::InputTag muonsInput;
  int maxNumObjt;

  //Declare some variables for storage
  std::ofstream myfile;
  int maxpart;
  std::ostringstream oss;
  std::string theHeader;


  //and declare variable that will go into the root tree
  int runno; //run number
  int evtno; //event number
  int nmu;//number of muons in the event
  std::string mu_partype; //type of particle
  std::vector<float> mu_e;
  std::vector<float> mu_pt;
  std::vector<float> mu_pth;
  std::vector<float> mu_ptm;
  std::vector<float> mu_ptl;
  std::vector<float> mu_px;
  std::vector<float> mu_py;
  std::vector<float> mu_pz;
  std::vector<float> mu_eta;
  std::vector<float> mu_phi;
  std::vector<float> mu_ch;
  std::string h1="Mu1";
  std::string h2="Mu20";
  std::string h3="Mu24";
  std::string h4="Mu23";
  std::string h5="Mu30";
  std::string h6="Mu40";
  std::string h7="Mu45";
  std::string h8="Mu60";
  std::string m1="Mu6";
  std::string m2="Mu7";
  std::string m3="Mu8";
  std::string l1="Mu3";
  std::string l2="Mu2";
  std::string l3="Mu4";
  std::string l4="Mu5";
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonObjectInfoExtractorToCsv::MuonObjectInfoExtractorToCsv(const edm::ParameterSet& iConfig):
  processName_(iConfig.getParameter<std::string>("processName")),
  triggerName_(iConfig.getParameter<std::string>("triggerName")),
  datasetName_(iConfig.getParameter<std::string>("datasetName")),
  triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(iConfig.getParameter<edm::InputTag>("triggerEvent"))
  
  {
   //now do what ever initialization is needed
  using namespace std;
  using namespace edm;
  muonsInput = iConfig.getParameter<edm::InputTag>("InputCollection");
  maxNumObjt = iConfig.getUntrackedParameter<int>("maxNumberMuons",10);
  //Print the configuration just to check
  cout << "Here is the information passed to the constructor:" <<endl;
  cout << "HLTEventAnalyzerAOD configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerName = " << triggerName_ << endl
       << "   DataSetName = " << datasetName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventTag = " << triggerEventTag_.encode() << endl;

}




MuonObjectInfoExtractorToCsv::~MuonObjectInfoExtractorToCsv()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// ------------ method called when starting to processes a run  ------------

void MuonObjectInfoExtractorToCsv::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)

//--------------------------------------------------------------------------

{
using namespace std;

  using namespace edm;


  bool changed(true);

  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {

    if (changed) {

      // check if trigger name in (new) config

      if (triggerName_!="@") { // "@" means: analyze all triggers in config

	const unsigned int n(hltConfig_.size());

	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));

	if (triggerIndex>=n) {

	  cout << "HLTEventAnalyzerAOD::analyze:"

	       << " TriggerName " << triggerName_ 

	       << " not available in (new) config!" << endl;

	  cout << "Available TriggerNames are: " << endl;

	  hltConfig_.dump("Triggers");

	}

      }


  } else {

    cout << "HLTEventAnalyzerAOD::analyze:"

	 << " config extraction failure with process name "

	 << processName_ << endl;

  }

  }

}//------------------- beginRun()
// ------------ method called for each event  ------------
void
MuonObjectInfoExtractorToCsv::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  using namespace std;
   using namespace edm;
    using namespace reco;
  using namespace trigger;

   //get the global information first
   runno = iEvent.id().run();
   evtno  = iEvent.id().event();
   
   
    Handle<reco::MuonCollection> mymuons;
   iEvent.getByLabel(muonsInput, mymuons); 
   iEvent.getByLabel(triggerResultsTag_, triggerResultsHandle_);

   analyzeMuons(iEvent,mymuons);
   dumpMuonsToCsv();
     
  
}


// ------------ function to analyze muons
void 
MuonObjectInfoExtractorToCsv::analyzeMuons(const edm::Event& iEvent, const edm::Handle<reco::MuonCollection> &muons)
{

 const unsigned int n(hltConfig_.size());

  //clear the storage containers for this objects in this event
nmu=0;
  mu_e.clear();
  mu_pt.clear();   //all
  mu_ptl.clear(); //low threshold
  mu_ptm.clear(); //medium threshold
  mu_pth.clear(); //high threshold
  mu_px.clear();
  mu_py.clear();
  mu_pz.clear();
  mu_eta.clear();
  mu_phi.clear();
  mu_ch.clear();
 

  //check if the collection is valid
  if(muons.isValid()){
    //get the number of muons in the event
	//loop over all the muons in this event
    int idx = 0;
    const edm::TriggerNames& triggerName = iEvent.triggerNames(*triggerResultsHandle_);
	for (reco::MuonCollection::const_iterator recoMu = muons->begin(); recoMu!=muons->end(); ++recoMu){
      //find only globlal muons for this specific example
      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis?rev=88
	  if(recoMu->isGlobalMuon()) {
	    mu_partype = "G"; 
	    mu_e.push_back(recoMu->energy());
	    mu_pt.push_back(recoMu->pt());
	 
	 
	 //Find PT whose threshold fits with a trigger	
	 for (unsigned int i=0; i!=n; ++i) { 
	 
	
	 if (((hltConfig_.triggerName(i).find(h1) != std::string::npos) or (hltConfig_.triggerName(i).find(h2)!= std::string::npos) or (hltConfig_.triggerName(i).find(h3) != std::string::npos) or (hltConfig_.triggerName(i).find(h4)!= std::string::npos) or (hltConfig_.triggerName(i).find(h5)!= std::string::npos) or (hltConfig_.triggerName(i).find(h6)!= std::string::npos) or (hltConfig_.triggerName(i).find(h7)!= std::string::npos) or (hltConfig_.triggerName(i).find(h8)!= std::string::npos) ) and (triggerResultsHandle_->accept(triggerName.triggerIndex(hltConfig_.triggerName(i)))==1)) //High threshold
          {
	  std::cout<<"Currently analyzing trigger "<< hltConfig_.triggerName(i)<< "Status"<<(triggerResultsHandle_->accept(triggerName.triggerIndex(hltConfig_.triggerName(i))))<<"H"<<std::endl;
          mu_pth.push_back(recoMu->pt());
	  mu_ptm.push_back(-999);
	  mu_ptl.push_back(-999);
	  goto Continue;
	  
          }
          else
          {
           if (((hltConfig_.triggerName(i).find(m1) != std::string::npos) or (hltConfig_.triggerName(i).find(m2)!= std::string::npos) or (hltConfig_.triggerName(i).find(m3) != std::string::npos)) and (triggerResultsHandle_->accept(triggerName.triggerIndex(hltConfig_.triggerName(i)))==1))  //Medium threshold
            { 
	    std::cout<<"Currently analyzing trigger "<< hltConfig_.triggerName(i)<< " Status "<<(triggerResultsHandle_->accept(triggerName.triggerIndex(hltConfig_.triggerName(i))))<<" M "<<std::endl;
	    
             mu_ptm.push_back(recoMu->pt());
	     mu_pth.push_back(-999);
	     mu_ptl.push_back(-999);
	     goto Continue;
            }
            else
            {
	    
             if (((hltConfig_.triggerName(i).find(l1) != std::string::npos) or (hltConfig_.triggerName(i).find(l2)!= std::string::npos) or (hltConfig_.triggerName(i).find(l3) != std::string::npos) or (hltConfig_.triggerName(i).find(l4)!= std::string::npos) ) and (triggerResultsHandle_->accept(triggerName.triggerIndex(hltConfig_.triggerName(i)))==1))  //low threshold
             {
	      std::cout<<"Currently analyzing trigger "<< hltConfig_.triggerName(i)<< "  Status "<<(triggerResultsHandle_->accept(triggerName.triggerIndex(hltConfig_.triggerName(i))))<<" L "<<std::endl;
              mu_ptl.push_back(recoMu->pt());
	      mu_ptm.push_back(-999);
	      mu_pth.push_back(-999);
	      goto Continue;
             }     
            
            }

          }
	  }
	  
	 
	   mu_ptm.push_back(-999);
	   mu_pth.push_back(-999);
	   mu_ptl.push_back(-999);
	  
	  Continue:
	 
	    mu_px.push_back(recoMu->px());
	    mu_py.push_back(recoMu->py());
	    mu_pz.push_back(recoMu->pz());
	    mu_eta.push_back(recoMu->eta());
	    mu_phi.push_back(recoMu->phi());
	    mu_ch.push_back(recoMu->charge());
	    ++idx;
	  }
	}
	nmu = idx;
  }
  
}

// ------------ function to analyze muons
void MuonObjectInfoExtractorToCsv::dumpMuonsToCsv()
{
  unsigned int maxnumobjt = maxNumObjt;
  if(nmu>0){
  oss.str("");oss.clear();oss<<runno;
  myfile<<oss.str();
  oss.str("");oss.clear();oss<<evtno;
  myfile<<","<<oss.str();
    for (unsigned int j=0;j<maxnumobjt;j++){
      oss.str("");oss.clear();oss<<mu_partype;
      myfile<<","<<oss.str();
      oss.str("");oss.clear();oss<<mu_e[j];
      j<mu_e.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      //      std::cout<<maxnumobjt<<"\t"<<nmu<<"\t"<<mu_e.size()<<"\t"<<mu_e[j]<<"\t"<<oss.str()<<std::endl;
      oss.str("");oss.clear();oss<<mu_px[j];
      j<mu_px.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_py[j];
      j<mu_py.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_pz[j];
      j<mu_pz.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_pt[j];
      j<mu_pt.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_pth[j];
      j<mu_pth.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_ptm[j];
      j<mu_ptm.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_ptl[j];
      j<mu_ptl.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_eta[j];
      j<mu_eta.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_phi[j];
      j<mu_phi.size() ? myfile<<","<<oss.str():myfile<<",0.0";
      oss.str("");oss.clear();oss<<mu_ch[j];
      j<mu_ch.size() ? myfile<<","<<oss.str():myfile<<",0.0";
    }
  myfile<<"\n";
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonObjectInfoExtractorToCsv::beginJob()
{
  //Define storage
  myfile.open("MuonObjectInfo.csv");
  //Write the header.
  //create the header string accordingly
  theHeader = "Run,Event";
  for(int j =1;j<maxNumObjt+1;j++){
    oss.str(""); oss<<j;
    std::string idxstr = oss.str();
    theHeader += ",type"+idxstr+",E"+idxstr+",px"+idxstr+",py"+idxstr+",pz"+idxstr+",pt"+idxstr+",pth"+idxstr+",ptm"+idxstr+",ptl"+idxstr+",eta"+idxstr+",phi"+idxstr+",Q"+idxstr;
  }
  
  myfile<<theHeader<<"\n";

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonObjectInfoExtractorToCsv::endJob() 
{

  //save file
  myfile.close();

}


// ------------ method called when ending the processing of a run  ------------
void 
MuonObjectInfoExtractorToCsv::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuonObjectInfoExtractorToCsv::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuonObjectInfoExtractorToCsv::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonObjectInfoExtractorToCsv::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonObjectInfoExtractorToCsv);std::ostringstream oss;
