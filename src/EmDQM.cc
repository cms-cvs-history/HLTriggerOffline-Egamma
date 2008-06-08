#include "HLTriggerOffline/Egamma/interface/EmDQM.h"

// Collaborating Class Header
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include <Math/VectorUtil.h>
using namespace ROOT::Math::VectorUtil ;

/// Constructor
EmDQM::EmDQM(const edm::ParameterSet& pset)  
{

  //paramters for generator study
  reqNum = pset.getParameter<unsigned int>("reqNum");
  pdgGen =  pset.getParameter<int>("pdgGen");
  genEtaAcc = pset.getParameter<double>("genEtaAcc");
  genEtAcc = pset.getParameter<double>("genEtAcc");
  //plotting paramters
  thePtMin = pset.getUntrackedParameter<double>("PtMin",0.);
  thePtMax = pset.getUntrackedParameter<double>("PtMax",1000.);
  theNbins = pset.getUntrackedParameter<unsigned int>("Nbins",40);
  
  //info for each filter-step
  std::vector<edm::ParameterSet> filters = pset.getParameter<std::vector<edm::ParameterSet> >("filters");

  for(std::vector<edm::ParameterSet>::iterator filterconf = filters.begin() ; filterconf != filters.end() ; filterconf++){
    theHLTCollectionLabels.push_back(filterconf->getParameter<edm::InputTag>("HLTCollectionLabels"));
    theHLTOutputTypes.push_back(filterconf->getParameter<unsigned int>("theHLTOutputTypes"));
    std::vector<double> bounds = filterconf->getParameter<std::vector<double> >("PlotBounds");
    assert(bounds.size() == 2);
    plotBounds.push_back(std::pair<double,double>(bounds[0],bounds[1]));
    isoNames.push_back(filterconf->getParameter<std::vector<edm::InputTag> >("IsoCollections"));
    assert(isoNames.back().size()>0);
    if (isoNames.back().at(0).label()=="none")
      plotiso.push_back(false);
    else{
       plotiso.push_back(true);
       //std::cout << "plotting isolation for: " <<  isoNames.back().at(0).label() << std::endl;
    }


  }
}


void EmDQM::beginJob(const edm::EventSetup&){
  edm::Service<TFileService> fs;
  std::string histoname="total eff";

  total = fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theHLTCollectionLabels.size()+2,0,theHLTCollectionLabels.size()+2);
  total->GetXaxis()->SetBinLabel(theHLTCollectionLabels.size()+1,"Total");
  total->GetXaxis()->SetBinLabel(theHLTCollectionLabels.size()+2,"Gen");
  for (unsigned int u=0; u<theHLTCollectionLabels.size(); u++){total->GetXaxis()->SetBinLabel(u+1,theHLTCollectionLabels[u].label().c_str());}

  TH1F* tmphisto;
  TH2F* tmpiso;

  histoname = "gen et";
  etgen =  fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theNbins,thePtMin,thePtMax);
  histoname = "gen eta";
  etagen = fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theNbins,-2.7,2.7);
 
  for(unsigned int i = 0; i< theHLTCollectionLabels.size() ; i++){
    histoname = theHLTCollectionLabels[i].label()+"et";
    tmphisto =  fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theNbins,thePtMin,thePtMax);
    ethist.push_back(tmphisto);
    
    histoname = theHLTCollectionLabels[i].label()+"eta";
    tmphisto =  fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theNbins,-2.7,2.7);
    etahist.push_back(tmphisto);          

    histoname = theHLTCollectionLabels[i].label()+"et MC matched";
    tmphisto =  fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theNbins,thePtMin,thePtMax);
    ethistmatch.push_back(tmphisto);
    
    histoname = theHLTCollectionLabels[i].label()+"eta MC matched";
    tmphisto =  fs->make<TH1F>(histoname.c_str(),histoname.c_str(),theNbins,-2.7,2.7);
    etahistmatch.push_back(tmphisto);          
    
    if(plotiso[i]){
      histoname = theHLTCollectionLabels[i].label()+"eta isolation";
      tmpiso = fs->make<TH2F>(histoname.c_str(),histoname.c_str(),theNbins,-2.7,2.7,theNbins,plotBounds[i].first,plotBounds[i].second);
    }
    else{
      tmpiso = NULL;
    }
    etahistiso.push_back(tmpiso);

    if(plotiso[i]){
      histoname = theHLTCollectionLabels[i].label()+"et isolation";
      tmpiso = fs->make<TH2F>(histoname.c_str(),histoname.c_str(),theNbins,thePtMin,thePtMax,theNbins,plotBounds[i].first,plotBounds[i].second);
    }
    else{
      tmpiso = NULL;
    }
    ethistiso.push_back(tmpiso);

  }
}


/// Destructor
EmDQM::~EmDQM(){
}

void EmDQM::analyze(const edm::Event & event , const edm::EventSetup& setup){



  // fill L1 and HLT info
  // get objects possed by each filter
  edm::Handle<trigger::TriggerEventWithRefs> triggerObj;
  event.getByLabel("hltTriggerSummaryRAW",triggerObj); 
  if(!triggerObj.isValid()) { 
    edm::LogWarning("EmDQM") << "RAW-type HLT results not found, skipping event";
    return;
  }
// throw(cms::Exception("Release Validation Error")<< "RAW-type HLT results not found" );



//  for(int i = 0; i<triggerObj->size() ;i++ ){
//    std::cout << triggerObj->filterTag(i) << std::endl;
//  }

  // total event number
  total->Fill(theHLTCollectionLabels.size()+0.5);

  // fill generator info
  edm::Handle<edm::HepMCProduct> genEvt;
  event.getByLabel("source", genEvt);
  
  std::vector<HepMC::GenParticle> mcparts;
  const HepMC::GenEvent * myGenEvent = genEvt->GetEvent();
  unsigned int ncand = 0;
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
    if (  !( abs((*p)->pdg_id())==pdgGen  && (*p)->status()==1 )   )  continue;
    float eta   =(*p)->momentum().eta();
    float e     =(*p)->momentum().e();
    float theta =2*atan(exp(-eta));
    float Et    =e*sin(theta);
    if(fabs(eta)<genEtaAcc  &&  Et > genEtAcc) {
      ncand++;
      etgen->Fill(Et);
      etagen->Fill(eta);
      mcparts.push_back(*(*p));
    }
  }//end of loop over MC particles
  if (ncand >= reqNum) total->Fill(theHLTCollectionLabels.size()+1.5);
	  


  
  for(unsigned int n=0; n < theHLTCollectionLabels.size() ; n++) { //loop over filter modules
    switch(theHLTOutputTypes[n]){
    case 82: // non-iso L1
      fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,mcparts);break;
    case 83: // iso L1
      fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,mcparts);break;
    case 91: //photon 
      fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,mcparts);break;
    case 92: //electron 
      fillHistos<reco::ElectronCollection>(triggerObj,event,n,mcparts);break;
    case 100: // TriggerCluster
      fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,mcparts);break;
    default: throw(cms::Exception("Release Validation Error")<< "HLT output type not implemented: theHLTOutputTypes[n]" );
    }
  }
}

template <class T> void EmDQM::fillHistos(edm::Handle<trigger::TriggerEventWithRefs>& triggerObj,const edm::Event& iEvent ,unsigned int n,std::vector<HepMC::GenParticle>& mcparts){
  
  std::vector<edm::Ref<T> > recoecalcands;
  if (!( triggerObj->filterIndex(theHLTCollectionLabels[n])>=triggerObj->size() )){ // only process if availabel
  
    // retrieve saved filter objects
    triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),theHLTOutputTypes[n],recoecalcands);
    //Danger: special case, L1 non-isolated
    // needs to be merged with L1 iso
    if(theHLTOutputTypes[n]==82){
      std::vector<edm::Ref<T> > isocands;
      triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),83,isocands);
      if(isocands.size()>0)
	for(unsigned int i=0; i < isocands.size(); i++)
	  recoecalcands.push_back(isocands[i]);
    }


    //fill filter objects into histos
    if (recoecalcands.size()!=0){
      if(recoecalcands.size() >= reqNum ) 
	total->Fill(n+0.5);
      for (unsigned int i=0; i<recoecalcands.size(); i++) {
	//unmatched
	ethist[n]->Fill(recoecalcands[i]->et() );
	etahist[n]->Fill(recoecalcands[i]->eta() );
	//matched
	math::XYZVector candDir=recoecalcands[i]->momentum();
	std::vector<HepMC::GenParticle>::iterator closest = mcparts.end();
	double closestDr=1000. ;
	for(std::vector<HepMC::GenParticle>::iterator mc = mcparts.begin(); mc !=  mcparts.end() ; mc++){
	  math::XYZVector mcDir( mc->momentum().px(),
				 mc->momentum().py(),
				 mc->momentum().pz());	  
	  double dr = DeltaR(mcDir,candDir);
	  if(dr < closestDr){
	    closestDr = dr;
	    closest = mc;
	  }
	}
	if (closest == mcparts.end())  edm::LogWarning("EmDQM") << "no MC match, may skew efficieny";
	else{
	  float eta   =closest->momentum().eta();
	  float e     =closest->momentum().e();
	  float theta =2*atan(exp(-eta));
	  float Et    =e*sin(theta);
	  ethistmatch[n]->Fill( Et );
	  etahistmatch[n]->Fill( eta );
	}
	
	//plot isolation variables (show not yet cut  iso, i.e. associated to next filter)
	if(n+1 < theHLTCollectionLabels.size()){ // can't plot beyond last
	  if(plotiso[n+1] ){
	    for(unsigned int j =  0 ; j < isoNames[n+1].size() ;j++  ){
	      edm::Handle<edm::AssociationMap<edm::OneToValue< T , float > > > depMap; 
	      if(depMap.isValid()){ //Map may not exist if only one candidate apsses a double filter
		iEvent.getByLabel(isoNames[n+1].at(j).label(),depMap);
		typename edm::AssociationMap<edm::OneToValue< T , float > >::const_iterator mapi = depMap->find(recoecalcands[i]);
		if(mapi!=depMap->end()){  // found candidate in isolation map! 
		  etahistiso[n+1]->Fill(recoecalcands[i]->eta(),mapi->val);
		ethistiso[n+1]->Fill(recoecalcands[i]->et(),mapi->val);
		break; // to avoid multiple filling we only look until we found the candidate once.
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void EmDQM::endJob(){
  //  total->Scale(1./total->GetBinContent(1));
  //for(unsigned int n= theHLTCollectionLabels.size()-1 ; n>0;n--){
  //  ethist[n]->Divide(ethist[n-1]);
  //  etahist[n]->Divide(etahist[n-1]);
  //}
}

DEFINE_FWK_MODULE(EmDQM);
