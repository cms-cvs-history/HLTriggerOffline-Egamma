////////////////////////////////////////////////////////////////////////////////
//                    Header file for this                                    //
////////////////////////////////////////////////////////////////////////////////
#include "HLTriggerOffline/Egamma/interface/EmDQM.h"

////////////////////////////////////////////////////////////////////////////////
//                    Collaborating Class Header                              //
////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
//                           Root include files                               //
////////////////////////////////////////////////////////////////////////////////
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include <Math/VectorUtil.h>
using namespace ROOT::Math::VectorUtil ;


////////////////////////////////////////////////////////////////////////////////
//                             Constructor                                    //
////////////////////////////////////////////////////////////////////////////////
EmDQM::EmDQM(const edm::ParameterSet& pset)  
{



  dbe = edm::Service < DQMStore > ().operator->();
  dbe->setVerbose(0);


  ////////////////////////////////////////////////////////////
  //          Read from configuration file                  //
  ////////////////////////////////////////////////////////////
  dirname_="HLT/HLTEgammaValidation/"+pset.getParameter<std::string>("@module_label");
  dbe->setCurrentFolder(dirname_);

  // paramters for generator study
  reqNum    = pset.getParameter<unsigned int>("reqNum");
  pdgGen    = pset.getParameter<int>("pdgGen");
  genEtaAcc = pset.getParameter<double>("genEtaAcc");
  genEtAcc  = pset.getParameter<double>("genEtAcc");
  // plotting paramters (untracked because they don't affect the physics)
  plotEtaMax = pset.getUntrackedParameter<double>("EtaMax", 4.0);
  plotPtMin  = pset.getUntrackedParameter<double>("PtMin" , 0.);
  plotPtMax  = pset.getUntrackedParameter<double>("PtMax" , 1000.);
  plotBins   = pset.getUntrackedParameter<unsigned int>("Nbins", 40);

  // preselction cuts
  gencutCollection_= pset.getParameter<edm::InputTag>("cutcollection");
  gencut_          = pset.getParameter<int>("cutnum");

  ////////////////////////////////////////////////////////////
  //         Read in the Vector of Parameter Sets.          //
  //           Information for each filter-step             //
  ////////////////////////////////////////////////////////////
  std::vector<edm::ParameterSet> filters = 
       pset.getParameter<std::vector<edm::ParameterSet> >("filters");

  for(std::vector<edm::ParameterSet>::iterator filterconf = filters.begin() ; filterconf != filters.end() ; filterconf++)
  {

    theHLTCollectionLabels.push_back(filterconf->getParameter<edm::InputTag>("HLTCollectionLabels"));
    theHLTOutputTypes.push_back(filterconf->getParameter<unsigned int>("theHLTOutputTypes"));
    std::vector<double> bounds = filterconf->getParameter<std::vector<double> >("PlotBounds");
    // If the size of plot "bounds" vector != 2, abort
    assert(bounds.size() == 2);
    plotBounds.push_back(std::pair<double,double>(bounds[0],bounds[1]));

    // Grab "IsoCollections" from config file
    isoNames.push_back(filterconf->getParameter<std::vector<edm::InputTag> >("IsoCollections"));
    // If the size of the isoNames vector is not greater than zero, abort
    assert(isoNames.back().size()>0);
    if (isoNames.back().at(0).label()=="none") {
      plotiso.push_back(false);
    } else {
      plotiso.push_back(true);
    }

  } // END of loop over parameter sets
}


////////////////////////////////////////////////////////////////////////////////
//       method called once each job just before starting event loop          //
////////////////////////////////////////////////////////////////////////////////
void 
EmDQM::beginJob(const edm::EventSetup&)
{
  //edm::Service<TFileService> fs;
  dbe->setCurrentFolder(dirname_);

  std::string histName  = "total eff";
  std::string histTitle = "Total Efficiency";

  total = dbe->book1D(histName.c_str(),histTitle.c_str(),theHLTCollectionLabels.size()+2,0,theHLTCollectionLabels.size()+2);
  total->setBinLabel(theHLTCollectionLabels.size()+1,"Total");
  total->setBinLabel(theHLTCollectionLabels.size()+2,"Gen");
  for (unsigned int u=0; u<theHLTCollectionLabels.size(); u++){total->setBinLabel(u+1,theHLTCollectionLabels[u].label().c_str());}

  MonitorElement* tmphisto;
  MonitorElement* tmpiso;

  // Generator-level histograms
  std::string pdgIdString= Form("w/pdgId=%d",pdgGen);

  histName  = "gen et";
  histTitle = "Et of Gen Particles "+pdgIdString;
  etgen     =       dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, plotPtMin, plotPtMax);
  histName  = "gen eta";
  histTitle = "Eta of Gen Particles "+pdgIdString;
  etagen    =       dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, -plotEtaMax, plotEtaMax);
  histName  = "gen et highest";
  histTitle = "Et of Highest Et Gen Particle "+pdgIdString;
  etgenHighestEt  = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, plotPtMin, plotPtMax);
  histName  = "gen eta highest";
  histTitle = "Eta of Highest Et Gen Particle "+pdgIdString;
  etagenHighestEt = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, -plotEtaMax, plotEtaMax); 
  histName  = "gen et 2nd highest";
  histTitle = "Et of 2nd Highest Et Gen Particle "+pdgIdString;
  etgen2ndHighestEt = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, plotPtMin, plotPtMax);
  histName  = "gen eta 2nd highest";
  histTitle = "Eta of 2nd Highest Et Gen Particle "+pdgIdString;
  etagen2ndHighestEt = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, -plotEtaMax, plotEtaMax);

  // HLT histograms
  for (unsigned int i = 0 ; i < theHLTCollectionLabels.size() ; i++){

    // Et distribution
    histName  = theHLTCollectionLabels[i].label()+"et";
    histTitle = theHLTCollectionLabels[i].label()+" Et";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax);
    ethist.push_back(tmphisto);
    
    // eta distribution
    histName  = theHLTCollectionLabels[i].label()+"eta";
    histTitle = theHLTCollectionLabels[i].label()+" eta";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
    etahist.push_back(tmphisto);

    // Et distribution of Monte-Carlo Matched objects
    histName  = theHLTCollectionLabels[i].label()+"et MC matched";
    histTitle = theHLTCollectionLabels[i].label()+" Et MC matched";
    tmphisto  = dbe->book1D(histName.c_str(),histName.c_str(),plotBins,plotPtMin,plotPtMax);
    ethistmatch.push_back(tmphisto);
    
    // eta distribution of Monte-Carlo Matched objects
    histName  = theHLTCollectionLabels[i].label()+"eta MC matched";
    histTitle = theHLTCollectionLabels[i].label()+" eta MC matched";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
    etahistmatch.push_back(tmphisto);          
   
    // Isolation values vs eta
    if (plotiso[i]) {
      histName  = theHLTCollectionLabels[i].label()+"eta isolation";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs eta";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax,plotBins,plotBounds[i].first,plotBounds[i].second);
    } else {
      tmpiso = NULL;
    }
    etahistiso.push_back(tmpiso);

    // Isolation values vs et
    if (plotiso[i]) {
      histName  = theHLTCollectionLabels[i].label()+"et isolation";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs Et";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax,plotBins,plotBounds[i].first,plotBounds[i].second);
    } else {
      tmpiso = NULL;
    }
    ethistiso.push_back(tmpiso);

  } // END of HLT histograms
}


////////////////////////////////////////////////////////////////////////////////
//                                Destructor                                  //
////////////////////////////////////////////////////////////////////////////////
EmDQM::~EmDQM(){
}


////////////////////////////////////////////////////////////////////////////////
//                     method called to for each event                        //
////////////////////////////////////////////////////////////////////////////////
void 
EmDQM::analyze(const edm::Event & event , const edm::EventSetup& setup)
{
  
  ////////////////////////////////////////////////////////////
  //           Check if there's enough gen particles        //
  //             of interest                                //
  ////////////////////////////////////////////////////////////
  edm::Handle< edm::View<reco::Candidate> > cutCounter;
  event.getByLabel(gencutCollection_,cutCounter);
  if (cutCounter->size() < (unsigned int)gencut_) {
    //edm::LogWarning("EmDQM") << "Less than "<< reqNum <<" gen particles with pdgId=" << pdgGen;
    return;
  }


  // fill L1 and HLT info
  // get objects possed by each filter
  edm::Handle<trigger::TriggerEventWithRefs> triggerObj;
  event.getByLabel("hltTriggerSummaryRAW",triggerObj); 
  if(!triggerObj.isValid()) { 
    edm::LogWarning("EmDQM") << "RAW-type HLT results not found, skipping event";
    return;
  }

  // total event number
  total->Fill(theHLTCollectionLabels.size()+0.5);


  ////////////////////////////////////////////////////////////
  //               Fill generator info                      //
  ////////////////////////////////////////////////////////////
  edm::Handle<edm::HepMCProduct> genEvt;
  event.getByLabel("source", genEvt);
  
  std::vector<HepMC::GenParticle> mcparts;
  const HepMC::GenEvent * myGenEvent = genEvt->GetEvent();
  unsigned int ncand = 0;
  float highestEtFound  = -1.0;
  float highestEtFound2 = -2.0;
  float etaOfHighestEtFound = -20.0;
  float etaOfHighestEtFound2 = -20.0;
  
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {

    // If the ID number is not what we're looking for or 
    //  it's status is !=1, go to the next particle
    if (  !( abs((*p)->pdg_id())==pdgGen  && (*p)->status()==1 )  )  continue;

    // Grab particle information
    float eta   = (*p)->momentum().eta();
    float e     = (*p)->momentum().e();
    float theta = 2*atan(exp(-eta));
    float Et    = e*sin(theta);

    if ( fabs(eta)<genEtaAcc  &&  Et > genEtAcc ) {
      // Store particle info if eta and Et within acceptance
      ncand++;
      etgen->Fill(Et);
      etagen->Fill(eta);
      mcparts.push_back(*(*p));
    }

    if ( Et > highestEtFound ) {
      // Store particle info if it's the highest Et found
      highestEtFound2      = highestEtFound;
      etaOfHighestEtFound2 = etaOfHighestEtFound;

      highestEtFound = Et;
      etaOfHighestEtFound = eta;
    } else if ( Et > highestEtFound2 ) {
      highestEtFound2 = Et;
      etaOfHighestEtFound2 = eta;
    }
  } // END of loop over Generated particles

  if (highestEtFound > 0.0) {
    // If we found a highest Et Gen particle, fill these
    etgenHighestEt->Fill(highestEtFound);
    etagenHighestEt->Fill(etaOfHighestEtFound);
  }
  if (highestEtFound2 > 0.0) {
    etgen2ndHighestEt->Fill(highestEtFound2);
    etagen2ndHighestEt->Fill(etaOfHighestEtFound2);
  }
  if (ncand >= reqNum) total->Fill(theHLTCollectionLabels.size()+1.5);
  ////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////
  //            Loop over filter modules                    //
  ////////////////////////////////////////////////////////////
  for(unsigned int n=0; n < theHLTCollectionLabels.size() ; n++) {
    // These numbers are from the Parameter Set, such as:
    //   theHLTOutputTypes = cms.uint32(100)
    switch(theHLTOutputTypes[n]) 
    {
      case 82: // Non-isolated Level 1
        fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,mcparts);break;
      case 83: // Isolated Level 1
        fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,mcparts);break;
      case 91: // Photon 
        fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,mcparts);break;
      case 92: // Electron 
        fillHistos<reco::ElectronCollection>(triggerObj,event,n,mcparts);break;
      case 100: // TriggerCluster
        fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,mcparts);break;
      default: 
        throw(cms::Exception("Release Validation Error") << "HLT output type not implemented: theHLTOutputTypes[n]" );
    }
  } // END of loop over filter modules
}


////////////////////////////////////////////////////////////////////////////////
// fillHistos                                                                 //
//   Called by analyze method.                                                //
//   
////////////////////////////////////////////////////////////////////////////////
template <class T> void EmDQM::fillHistos(edm::Handle<trigger::TriggerEventWithRefs>& triggerObj,const edm::Event& iEvent ,unsigned int n,std::vector<HepMC::GenParticle>& mcparts)
{
  std::vector<edm::Ref<T> > recoecalcands;
  if (!( triggerObj->filterIndex(theHLTCollectionLabels[n])>=triggerObj->size() )){ // only process if available
  
    ////////////////////////////////////////////////////////////
    //      Retrieve saved filter objects                     //
    ////////////////////////////////////////////////////////////
    triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),theHLTOutputTypes[n],recoecalcands);
    //Danger: special case, L1 non-isolated
    // needs to be merged with L1 iso
    if (theHLTOutputTypes[n] == 82)
    {
      std::vector<edm::Ref<T> > isocands;
      triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),83,isocands);
      if (isocands.size()>0) 
      {
        for (unsigned int i=0; i < isocands.size(); i++)
	  recoecalcands.push_back(isocands[i]);
      }
    } // END of if theHLTOutputTypes == 82

    ////////////////////////////////////////////////////////////
    //        Fill filter objects into histograms             //
    ////////////////////////////////////////////////////////////
    if (recoecalcands.size() != 0)
    {
      if (recoecalcands.size() >= reqNum ) 
	total->Fill(n+0.5);
      for (unsigned int i=0; i<recoecalcands.size(); i++) {

	//unmatched
	ethist[n]->Fill(recoecalcands[i]->et() );
	etahist[n]->Fill(recoecalcands[i]->eta() );

        // Loop over all Generated Particles
        //  to find one closest in delta-R
	math::XYZVector candDir=recoecalcands[i]->momentum();
	std::vector<HepMC::GenParticle>::iterator closest = mcparts.end();
	double closestDr = 1000. ;
	for(std::vector<HepMC::GenParticle>::iterator mc = mcparts.begin(); mc !=  mcparts.end() ; mc++){
	  math::XYZVector mcDir( mc->momentum().px(),
				 mc->momentum().py(),
				 mc->momentum().pz());
	  double dr = DeltaR(mcDir,candDir);
	  if (dr < closestDr) {
	    closestDr = dr;
	    closest = mc;
	  }
	} // END of loop over Gen Particles
	if (closest == mcparts.end())
          edm::LogWarning("EmDQM") << "Efficiency may become skewed: No match to gen particle with pdgId="<< pdgGen << " for RecoEcalCand with eta=" << recoecalcands[i]->eta() << ".";
	else {
	  float eta   = closest->momentum().eta();
	  float e     = closest->momentum().e();
	  float theta = 2*atan(exp(-eta));
	  float Et    = e*sin(theta);
	  ethistmatch[n]->Fill( Et );
	  etahistmatch[n]->Fill( eta );
	}

	////////////////////////////////////////////////////////////
	//  Plot isolation variables (show the not-yet-cut        //
        //  isolation, i.e. associated to next filter)            //
	////////////////////////////////////////////////////////////
	if (n+1 < theHLTCollectionLabels.size()) // can't plot beyond last
        {
	  if (plotiso[n+1] ){
	    for (unsigned int j =  0 ; j < isoNames[n+1].size() ;j++  ){
	      edm::Handle<edm::AssociationMap<edm::OneToValue< T , float > > > depMap; 
	      iEvent.getByLabel(isoNames[n+1].at(j).label(),depMap);
	      if (depMap.isValid()){ //Map may not exist if only one candidate passes a double filter
		typename edm::AssociationMap<edm::OneToValue< T , float > >::const_iterator mapi = depMap->find(recoecalcands[i]);
		if (mapi!=depMap->end()){  // found candidate in isolation map! 
		  etahistiso[n+1]->Fill(recoecalcands[i]->eta(),mapi->val);
		  ethistiso[n+1]->Fill(recoecalcands[i]->et(),mapi->val);
		  break; // to avoid multiple filling we only look until we found the candidate once.
		}
	      }
	    }
	  }
	} // END of if n+1 < then the number of hlt collections
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////// 
//      method called once each job just after ending the event loop          //
//////////////////////////////////////////////////////////////////////////////// 
void EmDQM::endJob(){
  // Normalize the histograms
  //  total->Scale(1./total->GetBinContent(1));
  //for(unsigned int n= theHLTCollectionLabels.size()-1 ; n>0;n--){
  //  ethist[n]->Divide(ethist[n-1]);
  //  etahist[n]->Divide(etahist[n-1]);
  //}
}

DEFINE_FWK_MODULE(EmDQM);
