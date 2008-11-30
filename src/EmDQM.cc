///////////////////////////////////////////////////////////////////////////////
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
  plotPtMin  = pset.getUntrackedParameter<double>("PtMin",0.);
  plotPtMax  = pset.getUntrackedParameter<double>("PtMax",1000.);
  plotEtaMax = pset.getUntrackedParameter<double>("EtaMax", 2.7);
  plotBins  = pset.getUntrackedParameter<unsigned int>("Nbins",40);

  //preselction cuts 
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
    isoNames.push_back(filterconf->getParameter<std::vector<edm::InputTag> >("IsoCollections"));
    // If the size of the isoNames vector is not greater than zero, abort
    assert(isoNames.back().size()>0);
    if (isoNames.back().at(0).label()=="none")
      plotiso.push_back(false);
    else
      plotiso.push_back(true);

  } // END of loop over parameter sets

  // Record number of HLTCollectionLabels
  numOfHLTCollectionLabels = theHLTCollectionLabels.size();
  
}


////////////////////////////////////////////////////////////////////////////////
//       method called once each job just before starting event loop          //
////////////////////////////////////////////////////////////////////////////////
void 
EmDQM::beginJob(const edm::EventSetup&)
{
  //edm::Service<TFileService> fs;
  dbe->setCurrentFolder(dirname_);

  ////////////////////////////////////////////////////////////
  //  Set up Histogram of Effiency vs Step.                 //
  //   theHLTCollectionLabels is a vector of InputTags      //
  //    from the configuration file.                        //
  ////////////////////////////////////////////////////////////

  std::string histName="total eff";
  std::string histTitle = "total events passing";
  // This plot will have bins equal to 2+(number of
  //        HLTCollectionLabels in the config file)
  total = dbe->book1D(histName.c_str(),histTitle.c_str(),numOfHLTCollectionLabels+2,0,numOfHLTCollectionLabels+2);
  total->setBinLabel(numOfHLTCollectionLabels+1,"Total");
  total->setBinLabel(numOfHLTCollectionLabels+2,"Gen");
  for (unsigned int u=0; u<numOfHLTCollectionLabels; u++){total->setBinLabel(u+1,theHLTCollectionLabels[u].label().c_str());}

  histName="total eff MC matched";
  histTitle="total events passing (mc matched)";
  totalmatch = dbe->book1D(histName.c_str(),histTitle.c_str(),numOfHLTCollectionLabels+2,0,numOfHLTCollectionLabels+2);
  totalmatch->setBinLabel(numOfHLTCollectionLabels+1,"Total");
  totalmatch->setBinLabel(numOfHLTCollectionLabels+2,"Gen");
  for (unsigned int u=0; u<numOfHLTCollectionLabels; u++){totalmatch->setBinLabel(u+1,theHLTCollectionLabels[u].label().c_str());}

  MonitorElement* tmphisto;
  MonitorElement* tmpiso;

  ////////////////////////////////////////////////////////////
  // Set up generator-level histograms                      //
  ////////////////////////////////////////////////////////////
  std::string pdgIdString;
  switch(pdgGen) {
  case 11:
    pdgIdString="Electron";break;
  case 22:
    pdgIdString="Photon";break;
  default:
    pdgIdString="Particle";
  }

  histName = "gen et";
  histTitle= "E_{T} of " + pdgIdString + "s" ;
  etgen =  dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax);
  histName = "gen eta";
  histTitle= "E_{T} of "+ pdgIdString +"s " ;
  etagen = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
 
  ////////////////////////////////////////////////////////////
  //  Set up histograms of HLT objects                      //
  ////////////////////////////////////////////////////////////
 
  for(unsigned int i = 0; i< numOfHLTCollectionLabels ; i++){
    // Et distribution of HLT objects passing filter i
    histName = theHLTCollectionLabels[i].label()+"et";
    histTitle = theHLTCollectionLabels[i].label()+" Et";
    tmphisto =  dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax);
    ethist.push_back(tmphisto);
    
    // Eta distribution of HLT objects passing filter i
    histName = theHLTCollectionLabels[i].label()+"eta";
    histTitle = theHLTCollectionLabels[i].label()+" #eta";
    tmphisto =  dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
    etahist.push_back(tmphisto);          

    // Et distribution of gen object matching HLT object passing filter i
    histName = theHLTCollectionLabels[i].label()+"et MC matched";
    histTitle = theHLTCollectionLabels[i].label()+" Et MC matched";
    tmphisto =  dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax);
    ethistmatch.push_back(tmphisto);
    
    // Eta distribution of gen object matching HLT object passing filter i
    histName = theHLTCollectionLabels[i].label()+"eta MC matched";
    histTitle = theHLTCollectionLabels[i].label()+" #eta MC matched";
    tmphisto =  dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
    etahistmatch.push_back(tmphisto);          
    
  if (!plotiso[i]) {
      tmpiso = NULL;
      etahistiso.push_back(tmpiso);
      ethistiso.push_back(tmpiso);
      etahistisomatch.push_back(tmpiso);
      ethistisomatch.push_back(tmpiso);
    } else {
      // 2D plot: Isolation values vs eta for all objects
      histName  = theHLTCollectionLabels[i].label()+"eta isolation";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs #eta";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      etahistiso.push_back(tmpiso);

      // 2D plot: Isolation values vs et for all objects
      histName  = theHLTCollectionLabels[i].label()+"et isolation";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs Et for highest Et obj";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      ethistiso.push_back(tmpiso);

      // 2D plot: Isolation values vs eta for matched objects
      histName  = theHLTCollectionLabels[i].label()+"eta isolation MC matched";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs #eta (mc matched)";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      etahistisomatch.push_back(tmpiso);

      // 2D plot: Isolation values vs et for matched objects
      histName  = theHLTCollectionLabels[i].label()+"et isolation MC matched";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs Et for 2nd highest Et obj";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      ethistisomatch.push_back(tmpiso);
    } // END of HLT histograms

  }
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
    //edm::LogWarning("EmDQM") << "Less than "<< gencut_ <<" gen particles with pdgId=" << pdgGen;
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

  ////////////////////////////////////////////////////////////
  //  Fill the bin labeled "Total"                          //
  //   This will be the number of events looked at.         //
  //////////////////////////////////////////////////////////// 
  total->Fill(numOfHLTCollectionLabels+0.5);
  totalmatch->Fill(numOfHLTCollectionLabels+0.5);


  ////////////////////////////////////////////////////////////
  //               Fill generator info                      //
  ////////////////////////////////////////////////////////////
  // the gencut_ highest Et generator objects of the preselected type are our matches

  std::vector<reco::Particle> sortedGen;
  for(edm::View<reco::Candidate>::const_iterator genpart = cutCounter->begin(); genpart != cutCounter->end();genpart++){
    reco::Particle tmpcand( *(genpart) );
    sortedGen.push_back(tmpcand);
  }
  std::sort(sortedGen.begin(),sortedGen.end(),pTComparator_ );

  for (unsigned int i = 0 ; i < gencut_ ; i++ ) {
    etgen->Fill( sortedGen[i].et() ); //validity has been implicitily checked by the cut on gencut_ above
    etagen->Fill( sortedGen[i].eta() );
  } // END of loop over Generated particles
  if (gencut_ >= reqNum) total->Fill(numOfHLTCollectionLabels+1.5); // this isn't really needed anymore keep for backward comp.
  if (gencut_ >= reqNum) totalmatch->Fill(numOfHLTCollectionLabels+1.5); // this isn't really needed anymore keep for backward comp.
	  

  ////////////////////////////////////////////////////////////
  //            Loop over filter modules                    //
  ////////////////////////////////////////////////////////////
  for(unsigned int n=0; n < numOfHLTCollectionLabels ; n++) {
    // These numbers are from the Parameter Set, such as:
    //   theHLTOutputTypes = cms.uint32(100)
    switch(theHLTOutputTypes[n]) 
    {
      case 82: // Non-isolated Level 1
        fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,sortedGen);break;
      case 83: // Isolated Level 1
        fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,sortedGen);break;
      case 91: // Photon 
        fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,sortedGen);break;
      case 92: // Electron 
        fillHistos<reco::ElectronCollection>(triggerObj,event,n,sortedGen);break;
      case 100: // TriggerCluster
        fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,sortedGen);break;
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
template <class T> void EmDQM::fillHistos(edm::Handle<trigger::TriggerEventWithRefs>& triggerObj,const edm::Event& iEvent ,unsigned int n,std::vector<reco::Particle>& sortedGen)
{
  std::vector<edm::Ref<T> > recoecalcands;
  if ( ( triggerObj->filterIndex(theHLTCollectionLabels[n])>=triggerObj->size() )){ // only process if available
    return;
  }

  ////////////////////////////////////////////////////////////
  //      Retrieve saved filter objects                     //
  ////////////////////////////////////////////////////////////
  triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),theHLTOutputTypes[n],recoecalcands);
  //Danger: special case, L1 non-isolated
  // needs to be merged with L1 iso
  if (theHLTOutputTypes[n] == 82){
    std::vector<edm::Ref<T> > isocands;
    triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),83,isocands);
    if (isocands.size()>0) 
      {
	for (unsigned int i=0; i < isocands.size(); i++)
	  recoecalcands.push_back(isocands[i]);
      }
  } // END of if theHLTOutputTypes == 82
  

  if (recoecalcands.size() < 1){ // stop if no object passed the previous filter
    return;
  }


  ////////////////////////////////////////////////////////////
  //        Fill all filter objects into histograms         //
  ////////////////////////////////////////////////////////////
  if (recoecalcands.size() >= reqNum ) 
    total->Fill(n+0.5);
  for (unsigned int i=0; i<recoecalcands.size(); i++) {
    //unmatched
    ethist[n]->Fill(recoecalcands[i]->et() );
    etahist[n]->Fill(recoecalcands[i]->eta() );
    ////////////////////////////////////////////////////////////
    //  Plot isolation variables (show the not-yet-cut        //
    //  isolation, i.e. associated to next filter)            //
    ////////////////////////////////////////////////////////////
    if (n+1 < numOfHLTCollectionLabels){ // can't plot beyond last
      if (plotiso[n+1] ){  // only plot if requested in config
	for (unsigned int j =  0 ; j < isoNames[n+1].size() ;j++  ){
	  edm::Handle<edm::AssociationMap<edm::OneToValue< T , float > > > depMap; 
	  iEvent.getByLabel(isoNames[n+1].at(j).label(),depMap);
	  if (depMap.isValid()){ //Map may not exist if only one candidate passes a double filter
	    typename edm::AssociationMap<edm::OneToValue< T , float > >::const_iterator mapi = depMap->find(recoecalcands[i]);
	    if (mapi!=depMap->end()){  // found candidate in isolation map! 
	      etahistiso[n+1]->Fill(recoecalcands[i]->eta(),mapi->val);
	      ethistiso[n+1]->Fill(recoecalcands[i]->et(),mapi->val);
	    }
	  }
	}
      }
    } // END of if n+1 < then the number of hlt collections
  }


  ////////////////////////////////////////////////////////////
  //        Fill mc matched objects into histograms         //
  ////////////////////////////////////////////////////////////
  unsigned int mtachedMcParts = 0;
  float mindist=0.3;
  if(n==0) mindist=0.5; //low L1-resolution => allow wider matching 
  for(unsigned int i =0; i < gencut_; i++){
    //match generator candidate    
    bool matchThis= false;
    math::XYZVector candDir=sortedGen[i].momentum();
    unsigned int closest = 0;
    double closestDr = 1000.;
    for(unsigned int trigOb = 0 ; trigOb < recoecalcands.size(); trigOb++){
      double dr = DeltaR(recoecalcands[trigOb]->momentum(),candDir);
      if (dr < closestDr) {
	closestDr = dr;
	closest = trigOb;
      }
      if(closestDr > mindist){ // it's not really a "match" if it's that far away
	closest = -1;
      }
      else{
	mtachedMcParts++;
	matchThis = true;
      }
    }
    if(! matchThis ) break; // only plot matched candidates
    //fill coordinates of mc particle matching trigger object
    ethistmatch[n]->Fill(sortedGen[i].et() );
    etahistmatch[n]->Fill(sortedGen[i].eta() );
    ////////////////////////////////////////////////////////////
    //  Plot isolation variables (show the not-yet-cut        //
    //  isolation, i.e. associated to next filter)            //
    ////////////////////////////////////////////////////////////
    if (n+1 < numOfHLTCollectionLabels){ // can't plot beyond last
      if (plotiso[n+1] ){  // only plot if requested in config
	for (unsigned int j =  0 ; j < isoNames[n+1].size() ;j++  ){
	  edm::Handle<edm::AssociationMap<edm::OneToValue< T , float > > > depMap; 
	  iEvent.getByLabel(isoNames[n+1].at(j).label(),depMap);
	  if (depMap.isValid()){ //Map may not exist if only one candidate passes a double filter
	    typename edm::AssociationMap<edm::OneToValue< T , float > >::const_iterator mapi = depMap->find(recoecalcands[closest]);
	    if (mapi!=depMap->end()){  // found candidate in isolation map! 
	      etahistisomatch[n+1]->Fill(sortedGen[i].eta(),mapi->val);
	      ethistisomatch[n+1]->Fill(sortedGen[i].et(),mapi->val);
	    }
	  }
	}
      }
    } // END of if n+1 < then the number of hlt collections
  }
  // fill total mc matched efficiency
  if (mtachedMcParts >= reqNum ) 
    totalmatch->Fill(n+0.5);
  

}


//////////////////////////////////////////////////////////////////////////////// 
//      method called once each job just after ending the event loop          //
//////////////////////////////////////////////////////////////////////////////// 
void EmDQM::endJob(){

}

DEFINE_FWK_MODULE(EmDQM);
