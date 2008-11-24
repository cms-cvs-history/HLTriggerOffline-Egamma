////////////////////////////////////////////////////////////////////////////////
//   Authors:  Matthias Mozer,
//             Michael B. Anderson
////////////////////////////////////////////////////////////////////////////////

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
  plotEtaMax = pset.getUntrackedParameter<double>("EtaMax", 3.0);
  plotPtMin  = pset.getUntrackedParameter<double>("PtMin" , 0.);
  plotPtMax  = pset.getUntrackedParameter<double>("PtMax" , 1000.);
  plotBins   = pset.getUntrackedParameter<unsigned int>("Nbins", 48);

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

  std::string histName  = "total eff";
  std::string histTitle = "Highest Et Object";
  // This plot will have bins equal to 2+(number of 
  //        HLTCollectionLabels in the config file)
  histHighestEt = dbe->book1D(histName.c_str(),histTitle.c_str(),numOfHLTCollectionLabels+2,0,numOfHLTCollectionLabels+2);
  histHighestEt->setBinLabel(numOfHLTCollectionLabels+1,"Total");
  histHighestEt->setBinLabel(numOfHLTCollectionLabels+2,"Gen");
  // Set the bin labels
  for (unsigned int bin=0; bin<numOfHLTCollectionLabels; bin++){
    histHighestEt->setBinLabel(bin+1, theHLTCollectionLabels[bin].label().c_str());
  }

  histName  = "total eff 2";
  histTitle = "2^{nd} Highest Et Object";
  histHighestEt2 = dbe->book1D(histName.c_str(),histTitle.c_str(),numOfHLTCollectionLabels+2,0,numOfHLTCollectionLabels+2);
  histHighestEt2->setBinLabel(numOfHLTCollectionLabels+1,"Total");
  histHighestEt2->setBinLabel(numOfHLTCollectionLabels+2,"Gen");
  // Set the bin labels
  for (unsigned int bin=0; bin<numOfHLTCollectionLabels; bin++){
    histHighestEt2->setBinLabel(bin+1, theHLTCollectionLabels[bin].label().c_str());
  }


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

  histName  = "gen et";
  histTitle = "E_{T} of all Generated "+pdgIdString+"s";
  etgen     =       dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, plotPtMin, plotPtMax);
  histName  = "gen eta";
  histTitle = "#eta of all Generated "+pdgIdString+"s";
  etagen    =       dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, -plotEtaMax, plotEtaMax);
  histName  = "gen et highest";
  histTitle = "E_{T} of Highest Et Gen "+pdgIdString;
  etgenHighestEt  = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, plotPtMin, plotPtMax);
  histName  = "gen eta highest";
  histTitle = "#eta of Highest Et Gen "+pdgIdString;
  etagenHighestEt = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, -plotEtaMax, plotEtaMax); 
  histName  = "gen et 2nd highest";
  histTitle = "E_{T} of 2^{nd} Highest Et Gen "+pdgIdString;
  etgen2ndHighestEt = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, plotPtMin, plotPtMax);
  histName  = "gen eta 2nd highest";
  histTitle = "#eta of 2^{nd} Highest Et Gen "+pdgIdString;
  etagen2ndHighestEt = dbe->book1D(histName.c_str(), histTitle.c_str(), plotBins, -plotEtaMax, plotEtaMax);


  ////////////////////////////////////////////////////////////
  //  Set up histograms of HLT objects                      //
  ////////////////////////////////////////////////////////////
  MonitorElement* tmphisto;
  MonitorElement* tmpiso;

  for (unsigned int i = 0 ; i < numOfHLTCollectionLabels ; i++){

    // Et distribution of HLT object matched to highest Et gen particle
    histName  = theHLTCollectionLabels[i].label()+"et";
    histTitle = theHLTCollectionLabels[i].label()+" Et";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax);
    histEtOfHighEtMatch.push_back(tmphisto);
    
    // eta distribution of HLT object matched to highest Et gen particle
    histName  = theHLTCollectionLabels[i].label()+"eta";
    histTitle = theHLTCollectionLabels[i].label()+" eta";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
    histEtaOfHighEtMatch.push_back(tmphisto);

    // Et distribution of HLT object matched to 2nd highest Et gen particle
    histName  = theHLTCollectionLabels[i].label()+"et2";
    histTitle = theHLTCollectionLabels[i].label()+" Et";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax);
    histEtOfHighEtMatch2.push_back(tmphisto);

    // eta distribution of HLT object matched to 2nd highest Et gen particle
    histName  = theHLTCollectionLabels[i].label()+"eta2";
    histTitle = theHLTCollectionLabels[i].label()+" eta";
    tmphisto  = dbe->book1D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax);
    histEtaOfHighEtMatch2.push_back(tmphisto);

    // Plot Isolation vs et and eta?
    if (!plotiso[i]) {
      tmpiso = NULL;
      histIsoVsEtaOfHighEtMatch.push_back(tmpiso);
      histIsoVsEtOfHighEtMatch.push_back(tmpiso);
      histIsoVsEtaOfHighEtMatch2.push_back(tmpiso);
      histIsoVsEtOfHighEtMatch2.push_back(tmpiso);
    } else {
      // 2D plot: Isolation values vs eta for highest Et obj
      histName  = theHLTCollectionLabels[i].label()+"eta isolation";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs eta for highest Et obj";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      histIsoVsEtaOfHighEtMatch.push_back(tmpiso);

      // 2D plot: Isolation values vs et for highest Et obj
      histName  = theHLTCollectionLabels[i].label()+"et isolation";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs Et for highest Et obj";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      histIsoVsEtOfHighEtMatch.push_back(tmpiso);

      // 2D plot: Isolation values vs eta for 2nd highest Et obj
      histName  = theHLTCollectionLabels[i].label()+"eta isolation2";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs eta for 2nd highest Et obj";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,-plotEtaMax,plotEtaMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      histIsoVsEtaOfHighEtMatch2.push_back(tmpiso);

      // 2D plot: Isolation values vs et for 2nd highest Et obj
      histName  = theHLTCollectionLabels[i].label()+"et isolation2";
      histTitle = theHLTCollectionLabels[i].label()+" isolation vs Et for 2nd highest Et obj";
      tmpiso    = dbe->book2D(histName.c_str(),histTitle.c_str(),plotBins,plotPtMin,plotPtMax,plotBins,plotBounds[i].first,plotBounds[i].second);
      histIsoVsEtOfHighEtMatch2.push_back(tmpiso);
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
    //edm::LogWarning("EmDQM") << "Less than "<< reqNum <<" gen particles with pdgId=" << pdgGen;
    return;
  }


  // fill L1 and HLT info
  // get objects possed by each filter
  edm::Handle<trigger::TriggerEventWithRefs> triggerObj;
  event.getByLabel("hltTriggerSummaryRAW",triggerObj); 
  if (!triggerObj.isValid()) { 
    edm::LogWarning("EmDQM") << "RAW-type HLT results not found, skipping event";
    return;
  }

  ////////////////////////////////////////////////////////////
  //  Fill the bin labeled "Total"                          //
  //   This will be the number of events looked at.         //
  ////////////////////////////////////////////////////////////
  histHighestEt->Fill(numOfHLTCollectionLabels+0.5);
  histHighestEt2->Fill(numOfHLTCollectionLabels+0.5);


  ////////////////////////////////////////////////////////////
  //               Fill generator info                      //
  ////////////////////////////////////////////////////////////
  edm::Handle<edm::HepMCProduct> genEvt;
  event.getByLabel("source", genEvt);

  // Store some of the gnerated particles  
  std::vector<HepMC::GenParticle> genParticlesPdgIdMatch;
  HepMC::GenParticle genParticleHiEt;
  HepMC::GenParticle genParticleHiEt2;

  const HepMC::GenEvent * myGenEvent = genEvt->GetEvent();
  unsigned int numOfGenParticlesPdgIdMatch = 0;
  float highestEtFound  = -1.0;
  float highestEtFound2 = -2.0;
  float etaOfHighestEtFound = -20.0;
  float etaOfHighestEtFound2 = -20.0;
  
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {

    // Grab particle information
    float eta   = (*p)->momentum().eta();
    float e     = (*p)->momentum().e();
    float theta = 2*atan(exp(-eta));
    float Et    = e*sin(theta);

    // If the ID number is not what we're looking for or 
    //  it's status is !=1, or Et is too low, go to the next particle
    if (  !( abs((*p)->pdg_id())==pdgGen  && (*p)->status()==1 && Et > genEtAcc)  )  continue;

    if ( fabs(eta)<genEtaAcc ) {
      // Store particle info if eta within acceptance
      numOfGenParticlesPdgIdMatch++;
      etgen->Fill(Et);
      etagen->Fill(eta);
    }

    if ( Et > highestEtFound ) {
      // Store particle info if it's the highest Et found
      highestEtFound2      = highestEtFound;
      etaOfHighestEtFound2 = etaOfHighestEtFound;
      genParticleHiEt2 = genParticleHiEt;

      highestEtFound = Et;
      etaOfHighestEtFound = eta;
      genParticleHiEt = *(*p);
    } else if ( Et > highestEtFound2 ) {
      highestEtFound2 = Et;
      etaOfHighestEtFound2 = eta;
      genParticleHiEt2 = *(*p);
    }
  } // END of loop over Generated particles

  // If the highest Et gen particle went out of the
  // detector acceptance don't bother going on
  if ( !( fabs(etaOfHighestEtFound)<genEtaAcc  && highestEtFound > genEtAcc ) ) {
    edm::LogWarning("EmDQM") << "Not continuing with event. Highest Et gen particle with pdgGen="<<pdgGen << " had eta="<<etaOfHighestEtFound << " which is beyond acceptance.";
    return;
  }

  // Fill histograms
  if (highestEtFound > 0.0) {
    // If we found a highest Et Gen particle, fill these
    etgenHighestEt->Fill(highestEtFound);
    etagenHighestEt->Fill(etaOfHighestEtFound);
    genParticlesPdgIdMatch.push_back(genParticleHiEt);

    // Fill total plot, bin labled "gen"
    histHighestEt->Fill(numOfHLTCollectionLabels+1.5);
  }
  if (highestEtFound2 > 0.0) {
    // If we found a 2nd highest Et Gen particle, fill these
    etgen2ndHighestEt->Fill(highestEtFound2);
    etagen2ndHighestEt->Fill(etaOfHighestEtFound2);
    genParticlesPdgIdMatch.push_back(genParticleHiEt2);

    // Fill total plot, bin labeled "gen"
    histHighestEt2->Fill(numOfHLTCollectionLabels+1.5);
  }


  ////////////////////////////////////////////////////////////
  //            Loop over filter modules                    //
  ////////////////////////////////////////////////////////////
  for(unsigned int n=0; n < numOfHLTCollectionLabels ; n++) {
    // These numbers are from the Parameter Set, such as:
    //   theHLTOutputTypes = cms.uint32(100)
    switch(theHLTOutputTypes[n]) 
    {
      case 82: // Non-isolated Level 1
        fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,genParticlesPdgIdMatch);break;
      case 83: // Isolated Level 1
        fillHistos<l1extra::L1EmParticleCollection>(triggerObj,event,n,genParticlesPdgIdMatch);break;
      case 91: // Photon 
        fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,genParticlesPdgIdMatch);break;
      case 92: // Electron 
        fillHistos<reco::ElectronCollection>(triggerObj,event,n,genParticlesPdgIdMatch);break;
      case 100: // TriggerCluster
        fillHistos<reco::RecoEcalCandidateCollection>(triggerObj,event,n,genParticlesPdgIdMatch);break;
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
template <class T> void EmDQM::fillHistos(edm::Handle<trigger::TriggerEventWithRefs>& triggerObj,const edm::Event& iEvent ,unsigned int n,std::vector<HepMC::GenParticle>& genParticlesPdgIdMatch)
{

  std::vector<edm::Ref<T> > recoecalcands;

  // Only Process if available
  if ( triggerObj->filterIndex(theHLTCollectionLabels[n]) >= triggerObj->size() ) { 
    return;
  }
  
  ////////////////////////////////////////////////////////////
  //      Retrieve saved filter objects                     //
  ////////////////////////////////////////////////////////////
  triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),theHLTOutputTypes[n],recoecalcands);
  //Danger: special case, L1 non-isolated
  // needs to be merged with L1 iso
  if (theHLTOutputTypes[n] == 82) {
    std::vector<edm::Ref<T> > isocands;
    triggerObj->getObjects(triggerObj->filterIndex(theHLTCollectionLabels[n]),83,isocands);
    if (isocands.size()>0) {
      for (unsigned int i=0; i < isocands.size(); i++)
        recoecalcands.push_back(isocands[i]);
    }
  } // END of if theHLTOutputTypes == 82


  if (recoecalcands.size() < 1) {
    return;
  }

  ////////////////////////////////////////////////////////////
  // Loop over the Gen Particles collection, this should    //
  // contain highest and then 2nd highest Et gen particle.  //
  ////////////////////////////////////////////////////////////
  int highestEtMatchIndex  = -1;
  int highestEtMatchIndex2 = -1;
  int currentGenParticleCounter = 0;
  for (std::vector<HepMC::GenParticle>::iterator currentGenParticle = genParticlesPdgIdMatch.begin(); currentGenParticle !=  genParticlesPdgIdMatch.end() ; currentGenParticle++) {

    // Find direction of this gen particle
    math::XYZVector currentGenParticleDir( currentGenParticle->momentum().px(),
			   currentGenParticle->momentum().py(),
			   currentGenParticle->momentum().pz());

    // Loop over all HLT objects in this collection
    // to find the one closest in delta-R to currentGenParticle
    double closestDr = 0.5 ;
    int closestEcalCandIndex = -1;
    for (unsigned int i=0; i<recoecalcands.size(); i++) {

      // Find direction of this HLT object
      math::XYZVector candDir=recoecalcands[i]->momentum();

      // Calculated delta-r between gen particle and HLT object
      double dr = DeltaR(currentGenParticleDir,candDir);

      if (dr < closestDr) {
	closestEcalCandIndex = i;
	closestDr = dr;
      }
    } // END loop over HLT objects

    // If we found a delta-R matching HLT object
    if ( closestEcalCandIndex >= 0 ) {
      // First time through the loop, we
      // are looking to match highest et particle.
      if (currentGenParticleCounter == 0) {
        // Fill "Highest Et" histogram
        histHighestEt->Fill(n+0.5);
	// Fill histograms of the HLT Object
	histEtOfHighEtMatch[n]->Fill(  recoecalcands[closestEcalCandIndex]->et() );
	histEtaOfHighEtMatch[n]->Fill( recoecalcands[closestEcalCandIndex]->eta() );

	highestEtMatchIndex=closestEcalCandIndex;
      } else if (currentGenParticleCounter == 1 ) {
        // Fill the "2nd Highest Et" histogram
	histHighestEt2->Fill(n+0.5);
	histEtOfHighEtMatch2[n]->Fill(  recoecalcands[closestEcalCandIndex]->et() );
        histEtaOfHighEtMatch2[n]->Fill( recoecalcands[closestEcalCandIndex]->eta() );

	highestEtMatchIndex2 = closestEcalCandIndex;
      }
    }

    currentGenParticleCounter++;
    if (currentGenParticleCounter > 1) {
      break;
    } 
  } // END loop over gen particles



  ////////////////////////////////////////////////////////////
  //  Plot isolation variables (show the not-yet-cut        //
  //  isolation, i.e. associated to next filter)            //
  ////////////////////////////////////////////////////////////
  bool histIsoVsEtOfHighEtMatchFilled = false;
  bool histIsoVsEtOfHighEtMatchFilled2 = false;
  if (n+1 < numOfHLTCollectionLabels) { // can't plot beyond last
   if (plotiso[n+1] ){
    for (unsigned int j =  0 ; j < isoNames[n+1].size() ;j++  ){
     edm::Handle<edm::AssociationMap<edm::OneToValue< T , float > > > depMap; 
     iEvent.getByLabel(isoNames[n+1].at(j).label(),depMap);
     if (depMap.isValid()){ //Map may not exist if only one candidate passes a double filter

      if (highestEtMatchIndex >= 0) {
       typename edm::AssociationMap<edm::OneToValue< T , float > >::const_iterator mapi = depMap->find(recoecalcands[highestEtMatchIndex]);
       if (mapi!=depMap->end()) {  // found candidate in isolation map!                            
	if (highestEtMatchIndex >= 0 && !histIsoVsEtOfHighEtMatchFilled) {
	  histIsoVsEtOfHighEtMatch[n+1]->Fill(recoecalcands[highestEtMatchIndex]->et(),mapi->val);
	  histIsoVsEtaOfHighEtMatch[n+1]->Fill(recoecalcands[highestEtMatchIndex]->eta(),mapi->val);
	  histIsoVsEtOfHighEtMatchFilled=true;
	}
       }
      }
      if (highestEtMatchIndex2 >= 0 ) {
       typename edm::AssociationMap<edm::OneToValue< T , float > >::const_iterator mapi2 = depMap->find(recoecalcands[highestEtMatchIndex2]);
       if (mapi2!=depMap->end()) {  // found candidate in isolation map!
        if (highestEtMatchIndex2 >= 0 && !histIsoVsEtOfHighEtMatchFilled2) {
         histIsoVsEtOfHighEtMatch2[n+1]->Fill(recoecalcands[highestEtMatchIndex2]->et(),mapi2->val);
         histIsoVsEtaOfHighEtMatch2[n+1]->Fill(recoecalcands[highestEtMatchIndex2]->eta(),mapi2->val);
	 histIsoVsEtOfHighEtMatchFilled2=true;
	}
       }
      }
     }
    }
   }
  } // END of if n+1 < then the number of hlt collections
}


//////////////////////////////////////////////////////////////////////////////// 
//      method called once each job just after ending the event loop          //
//////////////////////////////////////////////////////////////////////////////// 
void EmDQM::endJob(){
}

DEFINE_FWK_MODULE(EmDQM);
