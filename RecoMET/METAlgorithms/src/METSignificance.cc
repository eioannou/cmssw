
// -*- C++ -*-
//
// Package:    METAlgorithms
// Class:      METSignificance
// 
/**\class METSignificance METSignificance.cc RecoMET/METAlgorithms/src/METSignificance.cc
Description: [one line class summary]
Implementation:
[Notes on implementation]
*/
//
// Original Author:  Nathan Mirman (Cornell University)
//         Created:  Thu May 30 16:39:52 CDT 2013
//
//

#include "RecoMET/METAlgorithms/interface/METSignificance.h"


metsig::METSignificance::METSignificance(const edm::ParameterSet& iConfig) {

  edm::ParameterSet cfgParams = iConfig.getParameter<edm::ParameterSet>("parameters");

  double dRmatch = cfgParams.getParameter<double>("dRMatch");
  dR2match_ = dRmatch*dRmatch;
  
  jetThreshold_ = cfgParams.getParameter<double>("jetThreshold");
  jetEtas_ = cfgParams.getParameter<std::vector<double> >("jeta");
  jetParams_ = cfgParams.getParameter<std::vector<double> >("jpar");
  pjetParams_ = cfgParams.getParameter<std::vector<double> >("pjpar");
  
}

metsig::METSignificance::~METSignificance() {
}


reco::METCovMatrix
metsig::METSignificance::getCovariance(const edm::View<reco::Jet>& jets,
				       const std::vector< edm::Handle<reco::CandidateView> >& leptons,
				       const edm::Handle<edm::View<reco::Candidate> >& pfCandidatesH,
				       double rho,
				       JME::JetResolution& resPtObj,
				       JME::JetResolution& resPhiObj,
				       JME::JetResolutionScaleFactor& resSFObj,
				       bool isRealData) {

  //pfcandidates
  const edm::View<reco::Candidate>* pfCandidates=pfCandidatesH.product();

   // metsig covariance
   double cov_xx = 0;
   double cov_xy = 0;
   double cov_yy = 0;
  
   // for lepton and jet subtraction
   //std::set<reco::CandidatePtr> footprint;
   std::vector<reco::CandidatePtr> footprint;
   // subtract leptons out of sumPt
   for( const auto& lep_i : leptons ){
     for( const auto& lep : *lep_i ){
       if( lep.pt() > 10. ){
	 
	 for( unsigned int n=0; n < lep.numberOfSourceCandidatePtrs(); n++ ){
	   
	   if( lep.sourceCandidatePtr(n).isNonnull() and lep.sourceCandidatePtr(n).isAvailable() ){
	     footprint.push_back(lep.sourceCandidatePtr(n));
	   
	     //footprint.insert(lep->sourceCandidatePtr(n));
	   } 
	 }
       }
     }
   }

   
   
   // subtract jets out of sumPt
   //for(edm::View<reco::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
   for( const auto& jet : jets ) {
     // disambiguate jets and leptons
     //if(!cleanJet(*jet, leptons) ) continue;
     if( !cleanJet(jet, leptons) ) continue;
   

     for( unsigned int n = 0; n < jet.numberOfSourceCandidatePtrs(); n++){
       if( jet.sourceCandidatePtr(n).isNonnull() and jet.sourceCandidatePtr(n).isAvailable() ){
	 
	 footprint.push_back(jet.sourceCandidatePtr(n));
	 //footprint.insert(jet->sourceCandidatePtr(n));
       }
     }
   }

   // calculate sumPt
   double sumPt = 0;  
   std::cout << "Size of Candidates : " << pfCandidates->size() << std::endl;
   std::cout << "Size of footprint : " << footprint.size() << std::endl;
   for(size_t i = 0; i< pfCandidates->size();  ++i) {
     
     // check if candidate exists in a lepton or jet
     bool cleancand = true;
     if(std::find(footprint.begin(), footprint.end(), pfCandidates->ptrAt(i)) != footprint.end()){
       
       for( const auto& it : footprint){
	 if( (it->p4()-(*pfCandidates)[i].p4()).Et2()<0.000025){
	   cleancand = false;
	   break;
	 }
       }
     }
     // if not, add to sumPt
     if( cleancand ){
       sumPt += (*pfCandidates)[i].pt();
     }
   }
   
   
   // add jets to metsig covariance matrix and subtract them from sumPt
   for( const auto& jet : jets ){
     //for(edm::View<reco::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
     
     // disambiguate jets and leptons
     //if(!cleanJet(*jet, leptons) ) continue;
     if( !cleanJet(jet, leptons) ) continue;
     
     
     double jpt  = jet.pt();
     double jeta = jet.eta();
     double feta = std::abs(jeta);
     double c = jet.px()/jet.pt();
     double s = jet.py()/jet.pt();
     
     JME::JetParameters parameters;
     parameters.setJetPt(jpt).setJetEta(jeta).setRho(rho);
     
     // jet energy resolutions
     double sigmapt = resPtObj.getResolution(parameters);
     double sigmaphi = resPhiObj.getResolution(parameters);
     //double sigmapt = 1.;
     //double sigmaphi = 1.;
     // std::cout << " Sigmapt = " << sigmapt << " Sigmaphi = " << sigmaphi << std::endl;
     
     // SF not needed since is already embedded in the sigma in the dataGlobalTag
     //      double sigmaSF = isRealData ? resSFObj.getScaleFactor(parameters) : 1.0;
     
     // split into high-pt and low-pt sector
     if( jpt > jetThreshold_ ){
       // high-pt jets enter into the covariance matrix via JER
       
       double scale = 0;
       if(feta<jetEtas_[0]) scale = jetParams_[0];
       else if(feta<jetEtas_[1]) scale = jetParams_[1];
       else if(feta<jetEtas_[2]) scale = jetParams_[2];
       else if(feta<jetEtas_[3]) scale = jetParams_[3];
       else scale = jetParams_[4];
       
       //         double dpt = sigmaSF*scale*jpt*sigmapt;
       double dpt = scale*jpt*sigmapt;
       double dph = jpt*sigmaphi;
       
       cov_xx += dpt*dpt*c*c + dph*dph*s*s;
       cov_xy += (dpt*dpt-dph*dph)*c*s;
       cov_yy += dph*dph*c*c + dpt*dpt*s*s;
       
     } else {
       
       // add the (corrected) jet to the sumPt
       sumPt += jpt;
       
     }
     
   }
   
   
   std::cout << " SUMPT = " << sumPt << std::endl;
   //protection against unphysical events
   if(sumPt<0) sumPt=0;
   
   // add pseudo-jet to metsig covariance matrix
   cov_xx += pjetParams_[0]*pjetParams_[0] + pjetParams_[1]*pjetParams_[1]*sumPt;
   cov_yy += pjetParams_[0]*pjetParams_[0] + pjetParams_[1]*pjetParams_[1]*sumPt;

   reco::METCovMatrix cov;
   cov(0,0) = cov_xx;
   cov(1,0) = cov_xy;
   cov(0,1) = cov_xy;
   cov(1,1) = cov_yy;

   std::cout << " CovXX = " << cov_xx
             << " CovXY = " << cov_xy
             << " CovYY = " << cov_yy
	     << " SumPt = " << sumPt
             << std::endl;
   
   
   return cov;
}

double
metsig::METSignificance::getSignificance(const reco::METCovMatrix& cov, const reco::MET& met) {

   // covariance matrix determinant
   double det = cov(0,0)*cov(1,1) - cov(0,1)*cov(1,0);

   // invert matrix
   double ncov_xx = cov(1,1) / det;
   double ncov_xy = -cov(0,1) / det;
   double ncov_yy = cov(0,0) / det;

   // product of met and inverse of covariance
   double sig = met.px()*met.px()*ncov_xx + 2*met.px()*met.py()*ncov_xy + met.py()*met.py()*ncov_yy;

   return sig;
}


bool
metsig::METSignificance::cleanJet(const reco::Jet& jet,
				  const std::vector<edm::Handle<reco::CandidateView> >& leptons){
  
  for( const auto& lep_i : leptons ){
    for( const auto& lep : *lep_i ){
      if( reco::deltaR2(lep, jet) < dR2match_ ){
	return false;
      }
    }
  }
  return true;
}
/*bool
metsig::METSignificance::cleanJet(const reco::Jet& jet, 
      const std::vector< edm::Handle<reco::CandidateView> >& leptons ){

  for ( std::vector< edm::Handle<reco::CandidateView> >::const_iterator lep_i = leptons.begin();
        lep_i != leptons.end(); ++lep_i ) {
     for( reco::CandidateView::const_iterator lep = (*lep_i)->begin(); lep != (*lep_i)->end(); lep++ ){
        if ( reco::deltaR2(*lep, jet) < dR2match_ ) {
           return false;
        }
     }
  }
  return true;
}
*/
