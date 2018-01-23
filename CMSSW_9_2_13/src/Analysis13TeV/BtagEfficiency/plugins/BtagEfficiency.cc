// -*- C++ -*-
//
// Package:    Analysis13TeV/BtagEfficiency
// Class:      BtagEfficiency
// 
/**\class BtagEfficiency BtagEfficiency.cc Analysis13TeV/BtagEfficiency/plugins/BtagEfficiency.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Muhammad Shoaib
//         Created:  Tue, 21 Nov 2017 17:37:45 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/Framework/interface/Run.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"


//#include "Analysis13TeV/BtagEfficiency/interface/BtagEfficiency.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat; 


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class BtagEfficiency : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit BtagEfficiency(const edm::ParameterSet&);
      ~BtagEfficiency();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
	

	
      // ----------member data ---------------------------

	edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackTags_;

	edm::EDGetTokenT<reco::VertexCollection> PvtxToken_;
	edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;

	edm::EDGetTokenT<edm::View<pat::PackedCandidate> > CandidateToken_;

	edm::EDGetTokenT<pat::MuonCollection> muonToken_;
	edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_; //MiniAod PAt Electrons
	edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
	std::unordered_map<std::string,TH1*> histContainer_;


	edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
	edm::EDGetTokenT<GenEventInfoProduct> generatorevtToken_;
	edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;

	edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> pseudoTopToken_;
	edm::EDGetTokenT<std::vector<reco::GenJet>  > genLeptonsToken_,   genJetsToken_;


	edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
	edm::EDGetTokenT<double> rhoToken_;
	edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

	std::vector<std::string> bDiscriminators_; // string for b-jet discriminators 
	
	// declare a map of b-tag discriminator histograms
	
	std::map<std::string, TH1F *> bDiscriminatorsMap1;
	std::map<std::string, TH2F *> bDiscriminatorsMap2;
	


	edm::Service<TFileService> fs;
	TH1D *Tracks_Pt;
	TH1D *Tracks_dz;
	TH1D *Tracks_dxy;
	TH1D *Tracks_mass;
	TH1D *Tracks_energy;
	TH1D *histo_dxyError;

	TH1D *histo_Hits;
	TH1D *histo_PixelHits;

	TTree *ttree_;


	double bDiscriminator = 0. ; 
	double NoPVertices    = 0. ;
	double NoSecVertices    = 0. ;
	

	int    pvNTracks   ;
	
	Int_t nPUTrue_;    // true pile-up
	Int_t nPU_;        // generated pile-up
	Int_t nPV_;        // number of reconsrtucted primary vertices

	Float_t rho_;      // the rho variable

	Float_t genWeight_;
	//Float_t rho_;      // the rho variable


	// values of vriables copied from Deep_Csv_Code
	
//	float min3DIPValue=0.005;
//	float min3DIPSignificance=1.2;
//	int max3DIPValue=9999.;
//	int max3DIPSignificance=9999.;


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
BtagEfficiency::BtagEfficiency(const edm::ParameterSet& iConfig):
PvtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("generator")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

//	nVtx(0), vertex_id(0),
//        priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
//        priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),


// Universal tokens for AOD and miniAOD
	CandidateToken_	= consumes<edm::View<pat::PackedCandidate> > 
	(iConfig.getParameter <edm::InputTag>
	("trackLabel"));
	


	generatorevtToken_ = consumes<GenEventInfoProduct>
	
	(iConfig.getParameter <edm::InputTag>
	("genEventInfoProduct"));

	svToken_ = consumes<reco::VertexCompositePtrCandidateCollection>
	(iConfig.getParameter<edm::InputTag>
	("secVertices"));

	pileupToken_ = consumes<edm::View<PileupSummaryInfo> >
	(iConfig.getParameter <edm::InputTag>
	("pileup"));

	rhoToken_    = consumes<double>
	(iConfig.getParameter <edm::InputTag>
	("rho"));

	beamSpotToken_    = consumes<reco::BeamSpot>
	(iConfig.getParameter <edm::InputTag>
	("beamSpot"));

	electronToken_     = mayConsume<edm::View<pat::Electron> >
	(iConfig.getParameter<edm::InputTag>
	("electrons"));

	trackTags_ 	= consumes<edm::View<pat::PackedCandidate>>
	(iConfig.getParameter<edm::InputTag>
	("tracks"));

	bDiscriminators_   = (iConfig.getParameter<std::vector<std::string> >("bDiscriminators"));
	
	  std::string bDiscr_flav = "";


	   // initialize b-tag discriminator histograms
	   for( const std::string &bDiscr : bDiscriminators_ )
		{

	if( bDiscr.find("Counting") != std::string::npos ) // track counting discriminator can be both positive and negative and covers a wider range then other discriminators
       bDiscriminatorsMap1[bDiscr] = fs->make<TH1F>(bDiscr.c_str(), (bDiscr + ";b-tag discriminator").c_str(), 1100, -15, 40);
     else
       bDiscriminatorsMap1[bDiscr] = fs->make<TH1F>(bDiscr.c_str(), (bDiscr + ";b-tag discriminator").c_str(), 440, -11, 11);

	// for-loop over flavor
     for( const std::string &flav : {"b","c","udsg"} )
     {
       bDiscr_flav = bDiscr + "_" + flav;
       if( bDiscr.find("Counting") != std::string::npos ) // track counting discriminator can be both positive and negative and covers a wider range then other discriminators
         bDiscriminatorsMap2[bDiscr_flav] = fs->make<TH2F>(bDiscr_flav.c_str(), (bDiscr_flav + ";Jet p_{T} [GeV];b-tag discriminator").c_str(), 20, 0, 200, 11000, -15, 40);
       else
         bDiscriminatorsMap2[bDiscr_flav] = fs->make<TH2F>(bDiscr_flav.c_str(), (bDiscr_flav + ";Jet p_{T} [GeV];b-tag discriminator").c_str(), 20, 0, 200, 4400, -11, 11);
     }
		}// end of for-loop

	
	Tracks_Pt = fs->make<TH1D>("Tracks pt" , "Pt" , 50 , 0 , 50 );
	Tracks_mass = fs->make<TH1D>("Tracks mass" , "mass" , 50 , 0 , 50 );
	Tracks_energy = fs->make<TH1D>("Tracks energy" , "energy" , 50 , 0 , 50 );
	Tracks_dz = fs->make<TH1D>("Tracks dz" , "dz" , 100 , -50 , 50 );
	Tracks_dxy = fs->make<TH1D>("Tracks dxy" , "dxy" , 50 , 0 , 10 );
	histo_dxyError = fs->make<TH1D>("dxyError" , "dxyError" , 50 , 0 , 1 );

	histo_Hits = fs->make<TH1D>("Hits" , "Hits" , 30 , 0 , 30 );
	histo_PixelHits = fs->make<TH1D>("PixelHits" , "PixelHits" , 10 , 0 , 10 );

	ttree_ = fs->make<TTree>("tree", "tree");
	ttree_->Branch("nPV", &pvNTracks, "nPV/I");


}



BtagEfficiency::~BtagEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BtagEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


using std::vector;


// good info check this: https://github.com/jmejiagu/MiniAODexamples/blob/master/DemoAnalyzer/plugins/miniAODmuons.cc


	// Set event info
	cout<<"Run #:  "<<iEvent.id().run()<<endl;
	cout<<"Event #:  "<<iEvent.id().event()<<endl;
	cout<<"timestamp: "<<iEvent.time().unixTime()<<endl;



// Primary vertices

	reco::Vertex bestVtx;

	 edm::Handle<reco::VertexCollection> Pvertices;
        iEvent.getByToken(PvtxToken_, Pvertices);
        if (Pvertices->empty()) return ; // skip the event if no PV found
	
	bestVtx = *(Pvertices->begin());	

        const reco::Vertex &primVtx = Pvertices->front();
        reco::VertexRef primVtxRef(Pvertices,0);


        GlobalPoint pvp(primVtx.x(),primVtx.y(),primVtx.z());

	// Check there is at least one good primary vertex
	
	 bool foundGoodVertex = false;
        for(unsigned int i = 0; i< Pvertices->size() ; i++){
	
	//double Hello_Test = (*Pvertices)[i].dxz();
	//cout<<"Hello_Test:   "<<Hello_Test<<endl;

    if(!(*Pvertices)[i].isFake() && (*Pvertices)[i].tracksSize()>2
       && fabs( (*Pvertices)[i].z()) <= 24 && fabs( (*Pvertices)[i].position().rho()) <= 2) foundGoodVertex = true;
  }
        cout<<"foundGoodVertex?: "<<foundGoodVertex<<endl;

        NoPVertices = Pvertices->size();
        cout<<"No of Primary Vertices: "<<NoPVertices<<endl;

        pvNTracks = primVtx.nTracks();
        cout<<"pvNTracks....: "<<pvNTracks<<endl;


// Get gen weight info
	edm::Handle< GenEventInfoProduct > genWeightH;
	iEvent.getByToken(generatorevtToken_,genWeightH);
	genWeight_ = genWeightH->GenEventInfoProduct::weight();

	cout<<"genWeight_ ...:  "<< genWeight_ <<endl;



//-- Tracks Information -----



	// track builder
	 edm::ESHandle<TransientTrackBuilder> builder;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);



//	 double d03D = IPTools::absoluteImpactParameter3D(tk, pv).second;

	
	edm::Handle<edm::View<pat::PackedCandidate> > tracks; 
	iEvent.getByToken(CandidateToken_, tracks);


	//vector<TransientTrack> t_tks = (*builder).build(tracks);


	
	for(View<pat::PackedCandidate>::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack )
	     {

	if(itTrack->pt()<=0.5)continue; 
	if(itTrack->charge()==0) continue;// NO neutral objects
if(fabs(itTrack->pdgId())!=211) continue;//Due to the lack of the particle ID all the tracks for cms are pions(ID==211)
	if(!(itTrack->trackHighPurity())) continue;


/*
	TrackRef glbTrackP;	  
	TrackRef glbTrackM;	

	if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
*/

       Tracks_Pt ->Fill(itTrack->pt());
       Tracks_dz ->Fill(itTrack->dz());
	Tracks_mass ->Fill(itTrack->mass());
	Tracks_energy ->Fill(itTrack->energy());


	//*****Accessing track impact parameters **********
        Tracks_dxy ->Fill(itTrack->dxy());
	//histo_dxyError->Fill(itTrack->dxyError());


	cout<<"Pt: "<<itTrack->pt()<<"  Position:  "<<itTrack->dz()<<"dxy(): "<<itTrack->dxy()<<endl;
	cout<<"phi(): "<<itTrack->phi()<<" , mass():  "<<itTrack->mass()<<", energy(): "<<itTrack->energy()<<endl;
	cout<<"px(): "<<itTrack->px()<<" , p4():  "<<itTrack->p4()<<", polarP4():  " <<itTrack->polarP4()<<endl;
	
	cout<<"itTrack->charge():  " <<itTrack->charge()<<endl;


       //****** Accessing the hit pattern of a track  ********

	// count the number of valid tracker *** hits ***
	std::cout << "number of valid tracker hits is "  << itTrack->numberOfHits() << std::endl;
	

	// count the number of valid pixel *** hits ***
	std::cout << "number of valid pixel hits is " << itTrack->numberOfPixelHits() << std::endl;

	 histo_Hits->Fill(itTrack->numberOfHits());
       histo_PixelHits->Fill(itTrack->numberOfPixelHits());

	     }
	

/*
	std::vector<reco::TransientTrack> TTracks;	
	std::vector<float> masses;


	int Tracks_size = 0;
	Tracks_size = tracks->size();
	cout<<"Tracks_size: " <<Tracks_size <<endl;
	
	double TracksPt = 0.;

	for(size_t k = 0; k<tracks->size(); ++k) 
	{

	 if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5 && std::fabs(pvp.z()-builder->build(tracks->ptrAt(k)).track().vz())<0.5) {
            selectedTracks.push_back(builder->build(tracks->ptrAt(k)));
            masses.push_back(tracks->ptrAt(k)->mass());
        }

	TracksPt = (*tracks)[k].pt();
	if (TracksPt < 0.5) continue;

	cout<<"TracksPt:  "<<TracksPt <<endl;
	}//end for loop

*/

//for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
//std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*it, primVtx);
//        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*it, primVtx);
//        std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, primVtx);
//        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),primVtx, direction,it->field());
//}
//is the track in the jet cone?
//float angular_distance=std::sqrt(std::pow(jet.eta()-it->track().eta(),2) + std::pow(jet.phi()-it->track().phi(),2) );




// ----- Secondary Vertices ------------------

	edm::Handle<std::vector<reco::VertexCompositePtrCandidate> > secvertices;
	iEvent.getByToken(svToken_, secvertices);
	
	NoSecVertices = secvertices->size();	
	cout<<"Secondary Vertices:  "<<NoSecVertices<<endl;
// ------------- Pileup info -------------
 Handle<edm::View<PileupSummaryInfo> > pileupHandle;
  iEvent.getByToken(pileupToken_, pileupHandle);
  for( auto & puInfoElement : *pileupHandle){
	if( puInfoElement.getBunchCrossing() != 0 ) continue;
    if( puInfoElement.getBunchCrossing() == 0 ){
      nPU_    = puInfoElement.getPU_NumInteractions();
      nPUTrue_= puInfoElement.getTrueNumInteractions();
    }
  }

cout<<"nPU_ ...:  "<<nPU_<<"    nPUTrue_ ...   "<<nPUTrue_<<endl;

// Get rho value
	edm::Handle< double > rhoH;
	iEvent.getByToken(rhoToken_,rhoH);
	rho_ = *rhoH;
	cout<<"rho_ ...:   "<<rho_ <<endl;
 
// Get the beam spot

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(beamSpotToken_,theBeamSpot);
 
//Generator info
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);


 // ELECTRON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	edm::Handle<edm::View<pat::Electron> > electrons;
	iEvent.getByToken(electronToken_, electrons);     

	Int_t nele(0);
  for (const pat::Electron &el : *electrons) 
    {        
      const auto e = electrons->ptrAt(nele); 
	 if( el.pt() < 10 ) // keep only electrons above 10 GeV
	    continue;
      nele++;
	//kinematics cuts
	cout<<"no of Electrons:  "<<nele<<endl;
    }


//Gen studies ........................... 
 int ngjets(0),ngbjets(0),JetpdgId(-99);
 double GenJetspT(0.),GenJetsEta(0.), genJetPhi(0.), genJetMass(0.);

	for(auto genJet : *genJets)
    		{
		JetpdgId = genJet.pdgId();
		GenJetspT = genJet.pt();
		GenJetsEta = genJet.eta();
		genJetPhi = genJet.phi();
		genJetMass = genJet.mass();
cout<<"pdgId: "<<JetpdgId<<"pt: "<<GenJetspT<<"Eta: "<<GenJetsEta<<"phi: " <<genJetPhi<<"mass: "<<genJetMass<<endl;
		if(genJet.pt()> 25 && fabs(genJet.eta()) < 2.5)
			{
          		ngjets++;
          		if(abs(genJet.pdgId())==5) ngbjets++;

			} 

		}//for loop ends here


//JETS Study
edm::Handle<edm::View<pat::Jet>> MyJets;
iEvent.getByToken(jetToken_, MyJets);


std::string bDiscr_flav = "";

for(auto j = MyJets->begin();  j != MyJets->end(); ++j)
    {

//kinematic cuts --
if( j->pt()<30. || std::abs(j->eta())>2.4 ) continue; // skip jets with low pT or outside the tracker acceptance

//Jets flavor
int flavor = std::abs( j->hadronFlavour() );

//bdiscriminant

// fill discriminator histograms
for( const std::string &bDiscr : bDiscriminators_ )
{
	if( flavor==5 ) // b jet
        bDiscr_flav = bDiscr + "_b";
      else if( flavor==4 ) // c jets
        bDiscr_flav = bDiscr + "_c";
      else // light-flavor jet
        bDiscr_flav = bDiscr + "_udsg";

	bDiscriminatorsMap1[bDiscr]->Fill( j->bDiscriminator(bDiscr) );

	bDiscriminatorsMap2[bDiscr_flav]->Fill( j->pt(), j->bDiscriminator(bDiscr) );

}


bDiscriminator = j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

cout<<"bJet_discriminant_all_value: " <<bDiscriminator<<endl;


	}//jets loop ends here


ttree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
BtagEfficiency::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BtagEfficiency::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BtagEfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BtagEfficiency);
