////////////////////////////////////////////////////////////////////////
// Class:       PDVDCheckMichel

// Plugin Type: analyzer (art v3_05_01)
// File:        PDVDCheckMichel_module.cc
//
// Generate simple summary tree from tracks and hits to check sim and reco
// Follow the same logic as in CosmicsdQdx
//
// Generated at 2/06/2023 by Thibaut Houdy
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <utility>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "art_root_io/TFileService.h"

// LArSoft includes
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

#include "TPolyMarker3D.h"
#include "TPolyMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLorentzVector.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

// ROOT
#include "TTree.h"
#include "TLegend.h"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pdvdana {
  class PDVDCheckMichel;
}

class pdvdana::PDVDCheckMichel : public art::EDAnalyzer {
public:
  explicit PDVDCheckMichel(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDVDCheckMichel(PDVDCheckMichel const&) = delete;
  PDVDCheckMichel(PDVDCheckMichel&&) = delete;
  PDVDCheckMichel& operator=(PDVDCheckMichel const&) = delete;
  PDVDCheckMichel& operator=(PDVDCheckMichel&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  int      fFlagBasics = 1;
  int      fFlagInfos = 2;
  int      fFlagWarning = 3;
  int      fFlagDetails = 4;
// int      fFlagDebug = 5;

  string   fHitModuleLabel = "gaushit"; // or "hitpdune"
  string   fTrackModuleLabel = "pandoraTrack"; 
  string   fSimuModuleLabel = "IonAndScint::G4Stage2"; //or "largeant:LArG4DetectorServicevolTPCActive"
  string   fClusterModuleLabel = "pandora"; // or "Reco3D"
  bool     fIsRecoBool = false;
  string   fGeneratorTag = "generator";
  float    fTrackMinLen;
  float    fMichelRadiusSphere;
  int      ev_num;

  int michel_counter=0, muon_michel_counter=0;

  unsigned fTotalTracks;
  unsigned fSelectedTracks;

  // summary tree
  TTree *ftrackTree;
  TTree *fhitTree;

  // Michel electron trees
  TTree *fmichelTree;
  TTree *fmuonTree;
  TTree *fdepoTree;
  
  protoana::ProtoDUNETruthUtils        truthUtil;
  protoana::ProtoDUNETrackUtils        trackUtil;
  protoana::ProtoDUNEPFParticleUtils    pfpUtil;

  //
  unsigned fEventNum;
  unsigned fTrackId;
  unsigned fNtpcs, fNplanes;
  unsigned ftrackStartTick;
  unsigned ftrackEndTick;

  unsigned fMichTrackId;
  unsigned fMichNum;
  float fMothertrackLen;
  float fMichX;
  float fMichY;
  float fMichZ;

  vector<float> fMichDepX;
  vector<float> fMichDepY;
  vector<float> fMichDepZ;
  vector<float> fMichDepEne;

  int fMichelPDG;
  float fMichelEne;
  bool fIsMichelInside;
  
  unsigned fMuonTrackId;
  unsigned fMuonNum;
  float fMuonTrackLen;
  float fMuontrackDx;
  float fMuontrackDy;
  float fMuontrackDz;
  float fMuonTrackTheta;
  float fMuonTrackPhi;
  float fMuontrackStartX;
  float fMuontrackStartY;
  float fMuontrackStartZ;
  bool fIsMuonInside;
  float fMuonVisi;
  float fMuonEnergy;

  unsigned fDepoTrackId;
  unsigned fDepoNum;
  float fDepoMothertrackLen;
  float fDepoX;
  float fDepoY;
  float fDepoZ;
  int fDepoPDG;
  float fDepoEne;
  bool fIsDepoInside;


  float ftrackLen;
  float ftrackDx;
  float ftrackDy;
  float ftrackDz;
  float ftrackStartX;
  float ftrackStartY;
  float ftrackStartZ;
  float ftrackEndX;
  float ftrackEndY;
  float ftrackEndZ;

  float ftracktheta;
  float ftrackphi;
  float ftracknorm;

  float fhitX;
  float fhitY;
  float fhitZ;
  float fhitcharge;
  float fhitenergy;
  // unsigned fhitPlane, fhitTPC, fhitWire, fhitChannel;
  //int fhitTime;

  vector<int> fOffsetWireID_u;
  vector<int> fOffsetWireID_v;
  vector<int> fOffsetWireID_z;
  
  float fgeoXmin = 1e6; 
  float fgeoXmax =-1e6; 
  float fgeoYmin = 1e6; 
  float fgeoYmax =-1e6; 
  float fgeoZmin = 1e6; 
  float fgeoZmax =-1e6; 
  bool ftrackCheatPDG ,ftrackCheatIsMichMuon,ftrackIsInside, ftrackCheatIsMichMuon_Al;

  // detector geometry
  const geo::Geometry* fGeom;

  void checkTrackCharacs(TVector3 l_track, vector<float> &track_char );
  float GetVisiTrack(TVector3 PosInit, TVector3 PosFinal);
  bool IsThisPointInsideTPC(TVector3 PosPoint);
  bool IsThisTrackEnteringTPC(TVector3 PosPointA, TVector3 PosPointB);
  bool IsMichelMichel(simb::MCParticle l_particle, simb::MCParticle l_mother_particle);
  bool IsMuonMichel(simb::MCParticle l_particle);

  void Drawing3D_HitsAndTracks(TCanvas* l_c, vector<vector<double>> l_t, vector<vector<double>> l_v, vector<vector<double>> l_d);
  void DrawCube(TPolyLine3D *Cube, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max );

  string FromPDGToName(int pdgcode);

};


pdvdana::PDVDCheckMichel::PDVDCheckMichel(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fHitModuleLabel( p.get< std::string >("HitModuleLabel") ),
  fTrackModuleLabel( p.get< std::string >("TrackModuleLabel") ),
  fSimuModuleLabel( p.get< std::string >("SimuModuleLabel") ),
  fClusterModuleLabel( p.get< std::string >("ClusterModuleLabel") ),
  fIsRecoBool(p.get<bool>("IsRecoBool") ),
  fGeneratorTag( p.get<std::string>("GeneratorTag") ),
  fTrackMinLen( p.get< float  >("TrackMinLen") ),
  fMichelRadiusSphere(p.get< float  >("MichelRadiusSphere") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
  }

//
void pdvdana::PDVDCheckMichel::analyze(art::Event const& ev)
{

  ev_num++;
  const string myname = "pdvdana::PDVDCheckMichel::analyze: ";
  std::cout << "Entering into the event: " << ev_num << std::endl;

  art::InputTag sim_tag(fSimuModuleLabel); //"IonAndScint::G4Stage2" or "largeant:LArG4DetectorServicevolTPCActive"
  
  art::Handle<std::vector<sim::SimEnergyDeposit>> simEnergyHandle;
  ev.getByLabel(fSimuModuleLabel, simEnergyHandle);

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  vector<vector<double>> TracksSpacePoints;
  vector<vector<double>> MuTrackStartEnd;
  vector<vector<double>> MuMichSpacePoints;

  vector<vector<double>> mich_candidate_interest;  

  vector<int> MichelElectronParticleTrackID;
  vector<int> MuonParticleTrackID;

  int event_muon_counter, event_michel_electron_counter;
  event_michel_electron_counter=0;
  event_muon_counter=0;
  if( !fIsRecoBool){
    //Looking at the generated particles with simb::MCTruth

    const sim::ParticleList & plist = pi_serv->ParticleList();
    if (fLogLevel>=fFlagInfos){
      std::cout << "Particle list size: " << plist.size() << std::endl;
    }
    TVector3 pos_init, pos_final;

    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
        simb::MCParticle *particle = ipar->second;
        if( fabs(particle->PdgCode()) == 13){
          event_muon_counter++;
        }

        double KE = (particle->E() - particle->Mass())*1e3 ; //Kinetic energy in MeV
        if( KE < 0.1 ) //100 keV threshold
            continue;           

        TString ParticleCreationProcess = (TString) particle->Process();
        TString ParticleDeadProcess = (TString) particle->EndProcess();
        pos_init = particle->Position().Vect();
        pos_final = particle->EndPosition().Vect();
    
        //Looking for the mother particle. The mother method returnn the TrackID of the mother particle.
        // If primary, the mother methode should be return 0.
        if ( fabs(particle->PdgCode()) == 13 and // muon or antimuon
            (ParticleDeadProcess.Contains("Decay")))                       //decaying 
            { 
              vector<float> muon_charac(6);
              checkTrackCharacs((pos_final-pos_init), muon_charac);
              fMuonVisi = GetVisiTrack(pos_init, pos_final);
              fMuonTrackId = particle->TrackId();
              fMuontrackStartX = pos_init.X();
              fMuontrackStartY = pos_init.Y();
              fMuontrackStartZ = pos_init.Z();
              fMuontrackDx = muon_charac[0];
              fMuontrackDy = muon_charac[1];
              fMuontrackDz = muon_charac[2];
              fMuonTrackLen = muon_charac[3];
              fMuonTrackTheta = muon_charac[4];
              fMuonTrackPhi = muon_charac[5];
              fMuonNum = ev_num;
              fIsMuonInside = IsThisTrackEnteringTPC(pos_init, pos_final);
              fMuonEnergy = KE;
              fmuonTree->Fill();

              MuonParticleTrackID.push_back(fMuonTrackId);
              if (fLogLevel>=fFlagDetails){
                std::cout << "!MUON! This particle is " << particle->PdgCode() <<"  "<<fMuonTrackId<<std::endl;
               // std::cout<<"From :(" <<particle->Position().X()<<", "<<particle->Position().Y()<<",  "<<particle->Position().Z()
               // <<") to : ("<<particle->EndPosition().X()<<", "<<particle->EndPosition().Y()<<",  "<<particle->EndPosition().Z()<<")"<<endl;
              }

              if( IsThisTrackEnteringTPC(pos_init, pos_final) and IsThisPointInsideTPC(pos_final))
              { //muon crossing the volume and dying within the TPC (for Michel lepton)
                muon_michel_counter++;
                vector<double> pos_part = {pos_init.X(), pos_init.Y(), pos_init.Z(), pos_final.X(), pos_final.Y(), pos_final.Z()};
                MuTrackStartEnd.push_back(pos_part);
            }
          }
          else{  

            if((fabs(particle->PdgCode())==11) and particle->Mother()>0)
              {          

                const simb::MCParticle *l_MotherParticle = pi_serv->TrackIdToParticle_P(particle->Mother());  
                if(! (ParticleCreationProcess.Contains("Decay") and fabs(l_MotherParticle->PdgCode())==13))
                  continue;
                //Here we must have electron/positron coming from a decay
                event_michel_electron_counter++;
                michel_counter++;

              
                //std::cout<<"We found a Michel! "<<michel_counter<<std::endl;
                fMichTrackId = particle->TrackId();
                fMothertrackLen = ((l_MotherParticle->Position()).Vect()-(l_MotherParticle->EndPosition()).Vect()).Mag();
                fMichX = particle->Position().X();
                fMichY = particle->Position().Y();
                fMichZ = particle->Position().Z();
                fMichNum = ev_num;

                
                //Here we wonder how many are coming from a muon/antimuon
                fMichelPDG = particle->PdgCode();
                fMichelEne = KE;
                fIsMichelInside = IsThisPointInsideTPC(particle->Position().Vect());
                fmichelTree->Fill();
                MichelElectronParticleTrackID.push_back( fMichTrackId);

                if (fLogLevel>=fFlagDetails)
                  std::cout << "!MICHEL! This particle is " << fMichelPDG<<"  "<<fMichTrackId<<std::endl;
                
                double totaldep=0;
                std::vector<const sim::IDE *> ide_P = bt_serv->TrackIdToSimIDEs_Ps(fMichTrackId, geo::View_t(2));
               
                // if(fLogLevel>=fFlagInfos)
                  // std::cout<<"Looking at IDEs: " << ide_P.size()<<std::endl;

                for (size_t i = 0; i < ide_P.size(); ++i) {
                  auto currentIde = ide_P[i];
                  if(fLogLevel>=fFlagInfos)
                    std::cout<<i<<" "<<currentIde->x<<"  "<<currentIde->energy <<"  "<<std::endl;
                  fMichDepX.push_back(currentIde->x);
                  fMichDepY.push_back(currentIde->y);
                  fMichDepZ.push_back(currentIde->z);
                  fMichDepEne.push_back(currentIde->energy);
                  totaldep += currentIde->energy;
                }
                // if(totaldep!=KE && fLogLevel>=fFlagInfos)
                  // std::cout<<"WARNING : "<<totaldep<<" different than "<<KE<<std::endl;


                vector<double> pos_part = {pos_init.X(), pos_init.Y(), pos_init.Z(), pos_final.X(), pos_final.Y(), pos_final.Z()};

            }
          }
    }

    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag);
    
    if(!depolist.isValid()) 
    { 
      std::string sWarningEdep = "\nUnable to find std::vector<sim::SimEnergyDeposit> with module label: " + fSimuModuleLabel;
    }
    if(depolist.isValid())
    {
      size_t ndepos = depolist->size();
      if (fLogLevel>=fFlagDetails){
        cout << "  # deposition: " << ndepos << endl;
      }
      // loop over list of depos in this event
      for (size_t depo_i = 0; depo_i < ndepos; depo_i++) {
        const sim::SimEnergyDeposit& depo = depolist->at(depo_i);
        // Retrieve position of energy deposit
        if (depo.Energy()<0.001)
          continue;
        geo::Point_t xyz = depo.MidPoint();
        fDepoX = xyz.X();
        fDepoY = xyz.Y();
        fDepoZ = xyz.Z();
        fDepoEne = depo.Energy();
        fDepoNum = ev_num;
        fDepoPDG = depo.PdgCode();
        fDepoTrackId = depo.TrackID();
        fDepoMothertrackLen = depo.OrigTrackID();
        fdepoTree->Fill();

        vector<double> pos_depo = {xyz.X(),xyz.Y(),xyz.Z()};
        if( abs(fDepoPDG)==13)
          MuMichSpacePoints.push_back(pos_depo);

        TracksSpacePoints.push_back(pos_depo);

      }
    }

  
  }

//----------------------------------------------------------------------------------------------------------
//------------------------------CHECKING THE RECONSTRUCTION---------------------------------------------------
  if(fIsRecoBool){


    // get tracks
    string showertag = "";
    art::InputTag shower_tagt(showertag);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(ev);


    art::InputTag hit_tag(fHitModuleLabel); //"gaushit" or "hitpdune"
    art::InputTag trck_tag(fTrackModuleLabel); //"pandoraTrack"
    art::InputTag pandora_tag("pandora");

    string calotag = "pandoraGnocalo";//"calo";
    art::InputTag calo_tag(calotag);

    // art::InputTag cluster_tag(fClusterModuleLabel); //"pandora" or "reco3d"

    auto Tracks = ev.getValidHandle<vector<recob::Track>>(trck_tag);

    fTotalTracks += Tracks->size();
    fEventNum = ev_num;

    vector<float> track_charac(6);



    auto const hitHandle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
    auto PFParticles = ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");
    std::vector<art::Ptr<recob::Hit>> all_hits;
    art::fill_ptr_vector(all_hits, hitHandle);
    
    if (fLogLevel>=fFlagInfos){
      cout<<myname<<all_hits.size()<<" hits found"<<endl;
      cout<<myname<<Tracks->size()<<" tracks found"<<endl;
      cout<<myname<<PFParticles->size()<<" PFParticules found"<<endl;
    }
    std::cout << "All primary pfParticles = " <<  pfpUtil.GetNumberPrimaryPFParticle(ev,"pandora") << std::endl;

    // Get the hits associated with the space points
    const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, ev, pandora_tag);

    if (!fmsph.isValid()) {
      throw cet::exception("LArPandoraShowerCheatingAlg")
        << "Spacepoint and hit association not valid. Stopping.";
    }


    // for(unsigned int p = 0; p < recoParticles->size(); ++p){
    //   const recob::PFParticle* f_pfparticle = &(recoParticles->at(p));

    //   //  Only the primary particles have the slice association
    //   std::cout<<"PFParticle : "<<f_pfparticle->PdgCode()<<"  "<< f_pfparticle->IsPrimary()<<"   "<< f_pfparticle->NumDaughters()<< std::endl;
      
    //   const recob::Track* thisTrack   = pfpUtil.GetPFParticleTrack(*f_pfparticle, ev, "pandora", fTrackModuleLabel);
    //   if(thisTrack != 0x0){
    //     // Get the true mc particle
    //     const simb::MCParticle* mcparticle = truthUtil.GetMCParticleFromRecoTrack(clockData, *thisTrack, ev, fTrackModuleLabel);
    //     if(mcparticle!=0x0){
    //       int ftruthpdg=mcparticle->PdgCode();
    //       int ftruthid=mcparticle->TrackId();
    //       if (fLogLevel>=fFlagDetails){
    //         std::cout<<"Associated MCParticle: "<<ftruthid<<" --- PDG: "<<ftruthpdg<<std::endl;
    //       }
    //     }
    //   }
    // }

    std::cout<<"Entering tracks loop"<<endl;
    //----------------------------------------------------------------------------------------------
    //Loops on tracks
    for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {

      int daughterid_being_michel=-999;

      const recob::Track& track = Tracks->at(itrk);
      ftrackLen    = track.Length();     

      if (ftrackLen<fTrackMinLen)
        continue;

      const simb::MCParticle *MCParticle_from_track = truthUtil.GetMCParticleFromRecoTrack(clockData, track, ev, fTrackModuleLabel);
      if(!MCParticle_from_track) 
        continue;

      bool HasElectron = false, HasNuMu = false, HasNuE = false;

      if(abs(MCParticle_from_track->PdgCode())==13)
        event_muon_counter++;

      if((abs(MCParticle_from_track->PdgCode())==13) && (MCParticle_from_track->NumberDaughters())>0){       
        for (int ii=0; ii<(MCParticle_from_track->NumberDaughters());++ii){
          const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P((MCParticle_from_track->Daughter(ii)));
          if(!daughter1) continue;
          int d_pdg1 = abs(daughter1->PdgCode());
          if (d_pdg1 == 11){ 
            HasElectron = true;
            daughterid_being_michel=ii;
          }
          else if (d_pdg1 == 14) HasNuMu = true;
          else if (d_pdg1 == 12) HasNuE = true;
        }
      }

      ftrackStartX  = track.Start().X();
      ftrackStartY  = track.Start().Y();
      ftrackStartZ  = track.Start().Z();
      ftrackStartTick  = 0;    
      ftrackEndX    = track.End().X();
      ftrackEndY    = track.End().Y();
      ftrackEndZ    = track.End().Z();

      TVector3 pos_final(ftrackEndX, ftrackEndY, ftrackEndZ);
      TVector3 pos_init(ftrackStartX, ftrackStartY, ftrackStartZ);

      ftrackEndTick = 0;
      ftrackCheatPDG = MCParticle_from_track->PdgCode();
      ftrackCheatIsMichMuon_Al = (HasElectron and HasNuMu) and HasNuE;
      ftrackCheatIsMichMuon = IsMuonMichel(*MCParticle_from_track);
      ftrackIsInside = IsThisTrackEnteringTPC(pos_init, pos_final) and (IsThisPointInsideTPC(pos_final) or IsThisPointInsideTPC(pos_final));

      // Filling the track tree
      checkTrackCharacs( pos_final-pos_init, track_charac );
      ftrackDx = track_charac[0];
      ftrackDy = track_charac[1];
      ftrackDz = track_charac[2];
      ftracknorm = track_charac[3];
      ftracktheta = track_charac[4];
      ftrackphi = track_charac[5];
      fTrackId     = 0;
      ftrackTree->Fill();
      if(ftrackCheatIsMichMuon && ftrackIsInside){
        vector<double> pos_track = {ftrackStartX, ftrackStartY, ftrackStartZ, ftrackEndX, ftrackEndY, ftrackEndZ};
        //vector<double> pos_particule = {MCParticle_from_track->Position().X(), MCParticle_from_track->Position().Y(), MCParticle_from_track->Position().Z()        , MCParticle_from_track->EndPosition().X(), MCParticle_from_track->EndPosition().Y(), MCParticle_from_track->EndPosition().Z()};
        MuTrackStartEnd.push_back(pos_track);

        if(daughterid_being_michel<0)
          continue;

        
        //Filling the Michel electron tree
        const simb::MCParticle* daughter1 = pi_serv->TrackIdToParticle_P((MCParticle_from_track->Daughter(daughterid_being_michel)));
        double KE = (daughter1->E() - daughter1->Mass())*1e3 ; //Kinetic energy in MeV
        fMichTrackId = daughter1->TrackId();
        fMothertrackLen = ((MCParticle_from_track->Position()).Vect()-(MCParticle_from_track->EndPosition()).Vect()).Mag();
        fMichX = daughter1->Position().X();
        fMichY = daughter1->Position().Y();
        fMichZ = daughter1->Position().Z();
        fMichNum = ev_num;
        fMichelPDG = daughter1->PdgCode();
        fMichelEne = KE;
        fmichelTree->Fill();
        michel_counter++;
        event_michel_electron_counter++;
        muon_michel_counter++;
        if(IsThisPointInsideTPC(pos_init))
          mich_candidate_interest.push_back({ftrackStartX, ftrackStartY, ftrackStartZ, KE});
        if(IsThisPointInsideTPC(pos_final)) 
          mich_candidate_interest.push_back({ftrackEndX, ftrackEndY, ftrackEndZ, KE});
      }

      std::vector<anab::Calorimetry> calovector = trackUtil.GetRecoTrackCalorimetry(track, ev, fTrackModuleLabel, calotag);
      for (auto & calo : calovector){

        if (calo.PlaneID().Plane != 2) //only collection plane
          continue;

        for (size_t ihit = 0; ihit < calo.dQdx().size(); ++ihit){ //loop over hits
          // calo.ResidualRange()[ihit];
          // calo.dEdx()[ihit];
          // calo.TrkPitchVec()[ihit];
          //calo.dQdx()[ihit];
          const auto &primtrk_pos=(calo.XYZ())[ihit];
          float local_hit_x = primtrk_pos.X();
          float local_hit_y = primtrk_pos.Y();
          float local_hit_z = primtrk_pos.Z();        
          TracksSpacePoints.push_back({local_hit_x, local_hit_y, local_hit_z});          
        } //End of hits loop from 1 calo
        
      }//End of calo loop coming from 1 track

    }// End of Tracks loop

    for(auto local_candidate:mich_candidate_interest){      
      unsigned fNumHits=0;

      TVector3 local_center(local_candidate[0],local_candidate[1], local_candidate[2]);
      std::cout<<"Local Center : "<<local_center(0)<<", "<<local_center(1)<<", "<<local_center(2)<<", "<<std::endl;

      for (auto hit : all_hits) {
        if(!(hit->WireID().Plane ==2))
          continue;

        // Get hit position
         std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());

        if(sps.size() < 1)
           continue; 

        art::Ptr<recob::SpacePoint> sp = sps.front();
        TVector3 local_hit_pos(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2]);
         double distance = (local_center-local_hit_pos).Mag();

        //Attention this is assuming we have a correct X (drift direction) reconstruction. This is NOT the case so far        
        if(distance < fMichelRadiusSphere){
          std::cout<<"Within the Sphere"<<endl;
          fhitX+=local_hit_pos[0];
          fhitY+=local_hit_pos[1];
          fhitZ+=local_hit_pos[2];
          fhitcharge+=hit->SummedADC();
          fNumHits++;
          MuMichSpacePoints.push_back({fhitX, fhitY, fhitZ});
        }


      }

      if(fNumHits>0){
          fhitX=fhitX/fNumHits;
          fhitY=fhitY/fNumHits;
          fhitZ=fhitZ/fNumHits;
          fhitenergy=local_candidate[3];
          fhitTree->Fill();

      }



    }



  }// End of Reconstruction only part



  art::ServiceHandle<art::TFileService> tfs;
  
  gStyle->SetOptStat(0);
  //Drawing 3D tracks
  TString canvasName_3D = Form("canvas3D_%i", ev_num);
  TCanvas* canvas_3D = tfs->make<TCanvas>(canvasName_3D, canvasName_3D);
  Drawing3D_HitsAndTracks(canvas_3D, TracksSpacePoints , MuMichSpacePoints, MuTrackStartEnd);
  canvas_3D->Write(canvasName_3D);

  if (fLogLevel>=fFlagInfos){
    std::cout << myname<<" REPORT : in this event :   # muons : " << event_muon_counter <<"  # michel electrons :" 
    << event_michel_electron_counter << " for a total of: "<< michel_counter <<" michel e"<<std::endl;
  }


} 
//
void pdvdana::PDVDCheckMichel::beginJob()
{

  ev_num = 0;
  fTotalTracks    = 0;
  fSelectedTracks = 0;
  fNplanes = fGeom->Nplanes();
  fNtpcs = fGeom->NTPC();

  fOffsetWireID_u.clear();
  fOffsetWireID_v.clear();
  fOffsetWireID_z.clear();

  unsigned int mem_nwires_u = 0;
  unsigned int mem_nwires_v = 0;
  unsigned int mem_nwires_z = 0;
  
  if(fLogLevel>=fFlagBasics){
    cout<<"Geometry : "<<endl;
    cout<<"number of planes :  "<<fNplanes<<endl;
    cout<<"number of tpc :  "<<fNtpcs<<endl;
  }
  
  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){

    geo::TPCID tpcid{0, t_tpc_id};
    geo::PlaneID const uplane_id{tpcid, geo::View_t::kU};
    geo::PlaneID const vplane_id{tpcid, geo::View_t::kV};
    geo::PlaneID const zplane_id{tpcid, geo::View_t::kZ};

    unsigned int nwires_u = fGeom->Nwires(uplane_id);
    unsigned int nwires_v = fGeom->Nwires(vplane_id);
    unsigned int nwires_z = fGeom->Nwires(zplane_id);

    fOffsetWireID_u.push_back(mem_nwires_u);
    fOffsetWireID_v.push_back(mem_nwires_v);
    fOffsetWireID_z.push_back(mem_nwires_z);

    mem_nwires_u += nwires_u;
    mem_nwires_v += nwires_v;
    mem_nwires_z += nwires_z;

    if(fLogLevel>=fFlagInfos){
      std::cout << "  TPC " << t_tpc_id << " center: ("<< fGeom->TPC(tpcid).GetCenter().X()      << "," << fGeom->TPC(tpcid).GetCenter().Y()      << ","<< fGeom->TPC(tpcid).GetCenter().Z() << ")"
                                 << " box:  ["  << fGeom->TPC(tpcid).BoundingBox().MinX() << "," << fGeom->TPC(tpcid).BoundingBox().MaxX() << "]" 
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinY() << "," << fGeom->TPC(tpcid).BoundingBox().MaxY() << "]"
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinZ() << "," << fGeom->TPC(tpcid).BoundingBox().MaxZ() << "]" << std::endl;
      cout<< "Number of wires : "<<nwires_u<<"  "<<nwires_v <<"  "<<nwires_z<<endl;
    }

    if(fgeoXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fgeoXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
    if(fgeoXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fgeoXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
    if(fgeoYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fgeoYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
    if(fgeoYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fgeoYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
    if(fgeoZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fgeoZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
    if(fgeoZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fgeoZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();

  }


  // init summary tree
  art::ServiceHandle<art::TFileService> tfs;

  fmichelTree = tfs->make<TTree>("michel","Michel candidates");
  fmichelTree->Branch("EventNum", &fMichNum);
  fmichelTree->Branch("TrackId", &fMichTrackId);
  fmichelTree->Branch("MotherTrackLen",   &fMothertrackLen);
  fmichelTree->Branch("PosX",  &fMichX);
  fmichelTree->Branch("PosY",  &fMichY);
  fmichelTree->Branch("PosZ",  &fMichZ);
  fmichelTree->Branch("PDG",  &fMichelPDG);
  fmichelTree->Branch("Energy", &fMichelEne);

  fmichelTree->Branch("DepEnergy",&fMichDepEne);
  fmichelTree->Branch("DepPosX",  &fMichDepX);
  fmichelTree->Branch("DepPosY",  &fMichDepY);
  fmichelTree->Branch("DepPosZ",  &fMichDepZ);
  fmichelTree->Branch("IsInside", &fIsMichelInside);

  fdepoTree = tfs->make<TTree>("depo","Depo candidates");
  fdepoTree->Branch("EventNum", &fDepoNum);
  fdepoTree->Branch("TrackId", &fDepoTrackId);
  fdepoTree->Branch("MotherTrackLen",   &fDepoMothertrackLen);
  fdepoTree->Branch("PosX",  &fDepoX);
  fdepoTree->Branch("PosY",  &fDepoY);
  fdepoTree->Branch("PosZ",  &fDepoZ);
  fdepoTree->Branch("PDG",  &fDepoPDG);
  fdepoTree->Branch("Energy", &fDepoEne);
  fdepoTree->Branch("IsInside", &fIsDepoInside);

  fmuonTree = tfs->make<TTree>("muons","muons candidates");
  fmuonTree->Branch("EventNum", &fMuonNum);
  fmuonTree->Branch("TrackId", &fMuonTrackId);
  fmuonTree->Branch("TrackLen",   &fMuonTrackLen);
  fmuonTree->Branch("Dx",  &fMuontrackDx);
  fmuonTree->Branch("Dy",  &fMuontrackDy);
  fmuonTree->Branch("Dz",  &fMuontrackDz);
  fmuonTree->Branch("StartX",  &fMuontrackStartX);
  fmuonTree->Branch("StartY",  &fMuontrackStartY);
  fmuonTree->Branch("StartZ",  &fMuontrackStartZ);
  fmuonTree->Branch("IsMuonI",  &fIsMuonInside);
  fmuonTree->Branch("MuonInitialEnergy",  &fMuonEnergy);
  fmuonTree->Branch("MuonVisi", &fMuonVisi);
  



  ftrackTree = tfs->make<TTree>("tracks","Check tracks");
  ftrackTree->Branch("EventNum", &fEventNum, "EventNum/i");
  ftrackTree->Branch("TrackId", &fTrackId, "TrackId/i");
  ftrackTree->Branch("TrackLen",   &ftrackLen,   "TrackLen/F");
  ftrackTree->Branch("Dx",  &ftrackDx,  "Dx/F");
  ftrackTree->Branch("Dy",  &ftrackDy,  "Dy/F");
  ftrackTree->Branch("Dz",  &ftrackDz,  "Dz/F");
  ftrackTree->Branch("StartX",  &ftrackStartX,  "StartX/F");
  ftrackTree->Branch("StartY",  &ftrackStartY,  "StartY/F");
  ftrackTree->Branch("StartZ",  &ftrackStartZ,  "StartZ/F");
  ftrackTree->Branch("StartTick",  &ftrackStartTick,  "StartTick/i");
  
  ftrackTree->Branch("EndX",  &ftrackEndX,  "EndX/F");
  ftrackTree->Branch("EndY",  &ftrackEndY,  "EndY/F");
  ftrackTree->Branch("EndZ",  &ftrackEndZ,  "EndZ/F");
  ftrackTree->Branch("EndTick",  &ftrackEndTick,  "EndTick/i");

  ftrackTree->Branch("theta",&ftracktheta,  "theta/F");
  ftrackTree->Branch("phi",  &ftrackphi,  "phi/F");
  ftrackTree->Branch("norm", &ftracknorm,  "norm/F");
  ftrackTree->Branch("trackCheatPDG", &ftrackCheatPDG);
  ftrackTree->Branch("trackCheatIsMichMuon_Al", &ftrackCheatIsMichMuon_Al);
  ftrackTree->Branch("trackCheatIsMichMuon", &ftrackCheatIsMichMuon);
  ftrackTree->Branch("trackIsInside", &ftrackIsInside);


  fhitTree  = tfs->make<TTree>("hits","Check reconstruction");
  fhitTree->Branch("X", &fhitX, "X/F");
  fhitTree->Branch("Y", &fhitY, "Y/F");
  fhitTree->Branch("Z", &fhitZ, "Z/F");
  fhitTree->Branch("charge", &fhitcharge);
  fhitTree->Branch("energy", &fhitenergy);
  fhitTree->Branch("EventNum", &fEventNum, "EventNum/i");
  fhitTree->Branch("TrackId",  &fTrackId, "TrackId/i");

  // fhitTree->Branch("Plane",  &fhitPlane, "Plane/i");
  // fhitTree->Branch("TPC", &fhitTPC, "TPC/i");
  // fhitTree->Branch("Wire",  &fhitWire, "Wire/i");
  // fhitTree->Branch("Channel", &fhitChannel, "Channel/i");
  // fhitTree->Branch("Time", &fhitTime, "Time/i");

}

//
void pdvdana::PDVDCheckMichel::endJob()
{

  fMichDepX.clear();
  fMichDepY.clear();
  fMichDepZ.clear();
  fMichDepEne.clear();

  const string myname = "pdvdana::PDVDCheckMichel::endJob: ";
  if(fLogLevel>=fFlagBasics){
    cout<<myname<<"tracks processed total     : "<<fTotalTracks<<endl;
    cout<<myname<<"tracks processed selected  : "<<fSelectedTracks<<endl;
  }
}
void pdvdana::PDVDCheckMichel::DrawCube(TPolyLine3D *Cube, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max ){
  Cube->SetPoint(0,x_min, y_min, z_min);
  Cube->SetPoint(1,x_max, y_min, z_min);
  Cube->SetPoint(2,x_max, y_max, z_min);
  Cube->SetPoint(3,x_min, y_max, z_min);
  Cube->SetPoint(4,x_min, y_min, z_min);
  Cube->SetPoint(5,x_min, y_min, z_max);
  Cube->SetPoint(6,x_max, y_min, z_max);
  Cube->SetPoint(7,x_max, y_min, z_min);
  Cube->SetPoint(8,x_max, y_max, z_min);
  Cube->SetPoint(9,x_max, y_max, z_max);
  Cube->SetPoint(10,x_min, y_max, z_max);
  Cube->SetPoint(11,x_min, y_min, z_max);
  Cube->SetPoint(12,x_min, y_max, z_max);    
  Cube->SetPoint(13,x_min, y_max, z_min);
  Cube->SetPoint(14,x_min, y_max, z_max);
  Cube->SetPoint(15,x_max, y_max, z_max);
  Cube->SetPoint(16,x_max, y_min, z_max);
  Cube->SetLineColor(12);
  Cube->SetLineWidth(3);
  return;
}

void pdvdana::PDVDCheckMichel::Drawing3D_HitsAndTracks(TCanvas* canvas_3D, vector<vector<double>> l_muMichSpacePoints, vector<vector<double>> l_muSpacePoints, vector<vector<double>> l_trackStartEnd){
  //Drawing 3D tracks

  TH3F *axes = new TH3F("axes", "3D Hits and tracks distribution; X; Y; Z", 1, -400, 400, 1, -400, 400, 1, -10, 400);
  axes->SetDirectory(0);
  canvas_3D->cd();
  axes->Draw();

  TPolyMarker3D* HitsSpacePoints_from_normal_track = new TPolyMarker3D(l_muSpacePoints.size());
  int actual_point=0;
  for (auto const& sp : l_muMichSpacePoints) {
    HitsSpacePoints_from_normal_track->SetPoint(actual_point, sp[0], sp[1], sp[2]);
    actual_point++;
  }
  HitsSpacePoints_from_normal_track->SetMarkerStyle(8);
  HitsSpacePoints_from_normal_track->SetMarkerSize(1);    
  HitsSpacePoints_from_normal_track->SetMarkerColor(kBlue);

  TPolyMarker3D* HitsSpacePoints_from_mumich_track = new TPolyMarker3D(l_muMichSpacePoints.size());
  actual_point=0;
  for (auto const& sp : l_muMichSpacePoints) {
    HitsSpacePoints_from_mumich_track->SetPoint(actual_point, sp[0], sp[1], sp[2]);
    actual_point++;
  }
  HitsSpacePoints_from_mumich_track->SetMarkerStyle(8);
  HitsSpacePoints_from_mumich_track->SetMarkerSize(1);    
  HitsSpacePoints_from_mumich_track->SetMarkerColor(kGreen);

  TPolyLine3D* Detector3D = new TPolyLine3D(17);
  DrawCube(Detector3D, fgeoXmin, fgeoYmin, fgeoZmin, fgeoXmax, fgeoYmax, fgeoZmax);

  canvas_3D->cd();
  Detector3D->Draw();
  HitsSpacePoints_from_mumich_track->Draw();  
   HitsSpacePoints_from_normal_track->Draw();  
  
  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){
    geo::TPCID tpcid{{0}, t_tpc_id};
    TPolyLine3D* TPC_3D = new TPolyLine3D(17);
    DrawCube(TPC_3D, fGeom->TPC(tpcid).BoundingBox().MinX(), fGeom->TPC(tpcid).BoundingBox().MinY(), fGeom->TPC(tpcid).BoundingBox().MinZ(), 
            fGeom->TPC(tpcid).BoundingBox().MaxX(), fGeom->TPC(tpcid).BoundingBox().MaxY(), fGeom->TPC(tpcid).BoundingBox().MaxZ());
    TPC_3D->SetLineWidth(1);
    TPC_3D->SetLineColor(29);
    canvas_3D->cd();
    TPC_3D->Draw();
  }

  for (auto const& sp : l_trackStartEnd) {
    TPolyLine3D* t_line = new TPolyLine3D(2);
    t_line->SetPoint(0, sp[0], sp[1], sp[2]);
    t_line->SetPoint(1, sp[3], sp[4], sp[5]);
    t_line->SetLineColor(2);
    t_line->SetLineWidth(2);
    t_line->Draw();
  }
  return;
}

bool pdvdana::PDVDCheckMichel::IsMichelMichel(simb::MCParticle l_particle, simb::MCParticle l_mother_particle){

  TString ParticleCreationProcess = (TString) l_particle.Process();
  
  if(!((fabs(l_particle.PdgCode())==11) and (l_particle.Mother()>0)))
    return false;

  if(ParticleCreationProcess.Contains("Decay") and fabs(l_mother_particle.PdgCode())==13)
    return true;

  return false;

}
bool pdvdana::PDVDCheckMichel::IsMuonMichel(simb::MCParticle l_particle){

  bool t_inside = false;
  bool t_ismuon = false;
  TString ParticleDeadProcess = (TString) l_particle.EndProcess();
  TVector3 pos_init = l_particle.Position().Vect();
  TVector3 pos_final = l_particle.EndPosition().Vect();

  if(IsThisTrackEnteringTPC(pos_init, pos_final) and IsThisPointInsideTPC(pos_final))
    t_inside = true;
  if(fabs(l_particle.PdgCode()) == 13 and (ParticleDeadProcess.Contains("Decay")))
    t_ismuon=true;
  
  return (t_ismuon and t_inside);
}

string pdvdana::PDVDCheckMichel::FromPDGToName(int pdgcode){

  if( pdgcode == 11)
    return "electron";
  else if( pdgcode == -11)
    return "positron";
  else if( pdgcode == -13)
    return "antimuon";  
  else if( pdgcode == 13)
    return "muon";
  else if( pdgcode == 2112)
    return "neutron";
  else if( pdgcode == 2212)
    return "proton";
  else if( pdgcode == 22)
    return "gamma";
  return string("inconnu ID: "+std::to_string(pdgcode));  
}

float pdvdana::PDVDCheckMichel::GetVisiTrack( TVector3 PosInit, TVector3 PosFinal ){

  TVector3 VisiInit, VisiFinal;
  VisiInit.SetXYZ(fmin(fmax(PosInit.X(),fgeoXmin),fgeoXmax),fmin(fmax(PosInit.Y(),fgeoYmin),fgeoYmax),fmin(fmax(PosInit.Z(),fgeoZmin), fgeoZmax));
  VisiFinal.SetXYZ(fmin(fmax(PosFinal.X(),fgeoXmin),fgeoXmax),fmin(fmax(PosFinal.Y(),fgeoYmin),fgeoYmax),fmin(fmax(PosFinal.Z(),fgeoZmin), fgeoZmax));
  return (VisiFinal-VisiInit).Mag();

}


// check the characteristics of a track -> length in x, y, z, total length and theta/phi angles
void pdvdana::PDVDCheckMichel::checkTrackCharacs( TVector3 l_track, vector<float> &TrackCharac ){
  
  double norm = l_track.Mag();
  float Dx = l_track.X();
  float Dy = l_track.Y();
  float Dz = l_track.Z();

  if (l_track.X()>0){
    Dx = -Dx;
    Dz = -Dz;
  }

  TrackCharac[0] = Dx;
  TrackCharac[1] = Dy;
  TrackCharac[2] = Dz;
  TrackCharac[3] = norm;
  TrackCharac[4] = 0;
  TrackCharac[5] = 0;

  if (norm<1)
    return;

  // Determination of the angle
  double theta = 180 - acos(abs(Dx)/norm)*180/3.1415;
  double phi = 0;

  if(Dy >= 0){
    if (Dz>=0)
      phi = atan(abs(Dz)/abs(Dy))*180/3.1415 ;
    else{
      phi = 270+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  else{
    if (Dz<0){
      phi= atan(abs(Dz)/abs(Dy))*180/3.1415+180;
    }
    else{
      phi= 90+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  TrackCharac[4] = theta;
  TrackCharac[5] = phi;
  return;
}


bool pdvdana::PDVDCheckMichel::IsThisTrackEnteringTPC(TVector3 PointA, TVector3 PointB){

  //From A and B, check if a point in the segment AB is inside the TPC
  double alpha = (PointB.Y()*PointA.Z()-PointA.Y()*PointB.Z())/(PointA.X()*PointB.Y()-PointB.X()*PointA.Y());
  double beta = (PointA.Z()*PointB.X()-PointB.Z()*PointA.X())/(PointA.Y()*PointB.X()-PointB.Y()*PointA.X());
  float var_z=0;
  for (float var_x=fgeoXmin; var_x<fgeoXmax; var_x=var_x+1){
    for (float var_y=fgeoXmin; var_y<fgeoXmax; var_y=var_y+1){
      var_z=alpha*var_x+beta*var_y;
      if((var_z<fgeoZmin) && (var_z>fgeoZmin)){
        return false;
      }
    }
  }
  return true;
}

//From A check if the point is inside the TPC
bool pdvdana::PDVDCheckMichel::IsThisPointInsideTPC(TVector3 PointA){

  if(PointA.Y() < fgeoYmin || PointA.Y()>fgeoYmax)
    return false;
  if(PointA.Z() < fgeoZmin || PointA.Z()>fgeoZmax)
    return false;

  return true;
}

DEFINE_ART_MODULE(pdvdana::PDVDCheckMichel)
