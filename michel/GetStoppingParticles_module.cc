////////////////////////////////////////////////////////////////////////
// Class:       GetStoppingParticles 
// Plugin Type: analyzer (Unknown Unknown)
// File:        GetStoppingParticles_module.cc
//
// Generated at Mon Jan 30 17:27:10 2023 by Yoann Kermaidic using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include <limits>  // std::numeric_limits<>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// ROOT
#include <TStyle.h>
#include "TTree.h"
#include "TH2.h"


using std::vector;
using std::string;

namespace {
  constexpr unsigned int int_max_as_unsigned_int{std::numeric_limits<int>::max()};
}

namespace pddpana {
  class GetStoppingParticles;
}


class pddpana::GetStoppingParticles : public art::EDAnalyzer {
public:
  explicit GetStoppingParticles(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GetStoppingParticles(GetStoppingParticles const&) = delete;
  GetStoppingParticles(GetStoppingParticles&&) = delete;
  GetStoppingParticles& operator=(GetStoppingParticles const&) = delete;
  GetStoppingParticles& operator=(GetStoppingParticles&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  float    degTOrad = 3.14159/180.; // rad / deg

private:

  // Parameters set in the input FHICL file (verbosity, module name, thresholds, ...)
  int      fLogLevel;
  string   fTrackModuleLabel;
  string   fCaloModuleLabel;
  float    fFiducialCut;
  float    fTrackCurvature;
  float    fTrackDelta;
  float    fTrackDeltaStop;
  float    fTrackLenMin;
  float    fTrackLenMax;
  float    fTrackThetaMin;
  float    fTrackThetaMax;

  unsigned fTotalTracks;
  unsigned int fNplanes;

  // summary tree
  TTree   *fTree;

  // objects stored in the output ROOT file
  unsigned fEventNum;
  unsigned fTrackId;
  float    fLength;
  float    fTheta;
  float    fPhi;
  vector<vector<float> > fX;
  vector<vector<float> > fY;
  vector<vector<float> > fZ;
  vector<vector<float> > fResRange;
  vector<vector<float> > fdQdx;

  // cm - custom active volume boundaries
  float fAVx_min;
  float fAVx_max;
  float fAVy_min;
  float fAVy_max;
  float fAVz_min;
  float fAVz_max;
  float fAVl_min = 0;
  float fAVl_max = sqrt(pow(fAVx_max-fAVx_min,2)+pow(fAVy_max-fAVy_min,2)+pow(fAVz_max-fAVz_min,2));
 
  // cm - Reco uncertainty margin
  float eps   = 5;

  // cm - hit spacing threshold (holes in mis-reconstructed tracks)
  float hs_thrs = 50;

  // Store valid track start and end
  // also directly accessible through track.Start().X() / track.End().X() 
  struct Edge{
    float x;
    float y;
  };

  // detector geometry
  const geo::Geometry* fGeom;

  // Store dQ/dx vs residual range profile into TH2F
  vector<TH2F*> hRange;

  void  GetAngles(const recob::Track& track, float &theta, float &phi);
};


pddpana::GetStoppingParticles::GetStoppingParticles(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fTrackModuleLabel( p.get< std::string  >("TrackModuleLabel") ),
  fCaloModuleLabel(  p.get< std::string  >("CaloModuleLabel") ),
  fFiducialCut(   p.get< float  >("FiducialCut") ),
  fTrackCurvature(p.get< float  >("TrackCurvature") ), // cm - hit to straight track spacing threshold
  fTrackDelta(    p.get< float  >("TrackDelta") ),
  fTrackDeltaStop(p.get< float  >("TrackDeltaStop") ),
  fTrackLenMin(   p.get< float  >("TrackLenMin") ),
  fTrackLenMax(   p.get< float  >("TrackLenMax") ),
  fTrackThetaMin( p.get< float  >("TrackThetaMin") ),
  fTrackThetaMax( p.get< float  >("TrackThetaMax") ),
  fAVx_min( p.get< float  >("AVx_min") ),
  fAVx_max( p.get< float  >("AVx_max") ),
  fAVy_min( p.get< float  >("AVy_min") ),
  fAVy_max( p.get< float  >("AVy_max") ),
  fAVz_min( p.get< float  >("AVz_min") ),
  fAVz_max( p.get< float  >("AVz_max") )
   
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
}

void pddpana::GetStoppingParticles::analyze(art::Event const& e)
{
  if( fLogLevel >= 1 ) std::cout << "Start analysing file ..." << std::endl;

  // retrieve number of planes (2 for DP and 3 elsewhere)
  fNplanes = fGeom->Nplanes();
 
  bool point_set = false;
  Edge bot_edge;
  Edge top_edge;

  int hitout = 0;         // # of hits outside specs
  int isgood = 0;         // Flag bad = 0 or good = 1 tracks

  // get tracks
  auto Tracks   = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  fTotalTracks += Tracks->size();  

  // get associated calorimetry information
  art::FindManyP<anab::Calorimetry> fmcal(Tracks, e, fCaloModuleLabel);
  
  if( fLogLevel >= 1 )
    std::cout << "The event contains "<< Tracks->size() << " tracks\n";
  
  // retrieve input event number for output
  fEventNum = e.id().event();

  // loop over tracks
  int ntrks[3] = {0,0,0};   // Number of [good,bad,out] tracks taken into account in the analysis

  for (unsigned trk = 0; trk < Tracks->size(); ++trk) {
    const recob::Track& track = Tracks->at(trk);

    bot_edge.x = 0;
    bot_edge.y = 0;
    top_edge.x = 0;
    top_edge.y = 0;

    isgood = 1;
    hitout = 0;
    point_set = false;

    fLength = track.Length(0);
    fTheta = 0, fPhi = 0;
    GetAngles(track,fTheta,fPhi);

    if( fLogLevel >= 1 ) std::cout << "  #" << trk << " length: "<< fLength <<" cm / theta: " << fTheta << " deg \n";

    // Filter out tracks that are shorter/longer than specified threshold
    if(fLength < fTrackLenMin   || fLength > fTrackLenMax)  {ntrks[2]++; continue;}
    // Filter out tracks that are too horizontal or too vertical to speed-up processing 
    if(fTheta  < fTrackThetaMin || fTheta  > fTrackThetaMax){ntrks[2]++; continue;}      
  
    bool is_out = false;
/*
    // Filter out tracks that do not start at the anode (unknown T0)
    if(track.Start().Y() - fAVy_min < fFiducialCut || track.Start().Y() - fAVy_max > 100 - fFiducialCut) is_out = true;
    if(track.Start().Z() - fAVz_min < fFiducialCut || track.Start().Z() - fAVz_max > 100 - fFiducialCut) is_out = true;

    if(!is_out){
      // Filter out tracks that e (unknown T0)
      if(track.End().Y() - fAVy_min < fFiducialCut || track.End().Y() - fAVy_max > 100 - fFiducialCut) is_out = true;
      if(track.End().Z() - fAVz_min < fFiducialCut || track.End().Z() - fAVz_max > 100 - fFiducialCut) is_out = true;
    }*/

    
    if(is_out){ ntrks[2]++; std::cout << "OUT :  " << track.Start().Y() << "," << track.End().Y() << " - " << track.Start().Z() << "," << track.End().Z()  << std::endl;  continue;}

    // remove tracks that too short or too long with respect to the active volume
    if(fLength < fAVl_min-eps || fLength  > fAVl_max+eps){
      if(fLogLevel >= 1) std::cout << "  BAD: Length [" << fAVl_min-eps << "," << fAVl_max+eps << "] -> " <<  fLength << std::endl;
      isgood = 1;
    }

    // below we want to remove space points that are within boundaries 
    // mis-reconstructed or too large gap in between consecutive space points

    vector<bool> hitskip(track.NPoints(),false);
      
    for(size_t hit=0; hit<track.NPoints(); ++hit){
      if(!isgood) break;
      if(track.LocationAtPoint(hit).X() == -999) {hitskip[hit] = true; continue;}

      if(!point_set){
        bot_edge.x = track.LocationAtPoint(hit).X() - fAVx_min;
        bot_edge.y = sqrt(pow(track.LocationAtPoint(hit).Y() - fAVy_min,2) + pow(track.LocationAtPoint(hit).Z() - fAVz_min,2));
        point_set = true;
      }

      top_edge.x = track.LocationAtPoint(hit).X() - fAVx_min;
      top_edge.y = sqrt(pow(track.LocationAtPoint(hit).Y() - fAVy_min,2) + pow(track.LocationAtPoint(hit).Z() - fAVz_min,2));

      float hit_spacing = 0;
      if(hit>0) hit_spacing = sqrt(pow(track.LocationAtPoint(hit).X()-track.LocationAtPoint(hit-1).X(),2)
                                 + pow(track.LocationAtPoint(hit).Y()-track.LocationAtPoint(hit-1).Y(),2)
                                 + pow(track.LocationAtPoint(hit).Z()-track.LocationAtPoint(hit-1).Z(),2));

      if(track.LocationAtPoint(hit).X() < fAVx_min ||
         track.LocationAtPoint(hit).X() > fAVx_max ||
         track.LocationAtPoint(hit).Y() < fAVy_min-eps ||
         track.LocationAtPoint(hit).Y() > fAVy_max+eps ||
         track.LocationAtPoint(hit).Z() < fAVz_min-eps ||
         track.LocationAtPoint(hit).Z() > fAVz_max+eps) {hitskip[hit] = true; hitout++;}

      if(hitout > 5 || hit_spacing > hs_thrs) {
        std::cout << "  BAD: " << hit << " " << hitout << " "
                          << track.LocationAtPoint(hit).X() << " "
                          << track.LocationAtPoint(hit).Y() << " "
                          << track.LocationAtPoint(hit).Z() << std::endl;
        isgood = 1;
        break;
      }
    } // end of space point loop

    // Compute average track slope from the two track extremum points
    // used to compute the minimal distance between space point and straight line from extremums 

    if(top_edge.x-bot_edge.x == 0) isgood = 0;
    else if(isgood){
      float a = 666.; a = (top_edge.y-bot_edge.y)/(top_edge.x-bot_edge.x);
      float b = 666.; b = bot_edge.y - a*bot_edge.x;

      for(size_t hit=0; hit<track.NPoints(); ++hit){
        if(hitskip[hit]) continue;

        float xpt = track.LocationAtPoint(hit).X() - fAVx_min;
        float yzpt= sqrt(pow(track.LocationAtPoint(hit).Y() - fAVy_min,2) + pow(track.LocationAtPoint(hit).Z() - fAVz_min,2));

        float xx  = 1./(a*a+1)*(xpt+a*yzpt-a*b);
        float yy  = 1./a*(xpt-xx) + yzpt;

        if(sqrt(pow(xx-xpt,2) + pow(yy-yzpt,2)) > fTrackCurvature) {
          std::cout << "  BAD: " << hit << " distance: " << sqrt(pow(xx-xpt,2) + pow(yy-yzpt,2)) << " " << xpt << " " << yzpt << std::endl;
          isgood = 0;
          break;
        }
      }
    }

    if(fLogLevel >= 1)
        std::cout << "Event: #" << e.id().event() << ": track length: " << fLength << " cm / theta: " << fTheta << " deg found -> " << isgood << " status" << std::endl;

    ntrks[isgood]++;

    // if the particle is not stopping in the detector, go to the next one
    if(!isgood) continue;

    // Now the stopping track selection is done, we can look for the calorimetric information
    std::vector<art::Ptr<anab::Calorimetry>> calos = fmcal.at(trk);

    if(fLogLevel >= 1) std::cout << "        - calo size: " << calos.size() << std::endl;

    // loop over calos entry
    for(size_t ical = 0; ical<calos.size(); ++ical){
      if(!calos[ical])                    continue;
      if(!calos[ical]->PlaneID().isValid) continue;
      int view = calos[ical]->PlaneID().Plane;
      if(view<0 || view >= int(fNplanes))   continue;

      // retrieve the number of space points
      const size_t NHits = calos[ical] -> dEdx().size();

      if(fLogLevel >= 2) std::cout << "  Calo: " << ical << " " << NHits << std::endl;

      // initialize output vectors
      fX[view].resize(NHits,-999);
      fY[view].resize(NHits,-999);
      fZ[view].resize(NHits,-999);
      fResRange[view].resize(NHits,-999);
      fdQdx[view].resize(NHits,-999);

      // loop over track space points
      for(size_t iHit = 0; iHit < NHits; ++iHit){
        if((calos[ical] -> dQdx())[iHit] <= 0) continue;

        // retrieve the space points coordinates
        const auto& TrkPos = (calos[ical] -> XYZ())[iHit];

        fX[view][iHit]        = TrkPos.X();
        fY[view][iHit]        = TrkPos.Y();
        fZ[view][iHit]        = TrkPos.Z();
        fdQdx[view][iHit]     = (calos[ical] -> dQdx())[iHit];
        fResRange[view][iHit] = (calos[ical] -> ResidualRange())[iHit];

        // Fill TH2F dQ/dx profile
        hRange[view]->Fill((calos[ical] -> ResidualRange())[iHit],(calos[ical] -> dQdx())[iHit]);
      } // end of space point loop
    } // end of calo loop

    fTree->Fill();
  } // end of tracks loop

  if(fLogLevel >= 1){
    std::cout << "Total # of tracks: " << ntrks[0]+ntrks[1] << std::endl;
    std::cout << "  - bad:  " << ntrks[0] << std::endl;
    std::cout << "  - good: " << ntrks[1] << std::endl;
    std::cout << "  - out:  " << ntrks[2] << std::endl;
    std::cout << " " << std::endl;
  }
}

void pddpana::GetStoppingParticles::beginJob()
{
  if(fLogLevel >= 1){
    std::cout << " #planes:       " << fGeom->Nplanes() << std::endl;
    std::cout << " TrackMinLen:   " << fTrackLenMin << std::endl;
    std::cout << " TrackMaxLen:   " << fTrackLenMax << std::endl;
    std::cout << " TrackThetaMin: " << fTrackThetaMin << std::endl;
    std::cout << " TrackThetaMax: " << fTrackThetaMax << std::endl;
  }

  hRange.resize(fGeom->Nplanes());  
  fResRange.resize(fGeom->Nplanes());
  fdQdx.resize(fGeom->Nplanes());
  fX.resize(fGeom->Nplanes());
  fY.resize(fGeom->Nplanes());
  fZ.resize(fGeom->Nplanes());

  art::ServiceHandle<art::TFileService> tfs;
  for (unsigned plane_i = 0; plane_i < fGeom->Nplanes(); plane_i++){
    hRange[plane_i]   = tfs->make<TH2F>(Form("hRange_%d",plane_i),Form("Plane %d dE/dx vs Residual range",plane_i),600,0,600,2000,0,20000);
  }

  fTree = tfs->make<TTree>("fTree","Store dE/dx profile");
  fTree->Branch("EventNum",  &fEventNum,"EventNum/i"  );
  fTree->Branch("TrackId",   &fTrackId, "TrackId/i"   );
  fTree->Branch("TrackLen",  &fLength,  "TrackLen/F"  );
  fTree->Branch("TrackTheta",&fTheta,   "TrackTheta/F");
  fTree->Branch("TrackPhi",  &fPhi,     "TrackPhi/F"  );
  fTree->Branch("TrackX",    &fX);
  fTree->Branch("TrackY",    &fY);
  fTree->Branch("TrackZ",    &fZ);
  fTree->Branch("TrackRange",&fResRange);
  fTree->Branch("TrackDqdx", &fdQdx);
}

void pddpana::GetStoppingParticles::endJob()
{
  std::cout << "Done " << std::endl;
  std::cout << " " << std::endl;
}

void pddpana::GetStoppingParticles::GetAngles(const recob::Track& track, float &theta, float &phi){

  float x_len = track.End().X()-track.Start().X();
  float y_len = track.End().Y()-track.Start().Y();
  float z_len = track.End().Z()-track.Start().Z();

  phi = 0;
  if(     y_len >= 0 && z_len > 0) phi =                atan(y_len/z_len) / degTOrad;
  else if(y_len >= 0 && z_len < 0) phi = TMath::Pi()   -atan(y_len/z_len) / degTOrad;
  else if(y_len < 0  && z_len < 0) phi = TMath::Pi()   +atan(y_len/z_len) / degTOrad;
  else if(y_len < 0  && z_len > 0) phi = TMath::TwoPi()-atan(y_len/z_len) / degTOrad;

  theta = atan((x_len)/sqrt(pow(y_len,2)+pow(z_len,2))) / degTOrad;
}

DEFINE_ART_MODULE(pddpana::GetStoppingParticles)
