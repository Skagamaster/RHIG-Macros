/// C++ headers
#include <iostream>

/// PicoDst headers
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"

/// ROOT headers
#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"
#include "TH1.h"
#include "TMath.h"


R__LOAD_LIBRARY(libStPicoDst)

//_________________
void SkipRunAnalysis00(const Char_t *inFile = 
  "st_physics_19140004_raw_5000011.picoDst.root") {

  gSystem->Load("libStPicoDst.so");

  std::cout << "Excelsior! Command your physics servant, Master!" << std::endl;

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  
  /// This accelerates IO
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;
  
  //Long64_t events2read = 10;
  Long64_t events2read = picoReader->chain()->GetEntriesFast();

  std::cout << "Number of events to read: " << events2read
      << std::endl;

  /// Histogramming
  TH1F *hRefMult = new TH1F("hRefMult","Reference multiplicity;refMult",
        500, -0.5, 499.5);
  TH1F *hTransvMomentum = new TH1F("hTransvMomentum",
        "Track transverse momentum;p_{T} (GeV/c)",200, 0., 2.);
  TH1F *hBTofTrayHit = new TH1F("hBTofTrayHit","BTof tray number with the hit",
        120, -0.5, 119.5);
  TH1F *hPhiAngles = new TH1F("hPhiAngles",
        "#phi_{T} for Each Track;#phi_{T} (rad)",200,-3.2,3.2);
  TH1F *hV2 = new TH1F("hV2","v_{2} per Event;v_{2}",200,-0.4,0.4);
  TH1F *hPsiRP = new TH1F("hPsiRP","Reaction Plane Angle;#Psi_{RP} (rad)",
        200,-1.0,1.0);
  TH2F *hV2Calc = new TH2F("hV2Calc",
        "Sums for v_{2} Calculation;cos(2(#phi-#Psi_{RP});/N",
        200,-100.0,100.0,200,500.0,2000.0);

  /// Here are a couple of global variables for 27GeV v2.
    Double_t GeVv2 = 0.0;
    Double_t GeVsumVar = 0.0;
    Double_t GeVnTracks = 0.0;
  
  /// Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    /// Let's initialise some variables.
    Double_t Qx = 0.0;
    Double_t Qy = 0.0;
    Double_t sumVar = 0.0;
    Double_t v2 = 0.0;
    Double_t PsiRP = 0.0;

    //std::cout << "Obey your orders, Master! Working on event #["
        //<< (iEvent+1) << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->ReadPicoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..."
    << std::endl;
      break;
    }

    /// Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    /// Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..."
    << std::endl;
      break;
    }
    hRefMult->Fill( event->refMult() );

    /// Track analysis
    Int_t nTracks = dst->numberOfTracks();
    //std::cout << "Number of tracks in event: " << nTracks << std::endl;
    /// Loop over tracks
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      
      StPicoTrack *picoTrack = dst->track(iTrk);
            //std::cout << "Transverse Momentum: " << picoTrack->gPt()
            //<< std::endl;
      if(!picoTrack) continue;
      //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  
      //<< std::endl;

      /// Single-track cut example
      if( !picoTrack->isPrimary() ||
    picoTrack->nHits() < 15 ||
    TMath::Abs( picoTrack->gMom().pseudoRapidity() ) > 0.5 ||
    // picoTrack->gMom().pseudoRapidity() > 0.5 ||
    0.2 > picoTrack->gPt() > 2.0 ||
    nTracks < 5) {
  continue;
      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
      
      hTransvMomentum->Fill( picoTrack->gPt() );
      hPhiAngles->Fill(picoTrack->gMom().phi() );
      
      /// These are our variables for the reaction plane angle.
      Qx += (picoTrack->gPt())*cos(2*picoTrack->gMom().phi());
      Qy += (picoTrack->gPt())*sin(2*picoTrack->gMom().phi());
    }
    /// This is to skip events where PsiRP would be undefined.
    if (Qx == 0.0) { continue; } 

    /// Hit analysis
    Int_t nBTofHits = dst->numberOfBTofHits();
    //std::cout << "Number of btofHits in event: " << nBTofHits << std::endl;
    for(Int_t iHit=0; iHit<nBTofHits; iHit++) {

      StPicoBTofHit *btofHit = dst->btofHit(iHit);
      if( !btofHit ) continue;
      hBTofTrayHit->Fill( btofHit->tray() );
    }

    /// Now we can calculate the reaction plane angle.
    if (Qx != 0.0)
    {
      PsiRP += atan(Qy/Qx)/2;
      hPsiRP->Fill(PsiRP);
    }

    /// Now we can calculate v2 for the event.
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      
      StPicoTrack *picoTrack = dst->track(iTrk);
      if(!picoTrack) continue;
      //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  
      //<< std::endl;

      /// Single-track cut example
      if( !picoTrack->isPrimary() ||
    picoTrack->nHits() < 15 ||
    TMath::Abs( picoTrack->gMom().pseudoRapidity() ) > 0.5 ||
    // picoTrack->gMom().pseudoRapidity() > 0.5 ||
    0.2 > picoTrack->gPt() > 2.0 ||
    nTracks < 5) {
  continue;
      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
      sumVar += cos(2*(picoTrack->gMom().phi()-PsiRP));

    }

    /// Calculating v2 for the event and saving stuff.
    v2 += sumVar/nTracks;
    hV2->Fill(v2);
    hV2Calc->Fill(sumVar,nTracks);
    GeVsumVar += sumVar;
    GeVnTracks += nTracks;
    if ( v2 == 0.0) std::cout << "           PsiRP  ="<<PsiRP << std::endl;
    {
      /* code */
    }

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  picoReader->Finish();

  GeVv2 += GeVsumVar/GeVnTracks;
  std::cout << "Average v2 for all events: " << GeVv2 << std::endl;

  std::cout << "Analysis complete; the Nobel shall be ours, Master!" 
  << std::endl;
  
  string Tacos = "00PROCESSED";
  Tacos.append(inFile);
  TString Tacos1 = Tacos;
  TString pathSave = "/home/zander/Documents/PicoHists/";
  TFile *MyFile = TFile::Open(pathSave+Tacos1,"RECREATE");

    hRefMult->Write();
    hTransvMomentum->Write();
    hBTofTrayHit->Write();
    hV2->Write();
    hV2Calc->Write();
    hPhiAngles->Write();
    hPsiRP->Write();
    MyFile->Close();

  string Nachos = "00PsiCALC";
  Nachos.append(inFile);
  TString Nachos1 = Nachos;
  TString pathSave1 = "/home/zander/Documents/PicoPsiRP/";
  TFile *MyFile1 = TFile::Open(pathSave1+Nachos1,"RECREATE");

    hPsiRP->Write();
    MyFile1->Close();
}
