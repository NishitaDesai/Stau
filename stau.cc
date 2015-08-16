// Modification of Py8 example main24
// For study of long-lived stau; pass to DELPHES via HepMC
// Nishita Desai

#include <Pythia8/Pythia.h>
#include <Pythia8/ToyDetector.h>

// #include "Pythia8/Pythia8toHepMC.h"
// #include "HepMC/GenEvent.h"
// #include "HepMC/IO_GenEvent.h"

#if USE_ROOT
// ROOT, for histogramming.
#include <root/TH1.h>
// ROOT, for saving file.
#include <root/TFile.h>
#endif

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // // Interface for conversion from Pythia8::Event to HepMC event. 
  // HepMC::Pythia8ToHepMC ToHepMC;

  // // Specify file where HepMC events will be stored.
  // HepMC::IO_GenEvent ascii_io("hepmcoutstau.dat", std::ios::out);


  // Read in commands from external file.
  pythia.readFile("stau.cmnd");    

  bool isBatch = false; // Turn off for histograms

  if(argc < 2){
    cout << "Usage: ./stau.exe SLHA_inputfilename <outputfilename>" << endl;
    return 1;
  }
  string SLHAfile(argv[1]);
  cout << "Reading from " << SLHAfile <<endl;
  pythia.readString("SLHA:file="+SLHAfile);

  string outputFile="";
  if(argc == 3){
    outputFile += string(argv[2]);
  }

  ofstream op, op2, op3;
  if(outputFile.length() > 0){
    string out1 = outputFile+"-stau.dat";
    string out2 = outputFile+"-xsec-stau.dat";
    if(!isBatch) {
      op.open(out1.c_str(),ios::out);
      op3.open((SLHAfile+".out").c_str(),ios::out);
    }
    op2.open(out2.c_str(),ios::app);
  }
  else{
    // Output file for plots, cross sections 
    if(!isBatch) {
      op.open("stau-output.dat",ios::out);
      op3.open("stau.out",ios::out);
    }
    op2.open("cross-sec.dat",ios::app);
  }

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  // double eCM   = pythia.parm("Beams:eCM");

  // Initialize.
  pythia.init();
  CoupSUSY* coupSUSY = (CoupSUSY*) pythia.couplingsPtr;
  double m12 = (coupSUSY->slhaPtr)->minpar(2);
  pythia.rndm.readState("rndm.dat");



  // Histograms.
  Hist pTj1("pT of 1st jet",100,0,200);
  Hist met("Missing ET", 100, 0, 500);
  Hist met2("Missing ET (Toy)", 100, 0, 500);
  Hist phijm("Delta phi(j,MET)", 100, 0, 6.28);  
  Hist pT_trk("pT of track",100,0,300);
  Hist r_trk("radial track length",100,0,10000.0);
  Hist isocone("isolation",100,0,0.20);


#if USE_ROOT
  // Root Histograms
  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("hist.root", "RECREATE");
  
  TH1F *pT_stau = new TH1F("pT_stau","dN/dpT", 100, 0., 500.);
  TH1F *eta_stau = new TH1F("eta_stau","dN/deta", 100, -2.5, 2.5);
  TH1F *beta_stau = new TH1F("beta_stau","dN/dbeta", 100, 0., 1.);
  TH1F *dl = new TH1F("DecayLength_stau","dN/dL (mm)", 100, 0., 10000.);
  TH1F *rdl = new TH1F("RadialDecayLength_stau","dN/dR (mm)", 100, 0., 10000.);
  
  TH1F *pT_stau_acc = new TH1F("pT_stau_acc","dN/dpT", 100, 0., 500.);
  TH1F *eta_stau_acc = new TH1F("eta_stau_acc","dN/deta", 100, -2.5, 2.5);
  TH1F *beta_stau_acc = new TH1F("beta_stau_acc","dN/dbeta", 100, 0., 1.);
  TH1F *dl_acc = new TH1F("DecayLength_stau_acc","dN/dL (mm)", 100, 0., 10000.);
  TH1F *rdl_acc = new TH1F("RadialDecayLength_stau_acc","dN/dR (mm)", 100, 0., 10000.);
#endif

  // Begin event loop.
  int iAbort = 0;
  double lifetime;
  int neffA = 0, neffB = 0, neffC = 0, nStau = 0, nEv = 0;
  int nDecy = 0, nStable = 0;
  int nEvStable = 0, nEvDiss = 0;
  int nEvStableCand = 0, nEvDissCand = 0;
  int nDiss75 = 0, nDiss100 = 0, nDiss150 = 0, nDiss200 = 0;

  int neffMET = 0, neffdPhi = 0;

  bool printFlag = false;

  ToyDetector detector;

  // For validation
  // int iCand = pythia.mode("SUSY:idA");

  int iCand = 1000015;
  double mSt = pythia.particleData.m0(iCand);
  double mNeut = pythia.particleData.m0(1000022); 

  vector <Vec4> jets, leptons, electrons, muons;
  Vec4 MET;
  double rdist[10000];
  int iLen = 0;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    if(!detector.getObjects(event)) cout << "No objects found" << endl;  
    else {
      jets = detector.jets;
      leptons = detector.leptons;
      electrons = detector.electrons;
      muons = detector.muons;
      MET = detector.MET;
      // detector.printObjects();
    }

    int stauArray[100];
    int stauArrLen = 0;

    // Lepton veto in the final state for disappearing tracks
    bool flagLep = false;

    for (int i = 0; i < event.size(); ++i) {
      if (abs(event[i].id()) == iCand) {
	if(event[i].status() == -22 || event[i].isFinal()) {
	  nStau++;
	  stauArray[stauArrLen] = i;
	  stauArrLen++;
 	}
      }
      else if (event[i].status() > 0 && (abs(event[i].id()) == 11 || abs(event[i].id()) == 13)) {
	if (event[i].pT() > 10.0 ) flagLep = true;
      }
    } // End stau search loop
      
    if(stauArrLen < 1 ) continue;

    nEv++;

    // Look at highest pT track within eta < 2.4
    int iStau = -1;
    double maxpT = 0.0;
    for(int i = 0; i < stauArrLen; ++i){
      if (abs(event[stauArray[i]].eta()) < 2.4) {
	if(event[stauArray[i]].status() < 0) nDecy++;
	else nStable++;
	if(event[stauArray[i]].pT() > maxpT){
	  maxpT = event[stauArray[i]].pT();
	  iStau = stauArray[i];
	}
      }
    }
    if(iStau < 0) continue; 
    neffA++;
     
    if(!printFlag){
      lifetime = 6.58e-25 * 1.0e9/ event[stauArray[0]].mWidth();
      cout<<"Lifetime = "<<scientific<<lifetime<<endl;
      printFlag = true;
    }

    bool flagEsc = true;

    if(event[iStau].status() < 0){

      double decayLength = sqrt(pow2(event[iStau].xDec() - event[iStau].xProd()) 
				+ pow2(event[iStau].yDec() - event[iStau].yProd()) 
				+ pow2(event[iStau].zDec() - event[iStau].zProd()));

      // Compare with detector dimensions
      double r = sqrt(pow2(event[iStau].xDec() - event[iStau].xProd()) 
		      + pow2(event[iStau].yDec() - event[iStau].yProd()));
      double z = abs(event[iStau].zDec() - event[iStau].zProd());

      if (!isBatch){
	r_trk.fill(r);
      }

      rdist[iLen] = r;
      iLen++;


      // See how many escape the CMS detector

      if (r < 7.5e3 && abs(z) < 6.5e3 ) flagEsc = false;
      if (r < 7.0e3 && abs(z) < 10.0e3 ) flagEsc = false;


      // Look for disappearing track; ATLAS search
      
      if(!flagEsc && !flagLep) {

	nEvDissCand++;

	// Cuts for disappearing track search
	if (jets.size() < 1) continue;

	double phijmet = abs(jets.at(0).phi() - MET.phi());
	if(phijmet > 3.1416) phijmet -= 2*M_PI;	

	if (!isBatch) {
	  pTj1.fill(jets.at(0).pT()); 
	  met2.fill(detector.MET.pT());
	  phijm.fill(abs(phijmet));
	  op << event[iStau].eta() << "  "
	     << event[iStau].pAbs()/event[iStau].e() << "  "
	     << event[iStau].pT() << "  "	  
	     << r << "  ";
	}

	double beta = event[iStau].pAbs()/event[iStau].e();
	double htsum = 0;

	for (int j=0; j < jets.size(); j++)
	  htsum += jets.at(j).pT();
	if (!isBatch){
	  op << jets.at(0).pT() << "  "
	     << MET.pT() << "  "
	     << phijmet << "  "
	     << htsum  << "  "
	     << htsum/jets.size() << "  "
	     << jets.size() << "  "
	     << endl;
	}
#if USE_ROOT
	if (!isBatch) {
	  pT_stau ->Fill(event[iStau].pT());
	  eta_stau ->Fill(event[iStau].eta());
	  beta_stau ->Fill(beta);
	  dl ->Fill(decayLength); 
	  rdl ->Fill(r);
	}
#endif
	// Cuts

	if(jets.at(0).pT() < 90.0) continue;
	neffB++;
	if (MET.pT() < 90.0) continue;
	neffMET++;
	if (abs(phijmet) < 1.5) continue;
	neffdPhi++;

	// Track pT and isolation

	pT_trk.fill(event[iStau].pT());

	if(event[iStau].pT() < 15.0) continue;
	double ptSum = 0.0 ;
	for (int j = 0; j < event.size(); ++j) 
	  if (event[j].isFinal() && event[j].isCharged()){
	    
	    double deltaphi = abs(event[j].phi() - event[iStau].phi());
	    if(deltaphi > 3.1416) deltaphi -= 3.1416;
	    double sep = sqrt(pow2(event[j].eta() - event[iStau].eta()) + pow2(deltaphi)); 
	    if(sep < 0.4) ptSum += event[j].pT();
	  }
      
	double iso40 = ptSum/event[iStau].pT();
	if(iso40 > 0.04) continue;

	neffC++;      
	if (!isBatch) {
	  isocone.fill(iso40);
#if USE_ROOT
	  pT_stau_acc ->Fill(event[iStau].pT());
	  eta_stau_acc ->Fill(event[iStau].eta());
	  beta_stau_acc ->Fill(beta);
	  dl_acc ->Fill(decayLength); 
	  rdl_acc ->Fill(r);
#endif
	}
      
	// SCT and TRT dimensions
	if (abs(event[iStau].eta()) < 0.1 || abs(event[iStau].eta()) > 1.9) continue; 
	if (r < 255 || r > 563 ) continue;
      
	nEvDiss++;

	double pTstau = event[iStau].pT();
	if(pTstau > 75.0) nDiss75++;
	if(pTstau > 100.0) nDiss100++;
	if(pTstau > 150.0) nDiss150++;
	if(pTstau > 200.0) nDiss200++;

      } // End unstable stau if-block
    }  // End disappearing track search

    // Stable tracks
    if (event[iStau].status() > 0 || flagEsc) {
      nEvStableCand++;
      // Needs to reach the tracker and fall within volume
      if (event[iStau].pT() < 10.0) continue;
      if (abs(event[iStau].eta()) > 2.4 ) continue;
      nEvStable++;

      double ETmiss = MET.pT();
      if(event[iStau].status < 0) ETmiss = (MET+event[iStau].p()).pT();
      op3 << ETmiss << endl; 
    }
   
  } // End of event loop.

  // Sort the array of decay lengths to get the median

  if (iLen > 0) {
    for (int j1 = 0; j1 < iLen-1; j1++){
      for (int j2 = j1+1; j2 < iLen; j2++){
	if(rdist[j2] > rdist[j1]) {
	  double temp = rdist[j1];
	  rdist[j1] = rdist[j2]; 
	  rdist[j2] = temp; 
	}
      }
    }
  }


  // Final statistics and histogram output.
  pythia.stat();
  pythia.rndm.dumpState("rndm.dat");

  if (!isBatch){
#if USE_ROOT
    pT_stau ->Write();
    eta_stau ->Write();
    beta_stau ->Write();
    dl ->Write();
    rdl ->Write();
    pT_stau_acc ->Write();
    eta_stau_acc ->Write();
    beta_stau_acc ->Write();
    dl_acc ->Write();
    rdl_acc ->Write();
#endif
  }

  if (nDecy + nStable == 0 )
    cout << "No stau events found"<<endl;
  else {

    if (!isBatch){
      cout << pTj1;
      cout << pT_trk;
      cout << r_trk;
      cout << met;
      cout << met2;
      cout << phijm;
      cout << isocone;
    }

    double sigma75 = nDiss75 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent;
    double sigma100 = nDiss100 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent;
    double sigma150 = nDiss150 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent;
    double sigma200 = nDiss200 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent;
    double sigmaStable = nEvStable * pythia.info.sigmaGen(0) * 1.0e12 / nEvent;
    
    // CMS TK+TOF study; interpolate limit using (200,1.0 fb),(300,0.3 fb)
    
    double CMSlimit = 0.3;
    if(mSt < 300.0) CMSlimit = -7.0e-3 * mSt + 2.4;

    op2 << scientific << m12 << "  " 
	<< mSt - mNeut << "  "
	<< lifetime << "  "
	<< nEv * 1.0/nEvent << "  "
	<< sigma75 << "  " << sigma75/1.76 << "  "
	<< sigma100 << "  " << sigma100/1.02 << "  "
	<< sigma150 << "  " << sigma150/0.62 << "  "
	<< sigma200 << "  " << sigma200/0.44 << "  "
	<< sigmaStable << "  " << sigmaStable/CMSlimit << "  "
	<< endl;
    if(!isBatch){

      op3<< scientific << mSt << "  " 
	  << mSt - mNeut << "  " 
	  << lifetime << "  "
	  << nEvStable*1.0/nEvent << endl;
      op3<< "--------" << endl;
      op3<< "Fraction of candidate events:"<< scientific<< nEv*1.0/nEvent<<endl;
      op3<< "Fraction of events with decaying staus:"<< scientific<< nDecy*1.0/nStau<<endl;
      op3<< "Fraction of events with Stable staus:"<< scientific<< nStable*1.0/nStau<<endl;
      op3<< "--------" << endl;
      op3<< "Efficiency for at least one disappearing stau track:"<< scientific<< nEvDissCand*1.0/nEvent<<endl;
      op3<< "Fraction of events with atleast one stable stau:"<< scientific<< nEvStableCand*1.0/nEvent<<endl;
      op3<< "--------" << endl;
      op3<< "Efficiency for |eta| < 2.4: "<< scientific<< neffA*1.0/nEvent<<endl;
      op3<< "Efficiency for |pt_j| > 90: "<< scientific<< neffB*1.0/nEvent<<endl;
      op3<< "Efficiency for MET > 90   : "<< scientific<< neffMET*1.0/nEvent<<endl;
      op3<< "Efficiency for dPhi > 1.5 : "<< scientific<< neffdPhi*1.0/nEvent<<endl;
      op3<< "Efficiency for isolation  : "<< scientific<< neffC*1.0/nEvent<<endl;
      op3<< "Efficiency for disappearing track : "<< scientific<< nEvDiss*1.0/nEvent<<endl;
      op3<< "--------" << endl;
      op3<<"pT > 75 Cross section (fb): " << scientific <<  sigma75 <<"  "<< sigma75/1.76 << endl;
      op3<<"pT > 100 Cross section (fb): " << scientific << sigma100 <<"  "<< sigma100/1.02  << endl;
      op3<<"pT > 150 Cross section (fb): " << scientific << sigma150 <<"  "<< sigma150/0.62  << endl;
      op3<<"pT > 200 Cross section (fb): " << scientific << sigma200 <<"  "<< sigma200/0.44  << endl;
      op3<< "--------" << endl;
      op3 << "Max radial decay Length = "<< rdist[0] <<endl;
      op3 << "75% Quartile radial decay Length = "<< rdist[iLen/4] <<endl;
      op3 << "Median radial decay Length = "<< rdist[iLen/2] <<endl;
      op3 << "25% Quartile radial decay Length = "<< rdist[3*iLen/4] <<endl;
      op3 << "Min radial decay Length = "<< rdist[iLen] <<endl;
      op3<< "--------" << endl;
      op3<<"Production Cross section (fb):" << scientific <<  pythia.info.sigmaGen(0) * 1.0e12 << endl;
      op3<<"Disappearing track Cross section (fb):" << scientific <<  nEvDiss * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;
      op3 << "CMSlimit = " << CMSlimit <<endl;
      op3<<"Stable Stau Cross section (fb):" << scientific << sigmaStable  << endl;
    }
  }

  if (!isBatch) op.close();
  op2.close();
#if USE_ROOT
  delete outFile;
#endif
  return 0;
}

double efficiency(double mass){

  double eff = (6.8/9.0) * (1.56e-4 * mass + 0.059) * (-1.53e-4 * mass + 1.25);
  return eff;

}
