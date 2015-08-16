// Modification of Py8 example main24
// For study of long-lived stau; pass to DELPHES via HepMC
// Nishita Desai

#include "Pythia8/Pythia.h"
#include "Pythia8/ToyDetector.h"

// #include "Pythia8/Pythia8toHepMC.h"
// #include "HepMC/GenEvent.h"
// #include "HepMC/IO_GenEvent.h"

using namespace Pythia8; 

int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // // Interface for conversion from Pythia8::Event to HepMC event. 
  // HepMC::Pythia8ToHepMC ToHepMC;

  // // Specify file where HepMC events will be stored.
  // HepMC::IO_GenEvent ascii_io("hepmcoutstau.dat", std::ios::out);


  // Read in commands from external file.
  pythia.readFile("stau.cmnd");    

  // Output file for plot
  fstream op;
  op.open("stau-output.dat",ios::out);

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  // double eCM   = pythia.parm("Beams:eCM");

  // Initialize.
  pythia.init();
  pythia.rndm.readState("rndm.dat");

  // Histograms.
  Hist pTj1("pT of 1st jet",100,0,200);
  Hist met("Missing ET", 100, 0, 500);
  Hist met2("Missing ET (Toy)", 100, 0, 500);
  Hist phijm("Delta phi(j,MET)", 100, 0, 6.28);  
  Hist pT_trk("pT of track",100,0,300);
  Hist r_trk("radial track length",100,0,10000.0);
  Hist isocone("isolation",100,0,0.20);


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

  SlowJet slowJet(-1.0, 0.4, 20.0, 4.0);

  ToyDetector detector;

  //int iCand = pythia.mode("SUSY:idA");
  int iCand = 1000015;
  double mSt = pythia.particleData.m0(iCand);
  double mNeut = pythia.particleData.m0(1000022); 


  double rdist[10000];
  int iLen = 0;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    if(!detector.getObjects(event)) {
      cout << "No objects found" << endl;
    } 
    else {
      detector.printObjects();
    }

    int stauArray[3];
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
      
    }
      
    if(stauArrLen < 1 )  {
      // event.list();
      continue;
    }

    nEv++;

    // Look at highest pT track
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
    if(iStau > 0) neffA++;
     
    if(!printFlag){
      lifetime = 6.58e-25 * 1.0e9/ event[stauArray[0]].mWidth();
      cout<<"Lifetime = "<<scientific<<lifetime<<endl;
      printFlag = true;
    }

    bool flagEsc = true;

    if(event[iStau].status() < 0){

      // Converted to meters
      // double decayLength = 1.0e-3 * sqrt(pow2(event[iStau].xDec() - event[iStau].xProd()) 
      // 				       + pow2(event[iStau].yDec() - event[iStau].yProd()) 
      // 				       + pow2(event[iStau].zDec() - event[iStau].zProd()));

      // Compare with detector dimensions
      double r = sqrt(pow2(event[iStau].xDec() - event[iStau].xProd()) 
		      + pow2(event[iStau].yDec() - event[iStau].yProd()));
      double z = abs(event[iStau].zDec() - event[iStau].zProd());

      // op << event[iStau].eta() << "  "
      //    << event[iStau].pAbs()/event[iStau].e() << "  "
      //    << r << "  "
      //    << endl;

      r_trk.fill(r);
      rdist[iLen] = r;
      iLen++;


      // See how many escape the CMS detector

      if (r < 7.5e3 && abs(z) < 6.5e3 ) flagEsc = false;
      if (r < 7.0e3 && abs(z) < 10.0e3 ) flagEsc = false;


      // Look for disappearing track; ATLAS search
      
      if(!flagEsc && !flagLep) {

	nEvDissCand++;

	// Cuts for disappearing track search
	slowJet.analyze(event);
	if(slowJet.sizeJet() < 1) continue;
      
	int iJet;
	for(iJet = 0; iJet < slowJet.sizeJet(); iJet++)
	  if(abs(slowJet.y(iJet)) < 2.8) break;

	pTj1.fill(slowJet.pT(iJet));      

	// Missing energy
	double pxsum = 0, pysum = 0;
	for (int j = 0; j < event.size(); ++j) {
	  if (event[j].isFinal() && event[j].eta() < 4.5 && event[j].isVisible()){
	    pxsum += event[j].px();
	    pysum += event[j].py();
	  }
	}

	double pymiss = -pysum;
	double pxmiss = -pxsum;

	double etmiss = sqrt(pxsum*pxsum + pysum*pysum);
	met.fill(etmiss);
	met2.fill(detector.MET.pT());

	
	double phiMet = atan(abs(pymiss/pxmiss));

	if(pymiss < 0.0 && pxmiss > 0.0) phiMet = -phiMet;
	if(pymiss < 0.0 && pxmiss < 0.0) phiMet = -M_PI + phiMet;
	if(pymiss > 0.0 && pxmiss < 0.0) phiMet = M_PI - phiMet;

	double phiJet = slowJet.phi(iJet);
	//if (phiJet < 0) phiJet += 2* M_PI;

	double phijmet = abs(phiJet - phiMet);

	if(phijmet > 3.1416) phijmet -= 2*M_PI;

	phijm.fill(abs(phijmet));

	// op << event[iStau].eta() << "  "
	//    << event[iStau].pAbs()/event[iStau].e() << "  "
	//    << event[iStau].pT() << "  "	  
	//    << r << "  ";


	double htsum = 0;

	for (int j=0; j<slowJet.sizeJet(); j++)
	  htsum += slowJet.pT(j);

	op << slowJet.pT(iJet) << "  "
	   << etmiss << "  "
	   << detector.MET.pT() << "  "
	   << phijmet << "  "
	   << phiJet << "  "
	   << phiMet << "  "
	   << htsum  << "  "
	   << htsum/slowJet.sizeJet() << "  "
	   << slowJet.sizeJet() << "  "
	   << endl;

	// Cuts

	if(slowJet.pT(iJet) < 90.0) continue;
	neffB++;
	if (etmiss < 90.0) continue;
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
      
	isocone.fill(iso40);
      
	double sep = 10.0;
	for (int i = 0; i < slowJet.sizeJet(); i++){ 
	  if(slowJet.pT(i) > 45.0){
	    double deltaphi = abs(slowJet.phi(i) - event[iStau].phi());
	    if(deltaphi > 3.1416) deltaphi -= 3.1416;
	    sep = min(sep, sqrt(pow2(slowJet.y(i) - event[iStau].eta()) + pow2(deltaphi)));
	  }
	}
      
	if(sep < 0.4) continue;

	neffC++;
      
	// SCT and TRT dimensions
	if (abs(event[iStau].eta()) < 0.1 || abs(event[iStau].eta()) > 1.9) continue; 
	if (r < 255 || r > 563 ) continue;
      
	nEvDiss++;

	double pTstau = event[iStau].pT();
	if(pTstau > 75.0) nDiss75++;
	if(pTstau > 100.0) nDiss100++;
	if(pTstau > 150.0) nDiss150++;
	if(pTstau > 200.0) nDiss200++;

      }
  
    }  // End disappearing track search


    // Stable tracks
    if (event[iStau].status() > 0 || flagEsc) {
      nEvStableCand++;
      if (event[iStau].pT() < 10.0) continue;
      nEvStable++;
    }

    // HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    // ToHepMC.fill_next_event( pythia, hepmcevt );

    
    // End of event loop.
  }

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
  op.close();

  // // Write the HepMC event to file. Done with it.
  // ascii_io << hepmcevt;
  // delete hepmcevt;

  if (nDecy + nStable == 0 )
    cout << "No stau events found"<<endl;
  else {

    cout << pTj1;
    cout << pT_trk;
    cout << r_trk;
    cout << met;
    cout << met2;
    cout << phijm;
    cout << isocone;

    cout<< scientific << mSt << "  " 
	<< mSt - mNeut << "  " 
	<< lifetime << "  "
	<< nEvStable*1.0/nEvent << endl;
    cout<< "--------" << endl;
    cout<< "Fraction of candidate events:"<< scientific<< nEv*1.0/nEvent<<endl;
    cout<< "Fraction of events with decaying staus:"<< scientific<< nDecy*1.0/nStau<<endl;
    cout<< "Fraction of events with Stable staus:"<< scientific<< nStable*1.0/nStau<<endl;
    cout<< "--------" << endl;
    cout<< "Efficiency for at least one disappearing stau track:"<< scientific<< nEvDissCand*1.0/nEvent<<endl;
    cout<< "Fraction of events with atleast one stable stau:"<< scientific<< nEvStableCand*1.0/nEvent<<endl;
    cout<< "--------" << endl;
    cout<< "Efficiency for |eta| < 2.4: "<< scientific<< neffA*1.0/nEvent<<endl;
    cout<< "Efficiency for |pt_j| > 90: "<< scientific<< neffB*1.0/nEvent<<endl;
    cout<< "Efficiency for MET > 90   : "<< scientific<< neffMET*1.0/nEvent<<endl;
    cout<< "Efficiency for dPhi > 1.5 : "<< scientific<< neffdPhi*1.0/nEvent<<endl;
    cout<< "Efficiency for isolation  : "<< scientific<< neffC*1.0/nEvent<<endl;
    cout<< "Efficiency for disappearing track : "<< scientific<< nEvDiss*1.0/nEvent<<endl;
    cout<< "--------" << endl;
    cout<<"pT > 75 Cross section (fb): " << scientific <<  nDiss75 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;
    cout<<"pT > 100 Cross section (fb): " << scientific <<  nDiss100 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;
    cout<<"pT > 150 Cross section (fb): " << scientific <<  nDiss150 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;
    cout<<"pT > 200 Cross section (fb): " << scientific <<  nDiss200 * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;
    cout<< "--------" << endl;
    cout << "Max radial decay Length = "<< rdist[0] <<endl;
    cout << "75% Quartile radial decay Length = "<< rdist[iLen/4] <<endl;
    cout << "Median radial decay Length = "<< rdist[iLen/2] <<endl;
    cout << "25% Quartile radial decay Length = "<< rdist[3*iLen/4] <<endl;
    cout << "Min radial decay Length = "<< rdist[iLen] <<endl;
    cout<< "--------" << endl;
    cout<<"Production Cross section (fb):" << scientific <<  pythia.info.sigmaGen(0) * 1.0e12 << endl;
    cout<<"Disappearing track Cross section (fb):" << scientific <<  nEvDiss * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;
    cout<<"Stable Stau Cross section (fb):" << scientific <<  nEvStable * pythia.info.sigmaGen(0) * 1.0e12 / nEvent << endl;

  }

  pythia.rndm.dumpState("rndm.dat");

  return 0;
}

double efficiency(double mass){

  double eff = (6.8/9.0) * (1.56e-4 * mass + 0.059) * (-1.53e-4 * mass + 1.25);
  return eff;

}
