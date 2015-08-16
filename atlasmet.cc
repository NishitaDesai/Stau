// Written by Nishita Desai based on main24.cc Py8 example
// All input is specified in the so10.cmnd file.

#include "Pythia8/Pythia.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/ToyDetector.h"
#include <fstream>

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  if(argc < 2){
    cout << "Usage: ./atlasmet.x SLHA_filename <output>" << endl;
    return 1;
  }

  string SLHAfile(argv[1]);
  string SLHA="SLHA:file = "+ string(argv[1]);

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("atlasmet.cmnd");    

  pythia.readString(SLHA);
  cout << "Reading from " << SLHA << endl;

  string outputFile="";
  if(argc == 3){
    outputFile += string(argv[2]);
  }

  ofstream op, op2, opMet;
  if(outputFile.length() > 0){
    string out1 = outputFile+"-ATLAS.dat";
    string out2 = outputFile+"-METlimits.dat";
    op.open(out1.c_str(),ios::app);
    op2.open(out2.c_str(),ios::app);
  }
  else{
    op.open("output.dat",ios::app);
    op2.open("limits.dat",ios::app);
    opMet.open("met-meff.dat",ios::app);
  }
  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 

  // Initialize.
  pythia.init();

  pythia.rndm.readState("rndm.dat");

  // Histograms.

  // Begin event loop.
  int iAbort = 0;

  int neffMET = 0, neffj = 0, neffdPhi = 0;
  int n2jL = 0, n2jM = 0, n3jM = 0, n3jT = 0;
  int n4jM = 0, n4jT = 0, n5j = 0;
  int n6jL = 0, n6jM = 0, n6jT = 0;

  ToyDetector detector;

  ofstream op3;
  op3.open("plots.dat",ios::out);

  double msq = min(pythia.particleData.m0(2000001),pythia.particleData.m0(1000001));
  double mglu = pythia.particleData.m0(1000021); 

  double mhalf = ((CoupSUSY *) pythia.couplingsPtr)->slhaPtr->minpar(2);
  double m0 = ((CoupSUSY *) pythia.couplingsPtr)->slhaPtr->minpar(1);

  int iCand = 1000015;
  double mSt = pythia.particleData.m0(iCand);
  double mNeut = pythia.particleData.m0(1000022); 
  double deltam = mSt - mNeut;
  int nStau = 0;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    if(iEvent == 0) cout << "Starting event generation" << endl;

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    if(!detector.getObjects(event)) {

      cout << "No objects found" << endl;
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;

    } 
    else {
      // cout << "Size of (jets, electrons, muons) = ("<< jets.size() 
      // 	   <<", "<< electrons.size() << ", " << muons.size()<<")"<<endl;
      // detector.printObjects();
    }



    // Only look at events withno isolated leptons
    if(detector.leptons.size() > 0) continue;
      
    // Cuts for jets (starts with 0)
    if(detector.jets.size() < 2) continue;
    if(detector.jets.at(1).pT() < 60.0 ) continue;


    // Check if stau decays outside the detector

    int stauArray[10];
    int stauArrLen = 0;

    bool stauFlag = false;


    for (int i = 0; i < event.size(); ++i) {
      if (abs(event[i].id()) == iCand) {
	stauFlag = true;
	if(event[i].status() == -22 || event[i].status() > 0) {
	  stauArray[stauArrLen] = i;
	  stauArrLen++;
	}
      }
    }
    
    if(stauFlag) nStau++;

    // Compare with detector dimensions, add back to MET if decaying outside

    Vec4 StauMom = Vec4(0.,0.,0.,0.);

    for(int iS = 0; iS < stauArrLen; iS++){
      int iStau = stauArray[iS];
      bool flagEsc = true;
      if(event[iS].status() < 0) {
	double r = sqrt(pow2(event[iStau].xDec() - event[iStau].xProd()) 
			+ pow2(event[iStau].yDec() - event[iStau].yProd()));
	double z = abs(event[iStau].zDec() - event[iStau].zProd());
	
	// See how many decay inside the CMS detector
	
	if (r < 7.5e3 && abs(z) < 6.5e3 ) flagEsc = false;
	if (r < 7.0e3 && abs(z) < 10.0e3 ) flagEsc = false;
      }

      if(flagEsc) StauMom += event[iStau].p();
      // cout << r << " " << z << " " << flagEsc <<endl;
    }
    
    Vec4 eTmiss = detector.MET - StauMom;
    double etmiss = eTmiss.pT();
    
    // Missing energy angle
    double dphijetmin = 10.0;
    double dphijetmin3 = 10.0;


    for (int j=0; j < detector.jets.size(); j++) {

      if(detector.jets.at(j).pT() < 40.0) break;
      double phiJet = detector.jets.at(j).phi();
      double phijmet = phiJet - eTmiss.phi();
      if(abs(phijmet) < abs(dphijetmin)) dphijetmin = phijmet;
      if (j <= 2) dphijetmin3 = dphijetmin;
    }

    double htsum = 0;
    
    for (int j=0; j < detector.jets.size(); j++)
      htsum += detector.jets.at(j).pT();
    htsum += etmiss;

    // Plots of MET etc. for comparison
    op3 << detector.MET.pT() << "  " 
	<< etmiss << "  " << htsum << "  "
    	<< dphijetmin << "  " 
    	<< dphijetmin3 << "  " 
    	<< StauMom.pT() << "  "
    	<< endl;

    // Basic Cuts
    if (etmiss < 90.0) continue;
    neffMET++;
    if(detector.jets.at(0).pT() < 130.0) continue;
    if(detector.jets.at(1).pT() < 60.0) continue;
    neffj++;
    if(abs(dphijetmin3) < 0.4) continue;
    neffdPhi++;

    // 2j Loose
    if (abs(dphijetmin3) > 0.4  && 
	etmiss/htsum > 0.20 &&
	htsum > 1000) n2jL++;

    // 2j Medium
    if (abs(dphijetmin3) > 0.4  && 
	etmiss/sqrt(htsum) > 15.0 &&
	htsum > 1600) n2jM++;

    if(detector.jets.size() < 3) continue;

    // 3j Medium
    if (detector.jets.at(2).pT() > 60.0 &&
	etmiss/htsum > 0.3 &&
	htsum > 1800) n3jM++;

    // 3j Tight (modified to 1405.7875 3j cut)
    if (detector.jets.at(2).pT() > 60.0 &&
	etmiss/htsum > 0.3 &&
	htsum > 2200) n3jT++;

    if(detector.jets.size() < 4) continue;
    if (abs(dphijetmin) < 0.2) continue;

    if(detector.jets.at(3).pT() > 60.0) opMet << etmiss << "  "  << htsum <<endl;


    // 4j Medium
    if (detector.jets.at(3).pT() > 60.0 &&
	etmiss/htsum > 0.25 &&
	htsum > 1200) n4jM++;

    // 4j Tight
    if (detector.jets.at(3).pT() > 60.0 &&
	etmiss/htsum > 0.25 &&
	htsum > 2200) n4jT++;

    if(detector.jets.size() < 5) continue;

    // 5j 
    if (detector.jets.at(4).pT() > 60.0 &&
	etmiss/htsum > 0.2 &&
	htsum > 1600) n5j++;


    if(detector.jets.size() < 6) continue;

    // 6j Loose
    if (detector.jets.at(5).pT() > 60.0 &&
	etmiss/htsum > 0.15 &&
	htsum > 1000) n6jL++;

    // 6j Medium
    if (detector.jets.at(5).pT() > 60.0 &&
	etmiss/htsum > 0.2 &&
	htsum > 1200) n6jM++;

    // 6j Tight
    if (detector.jets.at(5).pT() > 60.0 &&
	etmiss/htsum > 0.25 &&
	htsum > 1500) n6jT++;

    
    // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.stat();
  double sigTot = pythia.info.sigmaGen(0) * 1.0e12;

  double eff_2jL = n2jL*1.0/nEvent;
  double eff_2jM = n2jM*1.0/nEvent;
  double eff_3jM = n3jM*1.0/nEvent;
  double eff_3jT = n3jT*1.0/nEvent;
  double eff_4jM = n4jM*1.0/nEvent;
  double eff_4jT = n4jT*1.0/nEvent;
  double eff_5j = n5j *1.0/nEvent;
  double eff_6jL = n6jL*1.0/nEvent;
  double eff_6jM = n6jM*1.0/nEvent;
  double eff_6jT = n6jT*1.0/nEvent;

  double lim_2jL = n2jL * sigTot / nEvent / 66.07;
  double lim_2jM = n2jM * sigTot / nEvent / 2.52;  
  double lim_3jM = n3jM * sigTot / nEvent / 0.73; 
  double lim_3jT = n3jT * sigTot / nEvent / 0.33;  
  double lim_4jM = n4jM * sigTot / nEvent / 4.00;  
  double lim_4jT = n4jT * sigTot / nEvent / 0.12;  
  double lim_5j = n5j * sigTot / nEvent / 0.77 ; 
  double lim_6jL = n6jL * sigTot / nEvent / 4.55;  
  double lim_6jM = n6jM * sigTot / nEvent / 1.41;  
  double lim_6jT = n6jT * sigTot / nEvent / 0.41;  

  op << scientific
     << mhalf << "  " << deltam << "  "
     << eff_2jL << "  " << sigTot * eff_2jL << "  " << lim_2jL << "  "
     << eff_2jM << "  " << sigTot * eff_2jM << "  " << lim_2jM << "  "
     << eff_3jM << "  " << sigTot * eff_3jM << "  " << lim_3jM << "  "
     << eff_3jT << "  " << sigTot * eff_3jT << "  " << lim_3jT << "  "
     << eff_4jM << "  " << sigTot * eff_4jM << "  " << lim_4jM << "  "
     << eff_4jT << "  " << sigTot * eff_4jT << "  " << lim_4jT << "  "
     << eff_5j << "  " << sigTot * eff_5j << "  " << lim_5j << "  "
     << eff_6jL << "  " << sigTot * eff_6jL << "  " << lim_6jL << "  "
     << eff_6jM << "  " << sigTot * eff_6jM << "  " << lim_6jM << "  "
     << eff_6jT << "  " << sigTot * eff_6jT << "  " << lim_6jT << "  "
     << endl;

  op2 << scientific
      << mhalf << "  " << deltam << "  "
      << nStau*1.0/nEvent << " "
      << lim_2jL << "  "
      << lim_2jM << "  "
      << lim_3jM << "  "
      << lim_3jT << "  "
      << lim_4jM << "  "
      << lim_4jT << "  "
      << lim_5j << "  "
      << lim_6jL << "  "
      << lim_6jM << "  "
      << lim_6jT << "  "
     << endl;

  // op.close();
  // op2.close();
  // op3.close();

  cout<< "--------" << endl;
  cout<< "Fraction of 2jL events:"<< scientific<< n2jL*1.0/nEvent<<endl;
  cout<< "Fraction of 2jM events:"<< scientific<< n2jM*1.0/nEvent<<endl;
  cout<< "Fraction of 3jM events:"<< scientific<< n3jM*1.0/nEvent<<endl;
  cout<< "Fraction of 3jT events:"<< scientific<< n3jT*1.0/nEvent<<endl;
  cout<< "Fraction of 4jM events:"<< scientific<< n4jM*1.0/nEvent<<endl;
  cout<< "Fraction of 4jT events:"<< scientific<< n4jT*1.0/nEvent<<endl;
  cout<< "Fraction of 5j events:"<< scientific<< n5j*1.0/nEvent<<endl;
  cout<< "Fraction of 6jL events:"<< scientific<< n6jL*1.0/nEvent<<endl;
  cout<< "Fraction of 6jM events:"<< scientific<< n6jM*1.0/nEvent<<endl;
  cout<< "Fraction of 6jT events:"<< scientific<< n6jT*1.0/nEvent<<endl;
  cout<< "--------" << endl;
  cout<<"Production Cross section (fb):" << scientific <<  sigTot << endl;
  cout<< "Cross Section of 2jL events:"<< scientific<< n2jL * sigTot / nEvent << "  " << n2jL * sigTot / nEvent / 66.07 << endl;
  cout<< "Cross Section of 2jM events:"<< scientific<< n2jM * sigTot / nEvent << "  " << n2jM * sigTot / nEvent / 2.52 << endl; 
  cout<< "Cross Section of 3jM events:"<< scientific<< n3jM * sigTot / nEvent << "  " << n3jM * sigTot / nEvent / 0.73 << endl; 
  cout<< "Cross Section of 3jT events:"<< scientific<< n3jT * sigTot / nEvent << "  " << n3jT * sigTot / nEvent / 0.33 << endl; 
  cout<< "Cross Section of 4jM events:"<< scientific<< n4jM * sigTot / nEvent << "  " << n4jM * sigTot / nEvent / 4.00 << endl; 
  cout<< "Cross Section of 4jT events:"<< scientific<< n4jT * sigTot / nEvent << "  " << n4jT * sigTot / nEvent / 0.12 << endl; 
  cout<< "Cross Section of 5j events:"<< scientific<< n5j * sigTot / nEvent << "  " << n5j * sigTot / nEvent / 0.77 << endl; 
  cout<< "Cross Section of 6jL events:"<< scientific<< n6jL * sigTot / nEvent << "  " << n6jL * sigTot / nEvent / 4.55 << endl; 
  cout<< "Cross Section of 6jM events:"<< scientific<< n6jM * sigTot / nEvent << "  " << n6jM * sigTot / nEvent / 1.41 << endl; 
  cout<< "Cross Section of 6jT events:"<< scientific<< n6jT * sigTot / nEvent << "  " << n6jT * sigTot / nEvent / 0.41 << endl; 
  cout<< "--------" << endl;

  // Limits
  if( n2jL * sigTot / nEvent > 66.07 ||
      n2jM * sigTot / nEvent > 2.52 ||
      n3jM * sigTot / nEvent > 0.73 ||
      n3jT * sigTot / nEvent > 0.33 ||
      n4jM * sigTot / nEvent > 4.00 ||
      n4jT * sigTot / nEvent > 0.12 ||
      n5j * sigTot / nEvent > 0.77 ||
      n6jL * sigTot / nEvent > 4.55 ||
      n6jM * sigTot / nEvent > 1.41 ||
      n6jT * sigTot / nEvent > 0.41 ) cout << "Point ruled out." << endl;
  else cout << "Point allowed." << 

  pythia.rndm.dumpState("rndm.dat");

  return 0;
}



