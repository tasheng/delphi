#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include "TTree.h"
#include "TVector3.h"
#include <cstdlib>

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

//DataProcessing dependencies
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventSelection.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/thrustTools.h"
#include "DataProcessing/include/sphericityTools.h"

using namespace std;

void convert(const char* inputFileName, const char* outputFileName,
             bool verbose=false)
{
    // Input and output files
    std::ifstream infile(inputFileName);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open input file " << inputFileName << std::endl;
        return;
    }

    TFile outFile(outputFileName, "RECREATE");
    if (outFile.IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFileName << std::endl;
        return;
    }

    // Create a TTree to store event-based data
    TTree * out_t= new TTree("t", "t");
    particleData    out_pData;
    eventData       out_eData;

    TTree * out_tgen= new TTree("tgen", "tgen");
    particleData    out_pData_gen;
    eventData       out_eData_gen;

    // Variables for event information and particle data
    float emf[particleData::nMaxPart];
    float hpc[particleData::nMaxPart];
    float hac[particleData::nMaxPart];
    float stic[particleData::nMaxPart];
    float lock[particleData::nMaxPart];
    //out_t->Branch("EMF", emf, "emf[nMaxPart]/F");
    //out_t->Branch("HPC", hpc, "hpc[nMaxPart]/F");
    //out_t->Branch("HAC", hac, "hac[nMaxPart]/F");
    //out_t->Branch("STIC", stic, "stic[nMaxPart]/F");
    //out_t->Branch("LOCK", lock, "lock[nMaxPart]/F");
    out_t->Branch("EMF", emf);
    out_t->Branch("HPC", hpc);
    out_t->Branch("HAC", hac);
    out_t->Branch("STIC", stic);
    out_t->Branch("LOCK", lock);

    // register branches
    do_chThrust           = false;
    do_neuThrust          = false;
    do_thrustCorr         = false;
    do_thrustCorrInverse  = false;
    do_thrustMissP        = true;
    out_pData.SetBranchWrite(out_t, 1);
    out_eData.SetBranchWrite(out_t, 0);

    out_pData_gen.SetBranchWrite(out_tgen, 1);
    out_eData_gen.SetBranchWrite(out_tgen, 0);

    TDatabasePDG* pdgDatabase = TDatabasePDG::Instance();

    // Process the log file for reco
    std::string line;
    bool doParticle=0;
    bool doGenParticle=0;
    int iEvent=0;
    int nParticleNoCut = 0; int nParticle = 0;
    int nGenParticleNoCut = 0; int nGenParticle = 0;
    int nParticleHP = 0; int nChargedParticle = 0; int nChargedParticleHP = 0;
    int nGenParticleHP = 0; int nGenChargedParticle = 0; int nGenChargedParticleHP = 0;
    float energy = 0;
    int RunNo;
    int EventNo;
    TVector3 netP(0, 0, 0);
    TVector3 netChargedP(0, 0, 0);
    TVector3 netPGen(0, 0, 0);
    TVector3 netChargedPGen(0, 0, 0);
    while (std::getline(infile, line)) {
        if (line.find("HAPPY CHECK: EVENT") != std::string::npos) {
            // reset
            nParticleNoCut = 0;
	  // Locate and extract the part after "HAPPY CHECK: EVENT   :"
	  if(iEvent%1000 == 0)  cout<<"\r event processed : "<<iEvent<<flush;
	  std::istringstream iss(line);
	  std::string dummy;
	  iss >> dummy >> dummy >> dummy >> RunNo >> EventNo;
	  out_pData.RunNo = RunNo;
	  out_pData.EventNo = EventNo;

	  out_pData_gen.RunNo = out_pData.RunNo;
	  out_pData_gen.EventNo = out_pData.EventNo;
	  
	  if (verbose) std::cout << "Run: " << out_pData.RunNo << ", Event: " << out_pData.EventNo << std::endl;
	  out_pData.year = 1998; // temp //
	  // out_pData.subDir = -999;
	  // out_pData.process = -999;
	  out_pData.source = 70; // temp //
	  out_pData.isMC = true;
	  out_pData.isOnres = false;
	  out_pData_gen.isMC = true;
	  out_pData_gen.isOnres = false;
	  // out_pData.uniqueID = 0;
	  // out_pData.Energy = 0;
	  // out_pData.bFlag = -999;
	  out_pData.particleWeight = 1;
	  out_pData_gen.particleWeight = 1;
	  // out_pData.bx = -999;
	  // out_pData.by = -999;
	  // out_pData.ebx = -999;
	  // out_pData.eby = -999;
	
	  nParticle = 0;
	  nParticleHP = 0;
	  nChargedParticle = 0;
	  nChargedParticleHP = 0;
	  netP = TVector3(0, 0, 0);
	  netChargedP = TVector3(0, 0, 0);
	  
	  nGenParticle = 0;
	  nGenParticleHP = 0;
	  nGenChargedParticle = 0;
	  nGenChargedParticleHP = 0;
	  netPGen = TVector3(0, 0, 0);
	  netChargedPGen = TVector3(0, 0, 0);
	  
	} else if (line.find("CHECK: ECM") != std::string::npos) {
	    std::istringstream iss(line);
	    std::string dummy;
            iss >> dummy >> dummy >> energy;
	    if (verbose) std::cout <<"ECM: "<< energy << std::endl;
	    out_pData.Energy = energy;
	    out_pData_gen.Energy = energy;
	} else if (line.find("CHECK: TRACKS") != std::string::npos) {
            std::istringstream iss(line);
            std::string dummy;
            iss >> dummy >> dummy >> nParticleNoCut >> dummy >> dummy >> dummy;
            if (verbose) std::cout <<"get "<< nParticleNoCut << " reco particles" << std::endl;
            doParticle=1;
            if (nParticleNoCut == 0) {
                std::cout << "Filling Evt " << EventNo << " with 0 reco particles" << std::endl;
                doParticle=0;
                out_t->Fill();
                continue;
            }
        } else if (line.find("CHECK: GEN") != std::string::npos) {
            std::istringstream iss(line);
            std::string dummy;
            iss >> dummy >> dummy >> nGenParticleNoCut;
            if (verbose) std::cout <<"get "<< nGenParticleNoCut << " gen particles" << std::endl;
            doGenParticle=1;
        } else if (!line.empty() && doParticle==1) {
            std::istringstream iss(line);
            int particleNumber;
            float q, pxVal, pyVal, pzVal, eVal, emfVal, hpcVal, hacVal, sticVal, lockVal, d0, z0, length;
            iss >> particleNumber >> q >> pxVal >> pyVal >> pzVal >> eVal >> emfVal >> hpcVal >> hacVal >> sticVal >> lockVal >> d0 >> z0 >> length;

            TLorentzVector temp(pxVal, pyVal, pzVal, eVal);

            ///// veto the beam background and other backgrounds /////
            bool pass(1);
            pass = !lockVal; // bad particle
            ///// veto the beam background and other backgrounds /////

            if (pass)
            {
                if (verbose) std::cout <<"do particle: " << nParticle << "-th (" << particleNumber << ")" <<std::endl;
                netP -= TVector3(pxVal, pyVal, pzVal);
		if (abs(q) > 0.5) {netChargedP -= TVector3(pxVal, pyVal, pzVal);}

                out_pData.px[nParticle] = pxVal;
                out_pData.py[nParticle] = pyVal;
                out_pData.pz[nParticle] = pzVal;
                emf[nParticle] = emfVal;
                hpc[nParticle] = hpcVal;
                hac[nParticle] = hacVal;
                stic[nParticle] = sticVal;
                lock[nParticle] = lockVal;
                out_pData.pt[nParticle] = temp.Pt();
                out_pData.pmag[nParticle]   = temp.P();
                out_pData.rap[nParticle]    = temp.Rapidity();
                out_pData.eta[nParticle]    = temp.Eta();
                out_pData.theta[nParticle]  = temp.Theta();
                out_pData.phi[nParticle]    = temp.Phi();
                out_pData.mass[nParticle]   = temp.M();
                out_pData.charge[nParticle] = q;
                out_pData.pwflag[nParticle] = (out_pData.charge[nParticle]!=0)? 0: 4;
		out_pData.d0[nParticle] = d0;
 		out_pData.z0[nParticle] = z0;
 		out_pData.ntpc[nParticle] = (out_pData.charge[nParticle]!=0)? 7: 0;
		out_pData.weight[nParticle] = length; // use this variable to store track length for DELPHI 

                // follow the same definition in eventSelection.h
		if (out_pData.pwflag[nParticle]<=2) {
		  out_pData.highPurity[nParticle]= out_pData.pwflag[nParticle]<=2 && temp.Pt() >= 0.2;
		} else if (out_pData.pwflag[nParticle]==4) {
		  out_pData.highPurity[nParticle]= out_pData.pwflag[nParticle]==4 && temp.Pt() >= 0.4;
		}

                if(out_pData.pwflag[nParticle]<=2) {
                    nChargedParticle++;
                    if (out_pData.highPurity[nParticle]) nChargedParticleHP++;
                }
                if (out_pData.highPurity[nParticle]) nParticleHP++;
                nParticle++;
            }

            if (particleNumber==nParticleNoCut)  {
	      out_pData.nParticle       = nParticle;
	      out_eData.nChargedParticle  = nChargedParticle;
	      out_eData.nParticleHP       = nParticleHP;
	      out_eData.nChargedParticleHP= nChargedParticleHP;
	      doParticle=0;
	      if (verbose) std::cout <<"End reading particles, recording " << out_pData.nParticle << " particles." <<std::endl;

	      out_eData.missP = netP.Mag();
	      out_eData.missPt = netP.Perp();
	      out_eData.missTheta = netP.Theta();
	      out_eData.missPhi = netP.Phi();

	      out_eData.missChargedP = netChargedP.Mag();
	      out_eData.missChargedPt = netChargedP.Perp();
	      out_eData.missChargedTheta = netChargedP.Theta();
	      out_eData.missChargedPhi = netChargedP.Phi();

	      TVector3 thrust             = getThrust(out_pData.nParticle, out_pData.px, out_pData.py, out_pData.pz, THRUST::OPTIMAL);
	      TVector3 thrustWithMissP    = getThrust(out_pData.nParticle, out_pData.px, out_pData.py, out_pData.pz, THRUST::OPTIMAL,false,false,NULL,true,out_pData.pwflag);

	      out_eData.Thrust        = thrust.Mag();
	      out_eData.TTheta        = thrust.Theta();
	      out_eData.TPhi          = thrust.Phi();
	      out_eData.ThrustWithMissP   = thrustWithMissP.Mag();
	      out_eData.TThetaWithMissP   = thrustWithMissP.Theta();
	      out_eData.TPhiWithMissP     = thrustWithMissP.Phi();
	      if ( do_thrustMissP ) {
		setThrustVariables(&out_pData, &out_eData, TVector3(),
				   TVector3(), TVector3(), TVector3(), TVector3(),
				   thrustWithMissP);
	      }
	      
	      Sphericity spher = Sphericity(out_pData.nParticle,
					    out_pData.px,
					    out_pData.py,
					    out_pData.pz,
					    out_pData.pwflag,
					    false);
	      spher.setTree(&out_eData);
	      
	      eventSelection eSelection;
	      eSelection.setEventSelection(&out_pData, &out_eData);
	      
	      if (verbose)
                {
		  printf("eSelection.getPassesNeuNch(): %o\n", eSelection.getPassesNeuNch());
		  printf("eSelection.getPassesSTheta(): %o\n", eSelection.getPassesSTheta());
		  printf("eSelection.getPassesTotalChgEnergyMin(): %o\n", eSelection.getPassesTotalChgEnergyMin());
		  printf("out_eData.nChargedParticleHP: %d\n", out_eData.nChargedParticleHP);
		  printf("out_eData.nParticleHP: %d\n", out_eData.nParticleHP);
		  printf("out_pData.nParticle: %d\n", out_pData.nParticle);
                }
	      out_eData.passesTotalChgEnergyMin = eSelection.getPassesTotalChgEnergyMin();
	      out_eData.passesNeuNch = eSelection.getPassesNeuNch();
	      out_eData.passesNTrkMin = eSelection.getPassesNTrkMin();
	      out_eData.passesSTheta = eSelection.getPassesSTheta();
	      out_eData.passesBELLE = eSelection.getPassesNeuNch() &&	\
		eSelection.getPassesSTheta() &&				\
		eSelection.getPassesTotalChgEnergyMin() &&		\
		eSelection.getPassesNTrkMin();
	      
	      if (verbose) {
		std::size_t arraySize = sizeof(out_pData.mass) / sizeof(out_pData.mass[0]);
		std::cout << "The particle array contains " << arraySize << " elements." << std::endl;
	      }
	      
	      out_eData.passesISR = eSelection.getPassesISR();
	      out_eData.passesWW = eSelection.getPassesWW();
	      out_eData.Mvis = eSelection.getMvis();
	      out_eData.sPrime = eSelection.getsPrime();
	      out_eData.d2 = eSelection.getd2();
	      out_eData.cW = eSelection.getcW();
	      
	      // if(!out_eData.passesBELLE || !out_eData.passesISR) continue;
	      out_t->Fill();
	      iEvent ++;
	      if (verbose) std::cout <<"fill"<<std::endl;
            }
	} else if (!line.empty() && doGenParticle==1) {	    
            std::istringstream iss(line);
            int particleNumber;
            float particleID, pxVal, pyVal, pzVal, eVal, mVal, status;
            iss >> particleNumber >> particleID >> pxVal >> pyVal >> pzVal >> eVal >> mVal >> status;

            TLorentzVector temp(pxVal, pyVal, pzVal, eVal);

            ///// veto the beam background and other backgrounds ///// 
            bool pass(1);
            pass = (status == 1); // final state particle
            ///// veto the beam background and other backgrounds /////
	    
            if (pass)
            {
	        // calculate charge
	        float q(0);
	        TParticlePDG* particle = pdgDatabase->GetParticle(particleID);
		q = particle->Charge() / 3.0;
                if (verbose) std::cout <<"do gen particle: " << nGenParticle << "-th (" << particleNumber << ")" <<std::endl;
                netPGen -= TVector3(pxVal, pyVal, pzVal);

		if (abs(q) > 0.5) {netChargedPGen -= TVector3(pxVal, pyVal, pzVal);}

                out_pData_gen.px[nGenParticle] = pxVal;
                out_pData_gen.py[nGenParticle] = pyVal;
                out_pData_gen.pz[nGenParticle] = pzVal;
                out_pData_gen.pt[nGenParticle] = temp.Pt();
                out_pData_gen.pmag[nGenParticle]   = temp.P();
                out_pData_gen.rap[nGenParticle]    = temp.Rapidity();
                out_pData_gen.eta[nGenParticle]    = temp.Eta();
                out_pData_gen.theta[nGenParticle]  = temp.Theta();
                out_pData_gen.phi[nGenParticle]    = temp.Phi();
                out_pData_gen.mass[nGenParticle]   = temp.M();
                out_pData_gen.charge[nGenParticle] = q;
		out_pData_gen.pid[nGenParticle] = particleID;
                out_pData_gen.pwflag[nGenParticle] = (out_pData_gen.charge[nGenParticle]!=0)? 0: 4; // charged = 0, neutral = 4
		
                // follow the same definition in eventSelection.h
		if (out_pData_gen.pwflag[nParticle]<=2) {
		  out_pData_gen.highPurity[nParticle]= out_pData_gen.pwflag[nParticle]<=2 && temp.Pt() >= 0.2;
		} else if (out_pData_gen.pwflag[nParticle]==4) {
		  out_pData_gen.highPurity[nParticle]= out_pData_gen.pwflag[nParticle]==4 && temp.Pt() >= 0.4;
		}

                if(out_pData_gen.pwflag[nGenParticle]<=2) {
                    nGenChargedParticle++;
                    if (out_pData_gen.highPurity[nGenParticle]) nGenChargedParticleHP++;
                }
                if (out_pData_gen.highPurity[nGenParticle]) nGenParticleHP++;
                nGenParticle++;
            }

            if (particleNumber==nGenParticleNoCut)  {
                out_pData_gen.nParticle       = nGenParticle;
                out_eData_gen.nChargedParticle  = nGenChargedParticle;
                out_eData_gen.nParticleHP       = nGenParticleHP;
                out_eData_gen.nChargedParticleHP= nGenChargedParticleHP;
                doGenParticle=0;
                if (verbose) std::cout <<"End reading particles, recording " << out_pData_gen.nParticle << " particles." <<std::endl;

                out_eData_gen.missP = netPGen.Mag();
                out_eData_gen.missPt = netPGen.Perp();
                out_eData_gen.missTheta = netPGen.Theta();
                out_eData_gen.missPhi = netPGen.Phi();

		out_eData_gen.missChargedP = netChargedPGen.Mag();
		out_eData_gen.missChargedPt = netChargedPGen.Perp();
		out_eData_gen.missChargedTheta = netChargedPGen.Theta();
		out_eData_gen.missChargedPhi = netChargedPGen.Phi();

                TVector3 thrust             = getThrust(out_pData_gen.nParticle, out_pData_gen.px, out_pData_gen.py, out_pData_gen.pz, THRUST::OPTIMAL);
                TVector3 thrustWithMissP    = getThrust(out_pData_gen.nParticle, out_pData_gen.px, out_pData_gen.py, out_pData_gen.pz, THRUST::OPTIMAL,false,false,NULL,true,out_pData_gen.pwflag);

                out_eData_gen.Thrust        = thrust.Mag();
                out_eData_gen.TTheta        = thrust.Theta();
                out_eData_gen.TPhi          = thrust.Phi();
                out_eData_gen.ThrustWithMissP   = thrustWithMissP.Mag();
                out_eData_gen.TThetaWithMissP   = thrustWithMissP.Theta();
                out_eData_gen.TPhiWithMissP     = thrustWithMissP.Phi();
                if ( do_thrustMissP ) {
                        setThrustVariables(&out_pData_gen, &out_eData_gen, TVector3(),
                                TVector3(), TVector3(), TVector3(), TVector3(),
                                thrustWithMissP);
                }

                Sphericity spher = Sphericity(out_pData_gen.nParticle,
                                  out_pData_gen.px,
                                  out_pData_gen.py,
                                  out_pData_gen.pz,
                                  out_pData_gen.pwflag,
                                  false);
                spher.setTree(&out_eData_gen);

                eventSelection eSelection;
                eSelection.setEventSelection(&out_pData_gen, &out_eData_gen);

                if (verbose)
                {
                    printf("eSelection.getPassesNeuNch(): %o\n", eSelection.getPassesNeuNch());
                    printf("eSelection.getPassesSTheta(): %o\n", eSelection.getPassesSTheta());
                    printf("eSelection.getPassesTotalChgEnergyMin(): %o\n", eSelection.getPassesTotalChgEnergyMin());
                    printf("out_eData_gen.nChargedParticleHP: %d\n", out_eData_gen.nChargedParticleHP);
                }
		out_eData_gen.passesTotalChgEnergyMin = eSelection.getPassesTotalChgEnergyMin();
		out_eData_gen.passesNeuNch = eSelection.getPassesNeuNch();
		out_eData_gen.passesNTrkMin = eSelection.getPassesNTrkMin();
		out_eData_gen.passesSTheta = eSelection.getPassesSTheta();
                out_eData_gen.passesBELLE = eSelection.getPassesNeuNch() && \
		  eSelection.getPassesSTheta() &&			\
		  eSelection.getPassesTotalChgEnergyMin() &&		\
		  (out_eData_gen.nChargedParticleHP>=5);
                out_eData_gen.passesISR = eSelection.getPassesISR();
                out_eData_gen.passesWW = eSelection.getPassesWW();
                // if(!out_eData_gen.passesBELLE || !out_eData_gen.passesISR) continue;

                out_tgen->Fill();
                iEvent ++;
                // if (iEvent==10) break;
                if (verbose) std::cout <<"fill"<<std::endl;
            }
        }
    }
    
    printf("Total %d events recorded", iEvent);

    // Write the tree to the output file
    out_t->Write();
    out_tgen->Write();
    outFile.Close();
    infile.close();

    std::cout << "Data successfully converted to " << outputFileName << std::endl;
}


int main(int argc, char const *argv[])
{
  convert(argv[1], argv[2]);
    return 0;
}
