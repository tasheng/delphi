#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>

void convert(const char* inputFileName, const char* outputFileName) {
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
    TTree tree("EventData", "Event-based particle data");

    // Variables for event information and particle data
    int eventNumber = 0;
    int nTracks = 0;
    std::vector<float> charge, px, py, pz, e, emf, hpc, hac, stic, lock;

    // Create branches
    tree.Branch("EventNumber", &eventNumber, "EventNumber/I");
    tree.Branch("NTracks", &nTracks, "NTracks/I");
    tree.Branch("Charge", &charge);
    tree.Branch("Px", &px);
    tree.Branch("Py", &py);
    tree.Branch("Pz", &pz);
    tree.Branch("E", &e);
    tree.Branch("EMF", &emf);
    tree.Branch("HPC", &hpc);
    tree.Branch("HAC", &hac);
    tree.Branch("STIC", &stic);
    tree.Branch("LOCK", &lock);

    // Process the log file
    std::string line;
    bool doParticle=0;
    while (std::getline(infile, line)) {
        if (line.find("HAPPY CHECK: EVENT") != std::string::npos) {
            // Extract event number
            std::istringstream iss(line);
            std::string dummy;
            char colon;
            int run, lumi;
            iss >> dummy >> dummy >> dummy >> run >> colon >> lumi >> colon >> eventNumber;

            // Clear the vectors for a new event
            charge.clear();
            px.clear();
            py.clear();
            pz.clear();
            e.clear();
            emf.clear();
            hpc.clear();
            hac.clear();
            stic.clear();
            lock.clear();
        } else if (line.find("CHECK: TRACKS") != std::string::npos) {
            // Extract the number of tracks
            std::istringstream iss(line);
            std::string dummy;
            iss >> dummy >> dummy >> dummy >> nTracks;
	    doParticle=1;
        } else if (!line.empty() && doParticle==1) {
            // Extract particle data
	    //std::cout <<"do particle"<<std::endl;
            std::istringstream iss(line);
            int particleNumber;
            float q, pxVal, pyVal, pzVal, eVal, emfVal, hpcVal, hacVal, sticVal, lockVal;
            iss >> particleNumber >> q >> pxVal >> pyVal >> pzVal >> eVal >> emfVal >> hpcVal >> hacVal >> sticVal >> lockVal;

            charge.push_back(q);
            px.push_back(pxVal);
            py.push_back(pyVal);
            pz.push_back(pzVal);
            e.push_back(eVal);
            emf.push_back(emfVal);
            hpc.push_back(hpcVal);
            hac.push_back(hacVal);
            stic.push_back(sticVal);
            lock.push_back(lockVal);
	    if (charge.size()==nTracks) doParticle=0; {
               if (!charge.empty()) {
                   tree.Fill();
	 	  //std::cout <<"fill"<<std::endl;
               }
	    }
        }
    }

    // Write the tree to the output file
    tree.Write();
    outFile.Close();
    infile.close();

    std::cout << "Data successfully converted to " << outputFileName << std::endl;
}


int main(int argc, char const *argv[])
{
    convert(argv[1], argv[2]);
    return 0;
}