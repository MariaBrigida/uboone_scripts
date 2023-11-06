//Usage example, where 2g1p.list is a text file with a list of sample paths, the second argument is the number of protons (0 or 1), and the third argument is for using implied direction from vertex (1) or fitted (0):
//   root -l
//   .L ../makeInvMassPlots.cc
//   makeInvMassPlots("2g1p.list","2g0p",1)

#include <algorithm>

typedef std::map<int, TH1D*> IntTypeHistoMap;
typedef std::map<int, std::string> IntTypeLabelMap;

void makePlot(IntTypeHistoMap &histos, const std::string &canvasName, const std::string &xLabel, const std::string &yLabel,
    const std::string &fileName, const std::vector<int> &interactionTypes, IntTypeLabelMap &labels) {
    TCanvas *c = new TCanvas(canvasName.c_str());
    //c->SetCanvasSize(1024, 768);
    std::string opt{""};
    int color{1};
    int style{1};
    double maxVal{0};
    for (int intType : interactionTypes)
    {
        histos[intType]->Draw(opt.c_str());
        maxVal = std::max(maxVal, histos[intType]->GetMaximum());
        opt = "SAME";
    }
    auto legend{new TLegend(0.7, 0.55, 0.99, 0.9)};
    for (int intType : interactionTypes)
    {
        histos[intType]->SetTitle("");
        histos[intType]->GetXaxis()->SetTitle(xLabel.c_str());
        histos[intType]->GetYaxis()->SetTitle(yLabel.c_str());
        histos[intType]->GetYaxis()->SetRangeUser(0, maxVal + std::ceil(std::sqrt(maxVal)));
        histos[intType]->SetStats(0);
        histos[intType]->SetLineWidth(2);
        if (intType == 1006)
        {
            color = 1;
            style = 1;
        }
        else if (intType == 1013)
        {
            color = 1;
            style = 2;
        }
        else if (intType == 1002)
        {
            style = 1;
        }
        histos[intType]->SetLineColor(color++);
        histos[intType]->SetLineStyle(style);
        legend->AddEntry(histos[intType], labels[intType].c_str(), "l");
    }
    legend->Draw();

    c->Print(fileName.c_str());
    delete c;
}

std::vector<char*> GetFileNames(const char* filenamelist) {
    std::ifstream inFile(filenamelist);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
    }

    std::vector<char*> lines;
    std::string line;
    while (getline(inFile, line)) {
        char* cstr = new char[line.size() + 1];
        std::strcpy(cstr, line.c_str());
        lines.push_back(cstr);
    }
    inFile.close();

    return lines;
}

void makeInvMassPlots(const char* inputFilesList, const char* channel, const int implDir) {

    std::vector<int> interactionTypes{1006, 1007, 1008, 1009, 1013, 1014, 1015, 1016, 1002, 1092, 1096};
    IntTypeLabelMap intTypeLabels{
        {1002, "NCQE"},
        {1006, "#nu + p #rightarrow #nu + p + #pi^{0}"},
        {1007, "#nu + p #rightarrow #nu + n + #pi^{+}"},
        {1008, "#nu + n #rightarrow #nu + n + #pi^{0}"},
        {1009, "#nu + n #rightarrow #nu + p + #pi^{-}"},
        {1013, "#bar{#nu} + p #rightarrow #bar{#nu} + p + #pi^{0}"},
        {1014, "#bar{#nu} + p #rightarrow #bar{#nu} + n + #pi^{+}"},
        {1015, "#bar{#nu} + n #rightarrow #bar{#nu} + n + #pi^{0}"},
        {1016, "#bar{#nu} + n #rightarrow #bar{#nu} + p + #pi^{-}"},
        {1092, "NCDIS"},
        {1096, "NCCOH"}
//        {1000, "1000"}
    };

    const char* treeName = "singlephoton/vertex_tree";
    TChain chain(treeName); // Assuming the tree in each file is named "tree"
    std::vector<char*> lines = GetFileNames(inputFilesList);
    for(size_t i = 0; i < lines.size(); ++i) chain.Add(lines[i]);

    //Set branch addresses
    int eventNumber{0};
    int trueInteractionType{0};
    //Truth
    double trueLeadingPhotonEnergy, trueSubleadingPhotonEnergy;
    std::vector<unsigned long> *trackIndices=0;
    std::vector<unsigned long> *showerIndices=0;
    std::vector<double> *trueLeadingPhotonMom=0;
    std::vector<double> *trueSubleadingPhotonMom=0;
    std::vector<double> *trueProtonPx=0;
    std::vector<double> *trueProtonPy=0;
    std::vector<double> *trueProtonPz=0;
    std::vector<double> *trueProtonEnergy=0;
    //Reco
    std::vector<double> *recoPhotonDirX=0;
    std::vector<double> *recoPhotonDirY=0;
    std::vector<double> *recoPhotonDirZ=0;
    std::vector<double> *recoPhotonImpliedDirX=0;
    std::vector<double> *recoPhotonImpliedDirY=0;
    std::vector<double> *recoPhotonImpliedDirZ=0;

    std::vector<double> *recoPhotonEnergy=0;
    std::vector<double> *recoProtonDirX=0;
    std::vector<double> *recoProtonDirY=0;
    std::vector<double> *recoProtonDirZ=0;
    std::vector<double> *recoProtonKineticEnergy=0;

    chain.SetBranchAddress("event_number", &eventNumber);
    chain.SetBranchAddress("mctruth_interaction_type", &trueInteractionType);
    chain.SetBranchAddress("i_trk", &trackIndices);
    chain.SetBranchAddress("i_shr", &showerIndices);
    chain.SetBranchAddress("mctruth_pi0_leading_photon_mom", &trueLeadingPhotonMom);
    chain.SetBranchAddress("mctruth_pi0_leading_photon_energy", &trueLeadingPhotonEnergy);
    chain.SetBranchAddress("mctruth_pi0_subleading_photon_mom", &trueSubleadingPhotonMom);
    chain.SetBranchAddress("mctruth_pi0_subleading_photon_energy", &trueSubleadingPhotonEnergy);
    chain.SetBranchAddress("mctruth_exiting_proton_px", &trueProtonPx);
    chain.SetBranchAddress("mctruth_exiting_proton_py", &trueProtonPy);
    chain.SetBranchAddress("mctruth_exiting_proton_pz", &trueProtonPz);
    chain.SetBranchAddress("mctruth_exiting_proton_energy", &trueProtonEnergy);
    chain.SetBranchAddress("reco_track_dirx",&recoProtonDirX);
    chain.SetBranchAddress("reco_track_diry",&recoProtonDirY);
    chain.SetBranchAddress("reco_track_dirz",&recoProtonDirZ);
    chain.SetBranchAddress("reco_track_proton_kinetic_energy",&recoProtonKineticEnergy);
    chain.SetBranchAddress("reco_shower_energy_max",&recoPhotonEnergy);
    chain.SetBranchAddress("reco_shower_implied_dirx",&recoPhotonImpliedDirX);
    chain.SetBranchAddress("reco_shower_implied_diry",&recoPhotonImpliedDirY);
    chain.SetBranchAddress("reco_shower_implied_dirz",&recoPhotonImpliedDirZ);
    chain.SetBranchAddress("reco_shower_dirx",&recoPhotonDirX);
    chain.SetBranchAddress("reco_shower_diry",&recoPhotonDirY);
    chain.SetBranchAddress("reco_shower_dirz",&recoPhotonDirZ);


    IntTypeHistoMap hTrueDeltaInvMass, hRecoDeltaInvMass, hDeltaInvMassRes, hTruePi0InvMass, hRecoPi0InvMass, hPi0InvMassRes, hLeadingPhotonEnergyRes, hSubleadingPhotonEnergyRes, hProtonEnergyRes, hRecoProtonPi0Angle, hTrueProtonPi0Angle, hProtonPi0AngleRes, hRecoLeadingPhotonEnergy, hTrueLeadingPhotonEnergy, hRecoSubleadingPhotonEnergy, hTrueSubleadingPhotonEnergy, hRecoOtherPhotonsEnergy;
    for (int intType : interactionTypes) {
        hTrueDeltaInvMass[intType] = new TH1D("hTrueDeltaInvMass","hTrueDeltaInvMass",20,0.8,1.8);
        hRecoDeltaInvMass[intType] = new TH1D("hRecoDeltaInvMass","hRecoDeltaInvMass",20,0.8,1.8);
        hDeltaInvMassRes[intType] = new TH1D("hRecoDeltaInvMassRes","hRecoDeltaInvMassRes",30,-0.6,0.6);
        hTruePi0InvMass[intType] = new TH1D("hTruePi0InvMass","hTruePi0InvMass",20,0,0.3);
        hRecoPi0InvMass[intType] = new TH1D("hRecoPi0InvMass","hRecoPi0InvMass",20,0,0.3);
        hPi0InvMassRes[intType] = new TH1D("hRecoPi0InvMassRes","hRecoPi0InvMassRes",20,-0.6,0.6);

        hLeadingPhotonEnergyRes[intType] = new TH1D("hLeadingPhotonEnergyRes","hLeadingPhotonEnergyRes",20,-0.6,0.6);
        hSubleadingPhotonEnergyRes[intType] = new TH1D("hSubleadingPhotonEnergyRes","hSubleadingPhotonEnergyRes",20,-0.6,0.6);
        hProtonEnergyRes[intType] = new TH1D("hProtonEnergyRes","hProtonEnergyRes",20,-0.4,0.4);
	    hRecoProtonPi0Angle[intType] = new TH1D("hRecoProtonPi0Angle","hRecoProtonPi0Angle",20,0,360);
	    hTrueProtonPi0Angle[intType] = new TH1D("hTrueProtonPi0Angle","hTrueProtonPi0Angle",20,0,360);
	    hProtonPi0AngleRes[intType] = new TH1D("hProtonPi0AngleRes","hProtonPi0AngleRes",20,-1,1);

	    hRecoLeadingPhotonEnergy[intType] = new TH1D("hRecoLeadingPhotonEnergy","hRecoLeadingPhotonEnergy",20,0,0.5);
	    hTrueLeadingPhotonEnergy[intType] = new TH1D("hTrueLeadingPhotonEnergy","hTrueLeadingPhotonEnergy",20,0,0.5);
	    hRecoSubleadingPhotonEnergy[intType] = new TH1D("hRecoSubleadingPhotonEnergy","hRecoSubleadingPhotonEnergy",20,0,0.4);
	    hTrueSubleadingPhotonEnergy[intType] = new TH1D("hTrueSubleadingPhotonEnergy","hTrueSubleadingPhotonEnergy",20,0,0.4);
	    hRecoOtherPhotonsEnergy[intType] = new TH1D("hRecoOtherPhotonsEnergy","hRecoOtherPhotonsEnergy",20,0,0.5);
    }
    Long64_t nEntries = chain.GetEntries();

    // Loop over entries in the chain
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        // Read data for the current entry
        chain.GetEntry(entry);
        if (trueInteractionType==1000) continue;
        if (std::find(interactionTypes.begin(), interactionTypes.end(), trueInteractionType) == interactionTypes.end())
            continue;
	    //Check that there is exactly 1 track, and at least 1 proton in event
	    if(trackIndices->size()!=1 || trueProtonEnergy->size()<1) continue;
	    int trackID=trackIndices->at(0);
	    double trueProtonE=trueProtonEnergy->at(trackID);
	    double trueProtonMomX=trueProtonPx->at(trackID);
	    double trueProtonMomY=trueProtonPy->at(trackID);
	    double trueProtonMomZ=trueProtonPz->at(trackID);

	    //Check that there is exactly 2 showers, and the true photon energies are > 0
	    if(showerIndices->size()<2) continue;
	    int showerID_0=showerIndices->at(0);
	    int showerID_1=showerIndices->at(1);
	    double trueLeadingPhotonE=trueLeadingPhotonEnergy;
	    double trueLeadingPhotonMomX=trueLeadingPhotonMom->at(0);
	    double trueLeadingPhotonMomY=trueLeadingPhotonMom->at(1);
	    double trueLeadingPhotonMomZ=trueLeadingPhotonMom->at(2);
	    double trueSubleadingPhotonE=trueSubleadingPhotonEnergy;
	    double trueSubleadingPhotonMomX=trueSubleadingPhotonMom->at(0);
	    double trueSubleadingPhotonMomY=trueSubleadingPhotonMom->at(1);
	    double trueSubleadingPhotonMomZ=trueSubleadingPhotonMom->at(2);

	    if(trueLeadingPhotonE<=0 || trueSubleadingPhotonE<=0) continue;

	    double Etot = trueProtonE+trueLeadingPhotonEnergy+trueSubleadingPhotonEnergy;
	    double Etot2 = Etot*Etot;
	    double Pxtot = trueLeadingPhotonMomX + trueSubleadingPhotonMomX + trueProtonMomX;
	    double Pytot = trueLeadingPhotonMomY + trueSubleadingPhotonMomY + trueProtonMomY;
	    double Pztot = trueLeadingPhotonMomZ + trueSubleadingPhotonMomZ + trueProtonMomZ;
	    double Ptot2 = Pxtot*Pxtot + Pytot*Pytot + Pztot*Pztot;
	    double trueDeltaInvMass = TMath::Sqrt(Etot2-Ptot2);
	    hTrueDeltaInvMass[trueInteractionType]->Fill(trueDeltaInvMass);
	    double EtotPi0 = trueLeadingPhotonEnergy+trueSubleadingPhotonEnergy;
	    double EtotPi02 = EtotPi0*EtotPi0;
	    double PxtotPi0 = trueLeadingPhotonMomX + trueSubleadingPhotonMomX;
	    double PytotPi0 = trueLeadingPhotonMomY + trueSubleadingPhotonMomY;
	    double PztotPi0 = trueLeadingPhotonMomZ + trueSubleadingPhotonMomZ;
	    double PtotPi02 = PxtotPi0*PxtotPi0 + PytotPi0*PytotPi0 + PztotPi0*PztotPi0;
	    double truePi0InvMass = TMath::Sqrt(EtotPi02-PtotPi02); 
	    hTruePi0InvMass[trueInteractionType]->Fill(truePi0InvMass);

            //Reco quantities
	    //Which is the most energetic photon?
	    int iRecoLeadingShower=-999;
	    int iRecoSubleadingShower=-999;
	    if(recoPhotonEnergy->at(showerIndices->at(0))>recoPhotonEnergy->at(showerIndices->at(1))) {iRecoLeadingShower = showerIndices->at(0);iRecoSubleadingShower=showerIndices->at(1);}
	    else {iRecoLeadingShower = showerIndices->at(1);iRecoSubleadingShower=showerIndices->at(0);}
	    double recoLeadingPhotonE = recoPhotonEnergy->at(iRecoLeadingShower);
	    double recoSubleadingPhotonE = recoPhotonEnergy->at(iRecoSubleadingShower);
		double recoLeadingPhotonDirX, recoLeadingPhotonDirY, recoLeadingPhotonDirZ, recoSubleadingPhotonDirX, recoSubleadingPhotonDirY, recoSubleadingPhotonDirZ;	
		if(implDir){
	        recoLeadingPhotonDirX = recoPhotonImpliedDirX->at(iRecoLeadingShower);	
	        recoLeadingPhotonDirY = recoPhotonImpliedDirY->at(iRecoLeadingShower);	
	        recoLeadingPhotonDirZ = recoPhotonImpliedDirZ->at(iRecoLeadingShower);	
	        recoSubleadingPhotonDirX = recoPhotonImpliedDirX->at(iRecoSubleadingShower);	
	        recoSubleadingPhotonDirY = recoPhotonImpliedDirY->at(iRecoSubleadingShower);	
	        recoSubleadingPhotonDirZ = recoPhotonImpliedDirZ->at(iRecoSubleadingShower);
        }
		else{
	        recoLeadingPhotonDirX = recoPhotonDirX->at(iRecoLeadingShower);	
	        recoLeadingPhotonDirY = recoPhotonDirY->at(iRecoLeadingShower);	
	        recoLeadingPhotonDirZ = recoPhotonDirZ->at(iRecoLeadingShower);	
	        recoSubleadingPhotonDirX = recoPhotonDirX->at(iRecoSubleadingShower);	
	        recoSubleadingPhotonDirY = recoPhotonDirY->at(iRecoSubleadingShower);	
	        recoSubleadingPhotonDirZ = recoPhotonDirZ->at(iRecoSubleadingShower);
		}	
	    double recoProtonKineticE = recoProtonKineticEnergy->at(trackIndices->at(0));	
	    double recoProtonDirx=recoProtonDirX->at(trackIndices->at(0));
	    double recoProtonDiry=recoProtonDirY->at(trackIndices->at(0));
	    double recoProtonDirz=recoProtonDirZ->at(trackIndices->at(0));

	    //Calibrate photon energies
	    recoLeadingPhotonE=1.21989*recoLeadingPhotonE + 8.50486;
	    recoSubleadingPhotonE=1.21989*recoSubleadingPhotonE + 8.50486;

	    double protonMassInGeV = 0.938272;
	    double recoProtonE = recoProtonKineticE + protonMassInGeV;
	    double recoProtonMom = TMath::Sqrt(recoProtonE*recoProtonE - protonMassInGeV*protonMassInGeV); 
	    double recoProtonMomX = recoProtonDirx*recoProtonMom;
	    double recoProtonMomY = recoProtonDiry*recoProtonMom;
	    double recoProtonMomZ = recoProtonDirz*recoProtonMom;
	    double recoLeadingPhotonMomX=recoLeadingPhotonE*recoLeadingPhotonDirX;
	    double recoLeadingPhotonMomY=recoLeadingPhotonE*recoLeadingPhotonDirY;
	    double recoLeadingPhotonMomZ=recoLeadingPhotonE*recoLeadingPhotonDirZ;
	    double recoSubleadingPhotonMomX=recoSubleadingPhotonE*recoSubleadingPhotonDirX;
	    double recoSubleadingPhotonMomY=recoSubleadingPhotonE*recoSubleadingPhotonDirY;
	    double recoSubleadingPhotonMomZ=recoSubleadingPhotonE*recoSubleadingPhotonDirZ;
	    double cosTheta=recoLeadingPhotonDirX*recoSubleadingPhotonDirX+recoLeadingPhotonDirY*recoSubleadingPhotonDirY+recoLeadingPhotonDirZ*recoSubleadingPhotonDirZ;
	    TVector3 recoProtonMomVect(recoProtonDirx,recoProtonDiry,recoProtonDirz);
	    TVector3 recoPi0MomVect(recoLeadingPhotonMomX+recoSubleadingPhotonMomX,recoLeadingPhotonMomY+recoSubleadingPhotonMomY,recoLeadingPhotonMomZ+recoSubleadingPhotonMomZ);
	    TVector3 trueProtonMomVect(trueProtonMomX,trueProtonMomY,trueProtonMomZ);
	    TVector3 truePi0MomVect(trueLeadingPhotonMomX+trueSubleadingPhotonMomX,trueLeadingPhotonMomY+trueSubleadingPhotonMomY,trueLeadingPhotonMomZ+trueSubleadingPhotonMomZ);

	    double recoProtonPi0Angle=recoProtonMomVect.Angle(recoPi0MomVect)/TMath::Pi()*360;
	    double trueProtonPi0Angle=trueProtonMomVect.Angle(truePi0MomVect)/TMath::Pi()*360;
	    double protonPi0AngleRes=(recoProtonPi0Angle-trueProtonPi0Angle)/recoProtonPi0Angle;

	    //std::cout << "DEBUG recoProtonE = " << recoProtonE << std::endl;
	    //std::cout << "DEBUG recoProtonMomX = " << recoProtonMomX << std::endl;
	    //std::cout << "DEBUG recoProtonMomY = " << recoProtonMomY << std::endl;
	    //std::cout << "DEBUG recoProtonMomZ = " << recoProtonMomZ << std::endl;
	    //std::cout << "DEBUG recoLeadingPhotonE = " << recoLeadingPhotonE << std::endl;
	    //std::cout << "DEBUG recoSubleadingPhotonE = " << recoSubleadingPhotonE << std::endl;
	    //std::cout << "DEBUG recoLeadingPhotonMomX = " << recoLeadingPhotonMomX << std::endl;
	    //std::cout << "DEBUG recoLeadingPhotonMomY = " << recoLeadingPhotonMomY << std::endl;
	    //std::cout << "DEBUG recoLeadingPhotonMomZ = " << recoLeadingPhotonMomZ << std::endl;
	    //std::cout << "DEBUG recoSubleadingPhotonMomX = " << recoSubleadingPhotonMomX << std::endl;
	    //std::cout << "DEBUG recoSubleadingPhotonMomY = " << recoSubleadingPhotonMomY << std::endl;
	    //std::cout << "DEBUG recoSubleadingPhotonMomZ = " << recoSubleadingPhotonMomZ << std::endl;
	    //std::cout << "DEBUG cosTheta = " << cosTheta << std::endl;

	    //Reco photon variables are in MeV...	
	    recoLeadingPhotonE=recoLeadingPhotonE/1000;
	    recoSubleadingPhotonE=recoSubleadingPhotonE/1000;
	    recoLeadingPhotonMomX=recoLeadingPhotonMomX/1000;
	    recoLeadingPhotonMomY=recoLeadingPhotonMomY/1000;
	    recoLeadingPhotonMomZ=recoLeadingPhotonMomZ/1000;
	    recoSubleadingPhotonMomX=recoSubleadingPhotonMomX/1000;
	    recoSubleadingPhotonMomY=recoSubleadingPhotonMomY/1000;
	    recoSubleadingPhotonMomZ=recoSubleadingPhotonMomZ/1000;


	    double recoDeltaInvMass = TMath::Sqrt(protonMassInGeV*protonMassInGeV + 2*recoLeadingPhotonE*recoSubleadingPhotonE*(1-cosTheta)+2*recoProtonE*(recoLeadingPhotonE+recoSubleadingPhotonE) -2*((recoLeadingPhotonMomX+recoSubleadingPhotonMomX)*recoProtonMomX+(recoLeadingPhotonMomY+recoSubleadingPhotonMomY)*recoProtonMomY+(recoLeadingPhotonMomZ+recoSubleadingPhotonMomZ)*recoProtonMomZ));
	    //std::cout << "trueDeltaInvMass = " << trueDeltaInvMass << std::endl;
	    //std::cout << "recoDeltaInvMass = " << recoDeltaInvMass << std::endl;
	    hRecoDeltaInvMass[trueInteractionType]->Fill(recoDeltaInvMass);
	    hDeltaInvMassRes[trueInteractionType]->Fill((recoDeltaInvMass-trueDeltaInvMass)/trueDeltaInvMass);

	    double recoPi0InvMass = TMath::Sqrt(2*recoLeadingPhotonE*recoSubleadingPhotonE*(1-cosTheta));
	    hRecoPi0InvMass[trueInteractionType]->Fill(recoPi0InvMass);
	    hPi0InvMassRes[trueInteractionType]->Fill((recoPi0InvMass-truePi0InvMass)/truePi0InvMass);
	    hLeadingPhotonEnergyRes[trueInteractionType]->Fill((recoLeadingPhotonE-trueLeadingPhotonE)/trueLeadingPhotonE);
	    hSubleadingPhotonEnergyRes[trueInteractionType]->Fill((recoSubleadingPhotonE-trueSubleadingPhotonE)/trueSubleadingPhotonE);
	    hProtonEnergyRes[trueInteractionType]->Fill((recoProtonE-trueProtonE)/trueProtonE);
	    hRecoProtonPi0Angle[trueInteractionType]->Fill(recoProtonPi0Angle);
	    hTrueProtonPi0Angle[trueInteractionType]->Fill(trueProtonPi0Angle);
	    hProtonPi0AngleRes[trueInteractionType]->Fill((recoProtonPi0Angle-trueProtonPi0Angle)/trueProtonPi0Angle);
	    hRecoLeadingPhotonEnergy[trueInteractionType]->Fill(recoLeadingPhotonE);
	    hTrueLeadingPhotonEnergy[trueInteractionType]->Fill(trueLeadingPhotonE);
	    hRecoSubleadingPhotonEnergy[trueInteractionType]->Fill(recoSubleadingPhotonE);
	    hTrueSubleadingPhotonEnergy[trueInteractionType]->Fill(trueSubleadingPhotonE);

	    //Fill other (i_shr > 2) photons energies plots  
        for(int iSh=0; iSh<showerIndices->size(); iSh++)
        {
            hRecoOtherPhotonsEnergy[trueInteractionType]->Fill(recoPhotonEnergy->at(showerIndices->at(iSh))/1000);
        }

    }
    TString implDirString="";
    if(implDir) implDirString="ImpliedShowerDirection_";
    else implDirString="FittedShowerDirection_";
    TString channelString=TString(channel)+"_";

    makePlot(hTrueDeltaInvMass, "trueDeltaInvMassCanvas", "Invariant Mass [GeV]", "N", (channelString+"trueDeltaInvMass.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hRecoDeltaInvMass, "recoDeltaInvMassCanvas", "Invariant Mass [GeV]", "N", (channelString+implDirString+"recoDeltaInvMass.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hDeltaInvMassRes, "deltaInvMassResCanvas", "Invariant Mass Resolution", "N", (channelString+implDirString+"deltaInvMassRes.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hTruePi0InvMass, "truePi0InvMassCanvas", "Invariant Mass [GeV]", "N", (channelString+"truePi0InvMass.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hRecoPi0InvMass, "recoPi0InvMassCanvas", "Invariant Mass [GeV]", "N", (channelString+implDirString+"recoPi0InvMass.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hPi0InvMassRes, "pi0InvMassResCanvas", "Invariant Mass Resolution", "N", (channelString+implDirString+"pi0InvMassRes.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hLeadingPhotonEnergyRes, "leadingPhotonEnergyResCanvas", "Energy Resolution", "N", (channelString+"leadingPhotonEnergyRes.pdf").Data(), interactionTypes,
        intTypeLabels);
    makePlot(hSubleadingPhotonEnergyRes, "subleadingPhotonEnergyResCanvas", "Energy Resolution", "N", (channelString+"subleadingPhotonEnergyRes.pdf").Data(), interactionTypes,
        intTypeLabels);
    makePlot(hRecoLeadingPhotonEnergy, "recoLeadingPhotonEnergy", "Energy [GeV]", "N", (channelString+"recoLeadingPhotonEnergy.pdf").Data(), interactionTypes,
        intTypeLabels);
    makePlot(hTrueLeadingPhotonEnergy, "trueLeadingPhotonEnergy", "Energy [GeV]", "N", (channelString+"trueLeadingPhotonEnergy.pdf").Data(), interactionTypes,
        intTypeLabels);
    makePlot(hRecoSubleadingPhotonEnergy, "recoSubleadingPhotonEnergy", "Energy [GeV]", "N", (channelString+"recoSubleadingPhotonEnergy.pdf").Data(), interactionTypes,
        intTypeLabels);
    makePlot(hTrueSubleadingPhotonEnergy, "trueSubleadingPhotonEnergy", "Energy [GeV]", "N", (channelString+"trueSubleadingPhotonEnergy.pdf").Data(), interactionTypes,
        intTypeLabels);
    makePlot(hRecoOtherPhotonsEnergy, "recoOtherPhotonsEnergy", "Energy [GeV]", "N", (channelString+"recoOtherPhotonsEnergy.pdf").Data(), interactionTypes,
        intTypeLabels);

    makePlot(hProtonEnergyRes, "protonEnergyResCanvas", "Invariant Mass Resolution", "N", (channelString+"protonEnergyRes.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hRecoProtonPi0Angle, "recoProtonPi0AngleCanvas", "Angle [deg]", "N", (channelString+implDirString+"recoProtonPi0Angle.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hTrueProtonPi0Angle, "trueProtonPi0AngleCanvas", "Angle [deg]", "N", (channelString+"trueProtonPi0Angle.pdf").Data(), interactionTypes, intTypeLabels);
    makePlot(hProtonPi0AngleRes, "protonPi0AngleResCanvas", "Angle Resolution", "N", (channelString+implDirString+"protonPi0AngleRes.pdf").Data(), interactionTypes, intTypeLabels);
}
