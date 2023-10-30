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

void makeInvMassPlots() {
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
    };

    const char* fileName = "sbnfit_2g0p_NextGen_v4_stage_3_NCPi0All.root";
    const char* treeName = "singlephoton/vertex_tree";

    // Open the ROOT file containing the TTree
    TFile* file = TFile::Open(fileName);

    if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open the input file!" << std::endl;
        return;
    }

    // Access the TTree from the file
    TTree* tree = dynamic_cast<TTree*>(file->Get(treeName));

    if (!tree) {
        std::cerr << "Error: Tree '" << treeName << "' not found in the file." << std::endl;
        file->Close();
        return;
    }

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
    std::vector<double> *recoPhotonEnergy=0;
    std::vector<double> *recoProtonDirX=0;
    std::vector<double> *recoProtonDirY=0;
    std::vector<double> *recoProtonDirZ=0;
    std::vector<double> *recoProtonKineticEnergy=0;

    tree->SetBranchAddress("event_number", &eventNumber);
    tree->SetBranchAddress("mctruth_interaction_type", &trueInteractionType);
    tree->SetBranchAddress("i_trk", &trackIndices);
    tree->SetBranchAddress("i_shr", &showerIndices);
    tree->SetBranchAddress("mctruth_pi0_leading_photon_mom", &trueLeadingPhotonMom);
    tree->SetBranchAddress("mctruth_pi0_leading_photon_energy", &trueLeadingPhotonEnergy);
    tree->SetBranchAddress("mctruth_pi0_subleading_photon_mom", &trueSubleadingPhotonMom);
    tree->SetBranchAddress("mctruth_pi0_subleading_photon_energy", &trueSubleadingPhotonEnergy);
    tree->SetBranchAddress("mctruth_exiting_proton_px", &trueProtonPx);
    tree->SetBranchAddress("mctruth_exiting_proton_py", &trueProtonPy);
    tree->SetBranchAddress("mctruth_exiting_proton_pz", &trueProtonPz);
    tree->SetBranchAddress("mctruth_exiting_proton_energy", &trueProtonEnergy);
    tree->SetBranchAddress("reco_track_dirx",&recoProtonDirX);
    tree->SetBranchAddress("reco_track_diry",&recoProtonDirY);
    tree->SetBranchAddress("reco_track_dirz",&recoProtonDirZ);
    tree->SetBranchAddress("reco_track_proton_kinetic_energy",&recoProtonKineticEnergy);
    tree->SetBranchAddress("reco_shower_energy_max",&recoPhotonEnergy);
    tree->SetBranchAddress("reco_shower_implied_dirx",&recoPhotonDirX);
    tree->SetBranchAddress("reco_shower_implied_diry",&recoPhotonDirY);
    tree->SetBranchAddress("reco_shower_implied_dirz",&recoPhotonDirZ);

    IntTypeHistoMap hTrueDeltaInvMass, hRecoDeltaInvMass, hDeltaInvMassRes, hTruePi0InvMass, hRecoPi0InvMass, hPi0InvMassRes, hLeadingPhotonEnergyRes, hSubleadingPhotonEnergyRes, hProtonEnergyRes, hRecoProtonPi0Angle, hTrueProtonPi0Angle, hProtonPi0AngleRes;
    for (int intType : interactionTypes) {
        hTrueDeltaInvMass[intType] = new TH1D("hTrueDeltaInvMass","hTrueDeltaInvMass",60,0.8,2);
        hRecoDeltaInvMass[intType] = new TH1D("hRecoDeltaInvMass","hRecoDeltaInvMass",60,0.8,2);
        hDeltaInvMassRes[intType] = new TH1D("hRecoDeltaInvMassRes","hRecoDeltaInvMassRes",60,-0.6,0.6);
        hTruePi0InvMass[intType] = new TH1D("hTruePi0InvMass","hTruePi0InvMass",60,0,0.3);
        hRecoPi0InvMass[intType] = new TH1D("hRecoPi0InvMass","hRecoPi0InvMass",60,0,0.3);
        hPi0InvMassRes[intType] = new TH1D("hRecoPi0InvMassRes","hRecoPi0InvMassRes",60,-0.6,0.6);

        hLeadingPhotonEnergyRes[intType] = new TH1D("hLeadingPhotonEnergyRes","hLeadingPhotonEnergyRes",30,-0.6,0.6);
        hSubleadingPhotonEnergyRes[intType] = new TH1D("hSubleadingPhotonEnergyRes","hSubleadingPhotonEnergyRes",30,-0.6,0.6);
        hProtonEnergyRes[intType] = new TH1D("hProtonEnergyRes","hProtonEnergyRes",30,-0.6,0.6);
	hRecoProtonPi0Angle[intType] = new TH1D("hRecoProtonPi0Angle","hRecoProtonPi0Angle",60,0,360);
	hTrueProtonPi0Angle[intType] = new TH1D("hTrueProtonPi0Angle","hTrueProtonPi0Angle",60,0,360);
	hProtonPi0AngleRes[intType] = new TH1D("hProtonPi0AngleRes","hProtonPi0AngleRes",30,-1,1);
    }

    Long64_t nEntries = tree->GetEntries();

    // Loop over entries in the tree
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        // Read data for the current entry
        tree->GetEntry(entry);
	//std::cout << "-------------------------------------------------------" << std::endl;
	//std::cout << "entry n. " << entry << " event number = " << eventNumber << std::endl;

    if (std::find(interactionTypes.begin(), interactionTypes.end(), trueInteractionType) == interactionTypes.end())
        continue;
	//Check that there is exactly 1 track, and at least 1 proton in event
	if(trackIndices->size()!=1 || trueProtonEnergy->size()<1) continue;
	//std::cout << "trueProtonEnergy->size() = " << trueProtonEnergy->size() << " trueProtonPx->size() = " << trueProtonPx->size() << " trackIndices->size() = " << trackIndices->size() << std::endl;
	int trackID=trackIndices->at(0);
	double trueProtonE=trueProtonEnergy->at(trackID);
	double trueProtonMomX=trueProtonPx->at(trackID);
	double trueProtonMomY=trueProtonPy->at(trackID);
	double trueProtonMomZ=trueProtonPz->at(trackID);

	//Check that there is exactly 2 showers, and the true photon energies are > 0
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

	if(trueLeadingPhotonE<=0 || trueSubleadingPhotonE<=0 || showerIndices->size()<2) continue;

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
	std::cout << "iRecoLeadingShower = " << iRecoLeadingShower << " recoPhotonEnergy->at(0) = " << recoPhotonEnergy->at(0) << " recoPhotonEnergy->at(1) = " << recoPhotonEnergy->at(1) << std::endl;

	double recoLeadingPhotonE = recoPhotonEnergy->at(iRecoLeadingShower);	
	double recoLeadingPhotonDirX = recoPhotonDirX->at(iRecoLeadingShower);	
	double recoLeadingPhotonDirY = recoPhotonDirY->at(iRecoLeadingShower);	
	double recoLeadingPhotonDirZ = recoPhotonDirZ->at(iRecoLeadingShower);	
	double recoSubleadingPhotonE = recoPhotonEnergy->at(iRecoSubleadingShower);	
	double recoSubleadingPhotonDirX = recoPhotonDirX->at(iRecoSubleadingShower);	
	double recoSubleadingPhotonDirY = recoPhotonDirY->at(iRecoSubleadingShower);	
	double recoSubleadingPhotonDirZ = recoPhotonDirZ->at(iRecoSubleadingShower);
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

	std::cout << "DEBUG recoProtonE = " << recoProtonE << std::endl;
	std::cout << "DEBUG recoProtonMomX = " << recoProtonMomX << std::endl;
	std::cout << "DEBUG recoProtonMomY = " << recoProtonMomY << std::endl;
	std::cout << "DEBUG recoProtonMomZ = " << recoProtonMomZ << std::endl;
	std::cout << "DEBUG recoLeadingPhotonE = " << recoLeadingPhotonE << std::endl;
	std::cout << "DEBUG recoSubleadingPhotonE = " << recoSubleadingPhotonE << std::endl;
	std::cout << "DEBUG recoLeadingPhotonMomX = " << recoLeadingPhotonMomX << std::endl;
	std::cout << "DEBUG recoLeadingPhotonMomY = " << recoLeadingPhotonMomY << std::endl;
	std::cout << "DEBUG recoLeadingPhotonMomZ = " << recoLeadingPhotonMomZ << std::endl;
	std::cout << "DEBUG recoSubleadingPhotonMomX = " << recoSubleadingPhotonMomX << std::endl;
	std::cout << "DEBUG recoSubleadingPhotonMomY = " << recoSubleadingPhotonMomY << std::endl;
	std::cout << "DEBUG recoSubleadingPhotonMomZ = " << recoSubleadingPhotonMomZ << std::endl;
	std::cout << "DEBUG cosTheta = " << cosTheta << std::endl;

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
	std::cout << "trueDeltaInvMass = " << trueDeltaInvMass << std::endl;
	std::cout << "recoDeltaInvMass = " << recoDeltaInvMass << std::endl;
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

//        double recoPhoton0DirX=recoPhotonDirX->at(0);
//        double recoPhoton0DirX=recoPhrecoPhotonDirY=0;
//        double recoPhoton0DirX=recoPhrecoPhotonDirZ=0;
//        double recoPhoton0DirX=recoPhrecoPhotonEnergy=0;
//        double recoPhoton0DirX=recoPhrecoProtonDirX=0;
//        double recoPhoton0DirX=recoPhrecoProtonDirY=0;
//        double recoPhoton0DirX=recoPhrecoProtonDirZ=0;
//        double recoPhoton0DirX=recoPhrecoProtonKineticEnergy=0;


    }

    makePlot(hTrueDeltaInvMass, "trueDeltaInvMassCanvas", "Invariant Mass [GeV]", "N", "trueDeltaInvMass.pdf", interactionTypes, intTypeLabels);
    makePlot(hRecoDeltaInvMass, "recoDeltaInvMassCanvas", "Invariant Mass [GeV]", "N", "recoDeltaInvMass.pdf", interactionTypes, intTypeLabels);
    makePlot(hDeltaInvMassRes, "deltaInvMassResCanvas", "Invariant Mass Resolution", "N", "deltaInvMassRes.pdf", interactionTypes, intTypeLabels);
    makePlot(hTruePi0InvMass, "truePi0InvMassCanvas", "Invariant Mass [GeV]", "N", "truePi0InvMass.pdf", interactionTypes, intTypeLabels);
    makePlot(hRecoPi0InvMass, "recoPi0InvMassCanvas", "Invariant Mass [GeV]", "N", "recoPi0InvMass.pdf", interactionTypes, intTypeLabels);
    makePlot(hPi0InvMassRes, "pi0InvMassResCanvas", "Invariant Mass Resolution", "N", "pi0InvMassRes.pdf", interactionTypes, intTypeLabels);
    makePlot(hLeadingPhotonEnergyRes, "leadingPhotonEnergyResCanvas", "Energy Resolution", "N", "leadingPhotonEnergyRes.pdf", interactionTypes,
        intTypeLabels);
    makePlot(hSubleadingPhotonEnergyRes, "subleadingPhotonEnergyResCanvas", "Energy Resolution", "N", "subleadingPhotonEnergyRes.pdf", interactionTypes,
        intTypeLabels);
    makePlot(hProtonEnergyRes, "protonEnergyResCanvas", "Invariant Mass Resolution", "N", "protonEnergyRes.pdf", interactionTypes, intTypeLabels);
    makePlot(hRecoProtonPi0Angle, "recoProtonPi0AngleCanvas", "Angle [deg]", "N", "recoProtonPi0Angle.pdf", interactionTypes, intTypeLabels);
    makePlot(hTrueProtonPi0Angle, "trueProtonPi0AngleCanvas", "Angle [deg]", "N", "trueProtonPi0Angle.pdf", interactionTypes, intTypeLabels);
    makePlot(hProtonPi0AngleRes, "protonPi0AngleResCanvas", "Angle Resolution", "N", "protonPi0AngleRes.pdf", interactionTypes, intTypeLabels);
}
