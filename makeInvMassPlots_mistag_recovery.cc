//Usage example, where 2g1p.list is a text file with a list of sample paths, the second argument is the number of protons (0 or 1), and the third argument is for using implied direction from vertex (1) or fitted (0):
//   root -l
//   .L ../makeInvMassPlots_mistage_recovery.cc
//   makeInvMassPlots("2g0p_stage-1.list","2g0p",1)

#include <algorithm>

typedef std::map<int, TH1D*> IntTypeHistoMap;
typedef std::map<int, TH2D*> IntTypeHistoMap2D;
typedef std::map<int, std::string> IntTypeLabelMap;


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

void makeInvMassPlots_mistag_recovery(const char* inputFilesList, const char* channel, const int implDir) {

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
    int runNumber{0};
    int subrunNumber{0};
    int trueInteractionType{0};
    int trueCCorNC{0};
    int trueNeutrinoPdg{0};
	int numExitingPhotons{0};
    double trueNeutrinoEnergy{0};
    double recoNuVertexX{0};
    double recoNuVertexY{0};
    double recoNuVertexZ{0};
    double trueNuVertexX{0};
    double trueNuVertexY{0};
    double trueNuVertexZ{0};

    std::vector<unsigned long> *trackIndices=0;
    std::vector<unsigned long> *showerIndices=0;

    //Truth
    double trueLeadingPhotonEnergy, trueSubleadingPhotonEnergy;
	std::vector<double> *trueLeadingPhotonMom=0;
	std::vector<double> *trueSubleadingPhotonMom=0;
    std::vector<unsigned long> *trueDaughtersPdg=0;
    std::vector<unsigned long> *trueDaughtersTrackID=0;		//I think here "trackID" is the index of true MC particle rather than a reco track ID
    std::vector<string> *trueDaughtersProcess=0;		//I think here "trackID" is the index of true MC particle rather than a reco track ID
    std::vector<unsigned long> *trueDaughtersMotherTrackID=0;		//I think here "trackID" is the index of true MC particle's mother
    std::vector<double> *trueDaughtersEnergy=0;
    std::vector<int> *trueDaughtersStatusCode=0;

	std::vector<int> *exitingPhotonTrackID=0;
	std::vector<int> *exitingPhotonMotherTrackID=0;
	std::vector<double> *exitingPhotonEnergy=0;
	std::vector<double> *exitingPhotonPx=0;
	std::vector<double> *exitingPhotonPy=0;
	std::vector<double> *exitingPhotonPz=0;
	



    //Sim
    std::vector<unsigned long> *simShowerPdg=0;
    std::vector<unsigned long> *simShowerTrackID=0;
    std::vector<unsigned long> *simShowerParentPdg=0;
    std::vector<unsigned long> *simShowerParentTrackID=0;
    std::vector<unsigned long> *simShowerIsTrueShower=0;
    std::vector<double> *simShowerNuScore=0;
    std::vector<unsigned long> *simShowerIsNuSlice=0;
    std::vector<unsigned long> *simShowerIsClearCosmic=0;
    std::vector<unsigned long> *simShowerSliceId=0;
    std::vector<double> *simShowerEnergy=0;

    std::vector<unsigned long> *simTrackPdg=0;
    std::vector<unsigned long> *simTrackTrackID=0;
    std::vector<unsigned long> *simTrackParentPdg=0;
    std::vector<unsigned long> *simTrackParentTrackID=0;
    std::vector<double> *simTrackNuScore=0;
    std::vector<unsigned long> *simTrackIsNuSlice=0;
    std::vector<unsigned long> *simTrackIsClearCosmic=0;
    std::vector<unsigned long> *simTrackSliceId=0;

    std::vector<double> *simTrackEnergy=0;
    std::vector<double> *simTrackOverlayFraction=0;
    std::vector<double> *simShowerOverlayFraction=0;
	std::vector<double> *simTrackPx=0;
	std::vector<double> *simTrackPy=0;
	std::vector<double> *simTrackPz=0;

    //Reco
    std::vector<double> *recoPhotonDirX=0;
    std::vector<double> *recoPhotonDirY=0;
    std::vector<double> *recoPhotonDirZ=0;
    std::vector<double> *recoPhotonImpliedDirX=0;
    std::vector<double> *recoPhotonImpliedDirY=0;
    std::vector<double> *recoPhotonImpliedDirZ=0;

    std::vector<double> *recoPhotonEnergy=0;
    std::vector<double> *recoTrackDirX=0;
    std::vector<double> *recoTrackDirY=0;
    std::vector<double> *recoTrackDirZ=0;
    std::vector<double> *recoProtonKineticEnergy=0;

    std::vector<double> *recoTrackCaloEnergyMax=0;
	std::vector<int> *recoTrackGoodCaloBestPlane=0;
	std::vector<double> *recoTrackGoodCaloPlane0=0;
	std::vector<double> *recoTrackGoodCaloPlane1=0;
	std::vector<double> *recoTrackGoodCaloPlane2=0;

	std::vector<unsigned long> *recoTrackSliceId=0;
	std::vector<unsigned long> *recoTrackNuScore=0;
	std::vector<unsigned long> *recoTrackIsNuSlice=0;
	std::vector<double> *recoTrackTrackScore=0;
	std::vector<unsigned long> *recoTrackIsClearCosmic=0;
	std::vector<unsigned int> *recoTrackPfParticlePdg=0;
 	std::vector<unsigned long> *recoShowerSliceId=0;
	std::vector<unsigned long> *recoShowerNuScore=0;
	std::vector<unsigned long> *recoShowerIsNuSlice=0;
	std::vector<double> *recoShowerTrackScore=0;
	std::vector<unsigned long> *recoShowerIsClearCosmic=0;
	std::vector<unsigned int> *recoShowerPfParticlePdg=0;

    chain.SetBranchAddress("event_number", &eventNumber);
    chain.SetBranchAddress("run_number", &runNumber);
    chain.SetBranchAddress("subrun_number", &subrunNumber);
    chain.SetBranchAddress("mctruth_nu_pdg", &trueNeutrinoPdg);
    chain.SetBranchAddress("mctruth_nu_E", &trueNeutrinoEnergy);
    chain.SetBranchAddress("mctruth_cc_or_nc", &trueCCorNC);
    chain.SetBranchAddress("mctruth_interaction_type", &trueInteractionType);
    chain.SetBranchAddress("mctruth_daughters_pdg", &trueDaughtersPdg);
    chain.SetBranchAddress("mctruth_daughters_E", &trueDaughtersEnergy);
    chain.SetBranchAddress("mctruth_daughters_trackID", &trueDaughtersTrackID);
    chain.SetBranchAddress("mctruth_daughters_mother_trackID", &trueDaughtersMotherTrackID);
    chain.SetBranchAddress("mctruth_daughters_process", &trueDaughtersProcess);
    chain.SetBranchAddress("mctruth_daughters_status_code", &trueDaughtersStatusCode);
    chain.SetBranchAddress("mctruth_pi0_leading_photon_energy", &trueLeadingPhotonEnergy);
    chain.SetBranchAddress("mctruth_pi0_subleading_photon_energy", &trueSubleadingPhotonEnergy);
    chain.SetBranchAddress("mctruth_pi0_leading_photon_mom", &trueLeadingPhotonMom);
    chain.SetBranchAddress("mctruth_pi0_subleading_photon_mom", &trueSubleadingPhotonMom);

	chain.SetBranchAddress("mctruth_num_exiting_photons", &numExitingPhotons);
	chain.SetBranchAddress("mctruth_exiting_photon_trackID", &exitingPhotonTrackID);
	chain.SetBranchAddress("mctruth_exiting_photon_mother_trackID", &exitingPhotonMotherTrackID);
	chain.SetBranchAddress("mctruth_exiting_photon_energy", &exitingPhotonEnergy);
	chain.SetBranchAddress("mctruth_exiting_photon_px", &exitingPhotonPx);
	chain.SetBranchAddress("mctruth_exiting_photon_py", &exitingPhotonPy);
	chain.SetBranchAddress("mctruth_exiting_photon_pz", &exitingPhotonPz);

	chain.SetBranchAddress("sim_track_pdg",&simTrackPdg);
	chain.SetBranchAddress("sim_track_parent_pdg",&simTrackParentPdg);
	chain.SetBranchAddress("sim_track_trackID",&simTrackTrackID);
	chain.SetBranchAddress("sim_track_energy",&simTrackEnergy);
	chain.SetBranchAddress("sim_track_sliceId",&simTrackSliceId);
	chain.SetBranchAddress("sim_track_nuscore",&simTrackNuScore);
	chain.SetBranchAddress("sim_track_isclearcosmic",&simTrackIsClearCosmic);
	chain.SetBranchAddress("sim_track_overlay_fraction",&simTrackOverlayFraction);
	chain.SetBranchAddress("sim_track_px",&simTrackPx);
	chain.SetBranchAddress("sim_track_py",&simTrackPy);
	chain.SetBranchAddress("sim_track_pz",&simTrackPz);


	chain.SetBranchAddress("sim_shower_pdg",&simShowerPdg);
	chain.SetBranchAddress("sim_shower_parent_pdg",&simShowerParentPdg);
	chain.SetBranchAddress("sim_shower_trackID",&simShowerTrackID);
	chain.SetBranchAddress("sim_shower_parent_trackID",&simShowerParentTrackID);
	chain.SetBranchAddress("sim_shower_energy",&simShowerEnergy);
	chain.SetBranchAddress("sim_shower_is_true_shower",&simShowerIsTrueShower);
	chain.SetBranchAddress("sim_shower_sliceId",&simShowerSliceId);
	chain.SetBranchAddress("sim_shower_nuscore",&simShowerNuScore);
	chain.SetBranchAddress("sim_shower_is_nusclice",&simShowerIsNuSlice);
	chain.SetBranchAddress("sim_shower_isclearcosmic",&simShowerIsClearCosmic);
	chain.SetBranchAddress("sim_shower_overlay_fraction",&simShowerOverlayFraction);

    chain.SetBranchAddress("i_trk", &trackIndices);
    chain.SetBranchAddress("i_shr", &showerIndices);
    chain.SetBranchAddress("reco_shower_energy_max",&recoPhotonEnergy);
    chain.SetBranchAddress("reco_shower_implied_dirx",&recoPhotonImpliedDirX);
    chain.SetBranchAddress("reco_shower_implied_diry",&recoPhotonImpliedDirY);
    chain.SetBranchAddress("reco_shower_implied_dirz",&recoPhotonImpliedDirZ);
    chain.SetBranchAddress("reco_shower_dirx",&recoPhotonDirX);
    chain.SetBranchAddress("reco_shower_diry",&recoPhotonDirY);
    chain.SetBranchAddress("reco_shower_dirz",&recoPhotonDirZ);
    chain.SetBranchAddress("reco_track_dirx",&recoTrackDirX);
    chain.SetBranchAddress("reco_track_diry",&recoTrackDirY);
    chain.SetBranchAddress("reco_track_dirz",&recoTrackDirZ);
    chain.SetBranchAddress("reco_track_proton_kinetic_energy",&recoProtonKineticEnergy);
	chain.SetBranchAddress("reco_track_calo_energy_max", &recoTrackCaloEnergyMax);
    chain.SetBranchAddress("reco_track_good_calo_best_plane", &recoTrackGoodCaloBestPlane);
    chain.SetBranchAddress("reco_track_good_calo_plane0", &recoTrackGoodCaloPlane0);
    chain.SetBranchAddress("reco_track_good_calo_plane1", &recoTrackGoodCaloPlane1);
    chain.SetBranchAddress("reco_track_good_calo_plane2", &recoTrackGoodCaloPlane2);

	chain.SetBranchAddress("reco_track_sliceId",&recoTrackSliceId);
	chain.SetBranchAddress("reco_track_nuscore",&recoTrackNuScore);
	chain.SetBranchAddress("reco_track_trackscore",&recoTrackTrackScore);
	chain.SetBranchAddress("reco_track_isclearcosmic",&recoTrackIsClearCosmic);
	chain.SetBranchAddress("reco_track_pfparticle_pdg",&recoTrackPfParticlePdg);
	chain.SetBranchAddress("reco_track_is_nuslice",&recoTrackIsNuSlice);
	chain.SetBranchAddress("reco_shower_sliceId",&recoShowerSliceId);
	chain.SetBranchAddress("reco_shower_nuscore",&recoShowerNuScore);
	chain.SetBranchAddress("reco_shower_trackscore",&recoShowerTrackScore);
	chain.SetBranchAddress("reco_shower_isclearcosmic",&recoShowerIsClearCosmic);
	chain.SetBranchAddress("reco_shower_pfparticle_pdg",&recoShowerPfParticlePdg);
	chain.SetBranchAddress("reco_shower_is_nuslice",&recoShowerIsNuSlice);

    Long64_t nEntries = chain.GetEntries();

    // Loop over entries in the chain
	int nAtLeastTwoGoodShowersOrAtLeastOneGoodTrackOneGoodShower(0), nOneGoodTrackOneGoodShower(0), nOneGoodPhotonTrackOneGoodPhotonShower(0), nOneTrackOneShower(0), nTotalExitingPhotons(0);

	//I will fill this with reco energy vs (true pi0-sim shower energy) for events with exactly 1t 1s that are pi0 daughters
	TH1D* hTrackMaxEnergyRes = new TH1D("hTrackMaxEnergyRes","hTrackMaxEnergyRes",40,-1,1);
	TH1D* hTrackBestPlaneEnergyRes = new TH1D("hTrackBestPlaneEnergyRes","hTrackBestPlaneEnergyRes",40,-1,1);
	TH1D* hShowerUncalibratedEnergyRes = new TH1D("hShowerUncalibratedEnergyRes","hShowerUncalibratedEnergyRes",40,-1,1);
	TH1D* hShowerCalibratedEnergyRes = new TH1D("hShowerCalibratedEnergyRes","hShowerCalibratedEnergyRes",40,-1,1);
	TH2D* hShowerUncalibratedEnergyVsEnergyRes = new TH2D("hShowerUncalibratedEnergyVsEnergyRes","hShowerUncalibratedEnergyVsEnergyRes",40,0,2,40,-1,1);
	TH2D* hShowerCalibratedEnergyVsEnergyRes = new TH2D("hShowerCalibratedEnergyVsEnergyRes","hShowerCalibratedEnergyVsEnergyRes",40,0,2,40,-1,1);
	TH2D* hTrackMaxEnergyVsMaxEnergyRes = new TH2D("hTrackMaxEnergyVsMaxEnergyRes","hTrackMaxEnergyVsMaxEnergyRes",40,0,2,40,-1,1);
	TH2D* hTrackMaxEnergyVsTrueEnergy = new TH2D("hTrackMaxEnergyVsTrueEnergy","hTrackMaxEnergyVsTrueEnergy",60,0,1,60,0,1);
	TH1D* hTrackMaxEnergyToTrueEnergyRatio = new TH1D("hTrackMaxEnergyToTrueEnergyRatio","hTrackMaxEnergyToTrueEnergyRatio",60,0,10);
	TH2D* hTrackTrueEnergyVsMaxEnergyToTrueEnergyRatio = new TH2D("hTrackTrueEnergyVsMaxEnergyToTrueEnergyRatio","hTrackTrueEnergyVsMaxEnergyToTrueEnergyRatio",60,0,1,60,0,10);
	TH2D* hTrackBestPlaneEnergyVsBestPlaneEnergyRes = new TH2D("hTrackBestPlaneEnergyVsBestPlaneEnergyRes","hTrackBestPlaneEnergyVsBestPlaneEnergyRes",40,0,2,40,-1,1);
	TH1D* hTrackDirectionRes = new TH1D("hTrackDirectionRes","hTrackDirectionRes",36,0,180);	
	//TH1D* hPi0Mass = new TH1D("hPi0Mass","hPi0Mass",30,0,0.3);	
	TH1D* hTrackEnergyRatio = new TH1D("hTrackEnergyRatio","hTrackEnergyRatio",60,-1,5);
	TH1D* hIsSimShowerCorrectMatch = new TH1D("hIsSimShowerCorrectMatch","hIsSimShowerCorrectMatch",2,0,2);

    TH1D* hShowerUncalibratedEnergyRes_correctMatch = new TH1D("hShowerUncalibratedEnergyRes_correctMatch","hShowerUncalibratedEnergyRes_correctMatch",40,-1,1);
    TH1D* hShowerCalibratedEnergyRes_correctMatch = new TH1D("hShowerCalibratedEnergyRes_correctMatch","hShowerCalibratedEnergyRes_correctMatch",40,-1,1);
	TH2D* hShowerUncalibratedEnergyVsEnergyRes_correctMatch = new TH2D("hShowerUncalibratedEnergyVsEnergyRes_correctMatch","hShowerUncalibratedEnergyVsEnergyRes_correctMatch",40,0,2,40,-1,1);
	TH2D* hShowerCalibratedEnergyVsEnergyRes_correctMatch = new TH2D("hShowerCalibratedEnergyVsEnergyRes_correctMatch","hShowerCalibratedEnergyVsEnergyRes_correctMatch",40,0,2,40,-1,1);

	TH1D* hMissingShowerEnergy = new TH1D("hMissingShowerEnergy","hMissingShowerEnergy",40,0,2);
	TH1D* hRecoTrackEnergy = new TH1D("hRecoTrackEnergy","hRecoTrackEnergy",40,0,2);

	TH2D* hMissingShowerVsRecoTrackEnergy = new TH2D("hMissingShowerVsRecoTrackEnergy","hMissingShowerVsRecoTrackEnergy",60,0,0.3,60,0,0.3);
	TH1D* hPi0InvMass = new TH1D("hPi0InvMass","hPi0InvMass",100,0,0.3);
	TH1D* hPi0InvMass_trackTruePhoton = new TH1D("hPi0InvMass_trackTruePhoton","hPi0InvMass_trackTruePhoton",100,0,0.3);
	TH1D* hPi0InvMass_showerTruePhoton = new TH1D("hPi0InvMass_showerTruePhoton","hPi0InvMass_showerTruePhoton",100,0,0.3);
	TH1D* hPi0InvMass_true = new TH1D("hPi0InvMass_true","hPi0InvMass_true",100,0,0.3);
	TH1D* hPi0InvMass_showerTruePhoton_trackTrueEnergy = new TH1D("hPi0InvMass_showerTruePhoton_trackTrueEnergy","hPi0InvMass_showerTruePhoton_trackTrueEnergy",100,0,0.3);
	TH1D* hPi0InvMass_showerTruePhoton_trackTrueDirection = new TH1D("hPi0InvMass_showerTruePhoton_trackTrueDirection","hPi0InvMass_showerTruePhoton_trackTrueDirection",100,0,0.3);
	TH1D* hPi0InvMass_showerRecoPhoton_trackTrueEnergy = new TH1D("hPi0InvMass_showerRecoPhoton_trackTrueEnergy","hPi0InvMass_showerRecoPhoton_trackTrueEnergy",100,0,0.3);
	TH1D* hPi0InvMass_showerRecoPhoton_trackTrueDirection = new TH1D("hPi0InvMass_showerRecoPhoton_trackTrueDirection","hPi0InvMass_showerRecoPhoton_trackTrueDirection",100,0,0.3);

    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        chain.GetEntry(entry);
        if (trueInteractionType==1000) continue;

		//Ask for exactly 1 shower and 1 track
		if(simShowerTrackID->size()!=1 || simTrackTrackID->size()!=1)  continue;

		//SELECTION
			int nSimTracksSmallOverlay(0), nPhotonGoodTracks(0), nPhotonGoodShowers(0), nPhotonTracks(0), nPhotonShowers(0);
		for(unsigned int iSimTrk=0; iSimTrk<simTrackTrackID->size(); iSimTrk++)
        {
			if(simTrackOverlayFraction->at(iSimTrk)<0.6)nSimTracksSmallOverlay++;
			if(simTrackOverlayFraction->at(iSimTrk)<0.6 && simTrackPdg->at(iSimTrk)==22 && simTrackParentPdg->at(iSimTrk)==111) nPhotonGoodTracks++;
			if(simTrackPdg->at(iSimTrk)==22 && simTrackParentPdg->at(iSimTrk)==111) nPhotonTracks++;
        }

		int nSimShowersNotClearCosmic(0), nSimShowersSmallOverlay(0), nGoodSimShowers(0);
		for(unsigned int iSimShw=0; iSimShw<simShowerTrackID->size(); iSimShw++)
        {
			if(simShowerOverlayFraction->at(iSimShw)<0.6)nSimShowersSmallOverlay++;
			if(simShowerOverlayFraction->at(iSimShw)<0.6 && simShowerPdg->at(iSimShw)==22 && simShowerParentPdg->at(iSimShw)==111) nPhotonGoodShowers++;
			if(/*simShowerPdg->at(iSimShw)==22 && */simShowerParentPdg->at(iSimShw)==111) nPhotonShowers++;
        }

		if(simShowerTrackID->size()==1 && simTrackTrackID->size()==1) nOneTrackOneShower++;
		if(nSimTracksSmallOverlay==1 && nSimShowersSmallOverlay==1) nOneGoodTrackOneGoodShower++;
		if(nPhotonGoodShowers==1 && nPhotonGoodTracks==1) nOneGoodPhotonTrackOneGoodPhotonShower++;
		if(nSimShowersSmallOverlay>=2 || (nSimTracksSmallOverlay >=1 && nSimShowersSmallOverlay>=1)) nAtLeastTwoGoodShowersOrAtLeastOneGoodTrackOneGoodShower++;

		//Only select events with exactly 1 track and 1 shower that are pi0 daughters and overlay fraction < 60%
		if(nSimTracksSmallOverlay !=1 || nSimShowersSmallOverlay!=1) continue;

		double missingShowerEnergy(-999999);
		bool isSimShowerCorrectMatch(false);
		bool showerIsLeadingPhoton(false);
		for(unsigned int iSimShw=0; iSimShw<simShowerTrackID->size(); iSimShw++)
        {
			if(simShowerEnergy->at(iSimShw)==trueLeadingPhotonEnergy || simShowerEnergy->at(iSimShw)==trueSubleadingPhotonEnergy) {
				if(simShowerEnergy->at(iSimShw)==trueLeadingPhotonEnergy) {showerIsLeadingPhoton=true; missingShowerEnergy=trueSubleadingPhotonEnergy;}
				else missingShowerEnergy=trueLeadingPhotonEnergy;
				hIsSimShowerCorrectMatch->Fill(1);
				isSimShowerCorrectMatch=true;
			}
			else hIsSimShowerCorrectMatch->Fill(0);
        }

		double showerEnergy(simShowerEnergy->at(0));
		double showerRecoEnergy(recoPhotonEnergy->at(showerIndices->at(0)));
		double calibratedShowerRecoEnergy = 1.21989*showerRecoEnergy + 8.50486;
        calibratedShowerRecoEnergy=calibratedShowerRecoEnergy/1000;
        showerRecoEnergy = showerRecoEnergy/1000;

		double trackEnergyMax=recoTrackCaloEnergyMax->at(trackIndices->at(0))/1000;
		double bestPlaneEnergy=recoTrackGoodCaloBestPlane->at(trackIndices->at(0));

		TVector3 simTrackDir(simTrackPx->at(0),simTrackPy->at(0), simTrackPz->at(0));
		TVector3 recoTrackDir(recoTrackDirX->at(trackIndices->at(0)),recoTrackDirY->at(trackIndices->at(0)),recoTrackDirZ->at(trackIndices->at(0)));

		hShowerUncalibratedEnergyRes->Fill((showerRecoEnergy-showerEnergy)/showerEnergy);
		hShowerUncalibratedEnergyVsEnergyRes->Fill(showerEnergy,(showerRecoEnergy-showerEnergy)/showerEnergy);
		hShowerCalibratedEnergyRes->Fill((calibratedShowerRecoEnergy-showerEnergy)/showerEnergy);
		hShowerCalibratedEnergyVsEnergyRes->Fill(showerEnergy,(calibratedShowerRecoEnergy-showerEnergy)/showerEnergy);
		
		TVector3 trackTrueMom(0,0,0),showerTrueMom(0,0,0);
		if(isSimShowerCorrectMatch){
			if(showerIsLeadingPhoton){
				trackTrueMom.SetX(trueSubleadingPhotonMom->at(0));
				trackTrueMom.SetY(trueSubleadingPhotonMom->at(1));
				trackTrueMom.SetZ(trueSubleadingPhotonMom->at(2));
				showerTrueMom.SetX(trueLeadingPhotonMom->at(0));
				showerTrueMom.SetY(trueLeadingPhotonMom->at(1));
				showerTrueMom.SetZ(trueLeadingPhotonMom->at(2));
			}
			else{
				trackTrueMom.SetX(trueLeadingPhotonMom->at(0));
				trackTrueMom.SetY(trueLeadingPhotonMom->at(1));
				trackTrueMom.SetZ(trueLeadingPhotonMom->at(2));
				showerTrueMom.SetX(trueSubleadingPhotonMom->at(0));
				showerTrueMom.SetY(trueSubleadingPhotonMom->at(1));
				showerTrueMom.SetZ(trueSubleadingPhotonMom->at(2));
			}

			hShowerUncalibratedEnergyRes_correctMatch->Fill((showerRecoEnergy-showerEnergy)/showerEnergy);
			hShowerUncalibratedEnergyVsEnergyRes_correctMatch->Fill(showerEnergy,(showerRecoEnergy-showerEnergy)/showerEnergy);
			hShowerCalibratedEnergyRes_correctMatch->Fill((calibratedShowerRecoEnergy-showerEnergy)/showerEnergy);
			hShowerCalibratedEnergyVsEnergyRes_correctMatch->Fill(showerEnergy,(calibratedShowerRecoEnergy-showerEnergy)/showerEnergy);
	 		hTrackEnergyRatio->Fill(trackEnergyMax/missingShowerEnergy);
			hTrackMaxEnergyRes->Fill((trackEnergyMax-missingShowerEnergy)/missingShowerEnergy);
			hTrackBestPlaneEnergyRes->Fill((bestPlaneEnergy-missingShowerEnergy)/missingShowerEnergy);
			hTrackMaxEnergyVsMaxEnergyRes->Fill(missingShowerEnergy,(trackEnergyMax-missingShowerEnergy)/missingShowerEnergy);
			hTrackBestPlaneEnergyVsBestPlaneEnergyRes->Fill(missingShowerEnergy,(bestPlaneEnergy-missingShowerEnergy)/missingShowerEnergy);
			hMissingShowerEnergy->Fill(missingShowerEnergy);
			hRecoTrackEnergy->Fill(trackEnergyMax);
			hMissingShowerVsRecoTrackEnergy->Fill(missingShowerEnergy,trackEnergyMax);

			TVector3 showerDirection(recoPhotonDirX->at(showerIndices->at(0)),recoPhotonDirY->at(showerIndices->at(0)),recoPhotonDirZ->at(showerIndices->at(0)));

			double pi0InvMass=TMath::Sqrt(2*calibratedShowerRecoEnergy*trackEnergyMax*(1-TMath::Cos(showerDirection.Angle(recoTrackDir))));
			double pi0InvMass_true=TMath::Sqrt(2*showerTrueMom.Mag()*trackTrueMom.Mag()*(1-TMath::Cos(showerTrueMom.Angle(trackTrueMom))));
			double pi0InvMass_trackTruePhoton=TMath::Sqrt(2*calibratedShowerRecoEnergy*trackTrueMom.Mag()*(1-TMath::Cos(showerDirection.Angle(trackTrueMom))));
			double pi0InvMass_showerTruePhoton=TMath::Sqrt(2*showerTrueMom.Mag()*trackEnergyMax*(1-TMath::Cos(showerTrueMom.Angle(recoTrackDir))));
			double pi0InvMass_showerTruePhoton_trackTrueEnergy=TMath::Sqrt(2*showerTrueMom.Mag()*trackTrueMom.Mag()*(1-TMath::Cos(showerTrueMom.Angle(recoTrackDir))));
			double pi0InvMass_showerRecoPhoton_trackTrueEnergy=TMath::Sqrt(2*calibratedShowerRecoEnergy*trackTrueMom.Mag()*(1-TMath::Cos(showerDirection.Angle(recoTrackDir))));
			double pi0InvMass_showerTruePhoton_trackTrueDirection=TMath::Sqrt(2*showerTrueMom.Mag()*trackEnergyMax*(1-TMath::Cos(showerTrueMom.Angle(trackTrueMom))));
			double pi0InvMass_showerRecoPhoton_trackTrueDirection=TMath::Sqrt(2*calibratedShowerRecoEnergy*trackEnergyMax*(1-TMath::Cos(showerDirection.Angle(trackTrueMom))));

			hTrackMaxEnergyVsTrueEnergy->Fill(trackTrueMom.Mag(),trackEnergyMax);
			hTrackMaxEnergyToTrueEnergyRatio->Fill(trackEnergyMax/trackTrueMom.Mag());
			hTrackTrueEnergyVsMaxEnergyToTrueEnergyRatio->Fill(trackTrueMom.Mag(),trackEnergyMax/trackTrueMom.Mag());

       		hPi0InvMass->Fill(pi0InvMass);
       		hPi0InvMass_true->Fill(pi0InvMass_true);
			hPi0InvMass_trackTruePhoton->Fill(pi0InvMass_trackTruePhoton);
			hPi0InvMass_showerTruePhoton->Fill(pi0InvMass_showerTruePhoton);
			hPi0InvMass_showerTruePhoton_trackTrueDirection->Fill(pi0InvMass_showerTruePhoton_trackTrueDirection);
			hPi0InvMass_showerTruePhoton_trackTrueEnergy->Fill(pi0InvMass_showerTruePhoton_trackTrueEnergy);
			hPi0InvMass_showerRecoPhoton_trackTrueDirection->Fill(pi0InvMass_showerRecoPhoton_trackTrueDirection);
			hPi0InvMass_showerRecoPhoton_trackTrueEnergy->Fill(pi0InvMass_showerRecoPhoton_trackTrueEnergy);
		}
    }

	/*TCanvas *c3 = new TCanvas();
	hTrackMaxEnergyRes->Draw();
	TCanvas *c4 = new TCanvas();
	hTrackBestPlaneEnergyRes->Draw();
	TCanvas *c5 = new TCanvas();
	hTrackMaxEnergyVsMaxEnergyRes->Draw("colz");
	TCanvas *c6 = new TCanvas();
	hTrackBestPlaneEnergyVsBestPlaneEnergyRes->Draw("colz");
	TCanvas *c7 = new TCanvas();
	hShowerUncalibratedEnergyRes->Draw();	
	TCanvas *c8 = new TCanvas();
	hShowerCalibratedEnergyRes->Draw();	
	TCanvas *c9 = new TCanvas();
	hShowerUncalibratedEnergyVsEnergyRes->Draw("colz");
	TCanvas *c10 = new TCanvas();
	hShowerCalibratedEnergyVsEnergyRes->Draw("colz");
	TCanvas *c11 = new TCanvas();
	hTrackDirectionRes->Draw();
	TCanvas *c12 = new TCanvas();
	hTrackEnergyRatio->Draw();
	TCanvas *c13 = new TCanvas();
	hIsSimShowerCorrectMatch->Draw();
	TCanvas *c14 = new TCanvas();
	hShowerUncalibratedEnergyRes_correctMatch->Draw();
	TCanvas *c15 = new TCanvas();
	hShowerUncalibratedEnergyVsEnergyRes_correctMatch->Draw("colz");
	TCanvas *c16 = new TCanvas();
	hShowerCalibratedEnergyRes_correctMatch->Draw();
	TCanvas *c17 = new TCanvas();
	hShowerCalibratedEnergyVsEnergyRes_correctMatch->Draw("colz");
	TCanvas *c18 = new TCanvas();
	hMissingShowerEnergy->Draw();
	TCanvas *c19 = new TCanvas();
	hRecoTrackEnergy->Draw();
	TCanvas *c20 = new TCanvas();
	hMissingShowerVsRecoTrackEnergy->Draw("colz");
	TCanvas *c22 = new TCanvas();
	hTrackMaxEnergyVsTrueEnergy->Draw("colz");
	TCanvas *c23 = new TCanvas();
	hTrackMaxEnergyToTrueEnergyRatio->Draw();
	TCanvas *c24 = new TCanvas();
    TProfile *profileY = hTrackTrueEnergyVsMaxEnergyToTrueEnergyRatio->ProfileY("profileY");
	profileY->Draw();
*/
    auto legend1{new TLegend(0.7, 0.55, 0.99, 0.9)};
	TCanvas *cPi0Mass = new TCanvas();
	hPi0InvMass_trackTruePhoton->SetLineColor(kRed);
	hPi0InvMass_trackTruePhoton->Draw();
	hPi0InvMass_showerTruePhoton->SetLineColor(kBlue);
	hPi0InvMass_showerTruePhoton->Draw("same");
	hPi0InvMass_showerRecoPhoton_trackTrueEnergy->SetLineColor(kOrange);
	hPi0InvMass_showerRecoPhoton_trackTrueEnergy->Draw("same");
	hPi0InvMass_showerRecoPhoton_trackTrueDirection->SetLineColor(kGreen);
	hPi0InvMass_showerRecoPhoton_trackTrueDirection->Draw("same");
	hPi0InvMass->SetLineColor(kMagenta);
	hPi0InvMass->Draw("same");
    legend1->AddEntry(hPi0InvMass, "reco (shower), reco (track)", "l");
    legend1->AddEntry(hPi0InvMass_trackTruePhoton, "reco (shower), true (track)", "l");
    legend1->AddEntry(hPi0InvMass_showerTruePhoton, "true (shower), reco (track)", "l");
    legend1->AddEntry(hPi0InvMass_showerRecoPhoton_trackTrueEnergy, "reco (shower), reco dir true energy (track)", "l");
    legend1->AddEntry(hPi0InvMass_showerRecoPhoton_trackTrueDirection, "reco (shower), true dir reco energy (track)", "l");
	legend1->Draw("same");

    auto legend2{new TLegend(0.7, 0.55, 0.99, 0.9)};
	TCanvas *cPi0Mass_2 = new TCanvas();
	hPi0InvMass_trackTruePhoton->SetLineColor(kRed);
	hPi0InvMass_trackTruePhoton->Draw();
	hPi0InvMass_showerTruePhoton->SetLineColor(kBlue);
	hPi0InvMass_showerTruePhoton->Draw("same");
	hPi0InvMass_showerTruePhoton_trackTrueEnergy->SetLineColor(kOrange);
	hPi0InvMass_showerTruePhoton_trackTrueEnergy->Draw("same");
	hPi0InvMass_showerTruePhoton_trackTrueDirection->SetLineColor(kGreen);
	hPi0InvMass_showerTruePhoton_trackTrueDirection->Draw("same");
	hPi0InvMass->SetLineColor(kMagenta);
    legend2->AddEntry(hPi0InvMass, "reco (shower), reco (track)", "l");
    legend2->AddEntry(hPi0InvMass_trackTruePhoton, "reco (shower), true (track)", "l");
    legend2->AddEntry(hPi0InvMass_showerTruePhoton, "true (shower), reco (track)", "l");
    legend2->AddEntry(hPi0InvMass_showerTruePhoton_trackTrueEnergy, "true (shower), reco dir true energy (track)", "l");
    legend2->AddEntry(hPi0InvMass_showerTruePhoton_trackTrueDirection, "true (shower), true dir reco energy (track)", "l");
	legend2->Draw("same");

	hPi0InvMass->Draw("same");

}
