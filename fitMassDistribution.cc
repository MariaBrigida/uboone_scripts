void fitMassDistribution(){
    int nProtons=1;
    const char* inputFileName="";
    const char* directoryPath="";
    const char* histogramName="";
    const char* canvasName="";
    double fitRangeMin(0);
    double fitRangeMax(0);
    double drawRangeMin(0);
    double drawRangeMax(0);

    if(!nProtons){
        inputFileName="pi0massPlotsToFit/SBNfit_variation_spectra_v8_d15_08_23_exclusive_2g0p.root";
        directoryPath="v8_d15_08_23_exclusive_2g0p_CV_Dir";
        histogramName="Sys2g0p_numerator_reco_Signal";
	canvasName="Pi0InvMass_2g0p.png";
	fitRangeMin=0.07;
	fitRangeMax=0.18;
	drawRangeMin=0.03;
	drawRangeMax=0.25;

    }
    else{
	inputFileName="pi0massPlotsToFit/SBNfit_variation_spectra_v8_d15_08_23_exclusive_2g1p.root";
        directoryPath="v8_d15_08_23_exclusive_2g1p_CV_Dir";
        histogramName="Sys2g1p_numerator_reco_Signal";
	canvasName="Pi0InvMass_2g1p.png";
	fitRangeMin=0.07;
	fitRangeMax=0.18;
	drawRangeMin=0.03;
	drawRangeMax=0.25;
    }
    TFile* inputFile = TFile::Open(inputFileName);
    TDirectory* directory = inputFile->GetDirectory(directoryPath);

    TH1* hist = dynamic_cast<TH1*>(directory->Get(histogramName));
    hist->Rebin(2);

    TCanvas* canvas = new TCanvas("canvas", "Histogram Plot", 800, 600);

    //double fitRangeMin=hist->GetXaxis()->GetXmin();
    //double fitRangeMax=hist->GetXaxis()->GetXmax();

    hist->Draw("hist");
    canvas->Update();
    canvas->Modified();

    std::cout << "hist->GetMaximum() = " << hist->GetMaximum() << " hist->GetMean() = " << hist->GetMean() << " hist->GetRMS() = "<< hist->GetRMS() << std::endl;

    //CrystalBall fit
    TF1* crystalBallFit = new TF1("crystalBallFunc", "[0]*ROOT::Math::crystalball_pdf(x,[1],[2],[3],[4])", fitRangeMin, fitRangeMax);
    crystalBallFit->SetParNames("Normalization", "Alpha", "N", "Sigma","Mean");
    // Set initial parameter values for the Crystal Ball function
    crystalBallFit->SetParameters(1, 1, 10, hist->GetRMS(), hist->GetMean());
    hist->Fit(crystalBallFit, "0R");
    crystalBallFit->SetRange(drawRangeMin,drawRangeMax);
    crystalBallFit->Draw("same");

    canvas->Update();
    canvas->Modified();
    canvas->Print(canvasName);
    canvas->WaitPrimitive();
    inputFile->Close();
}



