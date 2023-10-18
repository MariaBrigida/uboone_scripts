void calculateInvMass(){

    pMass=0.938272;
    pKineticEnergy=reco_track_proton_kinetic_energy[i_trl[0]];
    Ep=pMass+pKineticEnergy;
    Eg1=(1.21989*reco_shower_energy_max[i_shr[0]] + 8.50486);
    Eg2=(1.21989*reco_shower_energy_max[i_shr[1]] + 8.50486);
    E2=(Ep+Eg1+Eg2)*(Ep+Eg1+Eg2);
    Pg1x=Eg1*reco_shower_implied_dirx[i_shr[0]];
    Pg1y=Eg1*reco_shower_implied_diry[i_shr[0]];
    Pg1z=Eg1*reco_shower_implied_dirz[i_shr[0]];
    Pg2x=Eg2*reco_shower_implied_dirx[i_shr[1]];
    Pg2y=Eg2*reco_shower_implied_diry[i_shr[1]];
    Pg2z=Eg2*reco_shower_implied_dirz[i_shr[1]];
    Ppx=Ep*reco_track_dirx[i_trk[0]];
    Ppy=Ep*reco_track_diry[i_trk[0]];
    Ppz=Ep*reco_track_dirz[i_trk[0]];
    Px=Pg1x+Pg2x+Ppx;
    Py=Pg1y+Pg2y+Ppy;
    Pz=Pg1z+Pg2z+Ppz;
    P2=Px*Px+Py*Py+Pz*Pz
    InvMass=TMath::Sqrt(E2-P2);










}
