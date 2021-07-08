//OUTDATED DUE TO CHANGED NAMES

void HokeyPlottingMacro(){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();

  gStyle->SetOptStat(0);
  
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(1,"Other");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(2,"");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(3,"n");           
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(4,"p");          
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(5,"#pi^{0}");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(6,"#pi^{+}");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(7,"#pi^{-}");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(8,"#gamma");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(9,"e");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetBinLabel(10,"#mu");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetLabelSize(0.06);
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetTitle("Blob Primary Parent");
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetTitleSize(0.045);
  hw_tracker_primary_parent_Tejin->GetXaxis()->SetTitleOffset(1.075);
  hw_tracker_primary_parent_Tejin->GetYaxis()->SetTitle("Blobs");
  hw_tracker_primary_parent_Tejin->GetYaxis()->SetTitleSize(0.045);
  hw_tracker_primary_parent_Tejin->GetYaxis()->SetTitleOffset(1.075);
  hw_tracker_primary_parent_Tejin->Draw();
  c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_tracker_primary_parent_Tejin.pdf");

  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(1,"Other");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(2,"");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(3,"n");           
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(4,"p");           
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(5,"#pi^{0}");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(6,"#pi^{+}");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(7,"#pi^{-}");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(8,"#gamma");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(9,"e");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(10,"#mu");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetLabelSize(0.06);
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitle("Blob Primary Parent");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitleSize(0.045);
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitleOffset(1.075);
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitle("Blobs");
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitleSize(0.045);
  hw_tracker_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitleOffset(1.075);
  hw_tracker_primary_parent_Tejin_TrackerONLY->Draw();
  c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_tracker_primary_parent_Tejin_TrackerONLY.pdf");

      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Recoil->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Recoil->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_Recoil->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_Recoil->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_Recoil->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_tracker_primary_parent_Recoil");

      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(1,"Other");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(2,"");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(3,"n");           
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(4,"p");           
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(9,"e");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetBinLabel(10,"#mu");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetLabelSize(0.06);
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_CCQE->GetXaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_CCQE->GetYaxis()->SetTitle("Blobs");
      hw_tracker_primary_parent_CCQE->GetYaxis()->SetTitleSize(0.045);
      hw_tracker_primary_parent_CCQE->GetYaxis()->SetTitleOffset(1.075);
      hw_tracker_primary_parent_CCQE->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_tracker_primary_parent_CCQE.pdf");

      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_Tejin->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_Tejin->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_Tejin->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_Tejin->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_Tejin->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_target_primary_parent_Tejin.pdf");

      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Tejin_TrackerONLY->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_target_primary_parent_Tejin_TrackerONLY.pdf");

      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_Recoil->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_Recoil->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_Recoil->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_Recoil->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Recoil->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Recoil->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_Recoil->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_Recoil->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_Recoil->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_target_primary_parent_Recoil.pdf");

      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(1,"Other");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(2,"");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(3,"n");           
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(4,"p");           
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(9,"e");
      hw_target_primary_parent_CCQE->GetXaxis()->SetBinLabel(10,"#mu");
      hw_target_primary_parent_CCQE->GetXaxis()->SetLabelSize(0.06);
      hw_target_primary_parent_CCQE->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_target_primary_parent_CCQE->GetXaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_CCQE->GetXaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_CCQE->GetYaxis()->SetTitle("Blobs");
      hw_target_primary_parent_CCQE->GetYaxis()->SetTitleSize(0.045);
      hw_target_primary_parent_CCQE->GetYaxis()->SetTitleOffset(1.075);
      hw_target_primary_parent_CCQE->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_target_primary_parent_CCQE.pdf");

      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_Tejin->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_ALL_primary_parent_Tejin.pdf");

      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Tejin_TrackerONLY->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_ALL_primary_parent_Tejin_TrackerONLY.pdf");

      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Recoil->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Recoil->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_Recoil->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_Recoil->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_Recoil->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_ALL_primary_parent_Recoil.pdf");

      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(1,"Other");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(2,"");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(3,"n");           
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(4,"p");           
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(5,"#pi^{0}");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(6,"#pi^{+}");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(7,"#pi^{-}");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(8,"#gamma");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(9,"e");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetBinLabel(10,"#mu");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetLabelSize(0.06);
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetTitle("Blob Primary Parent");
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_CCQE->GetXaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_CCQE->GetYaxis()->SetTitle("Blobs");
      hw_ALL_primary_parent_CCQE->GetYaxis()->SetTitleSize(0.045);
      hw_ALL_primary_parent_CCQE->GetYaxis()->SetTitleOffset(1.075);
      hw_ALL_primary_parent_CCQE->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_ALL_primary_parent_CCQE.pdf");

      //Drawing target vs. tracker on same plot instead of stacked.

      hw_target_primary_parent_CCQE->SetLineColor(kBlue);
      //hw_target_primary_parent_CCQE->Sumw2(kFALSE);

      double scaler=hw_target_primary_parent_CCQE->Integral();
      hw_target_primary_parent_CCQE->Scale(1.0/scaler);
      hw_target_primary_parent_CCQE->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_CCQE->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_CCQE->SetLineColor(kRed);
      //hw_tracker_primary_parent_CCQE->Sumw2(kFALSE);

      scaler=hw_tracker_primary_parent_CCQE->Integral();
      hw_tracker_primary_parent_CCQE->Scale(1.0/scaler);
      hw_tracker_primary_parent_CCQE->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_CCQE->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_Recoil->SetLineColor(kBlue);
      //hw_target_primary_parent_Recoil->Sumw2(kFALSE);

      scaler=hw_target_primary_parent_Recoil->Integral();
      hw_target_primary_parent_Recoil->Scale(1.0/scaler);
      hw_target_primary_parent_Recoil->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_Recoil->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_Recoil->SetLineColor(kRed);
      //hw_tracker_primary_parent_Recoil->Sumw2(kFALSE);

      scaler=hw_tracker_primary_parent_Recoil->Integral();
      hw_tracker_primary_parent_Recoil->Scale(1.0/scaler);
      hw_tracker_primary_parent_Recoil->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_Recoil->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_Tejin->SetLineColor(kBlue);
      //hw_target_primary_parent_Tejin->Sumw2(kFALSE);

      scaler=hw_target_primary_parent_Tejin->Integral();
      hw_target_primary_parent_Tejin->Scale(1.0/scaler);
      hw_target_primary_parent_Tejin->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_Tejin->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_Tejin->SetLineColor(kRed);
      //hw_tracker_primary_parent_Tejin->Sumw2(kFALSE);

      scaler=hw_tracker_primary_parent_Tejin->Integral();
      hw_tracker_primary_parent_Tejin->Scale(1.0/scaler);
      hw_tracker_primary_parent_Tejin->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_Tejin->GetYaxis()->SetRangeUser(0.0,0.5);

      hw_target_primary_parent_Tejin_TrackerONLY->SetLineColor(kBlue);
      //hw_target_primary_parent_Tejin_TrackerONLY->Sumw2(kFALSE);

      scaler=hw_target_primary_parent_Tejin_TrackerONLY->Integral();
      hw_target_primary_parent_Tejin_TrackerONLY->Scale(1.0/scaler);
      hw_target_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_target_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetRangeUser(0.0,0.5);
      hw_tracker_primary_parent_Tejin_TrackerONLY->SetLineColor(kRed);
      //hw_tracker_primary_parent_Tejin_TrackerONLY->Sumw2(kFALSE);

      scaler=hw_tracker_primary_parent_Tejin_TrackerONLY->Integral();
      hw_tracker_primary_parent_Tejin_TrackerONLY->Scale(1.0/scaler);
      hw_target_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetTitle("Fraction of Blobs");
      hw_tracker_primary_parent_Tejin_TrackerONLY->GetYaxis()->SetRangeUser(0.0,0.5);

      TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
      legend->SetHeader("Blob Location");
      legend->AddEntry(hw_tracker_primary_parent_Tejin_TrackerONLY,"Tracker Region","L");
      legend->AddEntry(hw_target_primary_parent_Tejin_TrackerONLY,"Target Region","L");

      hw_tracker_primary_parent_CCQE->SetStats(kFALSE);
      hw_target_primary_parent_CCQE->SetStats(kFALSE);
      hw_tracker_primary_parent_CCQE->Draw();
      hw_target_primary_parent_CCQE->Draw("same");
      legend->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_both_regions_primary_parent_CCQE.pdf");

      hw_tracker_primary_parent_Recoil->SetStats(kFALSE);
      hw_target_primary_parent_Recoil->SetStats(kFALSE);
      hw_tracker_primary_parent_Recoil->Draw();
      hw_target_primary_parent_Recoil->Draw("same");
      legend->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_both_regions_primary_parent_Recoil.pdf");

      hw_tracker_primary_parent_Tejin->SetStats(kFALSE);
      hw_target_primary_parent_Tejin->SetStats(kFALSE);
      hw_tracker_primary_parent_Tejin->Draw();
      hw_target_primary_parent_Tejin->Draw("same");
      legend->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_both_regions_primary_parent_Tejin.pdf");

      hw_tracker_primary_parent_Tejin_TrackerONLY->SetStats(kFALSE);
      hw_target_primary_parent_Tejin_TrackerONLY->SetStats(kFALSE);
      hw_tracker_primary_parent_Tejin_TrackerONLY->Draw();
      hw_target_primary_parent_Tejin_TrackerONLY->Draw("same");
      legend->Draw();
      c1->Print("/minerva/data/users/dlast/TargetNeutronsAna/JunkPlots/TestQEOnly/h_both_regions_primary_parent_Tejin_TrackerONLY.pdf");
									    }
