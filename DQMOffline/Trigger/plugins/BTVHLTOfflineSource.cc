/*
  BTVHLTOffline DQM code
*/
//
// Originally created by:  Anne-Catherine Le Bihan
//                         June 2015
//                         John Alison <johnalison@cmu.edu>
//                         June 2020
// Following the structure used in JetMetHLTOfflineSource

// system include files
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>
#include <cmath>
#include <iostream>

// user include files
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMOffline/Trigger/plugins/TriggerDQMBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TMath.h"
#include "TPRegexp.h"

class BTVHLTOfflineSource : public DQMEDAnalyzer {
public:
  explicit BTVHLTOfflineSource(const edm::ParameterSet&);
  ~BTVHLTOfflineSource() override;

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  std::vector<const reco::Track*> getOfflineBTagTracks(float hltJetEta,
                                                       float hltJetPhi,
                                                       edm::Handle<edm::View<reco::BaseTagInfo>> offlineIPTagHandle,
                                                       std::vector<float>& offlineIP3D,
                                                       std::vector<float>& offlineIP3DSig);

  typedef reco::TemplatedSecondaryVertexTagInfo<reco::CandIPTagInfo, reco::VertexCompositePtrCandidate> SVTagInfo;

  template <class Base>
  std::vector<const reco::Track*> getOnlineBTagTracks(float hltJetEta,
                                                      float hltJetPhi,
                                                      edm::Handle<std::vector<Base>> jetSVTagsColl,
                                                      std::vector<float>& onlineIP3D,
                                                      std::vector<float>& onlineIP3DSig);

  void bookHistograms(DQMStore::IBooker&, edm::Run const& run, edm::EventSetup const& c) override;

  void dqmBeginRun(edm::Run const& run, edm::EventSetup const& c) override;

  std::string dirname_;
  std::string processname_;
  bool verbose_;

  std::vector<std::pair<std::string, std::string>> custompathnamepairs_;

  edm::InputTag triggerSummaryLabel_;
  edm::InputTag triggerResultsLabel_;

  float turnon_threshold_loose_;
  float turnon_threshold_medium_;
  float turnon_threshold_tight_;
  float turnon_threshold_offline_loose_;
  float turnon_threshold_offline_medium_;
  float turnon_threshold_offline_tight_;

  edm::EDGetTokenT<reco::JetTagCollection> offlineDiscrTokenb_;
  edm::EDGetTokenT<edm::View<reco::BaseTagInfo>> offlineIPToken_;

  edm::EDGetTokenT<std::vector<reco::Vertex>> hltPFPVToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> offlinePVToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsFUToken;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryToken;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryFUToken;

  edm::EDGetTokenT<std::vector<reco::DeepBoostedJetTagInfo>> DeepBoostedJetTagInfosTokenPf_;

  edm::EDGetTokenT<std::vector<SVTagInfo>> SVTagInfosTokenPf_;

  edm::EDGetTokenT<reco::JetTagCollection> pfTagsToken_;
  edm::Handle<reco::JetTagCollection> pfTags;

  float minDecayLength_;
  float maxDecayLength_;
  float minJetDistance_;
  float maxJetDistance_;
  float dRTrackMatch_;

  HLTConfigProvider hltConfig_;

  class PathInfo : public TriggerDQMBase {
  public:
    PathInfo()
        : prescaleUsed_(-1),
          pathName_("unset"),
          filterName_("unset"),
          processName_("unset"),
          objectType_(-1),
          triggerType_("unset") {}

    ~PathInfo() override = default;

    PathInfo(const int prescaleUsed,
             const std::string& pathName,
             const std::string& filterName,
             const std::string& processName,
             const int type,
             const std::string& triggerType)
        : prescaleUsed_(prescaleUsed),
          pathName_(pathName),
          filterName_(filterName),
          processName_(processName),
          objectType_(type),
          triggerType_(triggerType) {}

    const std::string getLabel() const { return filterName_; }
    void setLabel(std::string labelName) { filterName_ = std::move(labelName); }
    const std::string getPath() const { return pathName_; }
    const int getprescaleUsed() const { return prescaleUsed_; }
    const std::string getProcess() const { return processName_; }
    const int getObjectType() const { return objectType_; }
    const std::string getTriggerType() const { return triggerType_; }
    const edm::InputTag getTag() const { return edm::InputTag(filterName_, "", processName_); }
    const bool operator==(const std::string& v) const { return v == pathName_; }

    MonitorElement* Discr = nullptr;
    MonitorElement* Pt = nullptr;
    MonitorElement* Eta = nullptr;
    MonitorElement* Discr_HLTvsRECO = nullptr;
    MonitorElement* Discr_HLTMinusRECO = nullptr;
    ObjME Discr_turnon_loose;
    ObjME Discr_turnon_medium;
    ObjME Discr_turnon_tight;
    ObjME Pt_turnon_loose;
    ObjME Pt_turnon_medium;
    ObjME Pt_turnon_tight;
    ObjME Eta_turnon_loose;
    ObjME Eta_turnon_medium;
    ObjME Eta_turnon_tight;
    ObjME Phi_turnon_loose;
    ObjME Phi_turnon_medium;
    ObjME Phi_turnon_tight;
    MonitorElement* PVz = nullptr;
    MonitorElement* fastPVz = nullptr;
    MonitorElement* PVz_HLTMinusRECO = nullptr;
    MonitorElement* fastPVz_HLTMinusRECO = nullptr;
    MonitorElement* n_vtx = nullptr;
    MonitorElement* vtx_mass = nullptr;
    MonitorElement* n_vtx_trks = nullptr;
    MonitorElement* n_sel_tracks = nullptr;
    MonitorElement* h_3d_ip_distance = nullptr;
    MonitorElement* h_3d_ip_error = nullptr;
    MonitorElement* h_3d_ip_sig = nullptr;

    MonitorElement* h_pfcand_puppiw = nullptr;
    MonitorElement* h_pfcand_hcalFrac = nullptr;
    MonitorElement* h_pfcand_VTX_ass = nullptr;
    MonitorElement* h_pfcand_lostInnerHits = nullptr;
    MonitorElement* h_pfcand_quality = nullptr;
    MonitorElement* h_pfcand_charge = nullptr;
    MonitorElement* h_pfcand_isEl = nullptr;
    MonitorElement* h_pfcand_isMu = nullptr;
    MonitorElement* h_pfcand_isChargedHad = nullptr;
    MonitorElement* h_pfcand_isGamma = nullptr;
    MonitorElement* h_pfcand_isNeutralHad = nullptr;
    MonitorElement* h_pfcand_phirel = nullptr;
    MonitorElement* h_pfcand_etarel = nullptr;
    MonitorElement* h_pfcand_deltaR = nullptr;
    MonitorElement* h_pfcand_abseta = nullptr;
    MonitorElement* h_pfcand_ptrel_log = nullptr;
    MonitorElement* h_pfcand_erel_log = nullptr;
    MonitorElement* h_pfcand_pt_log = nullptr;
    MonitorElement* h_pfcand_drminsv = nullptr;
    MonitorElement* h_pfcand_drsubjet1 = nullptr;
    MonitorElement* h_pfcand_drsubjet2 = nullptr;
    MonitorElement* h_pfcand_normchi2 = nullptr;
    MonitorElement* h_pfcand_dz = nullptr;
    MonitorElement* h_pfcand_dzsig = nullptr;
    MonitorElement* h_pfcand_dxy = nullptr;
    MonitorElement* h_pfcand_dxysig = nullptr;
    MonitorElement* h_pfcand_dptdpt = nullptr;
    MonitorElement* h_pfcand_detadeta = nullptr;
    MonitorElement* h_pfcand_dphidphi = nullptr;
    MonitorElement* h_pfcand_dxydxy = nullptr;
    MonitorElement* h_pfcand_dzdz = nullptr;
    MonitorElement* h_pfcand_dxydz = nullptr;
    MonitorElement* h_pfcand_dphidxy = nullptr;
    MonitorElement* h_pfcand_dlambdadz = nullptr;
    MonitorElement* h_pfcand_btagEtaRel = nullptr;
    MonitorElement* h_pfcand_btagPtRatio = nullptr;
    MonitorElement* h_pfcand_btagPParRatio = nullptr;
    MonitorElement* h_pfcand_btagSip2dVal = nullptr;
    MonitorElement* h_pfcand_btagSip2dSig = nullptr;
    MonitorElement* h_pfcand_btagSip3dVal = nullptr;
    MonitorElement* h_pfcand_btagSip3dSig = nullptr;
    MonitorElement* h_pfcand_btagJetDistVal = nullptr;
    MonitorElement* h_pfcand_mask = nullptr;
    MonitorElement* h_pfcand_pt_log_nopuppi = nullptr;
    MonitorElement* h_pfcand_e_log_nopuppi = nullptr;
    MonitorElement* h_pfcand_ptrel = nullptr;
    MonitorElement* h_pfcand_erel = nullptr;

    MonitorElement* h_jet_pfcand_pt_log = nullptr;
    MonitorElement* h_jet_pfcand_energy_log = nullptr;
    MonitorElement* h_jet_pfcand_deta = nullptr;
    MonitorElement* h_jet_pfcand_dphi = nullptr;
    MonitorElement* h_jet_pfcand_eta = nullptr;
    MonitorElement* h_jet_pfcand_charge = nullptr;
    MonitorElement* h_jet_pfcand_frompv = nullptr;
    MonitorElement* h_jet_pfcand_nlostinnerhits = nullptr;
    MonitorElement* h_jet_pfcand_track_chi2 = nullptr;
    MonitorElement* h_jet_pfcand_track_qual = nullptr;
    MonitorElement* h_jet_pfcand_dz = nullptr;
    MonitorElement* h_jet_pfcand_dzsig = nullptr;
    MonitorElement* h_jet_pfcand_dxy = nullptr;
    MonitorElement* h_jet_pfcand_dxysig = nullptr;
    MonitorElement* h_jet_pfcand_etarel = nullptr;
    MonitorElement* h_jet_pfcand_pperp_ratio = nullptr;
    MonitorElement* h_jet_pfcand_ppara_ratio = nullptr;
    MonitorElement* h_jet_pfcand_trackjet_d3d = nullptr;
    MonitorElement* h_jet_pfcand_trackjet_d3dsig = nullptr;
    MonitorElement* h_jet_pfcand_trackjet_dist = nullptr;
    MonitorElement* h_jet_pfcand_nhits = nullptr;
    MonitorElement* h_jet_pfcand_npixhits = nullptr;
    MonitorElement* h_jet_pfcand_nstriphits = nullptr;
    MonitorElement* h_jet_pfcand_trackjet_decayL = nullptr;
    MonitorElement* h_jet_pfcand_puppiw = nullptr;

    ObjME OnlineTrkEff_Pt;
    ObjME OnlineTrkEff_Eta;
    ObjME OnlineTrkEff_3d_ip_distance;
    ObjME OnlineTrkEff_3d_ip_sig;
    ObjME OnlineTrkFake_Pt;
    ObjME OnlineTrkFake_Eta;
    ObjME OnlineTrkFake_3d_ip_distance;
    ObjME OnlineTrkFake_3d_ip_sig;
    // MonitorElement*  n_pixel_hits_;
    // MonitorElement*  n_total_hits_;

  private:
    int prescaleUsed_;
    std::string pathName_;
    std::string filterName_;
    std::string processName_;
    int objectType_;
    std::string triggerType_;
  };

  class PathInfoCollection : public std::vector<PathInfo> {
  public:
    PathInfoCollection() : std::vector<PathInfo>() {};
    std::vector<PathInfo>::iterator find(const std::string& pathName) { return std::find(begin(), end(), pathName); }
  };

  PathInfoCollection hltPathsAll_;
};

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;

BTVHLTOfflineSource::BTVHLTOfflineSource(const edm::ParameterSet& iConfig)
    : dirname_(iConfig.getUntrackedParameter("dirname", std::string("HLT/BTV/"))),
      processname_(iConfig.getParameter<std::string>("processname")),
      verbose_(iConfig.getUntrackedParameter<bool>("verbose", false)),
      triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummaryLabel")),
      triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResultsLabel")),
      turnon_threshold_loose_(iConfig.getParameter<double>("turnon_threshold_loose")),
      turnon_threshold_medium_(iConfig.getParameter<double>("turnon_threshold_medium")),
      turnon_threshold_tight_(iConfig.getParameter<double>("turnon_threshold_tight")),
      turnon_threshold_offline_loose_(iConfig.getParameter<double>("turnon_threshold_offline_loose")),
      turnon_threshold_offline_medium_(iConfig.getParameter<double>("turnon_threshold_offline_medium")),
      turnon_threshold_offline_tight_(iConfig.getParameter<double>("turnon_threshold_offline_tight")),
      offlineDiscrTokenb_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("offlineDiscrLabelb"))),
      offlineIPToken_(consumes<View<BaseTagInfo>>(iConfig.getParameter<edm::InputTag>("offlineIPLabel"))),

      hltPFPVToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("hltPFPVLabel"))),
      offlinePVToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlinePVLabel"))),
      triggerResultsToken(consumes<edm::TriggerResults>(triggerResultsLabel_)),
      triggerResultsFUToken(consumes<edm::TriggerResults>(
          edm::InputTag(triggerResultsLabel_.label(), triggerResultsLabel_.instance(), std::string("FU")))),
      triggerSummaryToken(consumes<trigger::TriggerEvent>(triggerSummaryLabel_)),
      triggerSummaryFUToken(consumes<trigger::TriggerEvent>(
          edm::InputTag(triggerSummaryLabel_.label(), triggerSummaryLabel_.instance(), std::string("FU")))),
      DeepBoostedJetTagInfosTokenPf_(
          consumes<vector<reco::DeepBoostedJetTagInfo>>(edm::InputTag("hltParticleNetJetTagInfos"))),
      SVTagInfosTokenPf_(consumes<std::vector<SVTagInfo>>(edm::InputTag("hltDeepSecondaryVertexTagInfosPF"))),
      pfTagsToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("onlineDiscrLabelPF"))),
      minDecayLength_(iConfig.getParameter<double>("minDecayLength")),
      maxDecayLength_(iConfig.getParameter<double>("maxDecayLength")),
      minJetDistance_(iConfig.getParameter<double>("minJetDistance")),
      maxJetDistance_(iConfig.getParameter<double>("maxJetDistance")),
      dRTrackMatch_(iConfig.getParameter<double>("dRTrackMatch")) {
  std::vector<edm::ParameterSet> paths = iConfig.getParameter<std::vector<edm::ParameterSet>>("pathPairs");
  for (const auto& path : paths) {
    custompathnamepairs_.push_back(
        make_pair(path.getParameter<std::string>("pathname"), path.getParameter<std::string>("pathtype")));
  }
}

BTVHLTOfflineSource::~BTVHLTOfflineSource() = default;

void BTVHLTOfflineSource::dqmBeginRun(const edm::Run& run, const edm::EventSetup& c) {
  bool changed = true;
  if (!hltConfig_.init(run, c, processname_, changed)) {
    LogDebug("BTVHLTOfflineSource") << "HLTConfigProvider failed to initialize.";
  }

  for (unsigned int idx = 0; idx != hltConfig_.size(); ++idx) {
    const auto& pathname = hltConfig_.triggerName(idx);

    for (const auto& custompathnamepair : custompathnamepairs_) {
      if (pathname.find(custompathnamepair.first) != std::string::npos) {
        hltPathsAll_.push_back(PathInfo(1, pathname, "dummy", processname_, 0, custompathnamepair.second));
      }
    }
  }
}

void BTVHLTOfflineSource::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults);
  if (!triggerResults.isValid()) {
    iEvent.getByToken(triggerResultsFUToken, triggerResults);
    if (!triggerResults.isValid()) {
      edm::LogInfo("BTVHLTOfflineSource") << "TriggerResults not found, skipping event";
      return;
    }
  }

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);

  edm::Handle<trigger::TriggerEvent> triggerObj;
  iEvent.getByToken(triggerSummaryToken, triggerObj);
  if (!triggerObj.isValid()) {
    iEvent.getByToken(triggerSummaryFUToken, triggerObj);
    if (!triggerObj.isValid()) {
      edm::LogInfo("BTVHLTOfflineSource") << "TriggerEvent not found, skipping event";
      return;
    }
  }

  edm::Handle<reco::JetTagCollection> pfTags;
  iEvent.getByToken(pfTagsToken_, pfTags);

  Handle<reco::VertexCollection> VertexHandler;

  Handle<reco::JetTagCollection> offlineJetTagHandlerb;
  iEvent.getByToken(offlineDiscrTokenb_, offlineJetTagHandlerb);

  Handle<View<BaseTagInfo>> offlineIPTagHandle;
  iEvent.getByToken(offlineIPToken_, offlineIPTagHandle);

  Handle<reco::VertexCollection> offlineVertexHandler;
  iEvent.getByToken(offlinePVToken_, offlineVertexHandler);

  if (verbose_ && iEvent.id().event() % 10000 == 0)
    cout << "Run = " << iEvent.id().run() << ", LS = " << iEvent.luminosityBlock()
         << ", Event = " << iEvent.id().event() << endl;

  if (!triggerResults.isValid())
    return;

  edm::Handle<std::vector<SVTagInfo>> jetSVTagsCollPF;

  for (auto& v : hltPathsAll_) {
    unsigned index = triggerNames.triggerIndex(v.getPath());
    if (!(index < triggerNames.size())) {
      edm::LogInfo("BTVHLTOfflineSource") << "Path " << v.getPath() << " not in menu, skipping event";
      continue;
    }

    if (!triggerResults->accept(index)) {
      edm::LogInfo("BTVHLTOfflineSource") << "Path " << v.getPath() << " not accepted, skipping event";
      continue;
    }

    iEvent.getByToken(SVTagInfosTokenPf_, jetSVTagsCollPF);

    // PF btagging
    if (v.getTriggerType() == "PF" && pfTags.isValid()) {
      const auto& iter = pfTags->begin();

      float Discr_online = iter->second;
      if (Discr_online < 0)
        Discr_online = -0.05;

      v.Discr->Fill(Discr_online);
      v.Pt->Fill(iter->first->pt());
      v.Eta->Fill(iter->first->eta());

      if (offlineJetTagHandlerb.isValid()) {
        for (auto const& iterOffb : *offlineJetTagHandlerb) {
          float DR = reco::deltaR(iterOffb.first->eta(), iterOffb.first->phi(), iter->first->eta(), iter->first->phi());
          if (DR < 0.3) {
            float Discr_offline = iterOffb.second;
            float Pt_offline = iterOffb.first->pt();
            float Eta_offline = iterOffb.first->eta();
            float Phi_offline = iterOffb.first->phi();

            if (Discr_offline < 0)
              Discr_offline = -0.05;
            v.Discr_HLTvsRECO->Fill(Discr_online, Discr_offline);
            v.Discr_HLTMinusRECO->Fill(Discr_online - Discr_offline);

            v.Discr_turnon_loose.denominator->Fill(Discr_offline);
            v.Discr_turnon_medium.denominator->Fill(Discr_offline);
            v.Discr_turnon_tight.denominator->Fill(Discr_offline);

            if (Discr_online > turnon_threshold_loose_)
              v.Discr_turnon_loose.numerator->Fill(Discr_offline);
            if (Discr_online > turnon_threshold_medium_)
              v.Discr_turnon_medium.numerator->Fill(Discr_offline);
            if (Discr_online > turnon_threshold_tight_)
              v.Discr_turnon_tight.numerator->Fill(Discr_offline);

            if (Discr_offline > turnon_threshold_offline_loose_) {
              v.Pt_turnon_loose.denominator->Fill(Pt_offline);
              v.Eta_turnon_loose.denominator->Fill(Eta_offline);
              v.Phi_turnon_loose.denominator->Fill(Phi_offline);
              if (Discr_online > turnon_threshold_loose_) {
                v.Pt_turnon_loose.numerator->Fill(Pt_offline);
                v.Eta_turnon_loose.numerator->Fill(Eta_offline);
                v.Phi_turnon_loose.numerator->Fill(Phi_offline);
              }
            }
            if (Discr_offline > turnon_threshold_offline_medium_) {
              v.Pt_turnon_medium.denominator->Fill(Pt_offline);
              v.Eta_turnon_medium.denominator->Fill(Eta_offline);
              v.Phi_turnon_medium.denominator->Fill(Phi_offline);
              if (Discr_online > turnon_threshold_medium_) {
                v.Pt_turnon_medium.numerator->Fill(Pt_offline);
                v.Eta_turnon_medium.numerator->Fill(Eta_offline);
                v.Phi_turnon_medium.numerator->Fill(Phi_offline);
              }
            }
            if (Discr_offline > turnon_threshold_offline_tight_) {
              v.Pt_turnon_tight.denominator->Fill(Pt_offline);
              v.Eta_turnon_tight.denominator->Fill(Eta_offline);
              v.Phi_turnon_tight.denominator->Fill(Phi_offline);
              if (Discr_online > turnon_threshold_tight_) {
                v.Pt_turnon_tight.numerator->Fill(Pt_offline);
                v.Eta_turnon_tight.numerator->Fill(Eta_offline);
                v.Phi_turnon_tight.numerator->Fill(Phi_offline);
              }
            }

            break;
          }
        }
      }  ///offline

      bool pfSVTagCollValid = (v.getTriggerType() == "PF" && jetSVTagsCollPF.isValid());
      if (offlineIPTagHandle.isValid() && pfSVTagCollValid) {
        std::vector<float> offlineIP3D;
        std::vector<float> offlineIP3DSig;
        std::vector<const reco::Track*> offlineTracks = getOfflineBTagTracks(
            iter->first->eta(), iter->first->phi(), offlineIPTagHandle, offlineIP3D, offlineIP3DSig);
        std::vector<const reco::Track*> onlineTracks;
        std::vector<float> onlineIP3D;
        std::vector<float> onlineIP3DSig;
        if (pfSVTagCollValid)
          onlineTracks = getOnlineBTagTracks<SVTagInfo>(
              iter->first->eta(), iter->first->phi(), jetSVTagsCollPF, onlineIP3D, onlineIP3DSig);

        for (unsigned int iOffTrk = 0; iOffTrk < offlineTracks.size(); ++iOffTrk) {
          const reco::Track* offTrk = offlineTracks.at(iOffTrk);
          bool hasMatch = false;
          float offTrkEta = offTrk->eta();
          float offTrkPhi = offTrk->phi();

          for (const reco::Track* onTrk : onlineTracks) {
            float DR = reco::deltaR(offTrkEta, offTrkPhi, onTrk->eta(), onTrk->phi());
            if (DR < dRTrackMatch_) {
              hasMatch = true;
            }
          }

          float offTrkPt = offTrk->pt();
          v.OnlineTrkEff_Pt.denominator->Fill(offTrkPt);
          if (hasMatch)
            v.OnlineTrkEff_Pt.numerator->Fill(offTrkPt);

          v.OnlineTrkEff_Eta.denominator->Fill(offTrkEta);
          if (hasMatch)
            v.OnlineTrkEff_Eta.numerator->Fill(offTrkEta);

          v.OnlineTrkEff_3d_ip_distance.denominator->Fill(offlineIP3D.at(iOffTrk));
          if (hasMatch)
            v.OnlineTrkEff_3d_ip_distance.numerator->Fill(offlineIP3D.at(iOffTrk));

          v.OnlineTrkEff_3d_ip_sig.denominator->Fill(offlineIP3DSig.at(iOffTrk));
          if (hasMatch)
            v.OnlineTrkEff_3d_ip_sig.numerator->Fill(offlineIP3DSig.at(iOffTrk));
        }

        for (unsigned int iOnTrk = 0; iOnTrk < onlineTracks.size(); ++iOnTrk) {
          const reco::Track* onTrk = onlineTracks.at(iOnTrk);
          bool hasMatch = false;
          float onTrkEta = onTrk->eta();
          float onTrkPhi = onTrk->phi();

          for (const reco::Track* offTrk : offlineTracks) {
            float DR = reco::deltaR(onTrkEta, onTrkPhi, offTrk->eta(), offTrk->phi());
            if (DR < dRTrackMatch_) {
              hasMatch = true;
            }
          }

          float onTrkPt = onTrk->pt();
          v.OnlineTrkFake_Pt.denominator->Fill(onTrkPt);
          if (!hasMatch)
            v.OnlineTrkFake_Pt.numerator->Fill(onTrkPt);

          v.OnlineTrkFake_Eta.denominator->Fill(onTrkEta);
          if (!hasMatch)
            v.OnlineTrkFake_Eta.numerator->Fill(onTrkEta);

          v.OnlineTrkFake_3d_ip_distance.denominator->Fill(onlineIP3D.at(iOnTrk));
          if (!hasMatch)
            v.OnlineTrkFake_3d_ip_distance.numerator->Fill(onlineIP3D.at(iOnTrk));

          v.OnlineTrkFake_3d_ip_sig.denominator->Fill(onlineIP3DSig.at(iOnTrk));
          if (!hasMatch)
            v.OnlineTrkFake_3d_ip_sig.numerator->Fill(onlineIP3DSig.at(iOnTrk));
        }
      }

      iEvent.getByToken(hltPFPVToken_, VertexHandler);
      if (VertexHandler.isValid()) {
        v.PVz->Fill(VertexHandler->begin()->z());
        if (offlineVertexHandler.isValid()) {
          v.PVz_HLTMinusRECO->Fill(VertexHandler->begin()->z() - offlineVertexHandler->begin()->z());
        }
      }
    }  // PFTagsValid

    // additional plots from tag info collections
    /////////////////////////////////////////////

    edm::Handle<std::vector<reco::DeepBoostedJetTagInfo>> DeepBoostedJetTagInfosPf;
    iEvent.getByToken(DeepBoostedJetTagInfosTokenPf_, DeepBoostedJetTagInfosPf);

    //    edm::Handle<std::vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<edm::RefVector<std::vector<reco::Track>, reco::Track, edm::refhelper::FindUsingAdvance<std::vector<reco::Track>, reco::Track> >, reco::JTATagInfo>, reco::Vertex> > > caloTagInfos;
    //    iEvent.getByToken(caloTagInfosToken_, caloTagInfos);

    //    edm::Handle<std::vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<edm::RefVector<std::vector<reco::Track>, reco::Track, edm::refhelper::FindUsingAdvance<std::vector<reco::Track>, reco::Track> >, reco::JTATagInfo>, reco::Vertex> > > pfTagInfos;
    //    iEvent.getByToken(pfTagInfosToken_, pfTagInfos);

    // first try to get info from DeepBoostedJetTagInfos ...
    if (v.getTriggerType() == "PF" && DeepBoostedJetTagInfosPf.isValid()) {
      const auto& DeepBoostedTagInfoCollection = DeepBoostedJetTagInfosPf;
      for (const auto& DeepBoostedTagInfo : *DeepBoostedTagInfoCollection) {
        const auto& features = DeepBoostedTagInfo.features();
        
        // Regular particle features
        const auto& pfcand_puppiw_values = features.get("pfcand_puppiw");
        for (const auto& val : pfcand_puppiw_values) {
          v.h_pfcand_puppiw->Fill(val);
        }
        
        const auto& pfcand_hcalFrac_values = features.get("pfcand_hcalFrac");
        for (const auto& val : pfcand_hcalFrac_values) {
          v.h_pfcand_hcalFrac->Fill(val);
        }
        
        const auto& pfcand_VTX_ass_values = features.get("pfcand_VTX_ass");
        for (const auto& val : pfcand_VTX_ass_values) {
          v.h_pfcand_VTX_ass->Fill(val);
        }
        
        const auto& pfcand_lostInnerHits_values = features.get("pfcand_lostInnerHits");
        for (const auto& val : pfcand_lostInnerHits_values) {
          v.h_pfcand_lostInnerHits->Fill(val);
        }
        
        const auto& pfcand_quality_values = features.get("pfcand_quality");
        for (const auto& val : pfcand_quality_values) {
          v.h_pfcand_quality->Fill(val);
        }
        
        const auto& pfcand_charge_values = features.get("pfcand_charge");
        for (const auto& val : pfcand_charge_values) {
          v.h_pfcand_charge->Fill(val);
        }
        
        const auto& pfcand_isEl_values = features.get("pfcand_isEl");
        for (const auto& val : pfcand_isEl_values) {
          v.h_pfcand_isEl->Fill(val);
        }
        
        const auto& pfcand_isMu_values = features.get("pfcand_isMu");
        for (const auto& val : pfcand_isMu_values) {
          v.h_pfcand_isMu->Fill(val);
        }
        
        const auto& pfcand_isChargedHad_values = features.get("pfcand_isChargedHad");
        for (const auto& val : pfcand_isChargedHad_values) {
          v.h_pfcand_isChargedHad->Fill(val);
        }
        
        const auto& pfcand_isGamma_values = features.get("pfcand_isGamma");
        for (const auto& val : pfcand_isGamma_values) {
          v.h_pfcand_isGamma->Fill(val);
        }
        
        const auto& pfcand_isNeutralHad_values = features.get("pfcand_isNeutralHad");
        for (const auto& val : pfcand_isNeutralHad_values) {
          v.h_pfcand_isNeutralHad->Fill(val);
        }
        
        const auto& pfcand_phirel_values = features.get("pfcand_phirel");
        for (const auto& val : pfcand_phirel_values) {
          v.h_pfcand_phirel->Fill(val);
        }
        
        const auto& pfcand_etarel_values = features.get("pfcand_etarel");
        for (const auto& val : pfcand_etarel_values) {
          v.h_pfcand_etarel->Fill(val);
        }
        
        const auto& pfcand_deltaR_values = features.get("pfcand_deltaR");
        for (const auto& val : pfcand_deltaR_values) {
          v.h_pfcand_deltaR->Fill(val);
        }
        
        const auto& pfcand_abseta_values = features.get("pfcand_abseta");
        for (const auto& val : pfcand_abseta_values) {
          v.h_pfcand_abseta->Fill(val);
        }
        
        const auto& pfcand_ptrel_log_values = features.get("pfcand_ptrel_log");
        for (const auto& val : pfcand_ptrel_log_values) {
          v.h_pfcand_ptrel_log->Fill(val);
        }
        
        const auto& pfcand_erel_log_values = features.get("pfcand_erel_log");
        for (const auto& val : pfcand_erel_log_values) {
          v.h_pfcand_erel_log->Fill(val);
        }
        
        const auto& pfcand_pt_log_values = features.get("pfcand_pt_log");
        for (const auto& val : pfcand_pt_log_values) {
          v.h_pfcand_pt_log->Fill(val);
        }
        
        const auto& pfcand_drminsv_values = features.get("pfcand_drminsv");
        for (const auto& val : pfcand_drminsv_values) {
          v.h_pfcand_drminsv->Fill(val);
        }
        
        const auto& pfcand_drsubjet1_values = features.get("pfcand_drsubjet1");
        for (const auto& val : pfcand_drsubjet1_values) {
          v.h_pfcand_drsubjet1->Fill(val);
        }
        
        const auto& pfcand_drsubjet2_values = features.get("pfcand_drsubjet2");
        for (const auto& val : pfcand_drsubjet2_values) {
          v.h_pfcand_drsubjet2->Fill(val);
        }
        
        const auto& pfcand_normchi2_values = features.get("pfcand_normchi2");
        for (const auto& val : pfcand_normchi2_values) {
          v.h_pfcand_normchi2->Fill(val);
        }
        
        const auto& pfcand_dz_values = features.get("pfcand_dz");
        for (const auto& val : pfcand_dz_values) {
          v.h_pfcand_dz->Fill(val);
        }
        
        const auto& pfcand_dzsig_values = features.get("pfcand_dzsig");
        for (const auto& val : pfcand_dzsig_values) {
          v.h_pfcand_dzsig->Fill(val);
        }
        
        const auto& pfcand_dxy_values = features.get("pfcand_dxy");
        for (const auto& val : pfcand_dxy_values) {
          v.h_pfcand_dxy->Fill(val);
        }
        
        const auto& pfcand_dxysig_values = features.get("pfcand_dxysig");
        for (const auto& val : pfcand_dxysig_values) {
          v.h_pfcand_dxysig->Fill(val);
        }
        
        const auto& pfcand_dptdpt_values = features.get("pfcand_dptdpt");
        for (const auto& val : pfcand_dptdpt_values) {
          v.h_pfcand_dptdpt->Fill(val);
        }
        
        const auto& pfcand_detadeta_values = features.get("pfcand_detadeta");
        for (const auto& val : pfcand_detadeta_values) {
          v.h_pfcand_detadeta->Fill(val);
        }
        
        const auto& pfcand_dphidphi_values = features.get("pfcand_dphidphi");
        for (const auto& val : pfcand_dphidphi_values) {
          v.h_pfcand_dphidphi->Fill(val);
        }
        
        const auto& pfcand_dxydxy_values = features.get("pfcand_dxydxy");
        for (const auto& val : pfcand_dxydxy_values) {
          v.h_pfcand_dxydxy->Fill(val);
        }
        
        const auto& pfcand_dzdz_values = features.get("pfcand_dzdz");
        for (const auto& val : pfcand_dzdz_values) {
          v.h_pfcand_dzdz->Fill(val);
        }
        
        const auto& pfcand_dxydz_values = features.get("pfcand_dxydz");
        for (const auto& val : pfcand_dxydz_values) {
          v.h_pfcand_dxydz->Fill(val);
        }
        
        const auto& pfcand_dphidxy_values = features.get("pfcand_dphidxy");
        for (const auto& val : pfcand_dphidxy_values) {
          v.h_pfcand_dphidxy->Fill(val);
        }
        
        const auto& pfcand_dlambdadz_values = features.get("pfcand_dlambdadz");
        for (const auto& val : pfcand_dlambdadz_values) {
          v.h_pfcand_dlambdadz->Fill(val);
        }
        
        const auto& pfcand_btagEtaRel_values = features.get("pfcand_btagEtaRel");
        for (const auto& val : pfcand_btagEtaRel_values) {
          v.h_pfcand_btagEtaRel->Fill(val);
        }
        
        const auto& pfcand_btagPtRatio_values = features.get("pfcand_btagPtRatio");
        for (const auto& val : pfcand_btagPtRatio_values) {
          v.h_pfcand_btagPtRatio->Fill(val);
        }
        
        const auto& pfcand_btagPParRatio_values = features.get("pfcand_btagPParRatio");
        for (const auto& val : pfcand_btagPParRatio_values) {
          v.h_pfcand_btagPParRatio->Fill(val);
        }
        
        const auto& pfcand_btagSip2dVal_values = features.get("pfcand_btagSip2dVal");
        for (const auto& val : pfcand_btagSip2dVal_values) {
          v.h_pfcand_btagSip2dVal->Fill(val);
        }
        
        const auto& pfcand_btagSip2dSig_values = features.get("pfcand_btagSip2dSig");
        for (const auto& val : pfcand_btagSip2dSig_values) {
          v.h_pfcand_btagSip2dSig->Fill(val);
        }
        
        const auto& pfcand_btagSip3dVal_values = features.get("pfcand_btagSip3dVal");
        for (const auto& val : pfcand_btagSip3dVal_values) {
          v.h_pfcand_btagSip3dVal->Fill(val);
        }
        
        const auto& pfcand_btagSip3dSig_values = features.get("pfcand_btagSip3dSig");
        for (const auto& val : pfcand_btagSip3dSig_values) {
          v.h_pfcand_btagSip3dSig->Fill(val);
        }
        
        const auto& pfcand_btagJetDistVal_values = features.get("pfcand_btagJetDistVal");
        for (const auto& val : pfcand_btagJetDistVal_values) {
          v.h_pfcand_btagJetDistVal->Fill(val);
        }
        
        const auto& pfcand_mask_values = features.get("pfcand_mask");
        for (const auto& val : pfcand_mask_values) {
          v.h_pfcand_mask->Fill(val);
        }
        
        const auto& pfcand_pt_log_nopuppi_values = features.get("pfcand_pt_log_nopuppi");
        for (const auto& val : pfcand_pt_log_nopuppi_values) {
          v.h_pfcand_pt_log_nopuppi->Fill(val);
        }
        
        const auto& pfcand_e_log_nopuppi_values = features.get("pfcand_e_log_nopuppi");
        for (const auto& val : pfcand_e_log_nopuppi_values) {
          v.h_pfcand_e_log_nopuppi->Fill(val);
        }
        
        const auto& pfcand_ptrel_values = features.get("pfcand_ptrel");
        for (const auto& val : pfcand_ptrel_values) {
          v.h_pfcand_ptrel->Fill(val);
        }
        
        const auto& pfcand_erel_values = features.get("pfcand_erel");
        for (const auto& val : pfcand_erel_values) {
          v.h_pfcand_erel->Fill(val);
        }
        
        // HLT particle features
        const auto& jet_pfcand_pt_log_values = features.get("jet_pfcand_pt_log");
        for (const auto& val : jet_pfcand_pt_log_values) {
          v.h_jet_pfcand_pt_log->Fill(val);
        }
        
        const auto& jet_pfcand_energy_log_values = features.get("jet_pfcand_energy_log");
        for (const auto& val : jet_pfcand_energy_log_values) {
          v.h_jet_pfcand_energy_log->Fill(val);
        }
        
        const auto& jet_pfcand_deta_values = features.get("jet_pfcand_deta");
        for (const auto& val : jet_pfcand_deta_values) {
          v.h_jet_pfcand_deta->Fill(val);
        }
        
        const auto& jet_pfcand_dphi_values = features.get("jet_pfcand_dphi");
        for (const auto& val : jet_pfcand_dphi_values) {
          v.h_jet_pfcand_dphi->Fill(val);
        }
        
        const auto& jet_pfcand_eta_values = features.get("jet_pfcand_eta");
        for (const auto& val : jet_pfcand_eta_values) {
          v.h_jet_pfcand_eta->Fill(val);
        }
        
        const auto& jet_pfcand_charge_values = features.get("jet_pfcand_charge");
        for (const auto& val : jet_pfcand_charge_values) {
          v.h_jet_pfcand_charge->Fill(val);
        }
        
        const auto& jet_pfcand_frompv_values = features.get("jet_pfcand_frompv");
        for (const auto& val : jet_pfcand_frompv_values) {
          v.h_jet_pfcand_frompv->Fill(val);
        }
        
        const auto& jet_pfcand_nlostinnerhits_values = features.get("jet_pfcand_nlostinnerhits");
        for (const auto& val : jet_pfcand_nlostinnerhits_values) {
          v.h_jet_pfcand_nlostinnerhits->Fill(val);
        }
        
        const auto& jet_pfcand_track_chi2_values = features.get("jet_pfcand_track_chi2");
        for (const auto& val : jet_pfcand_track_chi2_values) {
          v.h_jet_pfcand_track_chi2->Fill(val);
        }
        
        const auto& jet_pfcand_track_qual_values = features.get("jet_pfcand_track_qual");
        for (const auto& val : jet_pfcand_track_qual_values) {
          v.h_jet_pfcand_track_qual->Fill(val);
        }
        
        const auto& jet_pfcand_dz_values = features.get("jet_pfcand_dz");
        for (const auto& val : jet_pfcand_dz_values) {
          v.h_jet_pfcand_dz->Fill(val);
        }
        
        const auto& jet_pfcand_dzsig_values = features.get("jet_pfcand_dzsig");
        for (const auto& val : jet_pfcand_dzsig_values) {
          v.h_jet_pfcand_dzsig->Fill(val);
        }
        
        const auto& jet_pfcand_dxy_values = features.get("jet_pfcand_dxy");
        for (const auto& val : jet_pfcand_dxy_values) {
          v.h_jet_pfcand_dxy->Fill(val);
        }
        
        const auto& jet_pfcand_dxysig_values = features.get("jet_pfcand_dxysig");
        for (const auto& val : jet_pfcand_dxysig_values) {
          v.h_jet_pfcand_dxysig->Fill(val);
        }
        
        const auto& jet_pfcand_etarel_values = features.get("jet_pfcand_etarel");
        for (const auto& val : jet_pfcand_etarel_values) {
          v.h_jet_pfcand_etarel->Fill(val);
        }
        
        const auto& jet_pfcand_pperp_ratio_values = features.get("jet_pfcand_pperp_ratio");
        for (const auto& val : jet_pfcand_pperp_ratio_values) {
          v.h_jet_pfcand_pperp_ratio->Fill(val);
        }
        
        const auto& jet_pfcand_ppara_ratio_values = features.get("jet_pfcand_ppara_ratio");
        for (const auto& val : jet_pfcand_ppara_ratio_values) {
          v.h_jet_pfcand_ppara_ratio->Fill(val);
        }
        
        const auto& jet_pfcand_trackjet_d3d_values = features.get("jet_pfcand_trackjet_d3d");
        for (const auto& val : jet_pfcand_trackjet_d3d_values) {
          v.h_jet_pfcand_trackjet_d3d->Fill(val);
        }
        
        const auto& jet_pfcand_trackjet_d3dsig_values = features.get("jet_pfcand_trackjet_d3dsig");
        for (const auto& val : jet_pfcand_trackjet_d3dsig_values) {
          v.h_jet_pfcand_trackjet_d3dsig->Fill(val);
        }
        
        const auto& jet_pfcand_trackjet_dist_values = features.get("jet_pfcand_trackjet_dist");
        for (const auto& val : jet_pfcand_trackjet_dist_values) {
          v.h_jet_pfcand_trackjet_dist->Fill(val);
        }
        
        const auto& jet_pfcand_nhits_values = features.get("jet_pfcand_nhits");
        for (const auto& val : jet_pfcand_nhits_values) {
          v.h_jet_pfcand_nhits->Fill(val);
        }
        
        const auto& jet_pfcand_npixhits_values = features.get("jet_pfcand_npixhits");
        for (const auto& val : jet_pfcand_npixhits_values) {
          v.h_jet_pfcand_npixhits->Fill(val);
        }
        
        const auto& jet_pfcand_nstriphits_values = features.get("jet_pfcand_nstriphits");
        for (const auto& val : jet_pfcand_nstriphits_values) {
          v.h_jet_pfcand_nstriphits->Fill(val);
        }
        
        const auto& jet_pfcand_trackjet_decayL_values = features.get("jet_pfcand_trackjet_decayL");
        for (const auto& val : jet_pfcand_trackjet_decayL_values) {
          v.h_jet_pfcand_trackjet_decayL->Fill(val);
        }
        
        const auto& jet_pfcand_puppiw_values = features.get("jet_pfcand_puppiw");
        for (const auto& val : jet_pfcand_puppiw_values) {
          v.h_jet_pfcand_puppiw->Fill(val);
        }
      }
    }

    // ... otherwise from usual tag infos.
    // else
    // if (   (v.getTriggerType() == "PF"   && pfTagInfos.isValid())
    //     || (v.getTriggerType() == "Calo" && caloTagInfos.isValid()) )
    // {
    //   const auto & DiscrTagInfoCollection = (v.getTriggerType() == "PF") ? pfTagInfos : caloTagInfos;

    //   // loop over secondary vertex tag infos
    //   for (const auto & DiscrTagInfo : *DiscrTagInfoCollection) {
    //     v.n_vtx->Fill(DiscrTagInfo.nVertexCandidates());
    //     v.n_sel_tracks->Fill(DiscrTagInfo.nSelectedTracks());

    //     // loop over selected tracks in each tag info
    //     for (unsigned i_trk=0; i_trk < DiscrTagInfo.nSelectedTracks(); i_trk++) {
    //       const auto & ip3d = DiscrTagInfo.trackIPData(i_trk).ip3d;
    //       v.h_3d_ip_distance->Fill(ip3d.value());
    //       v.h_3d_ip_error->Fill(ip3d.error());
    //       v.h_3d_ip_sig->Fill(ip3d.significance());
    //     }

    //     // loop over vertex candidates in each tag info
    //     for (unsigned i_sv=0; i_sv < DiscrTagInfo.nVertexCandidates(); i_sv++) {
    //       const auto & sv = DiscrTagInfo.secondaryVertex(i_sv);
    //       v.vtx_mass->Fill(sv.p4().mass());
    //       v.n_vtx_trks->Fill(sv.nTracks());

    //       // loop over tracks for number of pixel and total hits
    //       const auto & trkIPTagInfo = DiscrTagInfo.trackIPTagInfoRef().get();
    //       for (const auto & trk : trkIPTagInfo->selectedTracks()) {
    //         v.n_pixel_hits->Fill(trk.get()->hitPattern().numberOfValidPixelHits());
    //         v.n_total_hits->Fill(trk.get()->hitPattern().numberOfValidHits());
    //       }
    //     }
    //   }
    // }
  }  //end paths loop
}

std::vector<const reco::Track*> BTVHLTOfflineSource::getOfflineBTagTracks(float hltJetEta,
                                                                          float hltJetPhi,
                                                                          Handle<View<BaseTagInfo>> offlineIPTagHandle,
                                                                          std::vector<float>& offlineIP3D,
                                                                          std::vector<float>& offlineIP3DSig) {
  std::vector<const reco::Track*> offlineTracks;

  for (auto const& iterOffIP : *offlineIPTagHandle) {
    float DR = reco::deltaR(iterOffIP.jet()->eta(), iterOffIP.jet()->phi(), hltJetEta, hltJetPhi);

    if (DR > 0.3)
      continue;

    const reco::IPTagInfo<vector<reco::CandidatePtr>, reco::JetTagInfo>* tagInfo =
        dynamic_cast<const reco::IPTagInfo<vector<reco::CandidatePtr>, reco::JetTagInfo>*>(&iterOffIP);

    if (!tagInfo) {
      throw cms::Exception("Configuration")
          << "BTagPerformanceAnalyzer: Extended TagInfo not of type TrackIPTagInfo. " << std::endl;
    }

    const GlobalPoint pv(tagInfo->primaryVertex()->position().x(),
                         tagInfo->primaryVertex()->position().y(),
                         tagInfo->primaryVertex()->position().z());

    const std::vector<reco::btag::TrackIPData>& ip = tagInfo->impactParameterData();

    std::vector<std::size_t> sortedIndices = tagInfo->sortedIndexes(reco::btag::IP2DSig);
    std::vector<reco::CandidatePtr> sortedTracks = tagInfo->sortedTracks(sortedIndices);
    std::vector<std::size_t> selectedIndices;
    vector<reco::CandidatePtr> selectedTracks;
    for (unsigned int n = 0; n != sortedIndices.size(); ++n) {
      double decayLength = (ip[sortedIndices[n]].closestToJetAxis - pv).mag();
      double jetDistance = ip[sortedIndices[n]].distanceToJetAxis.value();
      if (decayLength > minDecayLength_ && decayLength < maxDecayLength_ && fabs(jetDistance) >= minJetDistance_ &&
          fabs(jetDistance) < maxJetDistance_) {
        selectedIndices.push_back(sortedIndices[n]);
        selectedTracks.push_back(sortedTracks[n]);
      }
    }

    for (unsigned int n = 0; n != selectedIndices.size(); ++n) {
      const reco::Track* track = reco::btag::toTrack(selectedTracks[n]);
      offlineTracks.push_back(track);
      offlineIP3D.push_back(ip[n].ip3d.value());
      offlineIP3DSig.push_back(ip[n].ip3d.significance());
    }
  }
  return offlineTracks;
}

template <class Base>
std::vector<const reco::Track*> BTVHLTOfflineSource::getOnlineBTagTracks(float hltJetEta,
                                                                         float hltJetPhi,
                                                                         edm::Handle<std::vector<Base>> jetSVTagsColl,
                                                                         std::vector<float>& onlineIP3D,
                                                                         std::vector<float>& onlineIP3DSig) {
  std::vector<const reco::Track*> onlineTracks;

  for (auto iterTI = jetSVTagsColl->begin(); iterTI != jetSVTagsColl->end(); ++iterTI) {
    float DR = reco::deltaR(iterTI->jet()->eta(), iterTI->jet()->phi(), hltJetEta, hltJetPhi);
    if (DR > 0.3)
      continue;

    const auto& ipInfo = *(iterTI->trackIPTagInfoRef().get());
    const std::vector<reco::btag::TrackIPData>& ip = ipInfo.impactParameterData();

    unsigned int trackSize = ipInfo.selectedTracks().size();
    for (unsigned int itt = 0; itt < trackSize; ++itt) {
      const auto ptrackRef = (ipInfo.selectedTracks()[itt]);

      if (ptrackRef.isAvailable()) {
        const reco::Track* ptrackPtr = reco::btag::toTrack(ptrackRef);
        onlineTracks.push_back(ptrackPtr);
        onlineIP3D.push_back(ip[itt].ip3d.value());
        onlineIP3DSig.push_back(ip[itt].ip3d.significance());
      }
    }
  }
  return onlineTracks;
}

void BTVHLTOfflineSource::bookHistograms(DQMStore::IBooker& iBooker, edm::Run const& run, edm::EventSetup const& c) {
  iBooker.setCurrentFolder(dirname_);
  for (auto& v : hltPathsAll_) {
    std::string trgPathName = HLTConfigProvider::removeVersion(v.getPath());
    std::string subdirName = dirname_ + "/" + trgPathName + v.getTriggerType();
    std::string trigPath = "(" + trgPathName + ")";
    iBooker.setCurrentFolder(subdirName);

    std::string labelname("HLT");
    std::string histoname(labelname + "");
    std::string title(labelname + "");

    histoname = labelname + "_Discr";
    title = labelname + "_Discr " + trigPath;
    v.Discr = iBooker.book1D(histoname.c_str(), title.c_str(), 110, -0.1, 1);

    histoname = labelname + "_Pt";
    title = labelname + "_Pt " + trigPath;
    v.Pt = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 400);

    histoname = labelname + "_Eta";
    title = labelname + "_Eta " + trigPath;
    v.Eta = iBooker.book1D(histoname.c_str(), title.c_str(), 60, -3.0, 3.0);

    histoname = "HLTvsRECO_Discr";
    title = "online discr vs offline discr " + trigPath;
    v.Discr_HLTvsRECO = iBooker.book2D(histoname.c_str(), title.c_str(), 110, -0.1, 1, 110, -0.1, 1);

    histoname = "HLTMinusRECO_Discr";
    title = "online discr minus offline discr " + trigPath;
    v.Discr_HLTMinusRECO = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -1, 1);

    histoname = "Turnon_loose_Discr";
    title = "turn-on with loose threshold " + trigPath;
    v.bookME(iBooker, v.Discr_turnon_loose, histoname, title, 22, -0.1, 1.);

    histoname = "Turnon_medium_Discr";
    title = "turn-on with medium threshold " + trigPath;
    v.bookME(iBooker, v.Discr_turnon_medium, histoname, title, 22, -0.1, 1.);

    histoname = "Turnon_tight_Discr";
    title = "turn-on with tight threshold " + trigPath;
    v.bookME(iBooker, v.Discr_turnon_tight, histoname, title, 22, -0.1, 1.);

    histoname = "Turnon_loose_Pt";
    title = "turn-on with loose threshold " + trigPath;
    v.bookME(iBooker, v.Pt_turnon_loose, histoname, title, 50, 0., 500.);

    histoname = "Turnon_medium_Pt";
    title = "turn-on with medium threshold " + trigPath;
    v.bookME(iBooker, v.Pt_turnon_medium, histoname, title, 50, 0., 500.);

    histoname = "Turnon_tight_Pt";
    title = "turn-on with tight threshold " + trigPath;
    v.bookME(iBooker, v.Pt_turnon_tight, histoname, title, 50, 0., 500.);

    histoname = "Turnon_loose_Eta";
    title = "turn-on with loose threshold " + trigPath;
    v.bookME(iBooker, v.Eta_turnon_loose, histoname, title, 60, -3., 3.);

    histoname = "Turnon_medium_Eta";
    title = "turn-on with medium threshold " + trigPath;
    v.bookME(iBooker, v.Eta_turnon_medium, histoname, title, 60, -3., 3.);

    histoname = "Turnon_tight_Eta";
    title = "turn-on with tight threshold " + trigPath;
    v.bookME(iBooker, v.Eta_turnon_tight, histoname, title, 60, -3., 3.);

    histoname = "Turnon_loose_Phi";
    title = "turn-on with loose threshold " + trigPath;
    v.bookME(iBooker, v.Phi_turnon_loose, histoname, title, 60, -3., 3.);

    histoname = "Turnon_medium_Phi";
    title = "turn-on with medium threshold " + trigPath;
    v.bookME(iBooker, v.Phi_turnon_medium, histoname, title, 60, -3., 3.);

    histoname = "Turnon_tight_Phi";
    title = "turn-on with tight threshold " + trigPath;
    v.bookME(iBooker, v.Phi_turnon_tight, histoname, title, 60, -3., 3.);

    histoname = labelname + "_PVz";
    title = "online z(PV) " + trigPath;
    v.PVz = iBooker.book1D(histoname.c_str(), title.c_str(), 80, -20, 20);

    histoname = labelname + "_fastPVz";
    title = "online z(fastPV) " + trigPath;
    v.fastPVz = iBooker.book1D(histoname.c_str(), title.c_str(), 80, -20, 20);

    histoname = "HLTMinusRECO_PVz";
    title = "online z(PV) - offline z(PV) " + trigPath;
    v.PVz_HLTMinusRECO = iBooker.book1D(histoname.c_str(), title.c_str(), 200, -0.5, 0.5);

    histoname = "HLTMinusRECO_fastPVz";
    title = "online z(fastPV) - offline z(PV) " + trigPath;
    v.fastPVz_HLTMinusRECO = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -2, 2);

    histoname = "n_vtx";
    title = "N vertex candidates " + trigPath;
    v.n_vtx = iBooker.book1D(histoname.c_str(), title.c_str(), 10, -0.5, 9.5);

    histoname = "vtx_mass";
    title = "secondary vertex mass (GeV)" + trigPath;
    v.vtx_mass = iBooker.book1D(histoname.c_str(), title.c_str(), 20, 0, 10);

    histoname = "n_vtx_trks";
    title = "N tracks associated to secondary vertex" + trigPath;
    v.n_vtx_trks = iBooker.book1D(histoname.c_str(), title.c_str(), 20, -0.5, 19.5);

    histoname = "n_sel_tracks";
    title = "N selected tracks" + trigPath;
    v.n_sel_tracks = iBooker.book1D(histoname.c_str(), title.c_str(), 25, -0.5, 24.5);

    histoname = "3d_ip_distance";
    title = "3D IP distance of tracks (cm)" + trigPath;
    v.h_3d_ip_distance = iBooker.book1D(histoname.c_str(), title.c_str(), 40, -0.1, 0.1);

    histoname = "3d_ip_error";
    title = "3D IP error of tracks (cm)" + trigPath;
    v.h_3d_ip_error = iBooker.book1D(histoname.c_str(), title.c_str(), 40, 0., 0.1);

    histoname = "3d_ip_sig";
    title = "3D IP significance of tracks (cm)" + trigPath;
    v.h_3d_ip_sig = iBooker.book1D(histoname.c_str(), title.c_str(), 40, -40, 40);

    //DeepBoostedJetTagInfo
    histoname = "pfcand_puppiw";
    title = "PF candidate puppi weight" + trigPath;
    v.h_pfcand_puppiw = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "pfcand_hcalFrac";
    title = "PF candidate HCAL fraction" + trigPath;
    v.h_pfcand_hcalFrac = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "pfcand_VTX_ass";
    title = "PF candidate vertex association" + trigPath;
    v.h_pfcand_VTX_ass = iBooker.book1D(histoname.c_str(), title.c_str(), 10, -0.5, 9.5);

    histoname = "pfcand_lostInnerHits";
    title = "PF candidate lost inner hits" + trigPath;
    v.h_pfcand_lostInnerHits = iBooker.book1D(histoname.c_str(), title.c_str(), 10, -0.5, 9.5);

    histoname = "pfcand_quality";
    title = "PF candidate quality" + trigPath;
    v.h_pfcand_quality = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "pfcand_charge";
    title = "PF candidate charge" + trigPath;
    v.h_pfcand_charge = iBooker.book1D(histoname.c_str(), title.c_str(), 5, -2.5, 2.5);

    histoname = "pfcand_isEl";
    title = "PF candidate is electron" + trigPath;
    v.h_pfcand_isEl = iBooker.book1D(histoname.c_str(), title.c_str(), 2, -0.5, 1.5);

    histoname = "pfcand_isMu";
    title = "PF candidate is muon" + trigPath;
    v.h_pfcand_isMu = iBooker.book1D(histoname.c_str(), title.c_str(), 2, -0.5, 1.5);

    histoname = "pfcand_isChargedHad";
    title = "PF candidate is charged hadron" + trigPath;
    v.h_pfcand_isChargedHad = iBooker.book1D(histoname.c_str(), title.c_str(), 2, -0.5, 1.5);

    histoname = "pfcand_isGamma";
    title = "PF candidate is gamma" + trigPath;
    v.h_pfcand_isGamma = iBooker.book1D(histoname.c_str(), title.c_str(), 2, -0.5, 1.5);

    histoname = "pfcand_isNeutralHad";
    title = "PF candidate is neutral hadron" + trigPath;
    v.h_pfcand_isNeutralHad = iBooker.book1D(histoname.c_str(), title.c_str(), 2, -0.5, 1.5);

    histoname = "pfcand_phirel";
    title = "PF candidate relative phi" + trigPath;
    v.h_pfcand_phirel = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -3.2, 3.2);

    histoname = "pfcand_etarel";
    title = "PF candidate relative eta" + trigPath;
    v.h_pfcand_etarel = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 5);

    histoname = "pfcand_deltaR";
    title = "PF candidate delta R" + trigPath;
    v.h_pfcand_deltaR = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.5);

    histoname = "pfcand_abseta";
    title = "PF candidate absolute eta" + trigPath;
    v.h_pfcand_abseta = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 3);

    histoname = "pfcand_ptrel_log";
    title = "PF candidate log(ptrel)" + trigPath;
    v.h_pfcand_ptrel_log = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -10, 5);

    histoname = "pfcand_erel_log";
    title = "PF candidate log(erel)" + trigPath;
    v.h_pfcand_erel_log = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -10, 5);

    histoname = "pfcand_pt_log";
    title = "PF candidate log(pt)" + trigPath;
    v.h_pfcand_pt_log = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 10);

    histoname = "pfcand_drminsv";
    title = "PF candidate min deltaR to SV" + trigPath;
    v.h_pfcand_drminsv = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.5);

    histoname = "pfcand_drsubjet1";
    title = "PF candidate deltaR to subjet1" + trigPath;
    v.h_pfcand_drsubjet1 = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.5);

    histoname = "pfcand_drsubjet2";
    title = "PF candidate deltaR to subjet2" + trigPath;
    v.h_pfcand_drsubjet2 = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.5);

    histoname = "pfcand_normchi2";
    title = "PF candidate normalized chi2" + trigPath;
    v.h_pfcand_normchi2 = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "pfcand_dz";
    title = "PF candidate dz" + trigPath;
    v.h_pfcand_dz = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -20, 20);

    histoname = "pfcand_dzsig";
    title = "PF candidate dz significance" + trigPath;
    v.h_pfcand_dzsig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "pfcand_dxy";
    title = "PF candidate dxy" + trigPath;
    v.h_pfcand_dxy = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -2, 2);

    histoname = "pfcand_dxysig";
    title = "PF candidate dxy significance" + trigPath;
    v.h_pfcand_dxysig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "pfcand_dptdpt";
    title = "PF candidate dptdpt" + trigPath;
    v.h_pfcand_dptdpt = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "pfcand_detadeta";
    title = "PF candidate detadeta" + trigPath;
    v.h_pfcand_detadeta = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.1);

    histoname = "pfcand_dphidphi";
    title = "PF candidate dphidphi" + trigPath;
    v.h_pfcand_dphidphi = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.1);

    histoname = "pfcand_dxydxy";
    title = "PF candidate dxydxy" + trigPath;
    v.h_pfcand_dxydxy = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 0.1);

    histoname = "pfcand_dzdz";
    title = "PF candidate dzdz" + trigPath;
    v.h_pfcand_dzdz = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "pfcand_dxydz";
    title = "PF candidate dxydz" + trigPath;
    v.h_pfcand_dxydz = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -1, 1);

    histoname = "pfcand_dphidxy";
    title = "PF candidate dphidxy" + trigPath;
    v.h_pfcand_dphidxy = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -0.1, 0.1);

    histoname = "pfcand_dlambdadz";
    title = "PF candidate dlambdadz" + trigPath;
    v.h_pfcand_dlambdadz = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -0.1, 0.1);

    histoname = "pfcand_btagEtaRel";
    title = "PF candidate btag EtaRel" + trigPath;
    v.h_pfcand_btagEtaRel = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "pfcand_btagPtRatio";
    title = "PF candidate btag PtRatio" + trigPath;
    v.h_pfcand_btagPtRatio = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "pfcand_btagPParRatio";
    title = "PF candidate btag PParRatio" + trigPath;
    v.h_pfcand_btagPParRatio = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "pfcand_btagSip2dVal";
    title = "PF candidate btag Sip2dVal" + trigPath;
    v.h_pfcand_btagSip2dVal = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -20, 20);

    histoname = "pfcand_btagSip2dSig";
    title = "PF candidate btag Sip2dSig" + trigPath;
    v.h_pfcand_btagSip2dSig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "pfcand_btagSip3dVal";
    title = "PF candidate btag Sip3dVal" + trigPath;
    v.h_pfcand_btagSip3dVal = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -20, 20);

    histoname = "pfcand_btagSip3dSig";
    title = "PF candidate btag Sip3dSig" + trigPath;
    v.h_pfcand_btagSip3dSig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "pfcand_btagJetDistVal";
    title = "PF candidate btag JetDistVal" + trigPath;
    v.h_pfcand_btagJetDistVal = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -0.1, 0.01);

    histoname = "pfcand_mask";
    title = "PF candidate mask" + trigPath;
    v.h_pfcand_mask = iBooker.book1D(histoname.c_str(), title.c_str(), 2, -0.5, 1.5);

    histoname = "pfcand_pt_log_nopuppi";
    title = "PF candidate log(pt) no puppi" + trigPath;
    v.h_pfcand_pt_log_nopuppi = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 10);

    histoname = "pfcand_e_log_nopuppi";
    title = "PF candidate log(e) no puppi" + trigPath;
    v.h_pfcand_e_log_nopuppi = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 10);

    histoname = "pfcand_ptrel";
    title = "PF candidate ptrel" + trigPath;
    v.h_pfcand_ptrel = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "pfcand_erel";
    title = "PF candidate erel" + trigPath;
    v.h_pfcand_erel = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "jet_pfcand_pt_log";
    title = "Jet PF candidate log(pt)" + trigPath;
    v.h_jet_pfcand_pt_log = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 10);

    histoname = "jet_pfcand_energy_log";
    title = "Jet PF candidate log(energy)" + trigPath;
    v.h_jet_pfcand_energy_log = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 10);

    histoname = "jet_pfcand_deta";
    title = "Jet PF candidate deta" + trigPath;
    v.h_jet_pfcand_deta = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 5);

    histoname = "jet_pfcand_dphi";
    title = "Jet PF candidate dphi" + trigPath;
    v.h_jet_pfcand_dphi = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -3.2, 3.2);

    histoname = "jet_pfcand_eta";
    title = "Jet PF candidate eta" + trigPath;
    v.h_jet_pfcand_eta = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -3, 3);

    histoname = "jet_pfcand_charge";
    title = "Jet PF candidate charge" + trigPath;
    v.h_jet_pfcand_charge = iBooker.book1D(histoname.c_str(), title.c_str(), 5, -2.5, 2.5);

    histoname = "jet_pfcand_frompv";
    title = "Jet PF candidate from PV" + trigPath;
    v.h_jet_pfcand_frompv = iBooker.book1D(histoname.c_str(), title.c_str(), 10, -0.5, 9.5);

    histoname = "jet_pfcand_nlostinnerhits";
    title = "Jet PF candidate n lost inner hits" + trigPath;
    v.h_jet_pfcand_nlostinnerhits = iBooker.book1D(histoname.c_str(), title.c_str(), 10, -0.5, 9.5);

    histoname = "jet_pfcand_track_chi2";
    title = "Jet PF candidate track chi2" + trigPath;
    v.h_jet_pfcand_track_chi2 = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 20);

    histoname = "jet_pfcand_track_qual";
    title = "Jet PF candidate track quality" + trigPath;
    v.h_jet_pfcand_track_qual = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "jet_pfcand_dz";
    title = "Jet PF candidate dz" + trigPath;
    v.h_jet_pfcand_dz = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -20, 20);

    histoname = "jet_pfcand_dzsig";
    title = "Jet PF candidate dz significance" + trigPath;
    v.h_jet_pfcand_dzsig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "jet_pfcand_dxy";
    title = "Jet PF candidate dxy" + trigPath;
    v.h_jet_pfcand_dxy = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -2, 2);

    histoname = "jet_pfcand_dxysig";
    title = "Jet PF candidate dxy significance" + trigPath;
    v.h_jet_pfcand_dxysig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "jet_pfcand_etarel";
    title = "Jet PF candidate etarel" + trigPath;
    v.h_jet_pfcand_etarel = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -5, 5);

    histoname = "jet_pfcand_pperp_ratio";
    title = "Jet PF candidate pperp ratio" + trigPath;
    v.h_jet_pfcand_pperp_ratio = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "jet_pfcand_ppara_ratio";
    title = "Jet PF candidate ppara ratio" + trigPath;
    v.h_jet_pfcand_ppara_ratio = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    histoname = "jet_pfcand_trackjet_d3d";
    title = "Jet PF candidate trackjet d3d" + trigPath;
    v.h_jet_pfcand_trackjet_d3d = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -10, 10);

    histoname = "jet_pfcand_trackjet_d3dsig";
    title = "Jet PF candidate trackjet d3dsig" + trigPath;
    v.h_jet_pfcand_trackjet_d3dsig = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -100, 100);

    histoname = "jet_pfcand_trackjet_dist";
    title = "Jet PF candidate trackjet dist" + trigPath;
    v.h_jet_pfcand_trackjet_dist = iBooker.book1D(histoname.c_str(), title.c_str(), 100, -0.5, 0.5);

    histoname = "jet_pfcand_nhits";
    title = "Jet PF candidate nhits" + trigPath;
    v.h_jet_pfcand_nhits = iBooker.book1D(histoname.c_str(), title.c_str(), 40, -0.5, 39.5);

    histoname = "jet_pfcand_npixhits";
    title = "Jet PF candidate npixhits" + trigPath;
    v.h_jet_pfcand_npixhits = iBooker.book1D(histoname.c_str(), title.c_str(), 15, -0.5, 14.5);

    histoname = "jet_pfcand_nstriphits";
    title = "Jet PF candidate nstriphits" + trigPath;
    v.h_jet_pfcand_nstriphits = iBooker.book1D(histoname.c_str(), title.c_str(), 30, -0.5, 29.5);

    histoname = "jet_pfcand_trackjet_decayL";
    title = "Jet PF candidate trackjet decayL" + trigPath;
    v.h_jet_pfcand_trackjet_decayL = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 10);

    histoname = "jet_pfcand_puppiw";
    title = "Jet PF candidate puppi weight" + trigPath;
    v.h_jet_pfcand_puppiw = iBooker.book1D(histoname.c_str(), title.c_str(), 100, 0, 1);

    //end DeepBoostedJetTagInfo

    histoname = "OnlineTrkEff_Pt";
    title = "Relative Online Trk Efficiency vs Pt " + trigPath;
    v.bookME(iBooker, v.OnlineTrkEff_Pt, histoname, title, 50, -0.5, 20.);

    histoname = "OnlineTrkEff_Eta";
    title = "Relative Online Trk Efficiency vs Eta " + trigPath;
    v.bookME(iBooker, v.OnlineTrkEff_Eta, histoname, title, 60, -3.0, 3.0);

    histoname = "OnlineTrkEff_3d_ip_distance";
    title = "Relative Online Trk Efficiency vs IP3D " + trigPath;
    v.bookME(iBooker, v.OnlineTrkEff_3d_ip_distance, histoname, title, 40, -0.1, 0.1);

    histoname = "OnlineTrkEff_3d_ip_sig";
    title = "Relative Online Trk Efficiency vs IP3D significance " + trigPath;
    v.bookME(iBooker, v.OnlineTrkEff_3d_ip_sig, histoname, title, 40, -40, 40);

    histoname = "OnlineTrkFake_Pt";
    title = "Relative Online Trk Fake Rate  vs Pt " + trigPath;
    v.bookME(iBooker, v.OnlineTrkFake_Pt, histoname, title, 50, -0.5, 20.);

    histoname = "OnlineTrkFake_Eta";
    title = "Relative Online Trk Fake Rate vs Eta " + trigPath;
    v.bookME(iBooker, v.OnlineTrkFake_Eta, histoname, title, 60, -3.0, 3.0);

    histoname = "OnlineTrkFake_3d_ip_distance";
    title = "Relative Online Trk Fake Rate vs IP3D " + trigPath;
    v.bookME(iBooker, v.OnlineTrkFake_3d_ip_distance, histoname, title, 40, -0.1, 0.1);

    histoname = "OnlineTrkFake_3d_ip_sig";
    title = "Relative Online Trk Fake Rate vs IP3D significance " + trigPath;
    v.bookME(iBooker, v.OnlineTrkFake_3d_ip_sig, histoname, title, 40, -40, 40);

    // histoname = "n_pixel_hits";
    // title = "N pixel hits"+trigPath;
    // v.n_pixel_hits = iBooker.book1D(histoname.c_str(), title.c_str(), 16, -0.5, 15.5);

    // histoname = "n_total_hits";
    // title = "N hits"+trigPath;
    // v.n_total_hits = iBooker.book1D(histoname.c_str(), title.c_str(), 40, -0.5, 39.5);
  }
}

// Define this as a plug-in
DEFINE_FWK_MODULE(BTVHLTOfflineSource);
