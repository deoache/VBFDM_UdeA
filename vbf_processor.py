import numpy as np
import awkward as ak
import hist as hist2
from coffea import processor
from coffea.nanoevents.methods import candidate, vector
from coffea.analysis_tools import Weights


def invariant_mass(jet1, jet2):
    """di-jet invariant mass"""
    return np.sqrt(2 * jet1.pt * jet2.pt * (np.cosh(jet1.eta - jet2.eta) - np.cos(jet1.phi - jet2.phi)))


class BackgroundProcessor(processor.ProcessorABC):
    def __init__(self):
        
        # here we define an output
        self.make_output = lambda: {
            'sumw': 0.,
            "leading_jet_kin": hist2.Hist(
                hist2.axis.Regular(50, 30, 1000, name="jet_pt", label=r"Jet $p_T$ [GeV]"),
                hist2.axis.Regular(20, -5, 5, name="jet_eta", label=r"Jet $\eta$"),
                hist2.axis.Regular(30, 0, 200, name="jet_mass", label=r"Jet $m$ [GeV]"),
                hist2.storage.Weight(),
            ),
            "met_kin": hist2.Hist(
                hist2.axis.Regular(50, 0, 1000, name="met_pt", label=r"$p_T^{miss}$ [GeV]"),
                hist2.storage.Weight(),
            ),
            'kin': hist2.Hist(
                hist2.axis.Regular(30, 0, 1000, name='H', label='$H_T$ [GeV]'),
                hist2.axis.Regular(30, 0, 4, name='leadingJetsDphi', label='$|\Delta \phi (J_1, J_2)|$'),
                hist2.axis.Regular(30, 0, 5000, name='invariantMass', label=r"$m_{jj}^{max}$ [GeV]"),
                hist2.axis.Regular(20, 0, 10, name='maxDeltaEta', label="max $|\Delta \eta (J_i, J_j)|$"),
                hist2.axis.Regular(20, 0, 4, name="minDeltaPhi", label="min $|\Delta \phi (p_T^{miss}, J_i)|$"),
                hist2.storage.Weight(),
            ),
        }
        
    def process(self, events):
        dataset = events.metadata['dataset']
        weights = Weights(len(events), storeIndividual=True)
        weights.add('genweight', events.genWeight)
        
        selection = PackedSelection()
        
        output = self.make_output()
        output['sumw'] = ak.sum(events.genWeight)
        
        
        # jets
        jets = events.Jet
        
        # select events with at least two jets
        two_jets = ak.sum(jets.pt>0, axis=1) > 1
        jets = jets[two_jets]

        # baseline selection
        # leading jets: pT > 30 & |eta| < 5
        # leading jets in opposite hemispheres
        # H > 200 
        baseline = (
            (jets.pt[:,0] > 30)
            & (jets.pt[:,1] > 30)
            & (abs(jets.eta[:,0]) < 5)
            & (abs(jets.eta[:,1]) < 5)
            & (jets.eta[:,0] * jets.eta[:,1] < 0)
            & (ak.sum(jets.pt, axis=1) > 200)
        )
        selection.add("baseline", baseline) 
        
        # met
        met = events.MET
        met = met[two_jets]
        
        # delta phi for leading jets
        ledjets_dphi = jets[:,0].delta_phi(jets[:,1])
        
        selection.add("ledjets_dphi", abs(ledjets_dphi) > 2.3)
        
        # max invariant mass         
        jet_events = ak.zip(
            {
                "pt": jets.pt,
                "eta": jets.eta,
                "phi": jets.phi,
                "mass": jets.mass
            }
        )
        jet_pairs = ak.combinations(jet_events, 2)
        jet_pairs1, jet_pairs2 = ak.unzip(jet_pairs) 
        
        max_mjj = ak.max(invariant_mass(jet_pairs1, jet_pairs2), axis=1)
        
        selection.add("max_mjj", max_mjj > 1000)
        
        # minimun phi difference between jets and met
        min_dphi_met_jet = abs(ak.min(jets.delta_phi(met), axis=1))
        
        selection.add("min_dphi_met_jet", min_dphi_met_jet > 0.5)
        
        # max delta eta between jets
        max_delta_eta = abs(ak.max(jet_pairs1.eta - jet_pairs2.eta, axis=1))
        
        selection.add("max_delta_eta_l", max_delta_eta < 2.5)
        selection.add("max_delta_eta_r", max_delta_eta > 2.5)
        
        # function to normalize arrays after a cut or selection
        def normalize(val, cut=None):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar
    
        # here we fill our histogram
        output["leading_jet_kin"].fill(
            jet_pt=normalize(ak.firsts(jets.pt)),
            jet_eta=normalize(ak.firsts(jets.eta)),
            jet_mass=normalize(ak.firsts(jets.mass)),
            weight=weights.weight()[two_jets],
        )
        output["met_kin"].fill(
            met_pt=normalize(met.pt),
            weight=weights.weight()[two_jets],
        )
        output['kin'].fill(
            H=normalize(ak.sum(jets.pt, axis=1)),
            leadingJetsDphi = normalize(ledjets_dphi),
            invariantMass = normalize(max_mjj),
            maxDeltaEta = normalize(max_delta_eta),
            minDeltaPhi = normalize(min_dphi_met_jet),
            weight=weights.weight()[two_jets],
        )
        return {dataset: output}
            
    def postprocess(self, accumulator):
        return accumulator
