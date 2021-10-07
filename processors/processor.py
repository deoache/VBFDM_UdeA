import numpy as np
import awkward as ak
import hist as hist2
from coffea import processor
from coffea.analysis_tools import Weights, PackedSelection


def invariant_mass(obj1, obj2):
    """invariant mass"""
    return np.sqrt(2 * obj1.pt * obj2.pt * (np.cosh(obj1.eta - obj2.eta) - np.cos(obj1.phi - obj2.phi)))

def normalize(val, cut=None):
    """normalize arrays after a cut or selection"""
    if cut is None:
        ar = ak.to_numpy(ak.fill_none(val, np.nan))
        return ar
    else:
        ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
        return ar
        
class VBFProcessor(processor.ProcessorABC):
    def __init__(self):
        
        # here we define an output
        self.make_output = lambda: {
            'sumw': 0.,
            "jet_kin": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.Regular(30, 0, 1000, name='H', label='$H_T$ [GeV]'),
                hist2.axis.Regular(50, 30, 1000, name="jet_pt", label=r"Jet $p_T$ [GeV]"),
                hist2.axis.Regular(20, -5, 5, name="jet_eta", label=r"Jet $\eta$"),
                hist2.axis.Regular(30, 0, 200, name="jet_mass", label=r"Jet $m$ [GeV]"),
                hist2.storage.Weight(),
            ),
            "met_kin": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.Regular(50, 0, 1000, name="met_pt", label=r"$p_T^{miss}$ [GeV]"),
                hist2.axis.Regular(20, 0, 4, name="minDeltaPhi", label="min $|\Delta \phi (p_T^{miss}, J_i)|$"),
                hist2.storage.Weight(),
            ),
            "dijet_kin": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.Regular(30, 0, 4, name='leadingJetsDphi', label='$|\Delta \phi_{jj}|$'),
                hist2.axis.Regular(30, 0, 5000, name='invariantMass', label=r"$m_{jj}$ [GeV]"),
                hist2.axis.Regular(20, 0, 10, name='DeltaEta', label="$|\Delta \eta_{jj}|$"),
                hist2.storage.Weight(),
            ),
            "cutflow": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.IntCategory([0, 1, 2, 3], name='cut', growth=True),
                hist2.axis.Regular(50, 0, 1000, name="met_pt", label=r"$p_T^{miss}$ [GeV]"),
                hist2.storage.Weight(),
            ),
            "cutflow_DeltaEta": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.IntCategory([0, 1, 2, 3], name='cut', label='Cut index', growth=True),
                hist2.axis.Regular(20, 0, 10, name='maxDeltaEta', label="$|\Delta \eta_{jj}|$"),
                hist2.storage.Weight(),
            ),
            "cutflow_invariantMass": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.IntCategory([0, 1, 2, 3], name='cut', label='Cut index', growth=True),
                hist2.axis.Regular(30, 0, 5000, name='invariantMass', label=r"$m_{jj}$ [GeV]"),
                hist2.storage.Weight(),
            ),
            "cutflow_leadingJetsDphi": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.IntCategory([0, 1, 2, 3], name='cut', label='Cut index', growth=True),
                hist2.axis.Regular(30, 0, 4, name='leadingJetsDphi', label='$|\Delta \phi_{jj}|$'),
                hist2.storage.Weight(),
            ),
            "cutflow_minDeltaPhi": hist2.Hist(
                hist2.axis.StrCategory([], name='region', growth=True),
                hist2.axis.IntCategory([0, 1, 2, 3], name='cut', label='Cut index', growth=True),
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
        
        # jets (only first 4)
        jets = ak.pad_none(events.Jet, 4, clip=True)
        
        # leading jets
        ledjets_dphi = jets[:,0].delta_phi(jets[:,1])
        ledjets_deta = abs(jets[:,0].eta - jets[:,1].eta)
        ledjets_mass = invariant_mass(jets[:,0], jets[:,1])
        
        # minimun phi difference between jets and met
        met = events.MET
        min_dphi_met_jet = ak.min(
            ak.pad_none(
                abs(events.Jet.delta_phi(met)), 
                target=4, 
                clip=True
            ), 
            axis=1,
        )
        
        # baseline selection:
        #   leading jets pT > 30
        #   leading jets |eta| < 5
        #   leading jets in opposite hemispheres
        #   HT > 200 
        baseline = (
            (jets.pt[:,0] > 30)
            & (jets.pt[:,1] > 30)
            & (abs(jets.eta[:,0]) < 5)
            & (abs(jets.eta[:,1]) < 5)
            & (jets.eta[:,0] * jets.eta[:,1] < 0)
            & (ak.sum(jets.pt, axis=1) > 200)
        ) 
        
        # regions
        selection.add("baseline", baseline)
        selection.add("ledjets_dphi", abs(ledjets_dphi) > 2.3)
        selection.add("invariantMass", ledjets_mass > 1000)
        selection.add("min_dphi_met_jet", min_dphi_met_jet > 0.5)
        selection.add("ledjets_deta_l", ledjets_deta < 2.5)
        selection.add("ledjets_deta_r", ledjets_deta > 2.5)
        
        regions = {
            "no_selection": [],
            "baseline": ["baseline"],
            "high_mass": ["baseline","ledjets_dphi","invariantMass","min_dphi_met_jet","ledjets_deta_l"],
            "low_mass": ["baseline","ledjets_dphi","invariantMass","min_dphi_met_jet","ledjets_deta_r"],
        }
        
        def fill(region):
            selections = regions[region]
            cut = selection.all(*selections)
            
            output["jet_kin"].fill(
                region=region,
                H=normalize(ak.sum(jets.pt, axis=1), cut),
                jet_pt=normalize(ak.firsts(jets.pt), cut),
                jet_eta=normalize(ak.firsts(jets.eta), cut),
                jet_mass=normalize(ak.firsts(jets.mass), cut),
                weight=weights.weight()[cut],
            )
            output["met_kin"].fill(
                region=region,
                met_pt=normalize(met.pt, cut),
                minDeltaPhi=normalize(min_dphi_met_jet, cut),
                weight=weights.weight()[cut],
            )
            output['dijet_kin'].fill(
                region=region,
                leadingJetsDphi=normalize(ledjets_dphi, cut),
                invariantMass=normalize(ledjets_mass, cut),
                DeltaEta=normalize(ledjets_deta, cut),
                weight=weights.weight()[cut],
            )
            
            # cutflow
            allcuts = set([])
            cut = selection.all(*allcuts)
            output["cutflow"].fill(
                region=region,
                cut=0,
                met_pt=normalize(met.pt, cut),
                weight=weights.weight()[cut],
            )
            output["cutflow_DeltaEta"].fill(
                region=region,
                cut=0,
                maxDeltaEta=normalize(ledjets_deta, cut),
                weight=weights.weight()[cut],
            )
            output["cutflow_invariantMass"].fill(
                region=region,
                cut=0,
                invariantMass=normalize(ledjets_mass, cut),
                weight=weights.weight()[cut],
            )
            output["cutflow_leadingJetsDphi"].fill(
                region=region,
                cut=0,
                leadingJetsDphi=normalize(ledjets_dphi, cut),
                weight=weights.weight()[cut],
            )
            output["cutflow_minDeltaPhi"].fill(
                region=region,
                cut=0,
                minDeltaPhi=normalize(min_dphi_met_jet, cut),
                weight=weights.weight()[cut],
            )
            for i, cut in enumerate(regions[region]):
                allcuts.add(cut)
                cut = selection.all(*allcuts)
                
                output["cutflow"].fill(
                    region=region,
                    cut=i + 1,
                    met_pt=normalize(met.pt, cut),
                    weight=weights.weight()[cut],
                )
                output["cutflow_DeltaEta"].fill(
                    region=region,
                    cut=i + 1,
                    maxDeltaEta=normalize(ledjets_deta, cut),
                    weight=weights.weight()[cut],
                )
                output["cutflow_invariantMass"].fill(
                    region=region,
                    cut=i + 1,
                    invariantMass=normalize(ledjets_mass, cut),
                    weight=weights.weight()[cut],
                )
                output["cutflow_leadingJetsDphi"].fill(
                    region=region,
                    cut=i + 1,
                    leadingJetsDphi=normalize(ledjets_dphi, cut),
                    weight=weights.weight()[cut],
                )
                output["cutflow_minDeltaPhi"].fill(
                    region=region,
                    cut=i + 1,
                    minDeltaPhi=normalize(min_dphi_met_jet, cut),
                    weight=weights.weight()[cut],
                )
            
        for region in regions:
            fill(region)
                
        return {dataset: output}
            
    def postprocess(self, accumulator):
        return accumulator
