import numpy as np
import awkward as ak

from coffea import processor
from coffea.analysis_tools import PackedSelection


def invariant_mass(obj1, obj2):
    """returns the invariant mass of the system made of 'obj1' and 'obj2'"""
    return np.sqrt(2 * obj1.pt * obj2.pt * (np.cosh(obj1.eta - obj2.eta) - np.cos(obj1.phi - obj2.phi)))

def normalize(val, cut=None):
    """normalizes arrays after a cut or selection"""
    if cut is None:
        ar = ak.to_numpy(ak.fill_none(val, np.nan, axis=0))
        return ar
    else:
        ar = ak.to_numpy(ak.fill_none(val[cut], np.nan, axis=0))
        return ar

    
class VBFProcessor(processor.ProcessorABC):
    def __init__(self, year: int = 2017):
        self.year = year
    
        # variables, cutflow and output dictionaries
        self.variables = {}
        self.cutflow = {}
        self.output = {}
        
    def add_var(self, name: str, var: np.ndarray):
        """adds variable to output dictonary"""
        self.variables[name] = var
            
    def process(self, events):
        dataset = events.metadata['dataset']
        isMC = hasattr(events, "genWeight")
        
        # sum of genweights for MC data
        sumgenweight = ak.sum(events.genWeight) if isMC else None
        if isMC: self.cutflow["sumw"] = sumgenweight
        
        # jets (only the first 4) and missing energy
        jets = ak.pad_none(events.Jet, 4, clip=True)
        met = events.MET
            
        # variables 
        self.add_var("jet_pt", jets.pt)
        self.add_var("jet_eta", jets.eta)
        self.add_var("jet_mass", jets.mass)
        self.add_var("jet_H", ak.sum(jets.pt, axis=1))
        self.add_var("leading_jets_delta_phi", jets[:, 0].delta_phi(jets[:, 1]))
        self.add_var("leading_jets_delta_eta", abs(jets[:, 0].eta - jets[:, 1].eta))
        self.add_var("leading_jets_invariant_mass", invariant_mass(jets[:, 0], jets[:, 1]))
        self.add_var("met_pt", met.pt)
        self.add_var(
            name="min_met_jet_delta_phi", 
            var=ak.min(ak.pad_none(abs(events.Jet.delta_phi(met)), target=4, clip=True), axis=1)
        )
            
        # selections
        selection = PackedSelection()
        baseline = (
            (self.variables["jet_pt"][:, 0] > 30)
            & (self.variables["jet_pt"][:, 1] > 30)
            & (abs(self.variables["jet_eta"][:, 0]) < 5)
            & (abs(self.variables["jet_eta"][:, 1]) < 5)
            & (self.variables["jet_eta"][:, 0] * self.variables["jet_eta"][:, 1] < 0)
            & (self.variables["jet_H"] > 200)
            & (self.variables["met_pt"] > 50)
        )
        selection.add("baseline", baseline)
        selection.add("leading_jets_delta_phi", abs(self.variables["leading_jets_delta_phi"]) > 2.3)
        selection.add("leading_jets_invariant_mass", self.variables["leading_jets_invariant_mass"] > 1000)
        selection.add("min_met_jet_delta_phi", self.variables["min_met_jet_delta_phi"] > 0.5)
        selection.add("leading_jets_delta_eta_highMass", self.variables["leading_jets_delta_eta"] < 2.5)
        selection.add("leading_jets_delta_eta_lowMass", self.variables["leading_jets_delta_eta"] > 2.5)
        
        regions = {
            "high_mass": [
                "baseline", 
                "leading_jets_delta_phi",
                "leading_jets_invariant_mass",
                "min_met_jet_delta_phi",
                "leading_jets_delta_eta_highMass"
            ],
            "low_mass": [
                "baseline", 
                "leading_jets_delta_phi",
                "leading_jets_invariant_mass",
                "min_met_jet_delta_phi",
                "leading_jets_delta_eta_lowMass"
            ],
        }
        
        # initializing the cutflow dictionary
        for region in regions:
            cuts = {f"{cut}": 0 for cut in regions[region]}
            self.cutflow[region] = dict({"all": 0}, **cuts)
        
        for region in regions:
            # output variables
            selections = regions[region]
            cut = selection.all(*selections)
            
            self.output[region] = {
                key: processor.column_accumulator(normalize(val, cut))
                for key, val in self.variables.items()
            }
            
            # cutflow
            allcuts = set([])
            cut = selection.all(*allcuts)
            self.cutflow[region]["all"] = np.sum(cut)
            
            for sel in regions[region]:
                allcuts.add(sel)
                cut = selection.all(*allcuts)
                self.cutflow[region][sel] = np.sum(cut)

        return {
            dataset: {
                "variables": self.output, 
                "cutflow": self.cutflow
            }
        }
    
    def postprocess(self, accumulator):
        for dataset, output in accumulator.items():
            for region in output["variables"]:
                for var in output["variables"][region]:
                    output["variables"][region][var] = ak.from_numpy(output["variables"][region][var].value)
        return accumulator