import numpy as np
import pandas as pd
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, DelphesSchema
from coffea.nanoevents.methods.base import NanoEventsArray


def pad_none(array: ak.Array, target: int = 4, clip: bool = True) -> ak.Array:
    return ak.pad_none(array, target=target, clip=clip)


def min_delta_phi_met_jet(events: NanoEventsArray) -> ak.Array:
    """return the minimum delta_phi between missing energy and jets"""
    return ak.min(pad_none(abs(events.Jet.delta_phi(events.MissingET))), axis=1)


def add_vbf_composite_vars(events: NanoEventsArray, df: pd.DataFrame) -> None:
    """add some composite variables to df"""
    jets = pad_none(events.Jet)
    leading_jets = jets[:, 0]
    subleading_jets = jets[:, 1]

    df["h"] = ak.to_pandas(ak.sum(jets.pt, axis=1))
    df["invariant_mass"] = ak.to_pandas((leading_jets + subleading_jets).mass)
    df["delta_phi"] = ak.to_pandas(np.abs(leading_jets.delta_phi(subleading_jets)))
    df["delta_eta"] = ak.to_pandas(np.abs(leading_jets.eta - subleading_jets.eta))
    df["min_delta_phi_met_jet"] = ak.to_pandas(min_delta_phi_met_jet(events))


def add_feature(
    df: pd.DataFrame,
    events: NanoEventsArray,
    branch: str,
    leaf: str,
    target: int = 4,
    jagged: bool = False,
):
    """
    add the events[branch][leaf] feature as a column to df. If the feature array is jagged,
    then an unstacked dataframe with <target> columns is added
    """
    if jagged:
        columns = [f"{branch}_{leaf}_{i}".lower() for i in range(target)]
        df[columns] = ak.to_pandas(pad_none(events[branch][leaf])).unstack()
    else:
        df[f"{branch}_{leaf}".lower()] = ak.to_pandas(events[branch][leaf])


def get_dataframe(fname: str, vars: dict, dataset: str = None) -> pd.DataFrame:
    """
    return a pandas DataFrame from a ROOT file

    fname:
      path to the file
    vars:
      branches (keys) and leaves (values) to be read
    dataset:
      dataset name (optional)
    """
    events = NanoEventsFactory.from_root(
        fname, treepath="Delphes", schemaclass=DelphesSchema
    ).events()

    output_df = pd.DataFrame()
    output_df["dataset"] = dataset

    for branch in vars:
        for leaf in vars[branch]:
            if events[branch].layout.purelist_isregular:
                add_feature(output_df, events, branch, leaf, jagged=False)
            else:
                add_feature(output_df, events, branch, leaf, jagged=True)

    add_vbf_composite_vars(events, output_df)

    return output_df