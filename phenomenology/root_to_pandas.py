import pandas as pd
from utils import get_dataframe

datasets = {
    "dmspin0": {
        "my100_mx1000": [
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_01/tag_2_delphes_events.root",
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_02/tag_2_delphes_events.root",
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_03/tag_2_delphes_events.root",
        ],
        "my100_mx100": [
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_04/tag_2_delphes_events.root",
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_05/tag_1_delphes_events.root",
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_06/tag_1_delphes_events.root",
        ],
        "my5000_mx1000": [
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_07/tag_1_delphes_events.root",
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_08/tag_1_delphes_events.root",
            "/cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/Events/run_09/tag_1_delphes_events.root",
        ],
    }
}

variables = {
    "MissingET": ["Phi", "MET"],
    "Jet": ["PT", "Eta", "Phi", "Mass"],
}

output_df = pd.DataFrame()
for dataset, fnames in datasets["dmspin0"].items():
    for fname in fnames:
        output_df = pd.concat(
            [output_df, get_dataframe(fname, variables, label=dataset)], axis=0
        )

output_df.to_csv("dmspin0.csv")
