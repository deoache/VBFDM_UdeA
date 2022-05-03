import json
import argparse
import pandas as pd
from utils import get_dataframe


def main(args):
    # samples to be converted
    with open(
        "/home/cms-user0/VBFDM_UdeA/phenomenology/dmspin0_samples.json", "r"
    ) as f:
        dmspin0_samples = json.load(f)

    # variables to read
    variables = {"MissingET": ["Phi", "MET"], "Jet": ["PT", "Eta", "Phi", "Mass"]}

    output_df = pd.DataFrame()
    for dataset, fnames in dmspin0_samples.items():
        for fname in fnames:
            print(f"processing {dataset} {fname}...")
            output_df = pd.concat(
                [output_df, get_dataframe(fname, variables, label=dataset)], axis=0
            )

    print("saving csv file...")
    output_df.to_csv(f"{args.fpath}/dmspin0.csv")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fpath",
        dest="fpath",
        default="/.",
        type=str,
        help="folder to save the csv file",
    )
    args = parser.parse_args()

    main(args)
