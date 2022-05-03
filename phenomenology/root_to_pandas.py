import json
import argparse
import pandas as pd
from utils import get_dataframe


def main(args):
    # samples to be converted
    with open(args.spath, "r") as f:
        samples = json.load(f)

    # variables to read
    variables = {"MissingET": ["Phi", "MET"], "Jet": ["PT", "Eta", "Phi", "Mass"]}

    output_df = pd.DataFrame()
    for dataset, fnames in samples.items():
        for fname in fnames:
            print(f"processing {dataset} {fname}...")
            output_df = pd.concat(
                [output_df, get_dataframe(fname, variables, label=dataset)], axis=0
            )

    # save processed data to a parquet file
    print("saving parquet file...")
    output_df.to_parquet(f"{args.fpath}/{args.name}.parquet")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--spath", dest="spath", default=None, type=str, help="ROOT samples path")
    parser.add_argument("--fpath", dest="fpath", default="/.", type=str, help="folder path to save the parquet file")
    parser.add_argument("--name",  dest="name",  default=None, type=str, help="name of the parquet file")
    args = parser.parse_args()

    main(args)