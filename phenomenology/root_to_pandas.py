import json
import argparse
import numpy as np
import pandas as pd
from utils import get_dataframe


def main(args):
    # samples to be converted
    with open(args.spath, "r") as f:
        samples = json.load(f)

    # variables to read
    variables = {"MissingET": ["MET", "Eta", "Phi"], "Jet": ["PT", "Eta", "Phi", "Mass"]}

    output_df = pd.DataFrame()
    for dataset, fnames in samples.items():
        for fname in fnames:
            print(f"processing {dataset} {fname}...")
            output_df = pd.concat(
                [output_df, get_dataframe(fname, variables, label=dataset)], axis=0
            )

    # split dataframe and save processed data to parquet files
    print("saving compressed parquet files...")
    dfs = np.array_split(output_df, args.nfiles)
    for i, df in enumerate(dfs):
        df.to_parquet(f"{args.fpath}/{args.name}_{i}.parquet.gz", index=False, compression="gzip")
            

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--spath",  dest="spath",  default=None, type=str, help="ROOT samples path")
    parser.add_argument("--fpath",  dest="fpath",  default="/.", type=str, help="folder path to save the parquet file")
    parser.add_argument("--name",   dest="name",   default=None, type=str, help="name of the parquet file")
    parser.add_argument("--nfiles", dest="nfiles", default=3,    type=int, help="number of output parquet files")
    args = parser.parse_args()

    main(args)