### Data fileset

The fileset json files that contain a dictionary of the files per sample are in the `fileset` directory.

<details><summary>Re-making the input dataset files with DAS</summary>
<p>

```bash
# connect to LPC with a port forward to access the jupyter notebook server
ssh USERNAME@cmslpc-sl7.fnal.gov -L8xxx:localhost:8xxx

# create a working directory and clone the repo (if you have not done yet)
# git clone git@github.com:deoache/VBFDM_UdeA.git

# enable the coffea environment, either the python environment
source coffeaenv/bin/activate

# or the conda environment
conda activate coffea-env

# then activate your proxy
voms-proxy-init --voms cms --valid 100:00

# activate cmsset
source /cvmfs/cms.cern.ch/cmsset_default.sh

# the json files are in the fileset directory
cd fileset/

here should be a link looking like `http://localhost:8xxx/?token=...`, displayed in the output at this point, paste that into your browser.
You should see a jupyter notebook with a directory listing.
Open `filesetDAS.ipynb`.

The .json files containing the datasets to be run should be saved in the same `fileset/` directory.

</p>
</details>

