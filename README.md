## Setting up a coffea environment with conda

#### Install miniconda (if you do not have it already)
Preferably, in your `nobackup` area (in LPC) or in your local computer:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the instructions to finish the installation

# Make sure to choose `yes` for the following one to let the installer initialize Miniconda3
# > Do you wish the installer to initialize Miniconda3
# > by running conda init? [yes|no]
```
Verify the installation is successful by running conda info and check if the paths are pointing to your Miniconda installation. 
If you cannot run conda command, check if you need to add the conda path to your PATH variable in your bashrc/zshrc file, e.g.,
```
export PATH="$HOME/miniconda3/bin:$PATH"
```
To disable auto activation of the base environment:
```
conda config --set auto_activate_base false
```

#### Set up a conda environment and install the required packages
```
# create a new conda environment
conda create -n coffea-env python=3.7

# activate the environment
conda activate coffea-env

# install packages
pip install numpy matplotlib pandas scikit-learn coffea correctionlib mplhep hist

# install xrootd
conda install -c conda-forge xrootd
```


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

