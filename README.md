# Vector Boson Fusion Dark Matter 

## Setting up a coffea environment with conda

#### Install miniconda (if you do not have it already)
In your lxplus area or in your local computer:
```
# download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# run and follow instructions  
bash Miniconda3-latest-Linux-x86_64.sh

# Make sure to choose `yes` for the following one to let the installer initialize Miniconda3
# > Do you wish the installer to initialize Miniconda3
# > by running conda init? [yes|no]
```
Verify the installation is successful by running conda info and check if the paths are pointing to your Miniconda installation. 
If you cannot run conda command, check if you need to add the conda path to your PATH variable in your bashrc/zshrc file, e.g.,
```
export PATH="$HOME/nobackup/miniconda3/bin:$PATH"
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
conda install -c conda-forge coffea
conda install -c conda-forge numpy    
conda install -c conda-forge pandas
conda install -c conda-forge matplotlib
conda install -c conda-forge awkward
conda install -c conda-forge uproot
conda install -c conda-forge mplhep
conda install -c conda-forge hist
conda install -c conda-forge correctionlib
conda install -c conda-forge jupyterlab  
conda install -c conda-forge xrootd
```


## Data fileset

The fileset json files that contain a dictionary of the files per sample are in the `fileset` directory.

#### Re-making the input dataset files with DAS

```
# in your local machine, add the following to the ~/.ssh/config file
Host lxplus*
  HostName lxplus7.cern.ch
  User <your_username>
  ForwardX11 yes
  ForwardAgent yes
  ForwardX11Trusted yes
Host *_f
  LocalForward localhost:8800 localhost:8800
  ExitOnForwardFailure yes
  
# connect to lxplus with a port forward to access the jupyter notebook server
ssh lxplus_f

# create a working directory and clone the repo (if you have not done yet)
git clone https://github.com/deoache/VBFDM_UdeA.git

# enable the coffea environment, either the python environment
source coffeaenv/bin/activate

# or the conda environment
conda activate coffea-env

# then activate your proxy
voms-proxy-init --voms cms --valid 100:00

# activate cmsset
source /cvmfs/cms.cern.ch/cmsset_default.sh

# open the jupyter notebook on a browser
cd fileset/
jupyter notebook --no-browser --port 8800
```

there should be a link looking like `http://localhost:8800/?token=...`, displayed in the output at this point, paste that into your browser.
You should see a jupyter notebook with a directory listing.

Open `fileset/filesetDAS.ipynb` and run it. The json files containing the datasets to be run should be saved in the same `fileset/` directory.
