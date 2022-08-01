# prokaryotic-ipcr
This program is intended to assess efficiency of primers against prokaryotic organisms.

## Usage

First you should probably install the conda environment required:
```
conda env create -f environment.yml
```

After that's done, you should be good to initialize the database
```
python3 bin/main.py --init <name> -n  1000
```

And then run it:
```
python3 bin/main.py --run <name> -f  <forward-seq> -r  <forward-seq>
```
