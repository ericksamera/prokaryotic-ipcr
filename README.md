# prokaryotic-ipcr
This program is intended to assess efficacy of primers against prokaryotic organisms. It creates a database of sequences from Entrez fetch and performs in-silico PCR with `exonerate ipcress`. 

## Installation
This is pretty easy. Just use conda to create an environment with the required packages:
```
conda env create -n ipcr -f environment.yml
```
And then remember to activate:
```
conda activate ipcr
```

## Usage
After making sure the conda environment is active, initialize a database like so:
```
python3 bin/main.py --init <name> -n <number>
```

And then run in-silico PCR against that database: 
```
python3 bin/main.py --run <name> \
  -f <forward-seq> \
  -r <reverse-seq> \
  -m <mismatches> -M <memory (MB) to use> \
  -j '<job name>'
```

## Use-case example
This is an example use-case which compares the efficacy of 16S primers that target either both bacteria and archaea, or archaea specific primers:
```
#!/bin/bash

MISMATCH=6
DATABASE=prokaryota
MEMORY=64000

python3 bin/main.py --init $DATABASE -n <number>

python3 bin/main.py --run $DATABASE \
  -f GTGCCAGCMGCCGCGGTAA \
  -r GGACTACHVGGGTWTCTAAT \
  -m $MISMATCH -M $MEMORY \
  -j 'ILLUMINA V4'

python3 bin/main.py --run $DATABASE \
  -f GYGCASCAGKCGMGAAW \
  -r GGTDTTACCGCGGCKGCTG \
  -m $MISMATCH -M $MEMORY \
  -j 'ARCHAEA-SPECIFIC'
```

## Planned features
These are a list of things that I would eventually like to include in this pipeline.
- [ ] BLAST-based taxonomy classifier for in-silico PCR amplicons
- [ ] Functionality for batching primer sets with given names
- [ ] Better reporting
