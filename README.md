# BRIQ2
De novo modeling of RNA tertiary structures.

## Installation
Currently, BRIQ2 is only available for Unix-based systems.

### Prerequisites
- *nix operating system
- gcc 12.3

### Installation steps
1. Clone the GitHub repo of BRIQ2 and change to the directory.
```bash
$ git clone git@github.com:xionglab2023/BRIQ2.git
$ cd BRIQ2
```

2. Download data from [Google Drive]() and decompress it in BRIQ2 directory.
```bash
$ tar -xzf briq_data.tar.gz
```

3. Modify variable of `BRIQX_DATAPATH` in setup.sh.
 
4. Add BRIQ2 binary to system path, temporary change:
```bash
$ export PATH=$PATH:$HOME/apps/BRIQX/bin
```
Permanent change: Modify `$HOME/.bashrc`
 
5. Run setup program.
```bash
$ bash ./setup.sh
```

## Usage
Please refer to the [manual](./BRIQ_Manual.pdf).


## Citation
If you use BRIQ in your work, please cite as follows.
```txt
@article{xiong2021pairing,
  title={Pairing a high-resolution statistical potential with a nucleobase-centric sampling algorithm for improving RNA model refinement},
  author={Xiong, Peng and Wu, Ruibo and Zhan, Jian and Zhou, Yaoqi},
  journal={Nature communications},
  volume={12},
  number={1},
  pages={2777},
  year={2021},
  publisher={Nature Publishing Group UK London}
}
```
