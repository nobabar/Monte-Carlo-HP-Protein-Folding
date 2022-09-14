# Monte Carlo algorithm for protein folding in the HP model

## Setup your environment

Clone the repository:

```bash
git clone git@github.com:nobabar/Monte-Carlo-HP-Protein-Folding.git
```

Move to the new directory:

```bash
cd Monte-Carlo-HP-Protein-Folding
```

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Create the `monte-carlo` conda environment:

```
conda env create -f env.yml
```

Load the `monte-carlo` conda environment:

```
conda activate monte-carlo
```

## Run the program

### Main command 'fold.py'

```bash
python fold.py [-h] (-p PROTEIN | -f FILE) [-i {linear,random}] {MC,REMC} ...
```

| positional arguments |                                               |
| -------------------- | --------------------------------------------- |
| {MC,REMC}            | The algorithm to use.                         |
|                      | MC: Monte Carlo algorithm.                    |
|                      | REMC: Replica Exchange Monte Carlo algorithm. |

| options                                    |                                            |
| ------------------------------------------ | ------------------------------------------ |
| -h, --help                                 | show this help message and exit            |
| -p PROTEIN, --protein PROTEIN              | input protein sequence                     |
| -f FILE, --file FILE                       | input file containing the protein sequence |
| -i {linear,random}, --init {linear,random} | initial configuration of the protein       |

### Sub-command 'MC'

```bash
... MC [-h] [-n N_STEPS] [-t TEMPERATURE]
```

| options                                   |                                    | default |
| ----------------------------------------- | ---------------------------------- | ------- |
| -h, --help                                | show this help message and exit    |         |
| -n N_STEPS, --n-steps N_STEPS             | number of iterations in the search | 1000    |
| -t TEMPERATURE, --temperature TEMPERATURE | temperature of the system          | 200     |

### Sub-command 'REMC'

```bash
... REMC [-h] [-n N_REPLICA] [-e ENERGY_CUTOFF] [-m MAX_STEPS] [-l LOCAL_STEPS] [-tmin TEMPERATURE_MIN] [-tmax TEMPERATURE_MAX]
```

| options                                         |                                                  | default |
| ----------------------------------------------- | ------------------------------------------------ | ------- |
| -h, --help                                      | show this help message and exit                  |         |
| -n N_REPLICA, --n-replica N_REPLICA             | number of replicas to use                        | 10      |
| -e ENERGY_CUTOFF, --energy-cutoff ENERGY_CUTOFF | optimal energy to reach                          | -10     |
| -m MAX_STEPS, --max-steps MAX_STEPS             | maximum number of steps if cutoff is not reached | 1000    |
| -l LOCAL_STEPS, --local-steps LOCAL_STEPS       | number of steps to perform for each MC search    | 100     |
| -tmin TEMPERATURE_MIN, --temperature-min        | temperature of the first replica                 | 160     |
| -tmax TEMPERATURE_MAX, --temperature-max        | temperature of the last replica                  | 220     |

## Usage examples

### Monte Carlo algorithm

```bash
python fold.py -f data/1AF5.fasta -i random MC -n 5000 200
```

### Replica Exchange Monte Carlo algorithm

```bash
python fold.py -p HPHPPHHPHPPHPHHPPHPH REMC -n 5 -e -9
```

## Benchmark proteins

| ID  | Len | E^\* | Protein Sequence                                                                                     |
| --- | --- | ---- | ---------------------------------------------------------------------------------------------------- |
| S1  | 20  | -9   | HPHPPHHPHPPHPHHPPHPH                                                                                 |
| S2  | 24  | -9   | HHPPHPPHPPHPPHPPHPPHPPHH                                                                             |
| S3  | 25  | -8   | PPHPPHHPPPPHHPPPPHHPPPPHH                                                                            |
| S4  | 36  | -14  | PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP                                                                 |
| S5  | 48  | -23  | PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPPPHHPPHHHPPHHHHH                                                     |
| S6  | 51  | -21  | HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHPHHHHHPHPHPHH                                                  |
| S7  | 60  | -36  | PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP                                         |
| S8  | 64  | -42  | HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH                                     |
| S9  | 85  | -53  | HHHHPPPPHHHHHHHHHHHHPPPPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHHHHHHHHHHHHPPPHPPHHPPHHPPHPH                |
| S10 | 100 | -50  | PPPHHPPHHHHPPHHHPHHPHHPHHHHPPPPPPPPHHHHHHPPHHHHHHPPPPPPPPPHPHHPHHHHHHHHHHHPPHHHPHHPHPPHPHHHPPPPPPHHH |
| S11 | 100 | -48  | PPPPPPHPHHPPPPPHHHPHHHHHPHHPPPPHHPPHHPHHHHHPHHHHHHHHHHPHHPHHHHHHHPPPPPPPPPPPHHHHHHHPPHPHHHPPPPPPHPHH |
