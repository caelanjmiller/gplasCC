<div align="center"><img src="https://gitlab.com/mmb-umcu/gplascc/-/raw/master/figures/gplasCClogo.png?ref_type=heads" alt="gplasCC" width="600"/></div>

gplasCC: binning plasmid-predicted contigs
================

GplasCC is a tool to bin plasmid-predicted contigs based on sequence
composition, coverage and assembly graph information. 
GplasCC is a new version of [gplas](https://gitlab.com/sirarredondo/gplas) that allows for plasmid classification of any binary plasmid classifier and extends the possibility of accurately binning predicted
plasmid contigs into several discrete plasmid components by also attempting to place unbinned and repeat contigs into plasmid bins.

# Table of Contents
- [gplasCC: binning plasmid-predicted contigs](#gplascc-binning-plasmid-predicted-contigs)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  - [Requirements](#requirements)
  - [Installation using pip](#installation-using-pip)
- [Usage](#usage)
    - [Using gplasCC with plasmidCC](#using-gplascc-with-plasmidcc)
    - [Using gplasCC with an external classification tool](#using-gplascc-with-an-external-classification-tool)
- [Output files](#output-files)
    - [Intermediary results files](#intermediary-results-files)
- [Complete usage](#complete-usage)
- [Issues and Bugs](#issues-and-bugs)
- [Contributions](#contributions)
- [Citation](#citation)

# Installation

## Requirements
An installation of [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) is required if you are using [plasmidCC](https://gitlab.com/mmb-umcu/plasmidCC) as binary classifier (default). We reccomend using a conda environment with the [centrifuge-core](https://bioconda.github.io/recipes/centrifuge-core/README.html) package installed.

``` bash
conda create --name gplasCC -c conda-forge -c bioconda centrifuge-core=1.0.4.1 pip
conda activate gplasCC
```
If you prefer to use a different binary classifier, you can use gplasCC without installing Centrifuge.

## Installation using pip

The prefered way of installing gplasCC is through pip:

``` bash
pip install gplas
```
When this has finished, test the installation using 
``` bash
gplas --help
```
This should should show the help page of gplasCC.

# Usage

### Using gplasCC with plasmidCC

GplasCC comes built in with [plasmidCC](https://gitlab.com/mmb-umcu/plasmidCC) as a binary classifier. When using plasmidCC, gplasCC only requires one input file:

* An assembly graph in **.gfa** format. Such an assembly graph can be obtained by assembling quality trimmed reads using [Unicycler](https://github.com/rrwick/Unicycler) (preferred) or with [SPAdes genome assembler](https://github.com/ablab/spades).

Provide the path to your assembly graph with the **-i** flag, and select which plasmidCC database to use with the **-s** flag. Optionally, provide a custom name for your output with the **-n** flag. See example below: 

``` bash
gplas -i test_ecoli.gfa -s Escherichia_coli -n my_isolate
```
For an overview of [plasmidCC](https://gitlab.com/mmb-umcu/plasmidCC) supported species, use the **--speciesopts** flag:
``` bash
gplas --speciesopts
```

### Using gplasCC with an external classification tool

If you wish to use a different binary classifier, it is possible to provide your own external plasmid prediction file. We've listed and reviewed several other classifier tools [here](https://www.mdpi.com/2076-2607/9/8/1613). Although they are all compatible with gplasCC, extra preprocessing steps are required:

1) Use gplasCC to convert the nodes from the assembly graph to FASTA format (most binary classifiers only accept FASTA files as input). To do this, provide your assembly graph (**.gfa**) and include the **--extract** flag.

``` bash
gplas -i test_ecoli.gfa --extract -n my_isolate
```

The output FASTA file will be located in: __gplas_input/__*my_isolate*_contigs.fasta. By default, this file will only contain contigs larger than 1000 bp, however, this can be controlled with the **-l** flag. 

2) Use this FASTA file as an input for the binary classification tool of your choice. 

3) Format the output file: 

The output from the binary classification tool has to be formatted as a **tab separated** file containing specific columns and headers (**case sensitive**). See an example below:

``` bash
head -n 4 test_ecoli_plasmid_prediction.tab
```

| Prob\_Chromosome | Prob\_Plasmid |  Prediction  | Contig\_name                             | Contig\_length|
|-----------------:|--------------:|:-------------|:-----------------------------------------|--------------:|
|       1.0       |      0.0     |  Chromosome  |  S1\_LN:i:374865\_dp:f:1.0749885035087077   |      374865     |
|       1.0       |      0.0     |  Chromosome  | S10\_LN:i:198295\_dp:f:0.8919341045340952  |     198295    |
|       0.0       |      1.0     |    Plasmid   |  S20\_LN:i:91233\_dp:f:0.5815421095375989   |      91233     |

For proper compatability with gplasCC, please make sure your prediction file is **tab-separated**, and uses the correct (**case sensitive**) column names and prediction labels (Plasmid/Chromosome).

Once you've formatted the prediction file as above, move to [Predict plasmids](#predict-plasmids).

#### Predict plasmids <a name="predict-plasmids"></a>
After pre-processing, we are now ready to predict individual plasmids. 

Provide the paths to your assembly graph, using the **-i** flag, and to your binary classification file, with the **-P** flag. Optionally, provide a custom name for your output with the **-n** flag. See example below: 

``` bash
gplas -i test_ecoli.gfa -P test_ecoli_plasmid_prediction.tab -n my_isolate
```

# Output files

GplasCC will create a folder called ‘results’ with the following files:

``` bash
ls results/my_isolate*
```

    ## results/my_isolate_bin_0.fasta
    ## results/my_isolate_bin_1.fasta
    ## results/my_isolate_bin_2.fasta
    ## results/my_isolate_bins.tab
    ## results/my_isolate_chromosome_repeats.tab
    ## results/my_isolate_plasmidome_network.png
    ## results/my_isolate_results.tab

##### results/\*.fasta

Fasta files with the contigs belonging to each predicted plasmid bin.

``` bash
grep '>' results/my_isolate*.fasta
```

``` bash
>S20_LN:i:91233_dp:f:0.5815421095375989
>S1_LN:i:374865_dp:f:1.0749885035087077
>S32_LN:i:42460_dp:f:0.6016122804021161
>S44_LN:i:21171_dp:f:0.5924640018897323
>S47_LN:i:17888_dp:f:0.5893320957724726
>S48_LN:i:11703_dp:f:1.1884320594277211
>S50_LN:i:11225_dp:f:0.6758514700227541
>S56_LN:i:6837_dp:f:0.5759570101860518
>S59_LN:i:5519_dp:f:0.5544497698217399
>S67_LN:i:2826_dp:f:0.6746421335091037
>S70_LN:i:2125_dp:f:9.215759397832965
>S76_LN:i:1486_dp:f:1.3509551203209675
>S84_LN:i:1063_dp:f:3.2697611578099566
```

##### results/\*bins.tab

Tab delimited file containing a short overview showing the contigs that got assigned to each plasmid bin.

| number | Bin |
| ------ | --- |
| 1      | 1   |
| 20     | 0   |
| 32     | 1   |
| 44     | 1   |
| 47     | 1   |
| 48     | 1   |
| 50     | 1   |
| 56     | 1   |
| 59     | 1   |
| 67     | 1   |
| 70     | 1   |
| 76     | 1   |
| 84     | 1   |

##### results/\*results.tab

Tab delimited file containing the classification given by plasmidCC (or other binary classification tool) together with the bin prediction from gplasCC. The file contains
the following information: contig number, contig name, probability of being
chromosome-derived, probability of being plasmid-derived, class
prediction, length, k-mer coverage, assigned bin.

| Prob\_Chromosome | Prob\_Plasmid | Prediction | Contig\_name                             | number | length | coverage | Bin |
| ---------------- | ------------- | ---------- | ---------------------------------------- | ------ | ------ | -------- | --- |
| 1.0              | 0.0           | Repeat     | S1_LN:i:374865_dp:f:1.0749885035087077   | 1      | 374865 | 1.07     | 1   |
| 1.0              | 0.0           | Repeat     | S48_LN:i:11703_dp:f:1.1884320594277211   | 48     | 11703  | 1.19     | 1   |
| 0.5              | 0.5           | Repeat     | S70_LN:i:2125_dp:f:9.215759397832965     | 70     | 2125   | 9.22     | 1   |
| 0.0              | 1.0           | Repeat     | S76_LN:i:1486_dp:f:1.3509551203209675    | 76     | 1486   | 1.35     | 1   |
| 0.78             | 0.22          | Repeat     | S84_LN:i:1063_dp:f:3.2697611578099566    | 84     | 1063   | 3.27     | 1   |
| 0.0              | 1.0           | Plasmid    | S20_LN:i:91233_dp:f:0.5815421095375989   | 20     | 91233  | 0.58     | 0   |
| 0.0              | 1.0           | Plasmid    | S32_LN:i:42460_dp:f:0.6016122804021161   | 32     | 42460  | 0.6      | 1   |
| 0.0              | 1.0           | Plasmid    | S44_LN:i:21171_dp:f:0.5924640018897323   | 44     | 21171  | 0.59     | 1   |
| 0.0              | 1.0           | Plasmid    | S47_LN:i:17888_dp:f:0.5893320957724726   | 47     | 17888  | 0.59     | 1   |
| 0.0              | 1.0           | Plasmid    | S50_LN:i:11225_dp:f:0.6758514700227541   | 50     | 11225  | 0.68     | 1   |
| 0.0              | 1.0           | Plasmid    | S56_LN:i:6837_dp:f:0.5759570101860518    | 56     | 6837   | 0.58     | 1   |
| 0.0              | 1.0           | Plasmid    | S59_LN:i:5519_dp:f:0.5544497698217399    | 59     | 5519   | 0.55     | 1   |
| 0.0              | 1.0           | Plasmid    | S67_LN:i:2826_dp:f:0.6746421335091037    | 67     | 2826   | 0.67     | 1   |

##### results/\*chromosome_repeats.tab

Tab delimited file showing which contigs got assigned as chromosomal repeats.
| number | Bin |
| ------ | --- |
| 1      | Chromosome |
| 48     | Chromosome |
| 55     | Chromosome |
| 66     | Chromosome |
| 68     | Chromosome |
| 70     | Chromosome |
| 74     | Chromosome |
| 79     | Chromosome |
| 81     | Chromosome |
| 84     | Chromosome |

##### results/\*plasmidome\_network.png

A visual representation of the plasmidome network generated by gplasCC. The network is created using an undirected graph with edges between plasmid unitigs co-existing in the random walks created by gplasCC.

<div align="center"><img src="https://gitlab.com/mmb-umcu/gplascc/-/raw/master/figures/my_isolate_plasmidome_network.png?ref_type=heads" alt="plasmidome_network" width="600"/></div>

### Intermediary results files

If the **-k** flag is selected, gplasCC will also **keep** all intermediary files needed to construct the plasmid predictions. For example:

##### walks/normal_mode/\*solutions.tab

gplasCC generates plasmid-like walks for each plasmid starting node. These paths are later used to generate the edges of the plasmidome network, but they can also be useful to observe all the different walks starting
from a single node (plasmid unitig). These walks can be directly given to Bandage to visualize and manually inspect a walk.

In the example below, we find different possible plasmid walks starting from the node 67-. These paths may contain inversions and rearrangements since repeats units, such as transposases, can be present several times within the same plasmid sequence. In these cases, gplasCC can traverse the
sequence in different ways generating different plasmid-like paths.

``` bash
tail -n 10 walks/normal_mode/my_isolate_solutions.tab
```

``` bash
67-,70-,50-,143-
67-,70-,50-,143-
67-,70-,50-,143-
67-,70-,47+,117-,84-,59+,70-,50-,143-
67-,70-,50-,143-
67-,70-,50-,143-
67-,70-,47+,117-,84-,59+,70-,50-,143-
67-,70-,47+,117-,84-,59+,70-,50-,143-
67-,70-,50-,143-
67-,70-,50-,143-
```

We can use Bandage to inspect the following path on the assembly graph:
67-,70-,47+,117-,84-,59+,70-,50-,143-

<div align="center"><img src="https://gitlab.com/mmb-umcu/gplascc/-/raw/master/figures/bandage_path.jpg?ref_type=heads" alt="bandage_path" width="600"/></div>

# Complete usage

``` bash
gplas --help
```
``` bash
usage: gplas -i INPUT [-n NAME]
             (-s SPECIES | -p CUSTOM_DB_PATH | -P PREDICTION | --extract)
             [-t THRESHOLD_PREDICTION] [-b BOLD_COVERAGE_SD]
             [-x NUMBER_ITERATIONS] [-f FILT_GPLAS] [-e EDGE_THRESHOLD]
             [-q MODULARITY_THRESHOLD] [-l LENGTH_FILTER] [-k]
             [--speciesopts] [-v] [-h]

gplasCC: A tool for binning plasmid-predicted contigs into individual
predictions

General:
  -i INPUT              Path to the graph file in GFA (.gfa) format, used
                        to extract nodes and links
  -n NAME               Name prefix for output files (default: input file
                        name)
  -s SPECIES            Choose a species database for plasmidCC
                        classification. Use --speciesopts for a list of
                        all supported species
  -p CUSTOM_DB_PATH     Path to a custom Centrifuge database (name without
                        file extensions)
  -P PREDICTION         If not using plasmidCC. Provide a path to an
                        independent binary classification file
  --extract             extract FASTA sequences from the assembly graph to
                        use with an external classifier

Parameters:
  -t THRESHOLD_PREDICTION
                        Prediction threshold for plasmid-derived sequences
                        (default: 0.5)
  -b BOLD_COVERAGE_SD   Coverage variance allowed for bold walks to
                        recover unbinned plasmid-predicted nodes (default:
                        5)
  -x NUMBER_ITERATIONS  Number of walk iterations per starting node
                        (default: 20)
  -f FILT_GPLAS         filtering threshold to reject outgoing edges
                        (default: 0.1)
  -e EDGE_THRESHOLD     Edge threshold (default: 0.1)
  -q MODULARITY_THRESHOLD
                        Modularity threshold to split components in the
                        plasmidome network (default: 0.2)
  -l LENGTH_FILTER      Filtering threshold for sequence length (default:
                        1000)

Other:
  -k, --keep            Keep intermediary files

Info:
  --speciesopts         Prints a list of all supported species for the -s
                        flag
  -v, --version         Prints gplas version
  -h, --help            Prints this message

```

# Issues and Bugs
You can report any issues or bugs that you find while installing/running
gplasCC using the [issue tracker](https://gitlab.com/mmb-umcu/gplascc/-/issues).

# Contributions
GplasCC has been developed with contributions from Oscar Jordan, Julian Paganini, Jesse Kerkvliet, Malbet Rogers, Sergio Arredondo and Anita Schürch.

# Citation
If you use gplasCC in your study, please cite:
https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
