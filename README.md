# CLas_Phylogenomics
Collection of scripts for a range of genomic data processing and analysis.
Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Download Genbank assemblies](#Download-Genbank-assemblies)

## Download Genbank assemblies
The script `GenBankScan.py` automatically downloads [Genbank](https://www.ncbi.nlm.nih.gov/genbank/) assemblies for a given search term (e.g. a species name).

#### Example command
`python GenBankScan.py --term 'Candidatus Liberibacter asiaticus' --format 'GenBank' --email 'other@example.com''`

`python GenBankScan.py -h` Will print a full list of command arguments.

#### Command arguments
| Name | Description |
| :--: | ----------- | 
| `term` | Search term, usually organism name (required) |
| `format`  | Accession format (RefSeq or GenBank) (required) |
| `email`  | Provide your own email address here (required) |

#### Notes
To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [Bio](https://biopython.org/), [urllib](https://docs.python.org/3/library/urllib.html), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

___
## Citation

Please cite this when using these scripts:
Pruvost O., Boyer K., Labbé F., Weishaar M., Vynisale A., Melot C., Hoareau C., Cellier G., and Ravigné V., 2024. Genetic signatures of contrasted outbreak histories of ‘*Candidatus* Liberibacter asiaticus’, the bacterium that causes citrus huanglongbing in three outermost regions of the European Union. (*in prep*)
