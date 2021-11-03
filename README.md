# Deduper

## Explanation
Anytime you will be using PCR amplification as a step in a library prep, you have the chance of creating PCR duplicates. This will inhibit downstream analyses by creating false positives/over representation of a certain transcript. This can inhibit the accuracy of how much a gene is expressed. This tool can be used to remove all PCR duplicates within a given SAM file. By default this will return the first instance read encountered and will not store duplicates. This means it does not take into account the duplicates quality score. If you would like quality to be taken into consideration please specify it with a flag, as described below.

## How to Use
1. Please sort the SAM file by chromosome before running it through the script. If you do not, then you will not get correct results.
2. I would create a conda environment using the command `conda create -f conda_env/deduper_env.yml`. This will create a conda environment that will have all of the necessary packages installed to run Deduper. Activate this environment.
3. To run the file please type `python python_scripts/powers_deduper.py` followed by flags for the options that you are interested in using.
4. Look into the `output` folder for your deduped file. It will have the postfix `_deduped.sam`

## Flags
```-f``` Specify the path to the sorted SAM file that you will want to dedupe.

```-p``` Tell the script if the SAM file has paired ends or single end data. This is set False by default, meaning that it is assuming the data is single end reads. (Not currently finished, you will get an error if you pass True to this argument)

```-u``` Specify the path the the known UMIs that were used when constructing the library. This UMI file should have every UMI separated by newlines. If unknown UMIs were used, then do not pass anything to this argument. It will automatically default to 'random'. This will build a set of UMI from the headers in the SAM file. It will still discard any UMIs/reads that have Ns in the unknown UMI.

```-ds``` If you would like the duplicates to be stored in a separate file set this to True.

```-do``` If you set ```-ds``` to True then you need to specify a file name for them. If you do not then you will get an error.

```-q``` Returns only the highest quality read instead of the first instance read. Set this to True if you want the best quality read returned.


## Example Default Run

```
python python_scripts/powers_deduper.py -f <file> -u <known UMI file>
```
