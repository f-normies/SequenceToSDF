# SequenceToSDF

This script make SDF file from CSV file with protein/DNA/RNA sequences and their descriptions.

## Usage

<code>
python SeqToSDF.py path_to_config_file
</code>

Config is JSON file with script parameters. Here is description.

```json
{
  "input": "path to csv", 
  "output": "path to output directory", 
  "column": "column name with sequences",
  "charged": true, 
  "alphabet": "protein", 
  "threads": 10, 
  "separator": ";", 
  "filename": "1000" 
}
```

## About parameters

1. Required parameters are "input", "output", "column".
2. "Charged" works only  for protein sequences, otherwise it is ignored.
3. Alphabet can take following values: "protein", "DNA", "RNA".
4. Threads is number between 1 and 2 * number of cores.
5. Separator is CSV separartor. 
6. Filename is custom name for output sdf, otherwise it will be named by csv file.