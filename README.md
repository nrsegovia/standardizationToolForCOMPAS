# standardizationToolForCOMPAS

Python tools to process COMPAS .h5 files and build a standardized output

## Instructions

Run the python script alongside your COMPAS' output path, where the `.h5` and `Run_Details` files should be stored. Note that this assumes that your COMPAS logfile definitions includes (at least) what is shown in the `sample-logfile.txt` file.

```bash
python buildStandardOutput.py path/to/compas/output/directory
```

The standardized output will be created in the same COMPAS' output path.

Note: this script automatically skips systems where the simulation ended due to an error. Those systems are shown in the terminal whenever COMPAS finds one. For example:

```text
ERROR:  in BaseBinaryStar::EVOLUTION_STATUS BaseBinaryStar::Evolve(): Initial attributes are not valid - evolution not possible, (ObjectId = 323)
```

`ObjectId` is directly related to the line number in the grid file, so it might be a good idea to store it in case you need it later.
