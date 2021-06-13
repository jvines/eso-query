# eso-query
A quick script that queries the ESO raw archive for spectroscopic observations made with ESPRESSO, HARPS or FEROS.

Usage:

eso_query.py ra dec radius n output

The ra and dec are expected to be in degrees, while the radius is expected to be in arcminutes, n is the max number of rows you want to save and output is the output directory where a dat file will be saved with the results.

For example if you wanted to query a 2 arcmin box around NGTS-6 and save the results in a file called ngts-6.dat:

```bash
python eso-query.py 75.7955654461 -30.399375656 2 2000 ngts-6
```
