# NDS - Nuclear Data Source

NDS is a command-line tool for easily querying nuclear data. It provides
isotope-specific information, including atomic/molar masses, radioactivity
data (half-life, decay mode, energy spectra), and cross-section data.

The tool permits direct querying of information or executing Python 
expressions/files with a querying API.

**nds n \<nuclide\>** - 
Queries general isotope/element data. Nuclide format is either the symbol 
(and optionally an isotope number, e.g. "Co" or "Co60") or an MCNP-style
definition with element number followed by isotope number ("27000" or 
"27060" for Cobalt and Cobalt-60, respectively).

### Data Sources

* NUBASE2020
* Generic Periodic Table Data
* MCNP Atomic Weight Ratio (AWR) Data
* ENDF VIII.0 Decay Data

```
/data/
    awr_data
    ENDF-B_VIII.0_decay.zip
    nubase2020
    periodic.csv
```