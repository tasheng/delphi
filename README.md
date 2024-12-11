# delphi
Converting the delphi open data to TPCNtuple root format

## Converting DST to ROOT format
```
### setup to the DELPHI environment
source setup.sh
# TPCNtupleFormat/DataProcessing is based on janicechen git@github.com:janice-cat/StudyMult.git --depth 1
# we will have to re-define things in the event selection and track selection classes
### converting
./execute
```

---
# From 
## Event picking

This job loops the particles in each event and prints out some information in plain text.
Application: bit wise validation of the software stack.

A maximum of 10 events is dumped.

## Bugs
The code may need some cleanup.