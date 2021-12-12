# Protein-aggregation-vs.-subcellular-location

Add the following raw data to the folder **Code>Data**:
- **human_proteome_df.RData**
-  **domains.RData**

Download all human AlphaFold structures:
- Go to [the AlphaFold download page](https://www.alphafold.ebi.ac.uk/download) and download the entire human dataset.
- Save it to the **Code>Data>AF** folder and untar it.

The following processed data files can be added to speed up the process:
-  **complete_data.RData** and **hashed_proteins.RData** (human_proteome_df with annotated subcellular locations and secreted attribute)
-  **domains_with_COandTango.RData** or the intermediate file domains_with_contact_order.RData (domains with annotated contact order and tango scores)
