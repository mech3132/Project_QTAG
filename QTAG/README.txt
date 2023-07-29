############ Quantification of Turnover Across Gradients (QTAG) ####################





####### Before running #######

Make sure you have all modules required by moving into the QTAG directory and running:

pip install -r requirements.txt





####### To run #######

Make sure you have python3 installed. 
To see all options for QTAG, run:

python QTAG.py --help

Here is an example (files should be found in the 'examples' directory):

# Make sure you are in the QTAG directory
python QTAG.py \
-t examples/otu-table.txt \
-m examples/sample-metadata.txt \
-M Salinity 





####### Input files required #######

OTU/ASV table
--------taxa x sample (row x col)
--------header row MUST start with #OTU ID
--------rows preceding #OTU ID starting with # will be ignored
--------tab delimited

Sample metadata table
--------sample x variable (row x col)
--------first column must be sample IDs that match column names in OTU/ASV table. Name of column does not matter.
--------must have at least one additional column which is your "gradient" variable name (metadata_name)





