# gesar-file-utility

Utility for finding where the dataset for GESAR experiment originates from.

* To start the process, add VT_unfiltered_Feb.txt file to the data_files directory.
* Add un-zipped files from 48hd cloud for Feb 2018 with VT, VV, VG initials to data_files/unzipped_files directory.
    * Ensure that the files have been converted to .txt format.
* Run clean_48_zip_files from helper_functions.py once to remove the top few rows from the 48hd files.
* Install the requirements from requirements.yml file (you can use conda to create a separate environment)
* Run the two jupyter notebooks (compare_aa.ipynb, compare_nuc.ipynb) to run the rest of the program and regenerate the results we obtained.