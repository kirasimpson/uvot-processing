# uvot-processing
Main code and supporting scripts related to the UVOT data reduction pipeline.

This doc lays out the process for producing photometry, plotting SEDs, and fitting SEDs with MCSED.

The project is organized such that each galaxy has a directory containing all the data from HEASARC,
as well as all data products created by the projects. 

To do photometry and make the MCSED input:

- run the pipeline using the lvls_pipeline.py. This will:
    - check for new data in swift uv filters and download new stuff
    - stack and align the images
    - build the deep image
    - make a rgb data cube
- surface_phot.py (Lea’s tool)
    - does surface photometry
- phot_plot.py (Lea’s tool)
    - makes plots of surface brightness photometry
- make_phot_table.py
    - stitches together all the photometry from archival sources/extracted from Swift data
    - sources for the tables are listed in the documentation
    - converts all values listed as magnitudes into flux densities
- extinction_correction.py
    - some changes were added to sort out the photometry from other information.
- make_sed_plot.py
    - good sanity check, makes an SED plot .png from the MW extinction-corrected data.

From there you’re mostly good to start running MCSED. There are a couple of changes that will have
to be made to the output table, either manually or by editing make_phot_table.py:
    - change NaNs into -99s for the final table to match MCSED's desired inputs.
    - change optical/galex UV names to match MCSED’s .res files.
        - i.e.: f_U -> f_johnson_U, f_u -> f_sloan_u
