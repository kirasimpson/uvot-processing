#!/usr/bin/env python
# coding: utf-8

"""This script corrects the fluxes in the final_phot_table file for 
MW foreground extinction.  

In order to properly correct for the effects of dust, the input spectrum
must be known. The input spectrum and the filter transmission curve effects
the amount of extinction to correct for. Since we do not know a priori what
the true spectrum looks like, this correction will use a flat input spectrum. 
A more proper way to correct for the effects of dust would be to create model
galaxy spectra from BC03 and calculate what A_lambda/A_V is as a function of 
A_V = R_V E(B-V), where the color excess is known from Schlafly et al 2016. 
"""




import numpy as np
import pandas as pd
import matplotlib
import glob 
import scipy.integrate as integrate
import astropy.units as u
import os
import argparse
# this is the package I use for the dust extinction correction
# see link: https://dust-extinction.readthedocs.io/en/latest/index.html
from dust_extinction.parameter_averages import F99

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


def run(galaxy):


    os.chdir('../'+galaxy)



    # this is a file with the uncorrected magnitudes in each filter
    uncor_file = 'galphot_table.txt'
    df = pd.read_csv(uncor_file, delim_whitespace=True, header=0)



    # first split data frame into photometry and non phot pieces 

    df_phot = df.iloc[:, 6:]
    df_nonphot = df.iloc[:, :6]


    #my filters are: FUV  NUV F1.25(J)

    #F1.65(H) F2.17(K_s) F3.6 F4.5 F5.8 F8.0 F24 F70 F160 U B V R u g r i z

    #separate out different columns into different info flags, errors, and magnitudes
    extra =  ['l_', 'F-N']
    lims = [col for col in df_phot.columns if any(ent in col for ent in extra)]
    df_phot = df_phot.drop(lims, axis=1)


    err_cols = [col for col in df_phot.columns if 'e_' in col]


    df_err = df_phot[err_cols]



    df_flux = df_phot.drop(err_cols, axis=1)

    #for index, row in df_flux.iterrows():
     #   print(index, row)
    #input()

    df_flux = df_flux.nan_to_num(np.nan, nan=-99)




    #=================

    #putting a placeholder over a NaN so I can use delimiter in later script to make plots
    #df_mag = df_flux.replace(np.nan, -99)
    df_flux_f = pd.DataFrame(df_flux)




    df_flux_f



    # get central wavelengths for each filter based on transmission curve and flat spectrum 
    #in microns
    lambda_cen = {'galex_fuv':0.154, 'galex_nuv':0.231, 'johnson_U':0.365, 'johnson_B':0.445, 'johnson_V':0.551, 'johnson_R':0.658, 'sloan_u':0.356, 'sloan_g':0.472, 'sloan_r':0.618, 'sloan_i':0.750, 'sloan_z':0.896, 
                    '2MASS_J': 1.250, '2MASS_H':1.650, '2MASS_Ks':2.170, 'iracch1':3.56, 'iracch2':4.5, 'iracch3':5.74, 'iracch4':7.92, 'MIPS24um':23.84, 'MIPS70um':72.56,
                    'MIPS160um':156.96, 'swift_uvw2':0.214, 'swift_uvm2':0.227, 'swift_uvw1': 0.267} #, 

    # lambda_cen in microns 
    # lambda_cen uses values from Dale et al paper and lambda mean from SVO for swift filters

    #array of central wavelengths
    lc_array = list(lambda_cen.values())

    lc_array = np.array(lc_array)


    #lambda_cen_a = np.array(lambda_cen)*1e4    # to angstroms 




    # for every other columns/filters
    #     can calculate A_lamb_A_V from chosen extinction law
    #     F99 or ccm 89
    #     for each galaxy
    #         get E(B-V) and calculate A_V
    #         calculate corrected flux from A_V and measured flux
    #         calculate new error
    #    save new fluxes and errors 




    # initialize the dust model F99 and convert to inverse wavelength (wavenumbers) (see docs)
    dust_model = F99()
    inv_wave = 1/(np.array(lc_array))


    # loop through and get the A_lambda/A_V at each wavelength
    # since F99 is only defined on a certain wavlength range 
    # (approx lyman alpha to ~3 microns) if the central wavelength
    # falls outside limit, assume A_lambda/A_V=0 since extinction is
    # negligible

    ext = []
    for i in inv_wave:
        try:
            dust_model(i)
            ext.append(dust_model(i))
        except ValueError as e:
            ext.append(0)

    # from read in photometry, 
    # make a dataframe of errors and a data frame of fluxes


    cor_flux = []
    cor_mag = []
    header = []


    # for each row
    #    get E(B-V) and calculate A_V
    #    from extincted flux and A_V, calculated corrected flux
    #    save row in new data frame
    for ind, row in df.iterrows():

        ebv = row['E(B-V)']                 

        av = 3.1*ebv
        ax = av*np.array(ext) #E(B-V) extinction * 3.1 (reddening coefficient) * array of A_lambda/A_V values

        col = df_flux.columns 
        f = df_flux.iloc[ind].values #flux values


        f_corr = f*10**(ax/2.5) 


        # now correct errors
        err = df_err.iloc[ind].values

        err_corr = err * f_corr/f



        # #=======
        
        data = np.concatenate((f_corr, err_corr), axis=None) #array


        cor_flux.append(data) #list with array in it


                                 

        if ind==0: #makes the header for a new table
            flux_name = df_flux.columns.values
            err_name = df_err.columns.values
            head = np.concatenate((flux_name, err_name), axis=None)
            header.append(head)






    df_corr = pd.DataFrame(cor_flux)



    df_corr.columns = header[0]




    df_corr




    df_final = pd.concat([df_nonphot, df_corr], sort=False, axis=1)


    #df_final = df_final.replace('--', -99, regex=True)

    # save dataframe 
    df_final.to_csv('flat_spec_corr.txt', sep=' ', index=False)

#=====================================


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("galaxy",nargs = "?", help = "Enter the name of a galaxy", default = '')

    args = parser.parse_args()

    run(args.galaxy)

if __name__ == '__main__':
    main()






