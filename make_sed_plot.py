#makes a spectral energy distribution from a given set of photometry

#imports
import matplotlib
import matplotlib.pyplot as mp
import pandas as pd
import os
import numpy as np


galaxy = 'NGC4656'

os.chdir('/Users/kzs1031/Desktop/lvls/'+galaxy) #use placeholder gal

mags_file = 'flat_spec_corr.txt'
df = pd.read_csv(mags_file, delim_whitespace=True, header=0)

#separate out photometry

df_phot = df.iloc[:, 6:]

#central wavelengths for each filter based on transmission curve and flat spectrum 
#in microns
lambda_cen = {'f_galex_fuv':0.154, 'f_galex_nuv':0.231, 'f_2MASS_J': 1.250, 'f_2MASS_H':1.650, 'f_2MASS_Ks':2.170, 'f_iracch1':3.56, 'f_iracch2':4.5, 'f_iracch3':5.74, 'f_iracch4':7.92, 'f_MIPS24um':23.84, 'f_MIPS70um':72.56,
                'f_MIPS160um':156.96, 'f_johnson_U':0.365, 'f_johnson_B':0.445, 'f_johnson_V':0.551, 'f_johnson_R':0.658, 'f_sloan_u':0.356, 'f_sloan_g':0.472, 'f_sloan_r':0.618, 'f_sloan_i':0.750, 'f_sloan_z':0.896, 'f_swift_uvw2':0.214, 
              'f_swift_uvm2':0.227, 'f_swift_uvw1': 0.267} #x vals , 

f_dens = df_phot.to_dict(orient='list')

fluxes = []
#mwcorr_mags = []
for key in list(lambda_cen.keys()):
    if f_dens[key][0] != -99:
        flux = f_dens[key] 
        fluxes.append(flux[0])

    else:
        lambda_cen.pop(key)


#plot SED
mp.scatter(list(lambda_cen.values()), fluxes)
mp.title('Spectral energy distribution of '+galaxy)
mp.xlabel('wavelength (microns)')
mp.ylabel('flux density')
mp.xscale('log')
mp.yscale('log')
mp.savefig('gal_sed.png')
mp.show()