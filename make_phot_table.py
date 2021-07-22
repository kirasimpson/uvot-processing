from astropy.table import Table, QTable, Column
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import os.path
import glob
import argparse

import pdb


'''
Combine all the photometry into a table and write to a .txt file


Arguments:
- object name (str)
	- example command: 'python gal_table_se.py NGC4242'

Tables used for archival photometry can be found here:

Optical photometry from Cook+8:
https://academic.oup.com/mnras/article/445/1/881/984908#92778246
Apertures and UV magnitudes from Lee+11:
https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t2_mrt.txt
Infrared fluxes from Dale+23:
https://iopscience.iop.org/0004-637X/703/1/517/suppdata/apj314922t2_mrt.txt


Output is dictionary containing flux densities of all filters given
for a particular galaxy, which writes as a table to a .txt file
'''



#in the dictionary created by this script, replace the filter names in the tables with the names that match the .res files used by MCSED
#keys are the table filter names, values are the .res file filter names
replace = {'Fapmag':'galex_fuv', 'Napmag':'galex_nuv', 'U':'johnson_U', 'B':'johnson_B', 'V':'johnson_V', 'R':'johnson_R', 'u':'sloan_u',
			'g':'sloan_g', 'r':'sloan_r', 'i':'sloan_i', 'z':'sloan_z'}



def run(galaxy): #takes one galaxy name as a string 
	#makes a dictionary for all filter flux densities

	#initialize dictionary
	gal_dict = {galaxy:{}}


	#==============================================================
	#table that reports apertures and NUV/FUV aperture mags
	aperture_uvmags_table = Table.read('uvmags_yesheader.dat', format='ascii.cds')

	#split into photometry and other info
	nonphot = aperture_uvmags_table.colnames[0:6]
	phot = aperture_uvmags_table.colnames[7:14]
	filters = []
	uv_errs = []
	for entry in phot:
		if entry[0] == 'e':
			uv_errs.append(entry) #got to handle the errors
		elif entry[0] == 'l':
			continue
		else:
			filters.append(entry)


	all_gals = list(aperture_uvmags_table['Name'][:])
	index = all_gals.index(galaxy)

	#print(uv_columns)
	uv_galnames = list(aperture_uvmags_table['Name'][:])


	nonphot_dict = {key:aperture_uvmags_table[key][index] for key in nonphot}
	uv_phot_dict = {key:aperture_uvmags_table[key][index] for key in filters}
	uv_errs_dict = {key:aperture_uvmags_table[key][index] for key in uv_errs}

	#change magnitudes to flux densities (in mJy)
	#f_den = 3631*10**(mag/-2.5)
	for key, mag in uv_phot_dict.items():
		if 'apmag' in key:
			f_den = 10e3*3631*10**(mag/-2.5) #in mJy
			f_den *= 2.754 #adding skelton factor
			f_err = 10e3*(np.log(10)/2.5)*uv_errs_dict['e_'+key]*10**(mag/-2.5)
			f_err *= 2.754
			uv_phot_dict[key] = f_den
			uv_errs_dict['e_'+key] = f_err

	#update the established dictionary with the table entries
	gal_dict[galaxy].update(nonphot_dict)
	gal_dict[galaxy].update(uv_phot_dict)
	gal_dict[galaxy].update(uv_errs_dict)



	#==============================================================

	'''
	pulls the relevant spitzer photometry, sorts it into photometry, errors, and extra stuff, and adds it to the dictionary
	'''
	spitzer_ir_fluxes = Table.read('spitzer_ir_fluxes.dat', format='ascii.cds')

	#cut out all column names except the filters
	inf_columns = spitzer_ir_fluxes.colnames
	only_filters = []
	ir_errs = []
	for entry in inf_columns:
		if entry[0] == 'F':
			only_filters.append(entry)
		elif entry[0] == 'e':
			ir_errs.append(entry)
	inf_galnames = list(spitzer_ir_fluxes['Name'][:])

	index = inf_galnames.index(galaxy)


	ir_filternames = ['2MASS_J', '2MASS_H', '2MASS_Ks', 'iracch1', 'iracch2', 'iracch3', 'iracch4', 'MIPS24um', 'MIPS70um', 'MIPS160um']


	ir_dict = {'f_'+key:2.754*10e3*spitzer_ir_fluxes[only_filters[i]][index] for i, key in enumerate(ir_filternames)} #in mJy*Skelton factor
	ir_errs_dict = {'e_'+key:2.754*10e3*spitzer_ir_fluxes[ir_errs[i]][index] for i, key in enumerate(ir_filternames)}



	gal_dict[galaxy].update(ir_dict)
	gal_dict[galaxy].update(ir_errs_dict)


	#==============================================================

	#upload and store optical magnitudes

	with open('../python/optical_lvls_mags.dat','r') as f: 
		tab1 = f.readlines()

	# parse it
	# need lines 36 to 295 (the rest is comments and such)
	optical_data = tab1[35:295]

	n_gal1 = len(optical_data) #number of galaxies

	gal_list = list()

	#takes Name column of the optical magnitudes table and turns it into a reference list to find the correct row in the table
	for i in range(n_gal1):
		opt_line = optical_data[i].split()[0:-1] #makes list of each value separated by a tab
		opt_name = opt_line[0]
		gal_list.append(opt_name)


	#=============================

	# save the filters in a dictionary with nans to placehold
	opt_filter_list = ['U', 'B', 'V', 'R', 'u', 'g', 'r', 'i', 'z']
	opt_mags_dict = {'f_'+key:np.nan for key in opt_filter_list}


	opt_errs = {'e_'+key:np.nan for key in opt_filter_list}


	#===============================================================

	for i,gal in enumerate(gal_list):
		if gal == galaxy:
			opt_line = optical_data[i].split()[0:] #makes a table row into list of column entries separated by a tab
			opt_mags = opt_line[4::2]
			opt_mags_err = opt_line[5::2]
		else:
			continue
		#append galaxy name to a list, then stick the list as the value to the key

	#assign photometric values to the correct filter keys in the dictionary
	for p,val in enumerate(opt_mags):

		#if statement to handle any upper limit flags (denoted by '>')
		if val[0] == '>': 
			opt_mags_dict['f_'+opt_filter_list[p]] = float(val[1:])

		else:
			opt_mags_dict['f_'+opt_filter_list[p]] = float(val)


	#write in the errors as well
	for p,err in enumerate(opt_mags_err):
		opt_errs['e_'+opt_filter_list[p]] = float(err)



	for key, mag in opt_mags_dict.items():
		f_den = 10e3*3631*10**(mag/-2.5) #in mJy
		f_den *= 2.754
		f_err = 10e3*(np.log(10)/2.5)*opt_errs['e_'+key[-1]]*10**(mag/-2.5)
		f_err *= 2.754
		opt_mags_dict[key] = f_den
		opt_errs['e_'+key[-1]] = f_err
	

	gal_dict[galaxy].update(opt_mags_dict)
	gal_dict[galaxy].update(opt_errs)

	#==========================================================

	#find the files that report Swift UV photometry in the galaxy directory itself
	os.chdir('../'+galaxy)

	swift_filter_list = ['swift_uvw1', 'swift_uvw2', 'swift_uvm2']

	phot_files = glob.glob('*_phot_totmag.dat')


	if len(phot_files) != 0:
		swift_mag_w1 = np.loadtxt(galaxy+'_w1_phot_totmag.dat', dtype='str', comments='#').tolist() #I know how to index lists better than tables
		swift_mag_w2 = np.loadtxt(galaxy+'_w2_phot_totmag.dat', dtype='str', comments='#').tolist()
		swift_mag_m2 = np.loadtxt(galaxy+'_m2_phot_totmag.dat', dtype='str', comments='#').tolist()

		swift_uv = {'f_swift_uvw1':swift_mag_w1[2], 'f_swift_uvw2':swift_mag_w2[2], 'f_swift_uvm2':swift_mag_m2[2]}

		swift_uv_errs = {'e_swift_uvw1':swift_mag_w1[3], 'e_swift_uvw2':swift_mag_w2[3], 'e_swift_uvm2':swift_mag_m2[3]}
		
		#convert to flux densities
		for i, (key, mag) in enumerate(swift_uv.items()):
			mag = float(mag)
			f_den = 10e3*3631*10**(mag/-2.5) #in mJy
			f_den *= 2.754
			f_err = 10e3*(np.log(10)/2.5)*float(swift_uv_errs['e_'+swift_filter_list[i]])*10**(mag/-2.5)
			f_err *= 2.754
			swift_uv[key] = f_den
			swift_uv_errs['e_'+swift_filter_list[i]] = f_err

		gal_dict[galaxy].update(swift_uv)
		gal_dict[galaxy].update(swift_uv_errs)



	#make the dictionary into an ascii table, then write that table to a text file.

	phot_vals_table = Table()
	for key, val in gal_dict[galaxy].items(): 

		phot_vals_table[key] = [gal_dict[galaxy].get(key)]

	ascii.write(phot_vals_table, 'galphot_table.txt', overwrite=True)

	#and then you're done!


def main():
	parser = argparse.ArgumentParser()

	parser.add_argument("galaxy",nargs = "?", help = "Enter the name of a galaxy", default = '')

	args = parser.parse_args()

	run(args.galaxy)

if __name__ == '__main__':
	main()




