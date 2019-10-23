from astropy.table import Table
import numpy as np
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import os.path
import glob

import pdb

#input will have to be a galaxy name, either here at the beginning or later.
#def table_time(): 
'''
Combine all the photometry into a useful/accessible format

Optical photometry from Cook+8:
https://academic.oup.com/mnras/article/445/1/881/984908#92778246

Output will be dictionary containing all the magnitudes
for one particular galaxy
'''

#to make this a separate useful module, make an option that someone can use this on one galaxy or a list of them
#for now, restructure this to take just one galaxy input so it fits into the pipeline and feeds one dictionary into the sed code

#stopping point for next week: make code such that I feed in one galaxy and get one galaxy out

def make_sed_table(gal): #takes gal as the input from the other script, then have to make this function accessible
	#apertures used in processing
	aperture_table = Table.read('https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t2_mrt.txt', format='ascii.cds')
	#nuv-fuv mags
	galex_uvmag_table = Table.read('https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t2_mrt.txt', format='ascii.cds')

	#==============================================================

	'''
	Block for converting the infrared fluxes into magnitudes
	'''
	spitzer_infmag_table = Table.read('https://iopscience.iop.org/0004-637X/703/1/517/suppdata/apj314922t2_mrt.txt', format='ascii.cds')
	inf_columns = spitzer_infmag_table.colnames
	#print(spitzer_infmag_table)

	zeropoint = {'F1.25':1594.0, 'F1.65':1024.0, 'F2.17':666.7, 'F3.6':280.9, #janskys for 3.6-8.0 filters, but...hm. numbers/units feel wrong, denoted as F_(0,v)
					'F4.5':179.7, 'F5.8':115.0, 'F8.0':64.9, 'F24':7.17, 'F70':0.778, 'F160':0.159}

	zeropoint_err = {'F1.25':27.8, 'F1.65':20.0, 'F2.17':12.6, 'F3.6':4.1, 'F4.5':2.6, #zeropoints are flux density per unit frequency (janskys)
						'F5.8':1.7, 'F8.0':0.9, 'F24':0.11, 'F70':0.012, 'F160':0.020} #where are these damn zeropoints

	for colname in inf_columns:
		if 'e_' not in colname:
			continue
		else:
			flux_err_list = list(spitzer_infmag_table[colname][:])
			print(flux_err_list)
			input()

			filtname = colname.replace('e_', '')
			inf_flux = list(spitzer_infmag_table[filtname][:])
			inf_mag = list()
			inf_mag_err_list = list()
			for i in range(len(flux_err_list)):

				mag_err = np.sqrt( ( 2.5/np.log(10) * flux_err_list[i]/inf_flux[i] )**2 )

				mag = -2.5 * np.log10(inf_flux[i]/zeropoint.get(filtname))  #take some values and plug in manually to see what's wrong, zeropoint flux might be the problem
				
				inf_mag.append(mag)
				inf_mag_err_list.append(mag_err)

			spitzer_infmag_table[colname][:] = inf_mag_err_list
			spitzer_infmag_table[filtname][:] = inf_mag
			continue

	#the 

	#reference them the same as I would any of the other tables I've uploaded

	#==============================================================

	#Block for uploading and storing optical magnitudes

	with open('../python/optical_lvls_mags.dat','r') as f:
		tab1 = f.readlines()

	# parse it
	# need lines 36 to 295 (the rest is comments and such)
	optical_data = tab1[35:295]
	n_gal1 = len(optical_data)

	gal_list = list()

	for i in range(n_gal1):
		opt_line = optical_data[i].split()[0:-1] #makes list of each value separated by a tab
		opt_name = opt_line[0]
		gal_list.append(opt_name)

	# save them
	opt_filter_list = ['U', 'B', 'V', 'R',
					   'u', 'g', 'r', 'i', 'z']

	gal_dicts = {key:{filt:np.nan for filt in opt_filter_list} for key in gal_list}
	mag_temp = {key+'_err':np.nan for key in opt_filter_list}
	for gal in gal_list:
		gal_dicts[gal].update(mag_temp)

	#==========================================================

	for i,gal in enumerate(gal_list):
		opt_line = optical_data[i].split()[0:-1] #makes list of each value separated by a tab
		opt_mags = opt_line[4::2]
		opt_mags_err = opt_line[5::2]
		#append galaxy name to a list, then stick the list as the value to the key

		for p,val in enumerate(opt_mags): #all the splitting is completely unnecessary, yeehaw

			if val[0] == '>': #first character in a string of letters
				gal_dicts[gal][opt_filter_list[p]] = float(val[1:])

			else:
				gal_dicts[gal][opt_filter_list[p]] = float(val)


		for p,err in enumerate(opt_mags_err):
			if err[0] == '-':
				gal_dicts[gal][opt_filter_list[p]+'_err'] = float(err) #problem is with the _err addition, perhaps?
			else:
				gal_dicts[gal][opt_filter_list[p]+'_err'] = float(err)

	#==========================================================

	galex_filter_list = ['NUV', 'FUV']

	# gal_dicts = {key:{filt:np.nan for filt in opt_filter_list} for key in gal_list}
	# mag_temp = {key+'_err':np.nan for key in opt_filter_list}
	# for gal in gal_list:
	# 	gal_dicts[gal].update(mag_temp)

	galex_mags_dict = {key:np.nan for key in galex_filter_list}
	galex_mags_errs = {key+'_err':np.nan for key in galex_filter_list}
	for gal in gal_list:
		gal_dicts[gal].update(galex_mags_dict)
		gal_dicts[gal].update(galex_mags_errs)


	#==========================================================

	#gal comes from lvls_list, the 258 gals
	#have to update table by table
	#have gal as an input, do work from the outside?
	spitzer_filter_list = ['2MASS_J', '2MASS_H', '2MASS_KS', 'Spitzer_J', '4.5', \
							'5.8', '8.0', '24', '70', '160',] 

	spitzer_col_list = ["F1.25", "F1.65", "F2.17", "F3.6", "F4.5", "F5.8", "F8.0", "F24", "F70", "F160"]

	counter1 = 0
	for gal in gal_list: #the below if/else might not be functional anymore, the galaxy list comes from the optical table, which matches the order of the inf table
	#but it might come in handy for the NUV/FUV measurements
		# if entry != gal:
		counter1 += 1
		# else:
		if counter1 == 258: #find a better solution
			continue

		for ind, filt in enumerate(spitzer_filter_list):

			update = {filt:spitzer_infmag_table[spitzer_col_list[ind]][counter1]}
			gal_dicts[gal].update(update)

			#"2MASS_H": spitzer_infmag_table["F1.65"][counter1],

			update_err = {filt+'_err':spitzer_infmag_table['e_'+spitzer_col_list[ind]][counter1]}
			gal_dicts[gal].update(update_err)


	#print(gal_dicts['WLM'])
	#input('???? BANG')

	#what was all this for in the past
	counter2 = 0
	for entry in aperture_table["Name"]:
		if entry != gal:
			counter2 += 1
			continue

		else: 
			aperture_dict = {
			"Major Axis": aperture_table["MajAxis"][counter2],
			"Minor Axis": aperture_table["MinAxis"][counter2],
			"Position Angle": aperture_table["PosAng"][counter2]
			}

	#==========================================================


	#swift values are in three separate files, use loadtxt
	#for filt in ['w1','m2','w2']: #need something new every time, a new name that isn't filt_mag

	swift_filter_list = ['UVw1', 'UVw2', 'UVm2']

	#gal being inserted here...where do I make it come from
	#I'll have a directory problem here
	phot_files = glob.glob(gal+'_*_totmag.dat')

	for gal in gal_list:
		if len(phot_files) == 0:
			continue
		else:
			swift_mag_w1 = np.loadtxt(gal+'_w1_phot_totmag.dat', dtype='str', comments='#').tolist() #working with a few lists instead of full tables
			swift_mag_w2 = np.loadtxt(gal+'_w2_phot_totmag.dat', dtype='str', comments='#').tolist()
			swift_mag_m2 = np.loadtxt(gal+'_m2_phot_totmag.dat', dtype='str', comments='#').tolist()

			update = {'UVw1':swift_mag_w1[2], 'UVw2':swift_mag_w2[2], 'UVm2':swift_mag_m2[2]}
			gal_dicts[gal].update(update)

			update_err = {'UVw1_err':swift_mag_w1[3], 'UVw2_err':swift_mag_w2[3], 'UVm2_err':swift_mag_m2[3]}
			gal_dicts[gal].update(update_err)

			#dictionaries complete?


	print(gal_dicts)
	input()
		#test how this is printed out

	#first I need to extract swift values from these tables
	#put copy in a directory that has swift photometry

	#swift phot values are in files gal_filt_phot_totmag.dat, need 3rd and 4th col values
	#extract list of spitzer subset

	#==========================================================

	indices = np.arange(len(gal_list)) #so how does this section need to be different?

	arr = np.arange(len(gal_list)*(51)).reshape(len(gal_list), 51)

	#the names here are incorrect to begin with
	allphot_table = Table(arr, names=('Name', 'Maj_Axis', 'Min_Axis', 'U', 'B', 'V', 'R', \
									'u', 'g', 'r', 'i', 'z', 'Spitzer Global error', 'NUV', 'NUV error', \
									'FUV', 'FUV error', '2MASS J', '2MASS J error', 'H', 'H error', 'KS', \
									'KS error','Spitzer J', 'Spitzer J error', '4.5μm', '4.5μm error', \
									'5.8μm', '5.8μm error', '8.0μm', '8.0μm error', '24μm', '24μm error', \
									'70μm', '70μm error', '160μm', '160μm error', 'UVw1', 'UVw1 error', \
									'UVm2', 'UVm2 error', 'UVw2', 'UVw2 error'))

if __name__ == '__main__':

    make_sed_table()