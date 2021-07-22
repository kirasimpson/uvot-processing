#! /user/bin/env python

import pidly
#import argparse 
import numpy as np
#import subprocess
import os
import glob
import shutil
import stat
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip, biweight_location
import aplpy
import uvot_deep
import config_uvot_mosaic
import offset_mosaic
import swift_db_query as sdq
import query_heasarc as qh
import download_heasarc as dh
import surface_phot as sp
import phot_plot as phot
import gal_table_sed as sed

import pdb

def leading_zeros(table):
	'''
	Insert leading zeros so that NGC/UGC/UGCA galaxy names match between LVLS list of 258 
	galaxies and aperture/position tables.
	Parameters
	----------
	table : ASCII table
		name of a table that follows modification rules of an astropy Table object.
		https://docs.astropy.org/en/stable/table/modify_table.html
	'''

	#reduce this just a little
	for i, entry in enumerate(table["Name"]):
		if 'NGC' in entry:
			pieces = entry.split(entry[2])
			num = pieces[1]
			zerofill = num.zfill(4)
			newstring = 'NGC' + zerofill
			#print(newstring)
			table["Name"][i] = newstring
		elif 'KKH' in entry:
			pieces = entry.split(entry[2])
			num = pieces[1]
			zerofill = num.zfill(3)
			newstring = 'KKH' + zerofill
			table["Name"][i] = newstring
		elif 'IC' in entry:
			pieces = entry.split(entry[1])
			num = pieces[1]
			zerofill = num.zfill(4)
			newstring = 'IC' + zerofill
			table["Name"][i] = newstring
		elif 'BK' in entry:
			pieces = entry.split(entry[1])
			num = pieces[1]
			zerofill = num.zfill(3)
			newstring = 'BK' + zerofill
			table["Name"][i] = newstring
		elif 'UGC' in entry:
			if 'A' in entry:
				pieces = entry.split(entry[3])
				num = pieces[1]
				zerofill = num.zfill(3)
				newstring = 'UGCA' + zerofill
				table["Name"][i] = newstring
			else:
				pieces = entry.split(entry[2])
				num = pieces[1]
				zerofill = num.zfill(5)
				newstring = 'UGC' + zerofill
				table["Name"][i] = newstring

	return table

#=============================================================================	

def ra_dec_converter(RA_hr, RA_min, RA_sec, degrees, minutes, seconds, sign):
	'''
	RA/Dec HMS/DMS to Decimal converter for position of object
	Inputs are all table values from position_table in photometry_time function.
	'''
	RA_decimal = 15*(float(RA_hr) + float(RA_min)/60 + float(RA_sec)/3600)

	Dec_decimal = float(degrees) + float(minutes)/60 + float(seconds)/3600

	if sign == "-":
		Dec_decimal *= -1

	coord_dict = {"RA":RA_decimal, "Dec":Dec_decimal}

	return coord_dict

#=============================================================================

def photometry_time():
	'''
	Do photometry and produce plots!
	Extracts positions and aperture parameters from Lee+11 
	http://adsabs.harvard.edu/abs/2011ApJS..192....6L
	surface_phot and phot_plot made by Lea Hagen
	'''

	position_table = Table.read('https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t1_mrt.txt', format='ascii.cds')
	aperture_table = Table.read('https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t2_mrt.txt', format='ascii.cds')


	#photometric zero points for each UV filter
	zeropoint_dict = {'w1': 18.95,'w2': 19.11,'m2': 18.54}

	leading_zeros(position_table)
	leading_zeros(aperture_table)

	#print(leading_zeros(position_table))
	#input('STOP BASTARD')

	counter1 = 0

	for entry in position_table["Name"]:

		if entry != gal:
			counter1 += 1
			continue

		else:
			position_dict = {
			"RA_hr": position_table["RAh"][counter1],
			"RA_min": position_table["RAm"][counter1],
			"RA_sec": position_table["RAs"][counter1],

			"Sign": position_table["DE-"][counter1],
			"Degrees": position_table["DEd"][counter1],
			"Minutes": position_table["DEm"][counter1],
			"Seconds": position_table["DEs"][counter1]
			}


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

			#input('inordinate screeching!!')

	coordinates = ra_dec_converter(position_dict["RA_hr"], position_dict["RA_min"],
									position_dict["RA_sec"],position_dict["Degrees"],
									position_dict["Minutes"],position_dict["Seconds"],
									position_dict["Sign"])


	for key in zeropoint_dict:

		sp.surface_phot(gal, coordinates['RA'], coordinates['Dec'], 
						aperture_dict['Major Axis'], aperture_dict['Minor Axis'], 
						aperture_dict['Position Angle'], 5.0, key, zeropoint_dict[key]) #requires inputs from table2 as well

		#how does the annprofile.dat file get created? why is it not created for NGC0024?
		#. . . beef up error reporting/print statements
		phot.phot_plot([gal+'_'+key+'_phot_annprofile.dat'], [gal+'_'+key+'_phot_totprofile.dat'], key, gal+key+'phot_plot.png')

#=============================================================================

def make_cube(process=True):
	'''
	Turn all the UV data into a cube!
	'''

	cr_img_exists = glob.glob('*cr.fits') #makes a list of count rate images that exist in the current directory, checks if uvot_deep done

	filter_list = ['w2','m2','w1']

	if len(cr_img_exists) == 0:

		id_list = glob.glob('000*')

		try:
			uvot_deep.uvot_deep(id_list, gal+'_', filter_list) #tool adapted by Lea Hagen
		except FileNotFoundError:
			print('* Required files not made out of uvot_deep for '+gal+', moving on...')
			unprocessed_gals.write(gal + '\n')
			why_bypass.write(gal+': required files not made out of uvot_deep'+'\n')
			process = False

		except OSError:
			print('Too many open files error thrown, whatever that means. Moving on...') #sass maybe not good, maybe chill
			unprocessed_gals.write(gal + '\n')
			why_bypass.write(gal+': too many open files error' + '\n')
			process = False #does this work here?


	else:
		print("* uvot_deep already completed successfully, moving on!")
		#input()

	#offset_mosaic running over and over again, need to skip over this right now
	#if images have been created already, don't run it again
	#later will have to take into consideration whether or not there is new data that requires a rerun of uvot_deep

	offset_done = glob.glob(gal + '_offset_*')

	if len(offset_done) == 0:
		try: #make sure this does run again
			if True:
				print('* running offset_mosaic')
				offset_mosaic.offset_mosaic(gal+'_', gal+'_offset_', filter_list,
											min_exp_w2=150, min_exp_m2=150, min_exp_w1=150, restack_id=True)
		except IndexError:
			print("* index 0 is out of bounds for axis 0 with size zero")
			unprocessed_gals.write(gal+'\n')
			why_bypass.write(gal+': offset mosaic problem, index 0 is out of bounds for axis 0 with size zero'+'\n')
			#continue #fix continue
			process = False

	else:
		print('* offset_mosaic has been run already, moving on!')

	

	#list of uv filters, can be changed if someone needs optical filter images as well.
	#filter_list = ['w1', 'm2', 'w2']

	idl = pidly.IDL('/bulk/pkg/local/bin/idl')

	for filt in filter_list:
	
		#taken from Lea's code stack_uvot
		with fits.open(gal + "_" + filt + '_cr.fits') as hdu_cr: #cr = count rate

			pix_clip = sigma_clip(hdu_cr[0].data, sigma=2.5, maxiters=3)
			cr_mode = biweight_location(pix_clip.data[~pix_clip.mask])
			#print('cr_mode: ', cr_mode)

			hdu_cr[0].data -= cr_mode

			hdu_cr.writeto(gal + '_offset_' + filt + '_cr_bgsub.fits', overwrite=True)

		idl.pro('image_shift', gal) #idl procedure that aligns all the image_shift

	idl.close()

		#replaces mkcolor idl script (credits to Lea)
		# make RGB image

	cube_done = glob.glob(gal + '_rgb_cube*')
	if len(cube_done) == 0:
		#if no data for the filter, gotta skip over the image that would be created
		im_r = gal + '_offset_w1_cr_bgsub.fits'
		im_g = gal + '_offset_m2_cr_bgsub.fits'
		im_b = gal + '_offset_w2_cr_bgsub.fits'

		vmid_list = 0.001
		vmin_list = []
		
		for i,im in enumerate([im_r, im_g, im_b]):
			hdu_list = fits.open(im)
			filt = sigma_clip(hdu_list[0].data, sigma=2, maxiters=5)
			vmin_list.append( np.mean(filt.data[~filt.mask]) + 1.5 * np.std(filt.data[~filt.mask]) )
			hdu_list.close()

		# - create fits images in same projection
		print('* creating rgb fits cube for '+gal)
		aplpy.make_rgb_cube([im_r, im_g, im_b],
								gal+'_rgb_cube.fits', north=True)
		# - make the rgb image
		print('* creating rgb png image of '+gal)
		aplpy.make_rgb_image(gal+'_rgb_cube.fits', gal+'_image'+'.png',
								 stretch_r='log', stretch_g='log', stretch_b='log',
								 vmin_r=vmin_list[0], vmin_g=vmin_list[1], vmin_b=vmin_list[2],
								 vmid_r=vmid_list, vmid_g=vmid_list, vmid_b=vmid_list,
								 pmax_r=99.95, pmax_g=99.9, pmax_b=99.9,
								 make_nans_transparent=True)
	else:
		print('rgb_cube FITS file already exists, moving on to .png creation and photometry!')

	return process

#=============================================================================

def filter_check():
	filt_list = ['w2','m2','w1']
	filter_list = list()
	for filt in filt_list:
		filt_data_exists = glob.glob('*/*/image/*'+filt+'*')
		if len(filt_data_exists) > 0:
			filter_list.append(filt)

	if len(filter_list) == 1:
		print('Data only found for filter '+filter_list[0]+", processing will not continue. Please check to make sure data for each filter exists.")
		print('Moving on to next object...')
		unprocessed_gals.write(gal + "\n")
		why_bypass.write(gal+': missing data for some filters'+'\n')
		#continue #fix continue
		return False
	elif len(filter_list) == 2:
		print('Data only found for filters '+filter_list[0]+" and "+filter_list[1]+", processing will not continue. Please check to make sure data for each filter exists.")
		print('Moving on to next object...')
		unprocessed_gals.write(gal+'\n')
		why_bypass.write(gal+': missing data for some filters'+'\n')
		#continue #fix continue
		return False
	else:
		return True

#=============================================================================

def process_progress():
	'''
	Checks the progress of creating the data cube and the .png file
	by looking for the .png file and comparing the number of observation
	directories to the number of obsIDs returned by the heasarc query.
	If the file exists and there is no new data, it will skip over that 
	part of processing.
	'''
	processing_done = glob.glob(gal+"_image.png")
	obs_folder_exists = glob.glob("00*")

	with open('heasarc_obs.dat', 'r') as fh:
		rows_list = fh.readlines()

	total_obs = len(rows_list) - 4 #minus first four lines, which are just formatting

	if len(processing_done) > 0 and len(obs_folder_exists) == total_obs:
		print("* Nothing new to add for "+gal+". Moving along...")
		return True

	elif len(processing_done) > 0 and len(obs_folder_exists) != total_obs:
		os.remove(gal+'_image.png')
		print('Rerunning image processing...')
		make_cube()
		return True
	else:
		return False

#=============================================================================

unprocessed_gals = open('unprocessed_gals', 'w+')
why_bypass = open('why_bypass', 'w+')

#=============================================================================

def phot_progress(): #check for photometry files
	phot_files = glob.glob(gal+'_*_totmag.dat')
	if len(phot_files) == 0:
		return True

#=============================================================================
#I have a problem where two different lists of galaxies are clashing with each other.
#the easiest solution would seem to be making the lists match up, but what problems would that cause? any?

#currently pulling the list of galaxies from the Spitzer LVL optical magnitudes table (Cook et al.)
with open('optical_lvls_mags.dat','r') as f:
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

#print(gal_list)

for gal in gal_list:
	#currently skipping a bunch of galaxies until I can further investigate
	if gal == "UGC2847" or gal == 'KKH098' or gal == 'LSBCD565-06' or gal == 'NGC3628' or gal == 'Circinus' or gal == "LEDA166101" or gal == 'HS98_117' or gal == 'SCULPTOR-DE1' or gal == 'NGC3031' or gal == 'NGC3034' or gal == 'UGC05336' or gal == 'UGC07699' or gal == 'NGC5055' or '[' in gal or gal == 'ESO410-G005' or gal == 'F8D1' or gal == 'NGC4594' or gal == "NGC5195" or gal == "NGC5236" or gal == 'NGC5457': #UGC2034, Maffei2, UGC2259
		continue

		#[FM2000]1: files are not gunzipping for some reason

	if 'dag' in gal: #galaxies with this tag were removed from the sample
		continue

	os.chdir('..') #cd into full lvls directory so it doesn't make everything in your python directory
	
	filter_ch = ''

	
	#check if photometry has been done already


	print('querying for '+gal)
	qh.query_heasarc(gal) 

	os.chdir(gal)

	photometry_done = glob.glob('*plot.png')
	if len(photometry_done) > 0:
		print("Photometry plots made, moving along!")
		continue

	#ch
	progress = process_progress()
	if progress == True:

	#find a better solution
	#current solution: two versions of the same set of commands, one slightly modified to have one other command

		#make sure data for all filters exists for a given object
		filter_ch = filter_check()

		if filter_ch == False:
			continue

		processing = make_cube()
		if processing == False:
			continue

		
		value_error_list = list()
			#processing progress function and filter checks still need to go somewhere in here, but they can be jumped to instead of being part of the bigger function
		try:
			print('Doing photometry for '+gal+'...')
			photometry_time()
		except ValueError:
			print('Exception occurred.')
			print('Values in Aperture or Position table missing, skipping entry.')
			value_error_list.append(gal)
			continue
		except TypeError:
			print('Exception occurred.')
			print('TypeError, likely non-empty vector case')
			continue
		
	else:

		dh.download_heasarc('heasarc_obs.dat')

		filter_ch = filter_check()
	#	print(filter_check)
		#input('STOP')
		if filter_ch == False:
			continue

		processing = make_cube() #does it just run then??
		if processing == False:
			continue

		
		#make_cube()
		value_error_list = list()
			#processing progress function and filter checks still need to go somewhere in here, but they can be jumped to instead of being part of the bigger function
		try:
			print('Doing photometry for '+gal+'...')
			photometry_time()
		except ValueError:
			print('Exception occurred.')
			print('Values in Aperture or Position table missing, skipping entry.')
			value_error_list.append(gal)
			continue
		except TypeError:
			print('Exception occurred.')
			print('TypeError, likely non-empty vector case')
			continue
		



	#this is not working, so next step is to find out where exactly the script is stopping and put the fix there
	phot_done = phot_progress()
	if phot_done == True:
		continue


	#now ingest photometry from archival data
	#sed.make_sed_table(gal) #incomplete, won't make anything out of WLM yet


unprocessed_gals.close()

why_bypass.close()
