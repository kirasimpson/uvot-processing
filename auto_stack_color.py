#! /user/bin/env python

import pIDLy
#import argparse 
import numpy as np
import subprocess
import os
import glob
import shutil
import stat
from astropy.io import fits
from astropy.stats import sigma_clip, biweight_location
import aplpy
import uvot_deep
import config_uvot_mosaic
import offset_mosaic
import swift_db_query
import query_heasarc
import download_heasarc


def image_processing():
		#test this out, include checking for new data
		#galaxies with new data...how do I make a new list?
	idl = pIDLy.IDL()


	cr_img_exists = glob.glob('*cr.fits') #turn this into an if/and where there's a second condition that concerns a second data download...hm. what could I check?
	if len(cr_img_list) == 0 or inter_list:
		try:
			uvot_deep.uvot_deep(id_list, Name_Prefix+'_', ['w2','m2','w1']) #tool made by Lea Hagen
		except FileNotFoundError:
			print('* Required files not made out of uvot_deep for '+gal+', moving on...')
			unprocessed_gals.write(gal + '\n')
			why_bypass.write(gal+': required files not made out of uvot_deep'+'\n')
			continue
		except OSError:
			print('Too many open files error thrown, whatever that means. Moving on...') #sass maybe not good, maybe chill
			unprocessed_gals.write(gal + '\n')
			why_bypass.write(gal+': too many open files error' + '\n')
			continue
	else:
		print("* uvot_deep already completed successfully, moving on!")
		input()


	try:
		if True:
			offset_mosaic.offset_mosaic(gal+'_', gal+'_offset_', filter_list,
										min_exp_w2=150, min_exp_m2=150, min_exp_w1=150, restack_id=True)
	except IndexError:
		print("* index 0 is out of bounds for axis 0 with size zero")
		unprocessed_gals.write(gal+'\n')
		why_bypass.write(gal+': offset mosaic problem, index 0 is out of bounds for axis 0 with size zero'+'\n')
		continue

	#list of uv filters, can be changed if someone needs optical filter images as well.
	#filter_list = ['w1', 'm2', 'w2']

	Name_Prefix = gal

	for filt in filter_list:

		#taken from Lea's code stack_uvot
		with fits.open(Name_Prefix + "_" + filt + '_cr.fits') as hdu_cr: #cr = count rate

			pix_clip = sigma_clip(hdu_cr[0].data, sigma=2.5, iters=3)
			cr_mode = biweight_location(pix_clip.data[~pix_clip.mask])
			#print('cr_mode: ', cr_mode)

			hdu_cr[0].data -= cr_mode

			hdu_cr.writeto(gal + '_offset_' + filt + '_cr_bgsub.fits', overwrite=True)

		idl.pro('image_shift', Name_Prefix) #idl procedure that aligns all the images

	idl.close()

		#replaces mkcolor idl script (credits to Lea)
		# make RGB image

	#if no data for the filter, gotta skip over the image that would be created
	im_r = gal + '_offset_w1_cr_bgsub.fits'
	im_g = gal + '_offset_m2_cr_bgsub.fits'
	im_b = gal + '_offset_w2_cr_bgsub.fits'

	vmid_list = 0.001
	vmin_list = []
	
	for i,im in enumerate([im_r, im_g, im_b]):
		hdu_list = fits.open(im)
		filt = sigma_clip(hdu_list[0].data, sigma=2, iters=5)
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
		continue
	elif len(filter_list) == 2:
		print('Data only found for filters '+filter_list[0]+" and "+filter_list[1]+", processing will not continue. Please check to make sure data for each filter exists.")
		print('Moving on to next object...')
		unprocessed_gals.write(gal+'\n')
		why_bypass.write(gal+': missing data for some filters'+'\n')
		continue



def process_progress():
	#this piece is still useful, but is going to have to change since obslist will be created differently
	processing_done = glob.glob(gal+"_image.png")
	if len(processing_done) > 0 and len(obslist) == len(current_obslist):
		print("* Nothing new to add for "+gal+". Moving along...")
		return True			 
	else:
		#remove the png? and some other files? what else will be frustrating to overwrite?
		os.remove(gal+'_image.png')
		image_processing()





def data_download():
	unprocessed_gals = open('unprocessed_gals', 'w+')
	why_bypass = open('why_bypass', 'w+')

	swift_db_query() #script written to query the swift database, search the sotbackend fillin subdatabase and return all the targets that are marked as done and have
	ready_to_process = np.loadtxt('ready_to_process_galaxies.txt', dtype = 'str').tolist()

	for gal in ready_to_process: #gal_list now comes from swift_db_query, so is it better to just put the code in here??
		os.chdir('..') #cd into full lvls directory so it doesn't make everything in your python directory

		query_heasarc(gal) #only needs to accept one string object since we're iterating through list right now

		obs_folder_exists = glob.glob("00*")
		with open('heasarc_obs.dat', 'r') as fh:
			rows_list = fh.readlines()

		total_obs = len(rows_list) - 4
		if total_obs != len(obs_folder_exists): #i want this to happen if there is a new download, but also if this object is new then download has to be run
												#so if this needs to be run for a new object too, where do I account for that too??


			download_heasarc('heasarc_obs.dat')
		else:
			print("* Download already done")
			print("* No new observations to download")
			print('* Moving on to processing...')
			#uh hang on can I just stick the processing function into this?? OH YOU FOOL THIS IS THE USE OF FUNCTIONIZING
			#YOU MADMAN
			#skip on here to next shenanigans


		progress = process_progress()
		if progress = True:
			continue



		#processing progress function and filter checks still need to go somewhere in here, but they can be jumped to instead of being part of the bigger function
		

		image_processing()

	unprocessed_gals.close()

	why_bypass.close()



auto_color()
