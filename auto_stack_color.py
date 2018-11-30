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

def auto_color():
	gal_list = np.loadtxt("modified_galaxy_list", dtype = 'str').tolist()
	unprocessed_gals = open('unprocessed_gals', 'w+')
	#up_gals = list()
	for gal in gal_list:
		os.chdir('..')
		#cd into existing directory, then make a new one for the stuff we're working on

		dir_exists = glob.glob(gal)
		if len(dir_exists) == 0:
			os.mkdir(gal)

		os.chdir(gal)

		#this may not stay for later versions of this tool, but it's useful to me for now
		processing_done = glob.glob(gal+"_image.png")
		if len(processing_done) > 0: 
			print("* " + gal + "_image.png already created, moving on to next item...")
			continue				 

		'''processing is done if the .png file was created, so if it's there, then just move on
		   this is temporary, in case it doesn't finish a list of galaxies and the script needs to be run again'''

		#make the download script and store it in the new directory here?
		os.system('browse_extract_wget.pl table=swiftmastr position=' + gal + ' radius=5 fields=obsid,start_time outfile=data.dat')

		idl = pIDLy.IDL()

		#run browse_extract with all of the parameters needed to make data.dat

		#filename.open() new file here that will hold all the wget commands
		download_scr = open('download.scr', 'w+')

		obslist = np.loadtxt('data.dat', dtype = 'str', delimiter = '|', comments = 'S', skiprows = 2, usecols = (1,2)).tolist()
		id_list = list()

		#condition that handles cases where HEASARC query returns only one row
		if len(obslist[0]) > 2:
			#print(obslist)
			obsid = obslist[0]
			starttime = obslist[1]
			start_month = starttime[0:7]
			start_month = start_month.replace('-','_')
			id_list.append(obsid)
			wget_uvot = "wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + start_month + '//' + obsid + "/uvot/"
			wget_auxil = "wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + start_month + '//' + obsid + "/auxil/"
			download_scr.write(wget_uvot)
			download_scr.write("\n")
			download_scr.write(wget_auxil)
			download_scr.write("\n")
		else:
			for i in range(len(obslist)):
				#print(obslist[i])
				[obsid, starttime] = obslist[i]
		
				start_month = starttime[0:7]
				start_month = start_month.replace('-','_')
				id_list.append(obsid)
		
	#	string addition to make wget commands for data download
				wget_uvot = "wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + start_month + '//' + obsid + "/uvot/"
				wget_auxil = "wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + start_month + '//' + obsid + "/auxil/"
				download_scr.write(wget_uvot)
				download_scr.write("\n")
				download_scr.write(wget_auxil)
				download_scr.write("\n")

		#check for me to know if download has been done already
		#ooh problem, it will stop processing for all galaxies, I want it to pass on to the data processing stage
		#gotta fix that
		#shit
		download_scr.close()

		obs_folder_exists = glob.glob("000*")
		if len(id_list) == len(obs_folder_exists):
			print("* Download already done, moving on to data processing...")
			#make data processing chunk into a function
		else:
			#make download script executable
			st = os.stat('download.scr')
			os.chmod('download.scr', st.st_mode | stat.S_IEXEC)
			print("* running download script for "+gal)
			os.system('download.scr')

			#unzip all the downloaded data
			os.system('gunzip */*/*.gz')
			os.system('gunzip */*/*/*.gz')


		Name_Prefix = gal

		filt_list = ['w2','m2','w1']
		filter_list = list()
		for filt in filt_list:
			filt_data_exists = glob.glob('*/*/image/*'+filt+'*')
			if len(filt_data_exists) > 0: #glob useful for wildcards here, could be used to check images before processing anyway
				filter_list.append(filt)

		if len(filter_list) == 1:
			print('Data only found for filter '+filter_list[0]+", processing will not continue. Please check to make sure data for each filter exists.")
			print('Moving on to next object...')
			unprocessed_gals.write(gal)
			unprocessed_gals.write('\n')
			continue
		elif len(filter_list) == 2:
			print('Data only found for filters '+filter_list[0]+" and "+filter_list[1]+", processing will not continue. Please check to make sure data for each filter exists.")
			print('Moving on to next object...')
			unprocessed_gals.write(gal)
			unprocessed_gals.write('\n')
			continue

		uvot_deep.uvot_deep(id_list, Name_Prefix+'_', ['w2','m2','w1']) #tool made by Lea Hagen

		#filter list in most cases will be the same, but this section will take care of cases where there's no data for one of the UV filter

		#now check to make sure all the deep images needed for processing were made using glob
		deep_image_list = list()
		for filt in filter_list:
			deep_image_exists = glob.glob('*'+filt+'_sk.fits')
			if len(deep_image_exists) > 0: #glob useful for wildcards here, could be used to check images before processing anyway
				deep_image_list.append(filt)

		if len(filter_list) == 1:
			print('Deep image only created for '+filter_list[0]+", processing will not continue. Please check to make sure data for each filter exists.")
			print('Moving on to next object...')
			unprocessed_gals.write(gal)
			unprocessed_gals.write('\n')
			continue
		elif len(filter_list) == 2:
			print('Deep images only created for '+filter_list[0]+" and "+filter_list[1]+", processing will not continue. Please check to make sure data for each filter exists.")
			print('Moving on to next object...')
			unprocessed_gals.write(gal)
			unprocessed_gals.write('\n')
			continue

		if True:
			offset_mosaic.offset_mosaic(gal+'_', gal+'_offset_', filter_list,
										min_exp_w2=150, min_exp_m2=150, min_exp_w1=150, restack_id=True)

		#list of uv filters, can be changed if someone needs optical filter images as well.
		#filter_list = ['w1', 'm2', 'w2']

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

		#another block to handle the lack of data for one filter
		# im_list = [im_r, im_g, im_b]
		# image_list = list()
		# for image in im_list:
		# 	if os.path.isfile(image):
		# 		image_list.append(image)

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

	unprocessed_gals.close()

auto_color()

'''
latest issues:
image is produced, but background levels are still very obviously visible
	- lea has something called offset_mosaic in her repository
		- takes background 

'''

#at this point the script should be able to access browse_extract, return the data file made, open a download script, open a data.list
#for the entries in the list, store the three columns to three diff variables
#write the wget_uvot and auxil statements, then write them into the download script and data.list
