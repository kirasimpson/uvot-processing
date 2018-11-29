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
from astropy.stats import biweight_location, sigma_clip
import aplpy
import uvot_deep
import config_uvot_mosaic

def auto_color():
	gal_list = np.loadtxt("modified_galaxy_list", dtype = 'str').tolist()
	for gal in gal_list:
		os.chdir('..')
		#cd into existing directory, then make a new one for the stuff we're working on

		#using glob to look for files, so I don't delete/redownload every single fuckin' time

		if not os.path.exists("/Users/kzs1031/lvls/"+gal):	
			os.mkdir(gal)
		else:
			shutil.rmtree("/Users/kzs1031/lvls/"+gal)
			os.mkdir(gal)
		os.chdir(gal)
		#make the download script and store it in the new directory here?
		os.system('browse_extract_wget.pl table=swiftmastr position=' + gal + ' radius=5 fields=obsid,start_time outfile=data.dat')

		idl = pIDLy.IDL()

		#run browse_extract with all of the parameters needed to make data.dat

		#filename.open() new file here that will hold all the wget commands
		download_scr = open('download.scr', 'w+')

		obslist = np.loadtxt('data.dat', dtype = 'str', delimiter = '|', comments = 'S', skiprows = 2, usecols = (1,2)).tolist()
		id_list = list()

		for i in range(1,len(obslist)-1):
			#for thing in obslist[i]:
			#	print(thing)
			[obsid, starttime] = obslist[i]
			start_month = starttime[0:7]
			start_month = start_month.replace('-','_')
			id_list.append(obsid)
			#print(start_month)
			#print(obsid)
			#input()
		#	string addition to make wget commands for data download
			wget_uvot = "wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + start_month + '//' + obsid + "/uvot/"
			wget_auxil = "wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + start_month + '//' + obsid + "/auxil/"
			download_scr.write(wget_uvot)
			download_scr.write("\n")
			download_scr.write(wget_auxil)
			download_scr.write("\n")

		#run the download script here and put the results in the directories created at the beginning of the list
		download_scr.close()

		#make download script executable
		st = os.stat('download.scr')
		os.chmod('download.scr', st.st_mode | stat.S_IEXEC)
		print("* running download script for "+gal)
		os.system('download.scr')

		#unzip all the downloaded data
		os.system('gunzip */*/*.gz')
		os.system('gunzip */*/*/*.gz')



		Name_Prefix = gal

		uvot_deep.uvot_deep(id_list, Name_Prefix+'_', ['w2','m2','w1'])

		#list of uv filters, can be changed if someone needs optical filter images as well.
		filter_list = ['w1', 'm2', 'w2']
		bkg_list = list()


		for filt in filter_list: #exp_correct run for every single filter

			idl.pro('exp_correct', Name_Prefix+'_'+filt+"_cr.fits", Name_Prefix+'_'+filt+"_ex.fits", Name_Prefix+filt+"r.fits")

			#taken from Lea's code stack_uvot
			with fits.open(Name_Prefix + filt + 'r.fits') as hdu_cr: #cr = count rate

				# find mode of image using a sigma clip
				pix_clip = sigma_clip(hdu_cr[0].data, sigma=2.5, iters=3) #removes all values more than 2.5sigma from array
pix_clip.mask
				cr_mode = biweight_location(pix_clip.data[~pix_clip.mask])
				#print('cr_mode: ', cr_mode)

				hdu_cr[0].data -= cr_mode

				bkg_list.append(cr_mode)

				hdu_cr.writeto(Name_Prefix + '_offset_' + filt + '_cr_bgsub.fits', overwrite=True)

		print(bkg_list)

		idl.pro('image_shift', Name_Prefix) #idl procedure that aligns all the images

		idl.close()

		#replaces mkcolor idl script (credits to Lea)
		# make RGB image
		im_r = gal + '_offset_w1_cr_bgsub.fits'
		im_g = gal + '_offset_m2_cr_bgsub.fits'
		im_b = gal + '_offset_w2_cr_bgsub.fits'
	
		# - create fits images in same projection
		print('* creating rgb fits cube')
		aplpy.make_rgb_cube([im_r, im_g, im_b],f"{gal}_rgb_cube.fits", north=True)
		# - make the rgb image
		print('* creating rgb png image')
		data = f'{gal}_rgb_cube.fits'
		output = f'{gal}_image.png'
		try: 
		#aplpy crashes upon exiting and returns an error "string expected, got bytes instead"
		#that's gross and annoying
		#make_rgb_image does complete making the png file, so try/except helps me bypass the crash on exiting
		#and I still get the image at the end

			aplpy.make_rgb_image(data, output,
								 stretch_r='log', stretch_g='log', stretch_b='log',
								 vmin_r=bkg_list[0], vmin_g=bkg_list[1], vmin_b=bkg_list[2],
								 pmax_r=99.75, pmax_g=99.5, pmax_b=99.5,
								 make_nans_transparent=True)
		except:
			pass

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

#parser = argparse.ArgumentParser()
#parser.add_argument(description='Choose between inputting a list or a single input', 'integers', metavar='      ', type=int, nargs='+',
#                    help='an integer for the accumulator')
#args = parser.parse_args()
#print(args.accumulate(args.integers))