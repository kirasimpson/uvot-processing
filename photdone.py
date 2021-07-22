'''
Test file to check which galaxies have complete photometry by looking for a photometry plot .png file
'''


import os
import glob

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

print('There are ' + str(len(gal_list)) + ' galaxies in this sample')

os.chdir('..')

done_list = list()

for gal in gal_list:
	if 'dag' in gal:
		continue

	if gal == "UGC2847" or gal == 'KKH098' or gal == 'LSBCD565-06' or gal == 'NGC3628' or gal == 'Circinus' or gal == "LEDA166101" or gal == 'HS98_117' or gal == 'SCULPTOR-DE1' or gal == 'NGC3031' or gal == 'NGC3034' or gal == 'UGC05336' or gal == 'UGC07699' or gal == 'NGC5055' or '[' in gal or gal == 'ESO410-G005' or gal == 'F8D1' or gal == 'NGC4594' or gal == "NGC5195" or gal == "NGC5236" or gal == 'NGC5457': #UGC2034, Maffei2, UGC2259
		continue

	os.chdir(gal)

	photometry_done = glob.glob('*plot.png')
	if len(photometry_done) > 0:
		done_list.append(gal)
		
	os.chdir('..')

print('Photometry has been done without errors on: ' + str(len(done_list)) + ' galaxies.')