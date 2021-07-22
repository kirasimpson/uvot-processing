from astropy.table import Table

'''
Test script to fill in zeros in UGC/NGC column names the same length of characters
i.e. NGC42 to NGC0042
'''

position_table = Table.read('https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t1_mrt.txt', format='ascii.cds')
aperture_table = Table.read('https://iopscience.iop.org/0067-0049/192/1/6/suppdata/apjs375285t2_mrt.txt', format='ascii.cds')


#photometric zero points for each UV filter
#zeropoint_dict = {'w1': 18.95,'w2': 19.11,'m2': 18.54}

#print(position_table["Name"])
for i, entry in enumerate(position_table["Name"]):
	if 'NGC' in entry:
		pieces = entry.split(entry[2])
		num = pieces[1]
		zerofill = num.zfill(4)
		newstring = 'NGC' + zerofill
		print(newstring)
		position_table["Name"][i] = newstring
	elif 'UGC' in entry:
		pieces = entry.split(entry[2])
		num = pieces[1]
		zerofill = num.zfill(5)
		newstring = 'UGC' + zerofill
		position_table["Name"][i] = newstring
print(position_table["Name"])
#input('HALT BASTARD')