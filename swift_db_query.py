import mysql
import numpy as np

'''
Query the Swift database for targets marked 'done' in the Local Volume Legacy Survey,
then compare the total UVOT exposure time to the time requested to determine if requesting 
more time is necessary or the target is ready for processing.
'''

def compare_times():

	lvls_reqtime = np.loadtxt('requested_time_lvls.list', dtype = 'str').tolist()
	lvls_obstime = np.loadtxt('query_times.list', dtype = 'str').tolist()
	req_more_time = open('req_more_time_galaxies.txt', 'w+')
	process_time = open('ready_to_process_galaxies.txt', 'w+')

	obstimes_dict = dict()
	for el in lvls_obstime:
		obstimes_dict[el[0]] = int(el[2])

	reqtimes_dict = dict()
	for el in lvls_reqtime:
		reqtimes_dict[el[0]] = int(el[1])

	ready_to_process = list()

	for key in obstimes_dict.keys():
		if key in reqtimes_dict.keys():
			if obstimes_dict[key] > reqtimes_dict[key]:
				#just make a list for now of the stuff that's done, and then I'll see what other checks I have in place that might warrant changes
				ready_to_process.append(key) #so this will make a list of stuff that's ready for processing!
											 #how do I access this list later in the code?? yeek
				process_time.write(key + '\n')

			else:
				time_needed = reqtimes_dict[key] - obstimes_dict[key]
				req_more_time.write(key + ' ' + str(time_needed) + '\n')
		else:

			continue #does this continue make sense??

	process_time.close()
	req_more_time.close()



def obs_progress_query():

	BackendDB = mysql.CalendarDB()

	command = "SELECT `targname`,`targid`, `uvotexp` FROM fillin_targets WHERE `comment` = 'Lea Hagen - LVLS' AND `done` = 1"

	results = BackendDB.cur.execute(command)

	Targ_names_IDs = BackendDB.cur.fetchall()


	file = open("query_times.list", "w+")
	for i in range(len(Targ_names_IDs)):
		file.write(Targ_names_IDs[i][0] + " " + str(Targ_names_IDs[i][1]) + " " + str(Targ_names_IDs[i][2]) + "\n")
	file.close()

	compare_times()

	

def main():
	obs_progress_query()

if __name__ =="__main__":
    main()
