import mysql
import numpy as np



def obs_progress_query(compare_times=False):
	'''
	Execute an SQL query to the Swift backend for finished targets in SGWGSurvey.
	
	Outputs
	-------
	Writes out a file sgwgs_done.list that contains the names and IDs of fully observed targets
	'''

	BackendDB = mysql.CalendarDB()

	command = "SELECT `targname`,`targid` FROM fillin_targets WHERE `comment` LIKE '% Aaron Tohuvavohu%' AND `done` = 1"

	results = BackendDB.cur.execute(command)
	
	Targ_names_IDs = BackendDB.cur.fetchall()

	
	file = open("sgwgs_done.list", "w+")
	for i in range(len(Targ_names_IDs)):
		file.write(str(Targ_names_IDs[i][1])  + " " + str(Targ_names_IDs[i][0]) + "\n")
	file.close()



def main():
	obs_progress_query()

if __name__ =="__main__":
    main()
