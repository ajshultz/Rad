#! /usr/bin/env python

import re, sys, os, itertools, sets, getopt        # Load standard modules I often use
import MySQLdb   # Gives the ability to interact with the mysql database.
import indiv_read_pair_module
import pickle





def main(argv):
	try:
		opts,args = getopt.getopt(argv,'hl:u:p:c:n:m:s:t:',)
	except getopt.GetOptError:
		print "catalog_read_pair.py -l <mysql host (default localhost)> -u <mysql user (default root)> -p <mysql password> -s <Stacks database name> -c <path to catalog dictionaries (output from catalog_read_pair.py)>  -t <table to update in database (default = catalog_index)>"
		sys.exit(2)
			
	mysqlhost = 'localhost'
	mysqluser = 'root'
	mysqlpasswd = ''
	stacksdb = ''
	tabletoupdate = "catalog_index"
	catalogdir = '.'

	for opt, arg in opts:
		if opt == "-l":
			mysqlhost = arg
		elif opt == "-h":
			print "catalog_read_pair.py -l <mysql host (default localhost)> -u <mysql user (default root)> -p <mysql password> -s <Stacks database name> -c <path to catalog dictionaries (output from catalog_read_pair.py)> -n <number of samples> -t <table to update in database (default = catalog_index)>"
			sys.exit(2)
		elif opt == "-u":
			mysqluser = arg
		elif opt == "-p":
			mysqlpasswd = arg
		elif opt == "-c":
			catalogdir = arg
		elif opt == "-t":
			tabletoupdate = arg
		elif opt == "-s":
			stacksdb = arg

			
# mysqlpasswd = "hofi"
# stacksdb = "HFdenovo_PairTesting2_radtags"
# numsamples=140
# tabletoupdate = "catalog_index"
	catalogreadfile = "%s/single_val_catalog.p"%catalogdir
	catalogpairfile = "%s/catalogsinglepairsdict.p"%catalogdir

	##########################################################################################
	"""
	This part of the script will create a new catalog, with a single value for read, 1, 2, or 3 if both 1 and 2.  It will also print a warning for that locus if it is found on both reads 1 and 2.  Finally, it adds a new column to the mysql 'catalog_index' table called 'readend', and goes through to add each read designation for all cat_ids in the catalog. Note that the single catalog is saved as a pickle dump "single_val_catalog.p"
	"""

	singlecatalog = pickle.load(open(catalogreadfile,"r"))
	singlecatalogpairs = pickle.load(open(catalogpairfile,"r"))


	MyConnection = MySQLdb.connect( host = mysqlhost, user = mysqluser, \
										passwd = mysqlpasswd, db = stacksdb)
	MyCursor = MyConnection.cursor()

	#Add columns to table to interest for read information (which end, and pair) if they do not yet exist.
	try:
		SQL = """ALTER TABLE `%s` ADD COLUMN `readend` INT(1) NULL  AFTER `ref_id` ;"""%tabletoupdate
		MyCursor.execute(SQL)
	except:
		pass
	
	try:
		SQL = """ALTER TABLE `%s` ADD COLUMN `readpair` INT(1) NULL  AFTER `readend` ;"""%tabletoupdate
		MyCursor.execute(SQL)
	except:
		pass

	#Add read end information to mysql table.

	for loc in range(len(singlecatalog)):
		if singlecatalog[loc][1] != 0:
			SQL2 = """UPDATE `%s` SET `readend`=%d WHERE `cat_id`=%d;"""%(tabletoupdate,singlecatalog[loc][0],singlecatalog[loc][1])
			MyCursor.execute(SQL2)
			MyConnection.commit()
		else:
			pass
		
	#Add read pair information to mysql table.

	for loc in singlecatalogpairs:
		if singlecatalogpairs[loc] != 0:
			SQL3 = """UPDATE `%s` SET `readpair`=%d WHERE `cat_id`=%d;"""%(tabletoupdate,singlecatalogpairs[loc],loc)
			MyCursor.execute(SQL3)
			MyConnection.commit()
		else:
			pass
		
	MyCursor.close()
	MyConnection.close()

if __name__ == "__main__":
		main(sys.argv[1:])