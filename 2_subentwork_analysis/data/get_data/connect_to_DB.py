# LCSB EPFL 2020, Homa Mohamamdi
# These functions are used to connect to database

import mysql.connector
from mysql.connector import errorcode


def connectToDatabase():
	try:
		connection = mysql.connector.connect(
		  host="lcsbsrv2.epfl.ch",
		  user="Omid",
		  passwd="$c9PwwdhH2hveE[Y",
		  database="LCSB2"
		)
	except mysql.connector.Error as err:
		if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
		    print('Invalid credential. Unable to access database.')
		elif err.errno == errorcode.ER_BAD_DB_ERROR:
		    print('Database does not exists')
		else:
		    print('Failed to connect to database')
		quit()
	else:
		return connection


	
def fetch_LCSB_ID(_connection):
        print("fetch LCSB ID of compounds")
        query = _connection.cursor()
        print("Execute query....")
        query.callproc('LCSB2.EX_ID_to_refV_Inchikey(')
        for result in query.stored_results():
                print("Query done")
                table = result.fetchall()
        print("Result size: ", len(table))
        query.close()
        return [i[0] for i in table]

def closeConnection(_connection):
	_connection.close()
	print("Connection to DB closed")
    
def fetchResults(_connection, comp_IDs,n_comps):
   print("fetch results")
   query = _connection.cursor()
   print("Execute query....")
   args = (comp_IDs, n_comps)
   query.callproc('LCSB2.compounds_refV_2_rxn_refV', args)
   for result in query.stored_results():
      print("Query done")
      table = result.fetchall()
   print("Result size: ", len(table))
   query.close()
   return table


def fetchResults_comps_2_rxn(_connection, comp_IDs,n_comps):
    #comp_IDs: a string contating LCSB comps IDs, comma seperated, example:"1469424929,1467866080,1467866192,1467879801"
    #n_comps: a string defining number of compounds in comp_IDs, example: "4"
    #output: ID of reaction contating input compounds in LCSB
    # to visualize reactions, use following url (example: 2592969210)
    # https://lcsb-databases.epfl.ch/Graph2/loadReactionList/2592969210
    query = _connection.cursor()
    print("Execute query...")
    query.execute("select rxn.refV as rxn_ID from "
                  "(select link.refV from (select distinct refV, refV_2_CID(refV) as cid ,refV_2_IID(refV) as iid from"
                  " LCSB2.main where refV  in ("+comp_IDs+")) comp "
                  "inner join LCSB2.main link on link.cid=comp.cid and link.iid=comp.iid "
                  "group by link.refV having count(*)="+n_comps+") link "
                  "inner join LCSB2.main rxn on rxn.uid=link.refV "
                  "where rxn.cid=8 and LCSB2.lookup_beingBalance(rxn.refV)=0 ;")

    table = query.fetchall()
    print("Result size: ", len(table))
    query.close()
    return table
