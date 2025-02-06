
import mysql.connector

class IDMapper:
    def __init__(self, database_params, db_name):
        """ Set up a connection to the user-specified MySQL database """
        try:
            cnx = mysql.connector.connect(
                user=database_params["user"],
                password=database_params["password"],
                host=database_params["host"],
                port=database_params["port"],
                database=db_name,
            )
        except mysql.connector.Error as err:
            if err.errno == mysql.connector.errorcode.ER_ACCESS_DENIED_ERROR:
                print(f"Something is wrong with your user name or password.\n{err}")
            elif err.errno == mysql.connector.errorcode.ER_BAD_DB_ERROR:
                print("Database {database_params.db_name} does not exist.\n{err}")
            else:
                print(err)
            cnx = ""
        
        self.dbh = cnx

    def close(self):
        """ Close the database connection """
        if self.dbh.is_connected():
            self.dbh.close()

    def reverse_lookup(self, foreign_ids, batch_size=1000):
        """
        """
        # make sure the foreign_ids list is not empty
        if not foreign_ids:
            return [], []
        
        # prep the sql string for the list of foreign_ids
        # this search is grabbing both foreign_id and uniprot_id from the table
        # if the row's foreign_id is in the list of foreign_ids (input).
        placeholder = ", ".join(["%s"]*len(foreign_ids))
        sql = f"SELECT foreign_id, uniprot_id FROM idmapping WHERE foreign_id IN ({placeholder})"
        
        # spin up a cursor on the database handle
        cursor = self.dbh.cursor(dictionary=True)
        # execute query: 
        # for each foreign_id found by the IN statement, this query will
        # produce a dict with keys "foreign_id" and "uniprot_id" w/ associated 
        # values. 
        # fetchone() will return this dict, fetchall() or fetchmany() will 
        # return a list of dicts associated with the successful queries.
        # if a foreign_id is not found in the table but is in the list of 
        # foreign_ids, no dict will be returned. 
        cursor.execute(sql, foreign_ids)
        
        # prep the output arrays. will add found uniprot_ids to that list and
        # will remove the associated foreign_id from the no_matches list.
        uniprot_ids = []
        no_matches = set(foreign_ids)

        # instead of fetching results one or all at a time, let's grab a batch
        # with controlled size. This prevents inefficiencies of one at a time
        # while avoiding out of memory concerns if the query returns a huge 
        # number of results
        while True:
            # fetch a set of query results, number of rows dependent on
            # batch_size
            batch = cursor.fetchmany(batch_size)
            
            # if batch is an empty list, we've hit the end of query results so
            # break out of the while loop
            if not batch:
                break
            
            # loop over rows in the query results
            for row in batch:
                # add the uniprot_id to the list
                uniprot_ids.append(row['uniprot_id'])
                # remove the foreign_id from the no_matches list
                no_matches.discard(row['foreign_id'])
        
        # end cursor instances to keep dbh connections minimal
        cursor.close()
        
        return uniprot_ids, no_matches


