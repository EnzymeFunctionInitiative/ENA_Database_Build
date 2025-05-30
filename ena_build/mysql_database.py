
import mysql.connector

class IDMapper:
    def __init__(self, database_params: dict, db_name: str):
        """ 
        Set up a connection to the user-specified MySQL database 
        
        Parameters
        ----------
            database_params 
                dict or configparser.ConfigParser obj. necessary keys or 
                attributes are "user", "password", "host", "port"
            db_name
                string, name of the EFI database to be used to perform
                queries. 

        Attributes
        ----------
            self.dbh
                mysql.connector.MySQLConnection object

        """
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
                print("Something is wrong with your user name or password." 
                        + f"\n{err}")
            elif err.errno == mysql.connector.errorcode.ER_BAD_DB_ERROR:
                print(f"Database {database_params.db_name} does not exist."
                        + f"\n{err}")
            else:
                print(err)
            cnx = ""
        
        self.dbh = cnx

    def close(self):
        """ Close the database connection """
        if self.dbh.is_connected():
            self.dbh.close()

    def reverse_mapping(self, foreign_ids: list, batch_size: int = 1000):
        """
        Query the connected database with foreign_ids and return a dict that
        maps proteinIds (foreign_ids) to uniprotIds. 

        Parameters
        ----------
            foreign_ids
                list of strings, each string is a non-UniProtKB accession ID;
            batch_size
                int, default 1000, max number of query results to be gathered
                at once to avoid OutOfMemory errors while still efficiently
                gathering query results. 

        Returns
        -------
            mapping
                dict, keys are foreign_ids and values are sets of their 
                associated UniProtKB accession IDs
            no_matches
                set of strings, foreign_ids from the input list that do not 
                map to a UniProtKB accession ID in the database.
        
        """
        # make sure the foreign_ids list is not empty
        if not foreign_ids:
            return {}, []
        
        # spin up a cursor on the database handle
        cursor = self.dbh.cursor(dictionary=True)
        
        # prep the output arrays. will add found uniprot_ids to that list and
        # will remove the associated foreign_id from the no_matches list.
        mapping = {}
        no_matches = set(foreign_ids)

        # NOTE: add a loop over foreign_ids to prevent the executed query from
        # being too large

        # prep the sql string for the list of foreign_ids
        # this search is grabbing both foreign_id and uniprot_id from the
        # table if the row's foreign_id is in the list of foreign_ids (input).
        placeholder = ", ".join(["%s"]*len(foreign_ids))
        sql = f"SELECT foreign_id, uniprot_id FROM idmapping WHERE foreign_id IN ({placeholder})"
        
        # execute query: 
        # for each foreign_id found by the IN statement, this query will
        # produce a dict with keys "foreign_id" and "uniprot_id" w/ associated
        # values. 
        # fetchone() will return this dict, fetchall() or fetchmany() will 
        # return a list of dicts associated with the successful queries.
        # if a foreign_id is not found in the table but is in the list of 
        # foreign_ids, no dict will be returned. 
        cursor.execute(sql, foreign_ids)
        
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
                # get the set associated with "foreign_id" if the key is in
                # the mapping dict, else return an empty set
                uniprot_ids = mapping.get(row["foreign_id"],set())
                # add the row's uniprot_id to the set
                uniprot_ids.add(row['uniprot_id'])
                # reassign the foreign_id key's value to the filled set
                mapping[row["foreign_id"]] = uniprot_ids
                # remove the foreign_id from the no_matches list
                no_matches.discard(row["foreign_id"])
        
        # end cursor instances to keep dbh connections minimal
        cursor.close()
        
        return mapping, no_matches


