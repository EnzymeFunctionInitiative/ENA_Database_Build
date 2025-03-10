
import sys
import csv

db_type = sys.argv[1]

# handle import of the user-specified db_type's library
# both libraries seem to share the same api so just assign the specified 
# db_type's python library as "sql_lib" 
if db_type.lower() == "mysql":
    import mysql.connector as sql_lib
else:
    import sqlite3 as sql_lib

###############################################################################
# READ TSV AND MAKE SQL DATABASE FILE
###############################################################################

def create_db(tsv: str, columns: list, output: str) -> int:
    """
    """
    nRows = 0
    # spin up the new database
    with sql_lib.connect(output) as conn:
        cursor = conn.cursor()
            
        # Create table dynamically based on columns
        create_table_query = f"CREATE TABLE IF NOT EXISTS ENA ("
        for col in columns:
            create_table_query += f"{col[0]} {col[1]}, "
        create_table_query = create_table_query.rstrip(', ') + ")"
        # create the columns
        cursor.execute(create_table_query)
            
        # Prepare data for insertion
        placeholders = ', '.join(['?'] * len(columns))
        insert_query = f"INSERT INTO ENA VALUES ({placeholders})"

        # read the tsv file
        with open(tsv, 'r') as file:
            tsv_reader = csv.reader(file, delimiter='\t')
            
            # Insert data into the table at the cursor
            for row in tsv_reader:
                # try to execute the INSERT
                try:
                    cursor.execute(insert_query, row)
                # but, if a failure occurs, print out the error and row that 
                # fails, and then move onto the next row.
                except sql_lib.Error as e:
                    print(type(e))
                    print(e)
                    print(row)
                    pass
        
        # commit the changes to the database
        conn.commit()
        # close cursor
        cursor.close()


###############################################################################
# MAIN
###############################################################################

if __name__ == "__main__":
    # grab input arguments
    tsv = sys.argv[2]
    output = sys.argv[3]
    # hard-coded column names
    columns = [
        ("enaId","TINYTEXT"), 
        ("uniprotId","TINYTEXT"), 
        ("seqCount","SMALLINT"), 
        ("chr","TINYINT"), 
        ("dir","TINYINT"), 
        ("start","INT"), 
        ("stop","INT")
    ]
    # parse the tsv file into the sql database file
    nRows = create_db(tsv, columns, output)
    # log it
    print(f"Wrote {output} with 'ENA' table with {nRows} rows and columns:\n'enaId'*\n'uniprotId'*,\n'seqCount',\n'chr',\n'dir',\n'start',\n'stop'")


