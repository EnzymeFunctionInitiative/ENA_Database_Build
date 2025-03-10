
import sys
import csv

db_type = sys.argv[1]

# handle import of the data type
if db_type.lower() == "mysql":
    import mysql.connector as sql_lib
else:
    import sqlite3 as sql_lib

###############################################################################
# FUNCTIONS
###############################################################################


###############################################################################
# MAIN
###############################################################################

if __name__ == "__main__":
    # grab input arguments
    db1 = sys.argvs[2]
    db2 = sys.argvs[3]

    # QUERIES:
    #   - Grab enaId column from both and find disagreements (left and right joins?)
    #   - Grab (enaId, uniprotId) column pairs and find disagreements (join too?); for agreements, check whether their other columns match exactly
    #   - Grab uniprotId column and check for enaId disagreements
    #   - ...



