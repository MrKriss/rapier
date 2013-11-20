"""
Created on 19 Mar 2013

class to interface with a sqlite database

@author: musselle
"""
import os
import sys
import cPickle as pkl
from subprocess import PIPE, Popen, STDOUT
import tempfile, csv
import glob

import gzip
import numpy as np
from Bio import SeqIO

import sqlite3


class SQLdatabase(object):
    """ Class to handle all python communication with a sqlite database file 
    
    Note on using with statement as a context manager:
    
    # Successful, con.commit() is called automatically afterwards
    with con:
        con.execute("insert into person(firstname) values (?)", ("Joe",))

    # con.rollback() is called after the with block finishes with an exception, the
    # exception is still raised and must be caught
     
    """
    
    def __init__(self, dbfile="default.db", recbyname=True):
        """ 
        recbyname - sets returned records by select to be callable by column names 
        """        
        if  os.path.exists(dbfile):
            print 'Database found with matching file name.'
            print 'Connecting to database {0}'.format(dbfile)
        else:
            print 'Creating new Database file: {0}'.format(dbfile) 
                
        # Stored Vars
        self.con = sqlite3.connect(dbfile)
        self.dbfilepath = os.path.abspath(dbfile) 
        self.recbyname = recbyname
        self.tables = []

        # Update list of tables in database
        results = self.con.execute("select name from sqlite_master where type = 'table';")
        # Get list of tables         
        x = results.fetchall()
        x = [i[0] for i in x]
        self.tables = x

        if recbyname: 
            self.con.row_factory = sqlite3.Row
            print 'Setting Row_factory to named Rows' 

    def connect(self):
        ''' Open connection to the database '''
        self.con = sqlite3.connect(self.dbfilepath)
        if self.recbyname: 
            self.con.row_factory = sqlite3.Row
            print 'Setting Row_factory to named Rows' 
    
    def close(self):
        ''' Close connection to database'''
        self.con.close()

    def display(self, cmd, num=0, *args):
        """ Print out the records returned from the database querry 
        
        num = max number to display. default 0 is all records returned.
        """    
        with self.con as con:     
            cur = con.cursor()
            # SELECT column_name(s) FROM table_name
            cur.execute("SELECT " + cmd, *args)
            if num: 
                records = cur.fetchmany(num)
            else:  
                records = cur.fetchall()      
        
        if self.con.row_factory == sqlite3.Row:
            row_names = records[0].keys()
            print '  '.join(row_names)
            for row in records:
                print row
        else: 
            for row in records:
                print row
        
    def select_as_tuples(self, cmd, *args):
        ''' Select records from the database returned as tuple of tuples. '''
        tuples = []
        with self.con as con:     
            cur = con.cursor()
            # SELECT column_name(s) FROM table_name
            cur.execute("SELECT " + cmd, *args)
            records = cur.fetchall()
            for row in records:
                tuples.append(row)
        return tuples
         
    def update_binary(self, obj, col, target, value, table):
        """ Store a binary object to the database 
        Add Obj to field 'col' where field 'target' = 'value'
        """

        with self.con as con:        
            cur = con.cursor()
            b = sqlite3.Binary(pkl.dumps(obj))             
            # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value
            cur.execute('UPDATE {0} SET {1}=? WHERE {2}=?'.format(table, col, target), (b, value))
        
    def insert_binary(self, obj, col, table):
        """ Store a binary object to the database 
        Insert Obj to field 'col' on a new row. Return Rowid 
        """

        with self.con as con:        
            cur = con.cursor()
            b = sqlite3.Binary(pkl.dumps(obj))             
            # INSERT INTO table_name(column_name) VALUES (?) 
            cur.execute('INSERT INTO {0}({1}) VALUES (?)'.format(table, col), (b,))
        
        return cur.lastrowid
        
    def get_binary(self, col, target, value, table):
        """ Retrieve binary object in the database 
        Get Obj from field 'col' where field 'target' = 'value'
        """
        
        with self.con as con:        
            cur = con.cursor()
            cur.execute('SELECT {0} FROM {1} WHERE {2}=?'.format(col, table, target), (value,))        
            pickled_data = str(cur.fetchone()[0])
            data = pkl.loads(pickled_data)
        
        return data   
    
    def fastinsert(self,data,table):
        """
         Function for inserting data into a SQLite database via the sqlite3 command line client using .import 
         
         # Roughly 5x faster at dumping data than using the Python SQLite bindings
         
         # Your data should be a list of equal length tuples
         data = [(x,) for x in range(100)]
         
         # Source: https://gist.github.com/MattOates/6195743
         
        """ 
 
        #Name of the pipe, use tempfile to create some random filename usually in /tmp
        data_pipe = tempfile.mktemp('datapump')
        #Create the actual pipe 'file' where the OS knows you wanted a pipe type thing
        os.mkfifo( data_pipe, 0644 )
 
        #Create a child process to run in parallel
        child = os.fork()
 
        #If child is 0 we are the child
        if child == 0:
            #Have the child send the command to sqlite3
            #Create a forked process running sqlite3 with a pipe that we can send commands to
            sqlite3 = Popen(['sqlite3', self.dbfilepath], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
 
            #Tell sqlite3 to import from our named pipe
            db = sqlite3.communicate(input='.separator "\\t"\n.import %s %s\n.quit\n' % (data_pipe,table))[0]
            
            #The child exits so we stop waiting
            sys.exit(1)
 
        else:
            #File handle to pipe to write data into table
            data_pipe_write_handle = open(data_pipe,'w')
            #Send data down the pipe
            writer = csv.writer(data_pipe_write_handle, delimiter="\t")
            writer.writerows(data)
            data_pipe_write_handle.close()
            #Wait for SQLite3 to finish importing, waiting on all child processes
            os.wait()
            #Remove the named pipe file we created because its junk and we dont want a clash
            os.unlink(data_pipe)

    
