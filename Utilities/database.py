"""Database Class for connecting to sqlachemy
"""
import threading
import os
from sqlalchemy import create_engine, Column, Integer,  Table, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, mapper
from sqlalchemy.exc import OperationalError
from sqlalchemy.ext.declarative import declarative_base

__all__ = ["Database"]

_schema_version = 1
verbose=False
Base = declarative_base()

class column_dict:
    ''' A class to wrap around a dictionary that stores attributes for creating arbitrary column and, if relevant, the foreign key and relationship information associated with 
        that column object. This is basically a poor mans typedef for the structure'''
    def __init__(self, colName, colType, foreignKey=False, foreignClass=None, foreignTableName = '', foreignColName = '', primaryJoinString='', relationshipName=''):
        # must define column Name and Column type. Defaults to no foreign key. 
        self.columnDict = {
        'colName': colName,              #  a string representing the name of the column/field in the table and the attribute name in the class
        'colType': colType,            # Integer, Float, String or pickleType
        'foreignKey': foreignKey,
        'foreignClass': foreignClass,       # the name of the foreign class if there is a relationship
        'foreignTableName': foreignTableName,     # the name of the foreign table
        'foreignColName': foreignColName,       # the name of the column in the foreign table  
        'primaryJoinString': primaryJoinString,    # the string describing the relationship between the two classes.
        'relationshipName': relationshipName      # the name of the attribute describing the relationship in the current class. 
        } 


class Database(object):
    """Database storage class

    The Database class handles the connection to the database. 
    
    It has functions to create new objects and tables. The objects are persistent in the database and exist as
    soon as the Database class in connected to the database. If any value in the objects is changed,
    the changes are automatically persistent in the database (TODO: be careful, check commit transactions, ...)

    Database uses SQLAlchemy to connect to the database. Check the web page for available connectors. Unless
    you know better, the standard sqlite should be used. The database can be generated in memory (default) or
    written to a file if db is specified when creating the class.
    
    Parameters
    ----------
    db : string, optional
        filename of new or existing database to connect to.  
        default creates new database in memory.
    connect_string : string, optional
        connection string, default is sqlite database
    createdb : boolean, optional
        create database if not exists, default is true

    Attributes
    ----------
    engine : sqlalchemy database engine
    session : sqlalchemy session

    Examples
    --------

    """

    engine = None
    session = None
    connection = None
    # metadata = None # declared in declarative Base 
        
    def __init__(self, db=":memory:", connect_string='sqlite:///%s', createdb=True):

        if not os.path.isfile(db) or db == ":memory:":
            newfile = True
            if not createdb:
                raise IOError("createdb is False, but database does not exist (%s)" % db)
        else:
            newfile = False

        # set up the engine which will manage the backend connection to the database
        self.engine = create_engine(connect_string % db, echo=verbose)
        
        if not newfile and not self._is_valid_database():
            raise IOError("existing file (%s) is not a valid database." % db)
        
        # set up the tables and check the schema version
        if newfile:
            self._set_schema_version()
        self._check_schema_version()
        self._update_schema()

        # set up the session which will manage the frontend connection to the database
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        
        self.lock = threading.Lock()
        self.connection = self.engine.connect()

    def _is_valid_database(self):
        conn = self.engine.connect()
        result = True
        # if (<test Condition for valid database>):
        #    result = False
        conn.close()
        return result

    def _set_schema_version(self):
        conn = self.engine.connect()
        conn.execute("PRAGMA user_version = %d;"%_schema_version)
        conn.close()

    def _update_schema(self):
        conn = self.engine.connect()
        # self.metadata = MetaData(bind=self.engine)
        Base.metadata.create_all()
        conn.close()

    def _check_schema_version(self):
        conn = self.engine.connect()
        result=conn.execute("PRAGMA user_version;")
        schema = result.fetchone()[0]
        result.close()
        conn.close()
        if _schema_version != schema:
            raise IOError("database schema outdated, current (newest) version: "
                          "%d (%d)"%(schema, _schema_version))

    def add_records_to_table(self, cls, fields, valuesArray, commit=True):
        ''' Creates an array of instances of single type of class cls, and then populates the member variables
        of each of those objects using the valuesArray. 
        
        The names of the member variables to populate are held in the fields list. 
       
        The list of objects are then added to the table to which the class refers.
        
        Assumes the user has set up the schema properly using the make_table member function. 
        which uses the cls definition and fields variables to generate the table.

        This structure allows completely arbitrary tables to be created at runtime, so the user
        can add a parameter to the database easily via the controlling text file. 
        Obviously it is up to the user to ensure that the names in the fields list match the 
        class definition. '''
        success=False
        
        # create a list to store the instances of the class objects
        objList = []
        
        # for every record in valuesArray
        for values in valuesArray:
            # create a record object from the empty class
            obj = cls()
            
            # loop throught the fields and values in the records
            for field, value in zip(fields, values):
                # set the object attributes
                setattr(obj, field, value)
            # add object to list
            objList.append(obj)
        
        # grab the lock
        self.lock.acquire()
    
        # add all objects to database
        self.session.add_all(objList)
    
        # commit if requested (default is True)
        try:
            if commit:
                self.session.commit()
            success=True
        except OperationalError:
            print("\nWarning: Sqlalchemy Operational Error. Unable to write data to table. Table may not exist.\n")
            success=False
        
        # release the lock
        self.lock.release()
        
        return success

    def make_table(self, cls, colDictList):
        ''' Adds column and relationship objects to a table class that inherits the base declarative system. Instantiation creates 
        the schema that is defined in a list of columnDict objects.  
        
        colDictList is a list of objects each describing a column and, if the column is a foreign key, a relationship. Each object is a wrapper 
        for a dictionary which has the following keyWords.
        
        {
        'colName':   a string representing the name of the column/field in the table and the attribute name in the class
        'colType':  Integer, Float, String or pickleType
        'foreignKey: a boolean stating whether the column is a foreign key
        'foreignClass':  the name of the foreign class if there is a relationship
        'foreignTableName': the name of the foreign table
        'foreignColName': the name of the column in the foreign table  
        'primaryJoinString': the string describing the relationship between the two classes.
        'relationshipName':  the name of the attribute describing the relationship in the current class. 
        } 
        
        '''
        # loop through colDict List and add appropriate information for each column. 
        # Using this kind of dictionary makes for much more readable code.
        # Can over load this function to use more than just primary join relationship between tables, as necessary. 
        # Assuming primary join will be enough in most cases. 
        for columnObj in colDictList:
            colDict = columnObj.columnDict
            if colDict['foreignKey']:
                # add foreign key object to the column object and a relationship object to the class 
                setattr(cls, colDict['colName'], Column(colDict['colType'], ForeignKey(colDict['foreignTableName'] + '.' + colDict['foreignColName'])) )
                setattr(cls, colDict['relationshipName'], relationship(colDict['foreignClass'], primaryjoin=colDict['primaryJoinString']))
            else:
                # No foreign class specified so just adding a regular column
                setattr(cls,colDict['colName'], Column(colDict['colType']) )
                
        # use declarative system to generate the schema                    
        Base.metadata.create_all(self.engine)
    
    def make_Table_without_Base_declarative(self, cls, name, primary_id, columnNames, columnTypes, foreignKeyInfo=None):
        ''' Creates a table object in the metadata which is then mapped to the class cls,
        allowing record objects to be created using that class.  The list of column Names and column Types must be
        the same length. The types can be either Integer, Float, String or pickleType. If foreignKeyInfo is supplied for a field 
        then the ForeignKey keyword is set for that field and the relationships are generated after the table has been instantiated.
        
        Each entry in foreignKeyInfo is either a blank string '', or a List = 
        [ foreignTableId, TableAttributeRelationshipName,  ForeignClassName, ForeignClassAttribute ]
         
        '''
        # assume that no relations exist.
        relationsExist = False
        
        # if no foreignKeyInfo is provided we are making a straight table with no joins, so just add dummy information into the foreignKeyInfo list.
        if not foreignKeyInfo:
            foreignKeyInfo = [''] * len(columnNames)

        # check that all the input information lengths match
        if len(columnNames)==len(columnTypes)==len(foreignKeyInfo):
            # create array for output
            columnCommands = []
            
            # loop through column definitions and add appropriate information for each column
            for colName, colType, forKey in zip(columnNames, columnTypes, foreignKeyInfo):
                
                # Generate column objects.
                # only add foreignKey information if there is data in that record.
                if forKey: # false if forKey==''
                    relationsExist = True # at least one entry in the foreignKey list has data in it. 
                    columnCommands.append( Column(colName, colType, ForeignKey(forKey[0] + '.' + forKey[2]) ) )
                    # forkey[0] = foreign table name
                    # forkey[2] = Foreigh class attribute/colName
                    #
                    # Example from some other code:
                    # _minimum_id=Column(Integer, ForeignKey('tbl_minima._id'))
                                        
                else:
                    columnCommands.append( Column( colName, colType ) )

            # create a table with the primary id, and the List of column objects converted to a tuple of arbitrary length. 
            t = Table( name,
                       self.metadata,
                       Column(primary_id, Integer, primary_key=True),
                       *tuple(columnCommands) )

            # establish the relationship attributes for the table if there are any
            if relationsExist:
                # loop through the column information and create relationships with other tables as defined in the foreignKeyInfo 
                for colName, forKey in zip(columnNames, foreignKeyInfo):
                    if forKey:
                        setattr(t, forKey[3], relationship(forKey[1].__name__, primaryjoin=forKey[1].__name__ + '.' + forKey[2] + '==' + cls.__name__ + "." + colName))
                        # example:  
                        # t.currentClassRelationshipName = relationship( ForeignClassName, primaryjoin="ForeignClassName.ForeignClassAttribute==TableClassName.colName")
                        # forkey[0] = foreign table name
                        # forkey[1] = Foreign class name
                        # forkey[2] = Foreigh class attribute/colName
                        # forkey[3] = t.current class relationship name 
                        #
                        # Minimal Example from some other code:
                        # Class symmetry():
                        #    _minimum_id = Column(Integer, ForeignKey('tbl_minima._id')) 
                        #    minimum=relationship("Minimum", primaryjoin="Minimum._id==symmetry._minimum_id")
            
            self.metadata.create_all()
            mapper(cls, t)
        else:
            raise Exception("Inconsistent name, type and foreignKeyInfo array lengths while creating database schema.")
