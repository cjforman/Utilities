import unittest
import numpy as np
import os
import copy as cp
import Utilities.cartesian as cart
import Utilities.coordSystems as coords
import Utilities.fileIO as fIO
import Utilities.database as dbase
import Utilities.database.column_dict as cd
from sqlalchemy import Column, Float, Integer
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class additionalAsserts(unittest.TestCase):

    def assertAlmostEqualListOfVectors(self, List1, List2, delta=None):
        self.assertEqual(len(List1), len(List2))
        for a, b in zip(List1, List2):                    
            self.assertAlmostEqualVectors(a, b, delta=delta)
    
    
    def assertAlmostEqualVectors(self, vec1, vec2, delta=None):
        for i, j in zip(vec1, vec2):                    
            self.assertAlmostEqual(i, j, delta=delta)    


class TableForTesting(Base):
    __tablename__ = 'tbl_testing'
    _id = Column(Integer, primary_key=True) 

    
class TableForTestingRelationship1(Base):
    __tablename__ = 'tbl_relationship1'
    _id = Column(Integer, primary_key=True) 


class TableForTestingRelationship2(Base):
    __tablename__ = 'tbl_relationship2'
    _id = Column(Integer, primary_key=True) 


class Test_database(additionalAsserts):
    
    def test_dbase_no_table(self):
        if os.path.isfile('test2.sqlite'):
            os.remove('test2.sqlite')
        # create database
        db = dbase.Database(db="test2.sqlite")

        numFields = 10
        numRecords = 20

        import string
        import random
        
        # set up the alphabet
        letters = string.ascii_lowercase
        
        # generate some random field names
        fieldNames = [ ''.join(random.choice(letters) for _ in range(10)) for _ in range(numFields) ]

        # generate the mapping between testTable class and the database using the randomly generated fieldnames
        #db.make_Table(TableForTesting, "tbl_test_table", 'id', fieldNames, [Float]*(numFields))
        
        # set up a value array to go into table
        valArray = np.random.rand(numRecords, numFields)

        # add the valArray to the table - should fail due to absence of table in the database
        success = db.add_records_to_table(TableForTesting, fieldNames, valArray, True) 
        
        if success==False:
            print("sqlAlchemy OperationalError detected. Test successful.\n")
        
        self.assertEqual(success, False)
    
    def test_dbase_basic_table(self):
        if os.path.isfile('test.sqlite'):
            os.remove('test.sqlite')
        # create database
        db = dbase.Database(db="test.sqlite")

        numFields = 10
        numRecords = 20

        import string
        import random
        
        # set up the alphabet
        letters = string.ascii_lowercase
        
        # create a list of column Defining objects
        colDictList = []  
        for _ in range(numFields):
            colDict = cd( ''.join(random.choice(letters)), Float, foreignKey=None, foreignClass=None, foreignTableName='', foreignColName='', primaryJoinString='', relationshipName='')
            colDictList.append( cp.copy(colDict) )

        fieldNames = [colDict.columnDict['colName'] for colDict in colDictList]

        # generate the mapping between testTable class and the database using the randomly generated fieldnames
        db.make_Table(TableForTesting, colDictList)
        
        # set up a value array to go into table
        valArray = np.random.rand(numRecords, numFields)
        
        # add the valArray to the table
        db.add_records_to_table(TableForTesting, fieldNames, valArray, True) 
        
        # query the database to retrieve the table
        retrievedValues = db.session.query(TableForTesting).all()

        # assert we have the right number
        self.assertAlmostEqual(len(retrievedValues), numRecords, delta=1e-11)

        # check the retrieved records perefectly match the input values and randomly chosen fieldnames            
        for inputRecord, retrievedRecord in zip(valArray, retrievedValues):
            for value, fieldName in zip(inputRecord, fieldNames):
                self.assertAlmostEqual(value, getattr(retrievedRecord, fieldName), delta=1e-11)

    def test_dbase_relation_between_tables(self):
        if os.path.isfile('test_relationship.sqlite'):
            os.remove('test_relationship.sqlite')
        # create database
        db = dbase.Database(db="test_relationship.sqlite")

        numFieldsTable1 = 5
        numFieldsTable2 = 5
        numRecords = 20

        import string
        import random
        
        # set up the alphabet
        letters = string.ascii_lowercase

        # create a list of column Defining objects
        colDictListTable1 = []  
        for _ in range(numFieldsTable1):
            colDict = cd( ''.join(random.choice(letters)), Float, foreignKey=None, foreignClass=None, foreignTableName='', foreignColName='', primaryJoinString='', relationshipName='')
            colDictListTable1.append( cp.copy(colDict) )

        fieldNamesTable1 = [colDict.columnDict['colName'] for colDict in colDictListTable1]

        # generate the mapping between testTable class and the database using the randomly generated fieldnames
        db.make_Table(TableForTesting, colDictListTable1)
        
        # set up a value array to go into table
        valArray = np.random.rand(numRecords, numFieldsTable1)
        
        # add the valArray to the table
        db.add_records_to_table(TableForTesting, fieldNamesTable1, valArray, True) 



        # ***** UP to Here *******
        # create a list of column Defining objects
        colDictListTable2 = []  
        for _ in range(numFieldsTable1):
            colDict = cd( ''.join(random.choice(letters)), Float, foreignKey=None, foreignClass=None, foreignTableName='', foreignColName='', primaryJoinString='', relationshipName='')
            colDictListTable1.append( cp.copy(colDict) )





        # pick the number of fields to join between the tables - must be at least 1 and no more than numfieldsTable1
        numJoinedFields = random.randint(1, numFieldsTable1 - 1 ) 

        # select the fields from table 1 to join to randomly in table 2. Put them in table2 in the first N columns. No repetitions 
        joinedFields = []
        while len(joinedFields) < numJoinedFields:
            fieldNum = random.randint(1, numFieldsTable1 - 1 )
            if not fieldNum in joinedFields:
                joinedFields.append(fieldNum)  

        # nucleate the joined name array with the field names from table1 
        fieldNames2_joined = [ fieldNames1[ joinedField ] for joinedField in joinedFields ]  
        
        # generate the info needed about the foreign table in the first fields of foreignKeyInfo    
        foreignKeyInfo = [ [table1Name, TableForTestingJoins1, joinedFieldname, 'rel_' + joinedFieldname] for joinedFieldname in fieldNames2_joined ] 

        # add the remaining fields from table 2. 
        fieldNames2_joined += fieldNames2

        # append the empty foreign key info markers for remainging fields.  
        foreignKeyInfo += [''] * numFieldsTable2
        
        # generate the mapping between testTable2 class and the database using the randomly generated fieldnames, and the foreignKeyInfo
        db.make_Table(TableForTestingJoins2, table2Name, 'id', fieldNames2_joined, [Float]*(len(fieldNames2_joined)), foreignKeyInfo=foreignKeyInfo )
                
        # set up a value array to go into table1
        valArray1 = np.random.rand(numRecords, numFieldsTable1)
        valArray2 = np.random.rand(numRecords, numFieldsTable2)

        # retrieve data        
        joinedValArray = [ list1 + list2 for list1, list2 in zip(valArray1, valArray2) ]
        
        # add the valArray to the table 1
        db.add_records_to_table(TableForTestingJoins1, fieldNames1, valArray1, True) 
        db.add_records_to_table(TableForTestingJoins2, fieldNames2, valArray2, True) 

        # query the database to retrieve the values from table2, which correspond to the quizzed values in Table 1 
        retrievedValuesArray = []
        for fieldName in fieldNames1:
            retrievedValuesArray.append(self.session.query(TableForTestingJoins2).filter(getattr(TableForTestingJoins2, fieldName)==getattr(TableForTestingJoins1, fieldName)).all())
        

        for retrievedValues in retrievedValuesArray: 
            # assert we have the right number
            self.assertAlmostEqual(len(retrievedValues), numRecords, delta=1e-11)

            # check the retrieved records perfectly match the input values and randomly chosen fieldnames            
            for inputRecord, retrievedRecord in zip(joinedValArray, retrievedValues):
                for value, fieldName in zip(inputRecord, fieldNames2_joined):
                    self.assertAlmostEqual(value, getattr(retrievedRecord, fieldName), delta=1e-11)



class Test_fileIO(additionalAsserts):

    def test_writeTextFile(self):

        lines = ["Test String 1\n", "Test String 2"]
        filename = "testWriteTextFile.txt"
        # make sure the test file doesn't already exist 
        try:
            os.remove(filename)
        except:
            pass
        
        # use the test function to write the data to file
        goodWrite = fIO.writeTextFile(lines, filename)
        
        # assume the file does exist
        noFile = False
        
        # attempt to read the textfile back in
        try:
            with open(filename, 'r') as fh:
                linesIn = fh.readlines()
                fh.close()
        except FileNotFoundError:
            # if file doesn't exist then write function likely failed 
            noFile = True
        
        # remove the file for the next time round
        try:
            os.remove(filename)
        except:
            pass

        # To pass the test the write function must return true, 
        # the read in data must match the test data
        # and the written file had to have been found
        self.assertEqual(linesIn, lines)
        self.assertEqual(False, noFile)
        self.assertEqual(True, goodWrite)

        
    def test_readTextFile(self):

        # generate the test data        
        lines = ["Test String 1\n", "Test String 2"]
        linesIn = []
        filename = "testReadTextFile.txt"

        # ensure file does not exist already
        try:
            os.remove(filename)
        except:
            pass

        # assume that an error will occur
        noError = False
        
        # assume file does exist
        noFile = False

        try:
            # open the file for write 
            with open(filename, 'w') as fh:
                
                # write the test data
                for line in lines:
                    fh.write(line)
                
                #close the file
                fh.close()

                # if we get to this point then the write function wrote something.                
                noError = True
        
                # attempt to read in the data        
                linesIn = fIO.readTextFile(filename)
                

        except FileNotFoundError:
            # if the file wasn't found by the read file raise a flag
            noFile = True

        # delete test file
        try:
            os.remove(filename)
        except:
            pass
        
        # to pass the test the data readin from the file must match the output data
        # and all the files exists flags must be correct. 
        self.assertEqual(linesIn, lines)
        self.assertEqual(True, noError)
        self.assertEqual(False, noFile)


class Test_CoordSystems(additionalAsserts):
    
    def setUp(self):
        self.maxRndTests = 10000
    
    
    def test_sphericalPolar2XYZ(self):
        r = 4
        theta = np.pi/4
        phi = np.pi/4
        x = 2
        y = 2
        z = 2 * np.sqrt(2.0)
        inp = np.array([r, theta,phi])
        out = np.array([x, y, z])
        self.assertAlmostEqualVectors(out, coords.sphericalPolar2XYZ(inp), 1e-12)

    def test_bondAngleDihedral2TNB(self):
        inp = np.array([2.0, 3 * np.pi/4, np.pi/4])
        out =  np.array([np.sqrt(2), 1, 1])
        self.assertAlmostEqualVectors(coords.bondAngleDihedral2TNB(inp), out, 1e-12)

    def test_polarToUnitSphereXYZ(self):
        theta = np.pi/4;
        phi = np.pi/4;
        self.assertAlmostEqualVectors(coords.polarToUnitSphereXYZ(theta, phi), np.array([0.5, 0.5, np.sqrt(2.0)/2.0]), 1e-12)
        
    def test_XYZ2SphericalPolar(self):
        r = 4
        theta = np.pi/4
        phi = np.pi/4
        x = 2
        y = 2
        z = 2 * np.sqrt(2.0)
        inp = np.array([x, y, z])
        out = np.array([r, theta,phi])
        self.assertAlmostEqualVectors(out, coords.XYZ2SphericalPolar(inp), 1e-12)
    
    def test_cyl2XYZ(self):
        x = np.sqrt(2.0)
        y = np.sqrt(2.0)
        z = 2.0
        r = 2.0
        phi = np.pi/4
        inp = np.array([r, phi, z])
        out = np.array([x, y, z])
        self.assertAlmostEqualVectors(out, coords.cyl2XYZ(inp), 1e-12)

    def test_XYZ2Cyl(self):
        x = np.sqrt(2.0)
        y = np.sqrt(2.0)
        z = 2.0
        r = 2.0
        phi = np.pi/4
        inp = np.array([x, y, z])
        out = np.array([r, phi, z])
        self.assertAlmostEqualVectors(out, coords.XYZ2Cyl(inp), 1e-12)

    def test_azToRhoCyl(self):
        az = np.pi/4;
        r_x = 3;
        r_y = 2;
        rho = r_x * r_y /np.sqrt( (r_y * np.cos(az))**2  +  (r_x * np.sin(az))**2  )
        self.assertAlmostEqual(rho, coords.azToRhoCyl(az, r_x, r_y), 1e-12)

    def test_ellipsoidPolarUVW(self):
        theta = np.pi/4
        phi = np.pi/4
        r_x =3 
        r_y = 4
        r_z =5 
        rho = 3.9714203539289703
        dirn = np.array([0.5, 0.5, np.sqrt(2)/2.0])
        rho_out, dirn_out = coords.ellipsoidalPolarUVW(theta, phi, r_x, r_y, r_z)
        self.assertAlmostEqual(rho, rho_out, 1e-12)
        self.assertAlmostEqualVectors(dirn, dirn_out, 1e-12)
     
    def test_ellipsoidPolarToXYZ(self):
        theta = np.pi/4
        phi = np.pi/4
        r_x =3 
        r_y = 4
        r_z =5 
        rho = 3.9714203539289703
        dirn = rho * np.array([0.5, 0.5, np.sqrt(2)/2.0])
        self.assertAlmostEqualVectors(dirn, coords.ellipsoidalPolarToXYZ(theta, phi, r_x, r_y, r_z), 1e-12)      
        
        
    def test_TNB2XYZ(self):
        s2 = np.sqrt(2.0)
        T = np.array([1.0, 1.0, 0])/s2
        N = np.array([0.0, 0.0, 1.0])
        B = np.array([1.0, -1.0, 0.0])/s2
        TNBFrame= np.array( [T, N, B] ) 
        posTNB = np.array([3.0, 4.0, 5.0])
        outVec = np.array([3.0/s2 + 0.0 + 5.0/s2, 3.0/s2 + 0.0 - 5.0/s2, 0.0 + 4.0 + 0.0])
        self.assertAlmostEqualVectors(outVec, coords.TNB2XYZ(TNBFrame, posTNB), 1e-12)

    def test_isColinear1(self):
        v1= np.array([1.0, -1.0, 3.0])
        self.assertEqual(True, coords.isColinear(v1, 3.0 * v1, 1e-12),'Two colinear vectors reported as not co-linear.')

    def test_isColinear2(self):
        v1= np.array([1.0, -1.0, 3.0])
        v2= np.array([1.0, 1.0, 3.0])
        self.assertEqual(False, coords.isColinear(v1, v2, 1e-12),'Two not colinear vectors reported as co-linear.')
    
    def test_isEqual1(self):
        v1= np.array([1.0, -1.0, 3.0])
        self.assertEqual(True, coords.isEqual(v1, v1, 1e-12),'Two identical vectors reported as not equal.')
    
    def test_isEqual2(self):
        v1 = np.array([1.0, -1.0, 3.0])
        v2 = np.array([1.0,  1.0, 3.0])
        self.assertEqual(False, coords.isEqual(v1, v2, 1e-12),'Two non-identical vectors reported as equal.')

    def test_isZero1(self):
        v1 = np. array([0.0, 0.0, 0.0])
        self.assertEqual(True, coords.isZero(v1, 1e-12),'Zero vector reported as non-zero.')

    def test_isZero2(self):
        v1 = np. array([1e-11, 0.0, 0.0])
        self.assertEqual(False, coords.isZero(v1, 1e-12),'Non-Zero vector reported as zero.')

    def test_axisFromHelix(self):
        # construct a perfect helix with a known axis along z.
        t = np.linspace(0.0, 2.0 * np.pi, 100 )
        r = 3.0
        axis = np.array([0.0, 0.0, 1.0])
        vecList = [ np.array([x, y, z]) for x, y, z in zip(r * np.cos(t), r * np.sin(t), t) ]
        self.assertAlmostEqualVectors(axis, 
                                      coords.axisFromHelix( vecList ),
                                      1e-12)

    def test_transformFromBlockFrameToLabFrame(self):
        # Start with helix radius 3 around z axis, with first point on x-axis. (3,0,0) 
        # Rotate helix to a new director at 45 degrees in x-z plane.
        # rotate helix 90 degrees about new director, then translate to new point.  
        # first point ends on y-axis + translate vector.
        
        t = np.linspace(0.0, 2.0 * np.pi, 100 )
        r = 3.0
        blockDirector = np.array([0.0, 0.0, 1.0])
        blockRefPoint = np.array([0.0, 0.0, 0.0])
        vecList = [ np.array([x, y, z]) for x, y, z in zip(r * np.cos(t), r * np.sin(t), t) ]
        
        labDirector = np.array([1.0, 0.0, -1.0])/np.sqrt(2) 
        labRefPoint = np.array([-3, 2, 1])
        labRotation = np.pi/2

        newVecs = coords.transformFromBlockFrameToLabFrame(labDirector, 
                                                           labRefPoint, 
                                                           labRotation, 
                                                           blockDirector, 
                                                           blockRefPoint, 
                                                           vecList)

        # make sure axis is the right direction
        self.assertAlmostEqualVectors(labDirector, coords.axisFromHelix(newVecs), 1e-12)
        
        # First point should end on y axis, and then be translated relative to new Ref Point
        self.assertAlmostEqualVectors(np.array([0.0, 3.0, 0.0]) + labRefPoint, newVecs[0], 1e-12)
        
        
    def test_transformFromLabFrameToBlockFrame(self):
        # Use Block to Lab transform and then Lab to Block.
        # should recover original structure.  
        
        t = np.linspace(0.0, 2.0 * np.pi, 100 )
        r = 3.0
        blockDirector = np.array([0.0, 0.0, 1.0])
        blockRefPoint = np.array([0.0, 0.0, 0.0])
        vecList = [ np.array([x, y, z]) for x, y, z in zip(r * np.cos(t), r * np.sin(t), t) ]
        
        labDirector = np.array([1.0, 0.0, -1.0])/np.sqrt(2) 
        labRefPoint = np.array([-3, 2, 1])
        labRotation = np.pi/2

        newVecs = coords.transformFromBlockFrameToLabFrame(labDirector, 
                                                           labRefPoint, 
                                                           labRotation, 
                                                           blockDirector, 
                                                           blockRefPoint, 
                                                           vecList)

        newNewVecs = coords.transformFromLabFrameToBlockFrame(labDirector, 
                                                           labRefPoint, 
                                                           labRotation, 
                                                           blockDirector, 
                                                           blockRefPoint, 
                                                           newVecs)

        # Block to Lab is tested separately. Assume it works. Should end up with original list
        # if it doesn't work then Lab to Block may be broken.         
        self.assertAlmostEqualListOfVectors(newNewVecs, vecList, 1e-12)
        
    def test_convertTriplesToEllipsoids(self):
        # construct triple of xyzpoints lying along x axis
        blocknames = ['O', 'P', 'Ca']
        xyzVals = [ np.array([-1.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]  
        aAxis = 3.0;
        cAxis = 0.5;
        eAtomNames, ePositions, eSizes, eRs, eRotVecs = coords.convertTriplesToEllipsoids(blocknames, xyzVals, aAxis, cAxis)
        self.assertEqual('P', eAtomNames[0])
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, 0.0]), ePositions[0])
        self.assertAlmostEqual(6.0, eSizes[0][0], 1e-12)
        self.assertAlmostEqual(2.0, eSizes[0][1], 1e-12)
        self.assertAlmostEqual(1.0, eSizes[0][2], 1e-12)
        self.assertAlmostEqualVectors(np.array([0.0, 1.0, 0.0]), eRs[0][0], 1e-12)
        self.assertAlmostEqualVectors(np.array([-1.0, 0.0, 0.0]), eRs[0][1], 1e-12)
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, 1.0]), eRs[0][2], 1e-12)
        self.assertAlmostEqualVectors(np.array([-1.0, 0.0, 0.0]), eRotVecs[0], 1e-12)
        
    def test_XYZPointsToOrientationMat(self):
        pos1 = np.array([-1.0, 0.0, 0.0])
        pos2 = np.array([ 1.0, 0.0, 0.0])
        rotMat = coords.XYZPointsToOrientationMat(pos1, pos2)
        self.assertAlmostEqualVectors(np.array([0.0, 1.0, 0.0]), rotMat[0], 1e-12)
        self.assertAlmostEqualVectors(np.array([-1.0, 0.0, 0.0]), rotMat[1], 1e-12)
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, 1.0]), rotMat[2], 1e-12)
        

    def test_constructTNBFrame1(self):
        p1 = np.array([-1.0, 0.0, -1.0])
        p2 = np.array([ 0.0, 0.0, 0.0])
        p3 = np.array([ 0.0, 0.0, 1.0])
        
        TNB = coords.constructTNBFrame(p1, p2, p3)
        
        self.assertAlmostEqualVectors(np.array([ 0.0,  0.0,  1.0]), TNB[0], 1e-12)
        self.assertAlmostEqualVectors(np.array([ 0.0, -1.0,  0.0]), TNB[1], 1e-12)
        self.assertAlmostEqualVectors(np.array([-1.0,  0.0,  0.0]), TNB[2], 1e-12)
        
    def test_constructTNBFrame2(self):
        p1 = np.array([ 0.0, 0.0, -1.0])
        p2 = np.array([ 0.0, 0.0, 0.0])
        p3 = np.array([ 0.0, 0.0, 1.0])
        
        TNB = coords.constructTNBFrame(p1, p2, p3)
        
        self.assertEqual(None, TNB)
        
    def test_bondAngle(self):
        p1 = np.array([-1.0, 0.0, -1.0])
        p2 = np.array([ 0.0, 0.0, 0.0])
        p3 = np.array([ 0.0, 0.0, 1.0])
        self.assertAlmostEqual(3.0 * np.pi/4.0, coords.bondAngle(p1, p2, p3), delta=1e-12 )

    def test_Dihedral(self):
        p1 = np.array([-1.0, 0.0, 0.0])
        p2 = np.array([ 0.0, 0.0, 0.0])
        p3 = np.array([ 0.0, 0.0, 1.0])
        p4 = np.array([ 0.0, 1.0, 1.0])
        self.assertAlmostEqual(-np.pi/2, coords.Dihedral(p1, p2, p3, p4), delta=1e-12)

    def test_measureAnglesAtConnection(self):
        s0 = np.array([-1.0, 0.0, 0.0])
        s1 = np.array([ 0.0, 0.0, 0.0])
        s2 = np.array([ 0.0, 0.0, 1.0])
        m2 = np.array([ 0.0, 1.0, 1.0])
        m1 = np.array([ -1.0, 1.0, 2.0])
        m0 = np.array([  0.0, 1.0, 3.0])

        d, alpha1, beta1, alpha2, beta2, alpha3 = coords.measureAnglesAtConnection(s0, s1, s2, m2, m1, m0)
        
        self.assertAlmostEqual(d, 1.0, delta=1e-12)
        self.assertAlmostEqual(alpha1, -np.pi/2.0, delta=1e-12)
        self.assertAlmostEqual(beta1, np.pi/2.0, delta=1e-12)
        self.assertAlmostEqual(alpha2, 3 * np.pi/4.0, delta=1e-12)
        self.assertAlmostEqual(beta2, np.pi/2, delta=1e-12)
        self.assertAlmostEqual(alpha3, np.pi/2, delta=1e-12)
        
    def test_computeDistBetweenM1AndS3(self):
        angle = np.pi/4
        # Used for the hairPinResidueToAtoms function.
        # Compute the difference between bondLength and the distance between M1 and S3 for a given 
        # value of dihedral around the S2 to M0 axis
        s2 = np.array([ 0.0, 0.0, 0.0])
        s3 = np.array([ 0.0, 0.0, 1.0])
        m3 = np.array([ 0.0, 1.0, 1.0])
        tnb = [s2,s3,m3]
        beta = np.pi/4
        bondLength = np.linalg.norm(s3 - s2)
        M1 = coords.generateTNBVecXYZ(tnb, beta, angle)
        self.assertAlmostEqual(np.abs(np.linalg.norm(s3 - M1) - bondLength), coords.computeDistBetweenM1AndS3(angle, tnb, beta, s3, bondLength), 1e-12) 

    def test_computeDistBetweenPoints(self):
        alphaA = np.pi/4
        alphaB = np.pi/3
        betaA = 165 * np.pi/180
        betaB = 150 * np.pi/180
        
        p1 = np.array([0.0, 0.0, 0.0])
        p2 = np.array([1.0, 0.0, 0.0])
        TNBA = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])]
        TNBB = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])]
        bondLength = 1.0
        targetDistance = 0.5
        
        # pointA = np.array([0.96592583, 0.1830127, 0.1830127 ]) # p1 + bondLength * coords.generateTNBVecXYZ(TNBA, betaA, alphaA)
        # pointB = np.array([1.8660254, 0.4330127, 0.25]) # p2 + bondLength * coords.generateTNBVecXYZ(TNBB, betaB, alphaB)
        
        DistBetweenPoints = 0.4365716990786799 # np.linalg.norm(pointA - pointB) - targetDistance
    
        self.assertAlmostEqual(DistBetweenPoints, 
                               coords.computeDistBetweenPoints([alphaA, alphaB],
                                                                p1, p2, TNBA, TNBB, betaA, betaB, 
                                                                bondLength, targetDistance),
                               delta=1e-12)


    def test_cpoa(self):

        p1List = [np.array([-1.0, 0.0, 0.0])]
        p2List = [np.array([ 1.0, 0.0, 0.0])]
        p3List = [np.array([ 0.0, -2.0, 1.0])]
        p4List = [np.array([ 0.0,  1.0, 1.0])]
        
        #muaList, paList, mubList, pbList
        [mua, pa, mub, pb] = coords.cpoa(p1List, p2List, p3List, p4List)
        
        self.assertAlmostEqual(1.0, mua[0], delta=1e-12)
        self.assertAlmostEqual(2.0, mub[0], delta=1e-12)
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, 0.0]), pa[0], delta=1e-12)
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, 1.0]), pb[0], delta=1e-12)

    def test_CylArea(self):
        Area = 2 * np.pi * 5.0 * 10.0
        self.assertAlmostEqual(Area, coords.CylArea(0, 2 * np.pi, 5.0, 10), delta=1e-12)

    def test_SphereArea(self):
        Area = 4.0 * np.pi * 5.0 * 5.0
        self.assertAlmostEqual(Area, coords.SphereArea(-np.pi, np.pi, -np.pi/2, np.pi/2, 5.0), delta=1e-12)

    def test_FrustumRadiusAtGivenZ(self):
        self.assertAlmostEqual(5.0, coords.FrustumRadiusAtGivenZ(5.0, 0, 10, 0, 10), delta=1e-12)

    def test_FrustumZAtZeroRadius(self):
        self.assertAlmostEqual(0.0, coords.FrustumZAtZeroRadius(5.0, 10, 5, 10), delta=1e-12)

    def test_checkPointInFrustum1(self):
        self.assertEqual(True, coords.checkPointInFrustum(np.array([0.0, 0.0, 7.5]), 5.0, 5.0, 10.0, 10.0, 0))
        
    def test_checkPointInFrustum2(self):
        self.assertEqual(False, coords.checkPointInFrustum(np.array([8.0, 0.0, 7.5]), 5.0, 5.0, 10.0, 10.0, 0))
    
    def test_rotatePointInChain(self):
        points = [ np.array([0.0, 0.0, 0.0]), np.array([0.0, 1.0, 1.0]), np.array([0.0, 0.0, 2.0])]
        self.assertAlmostEqualVectors(np.array([1.0, 0.0, 1.0]), coords.rotatePointInChain(points, -np.pi/2), delta=1e-12)
        
    def test_rotateDihedral(self):
        points = [ np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 2.0])]
        self.assertAlmostEqualVectors(np.array([1.0, 0.0, 2.0]), coords.rotateDihedral(points, np.pi/2), delta=1e-12)
    
    def test_pickRandomPointInFrustumAtGivenZHeightXYZ(self):
        point = coords.pickRandomPointInFrustumAtGivenZHeightXYZ(7.5, 5, 10, 5, 10)
        self.assertEqual(True, coords.checkPointInFrustum(point, 5, 5, 10, 10, 0))
        self.assertAlmostEqual(7.5, point[2], delta=1e-12)
        
    def test_pickRandomPointInFrustumXYZ(self):
        point = coords.pickRandomPointInFrustumXYZ(5, 10, 5, 10)
        self.assertEqual(True, coords.checkPointInFrustum(point, 5, 5, 10, 10, 0))

    def test_generateTNBVecXYZ(self):
        TNB = [ np.array([0.0, 0.0, 1.0]), np.array([1.0, 1.0, 0.0])/np.sqrt(2), np.array([-1.0, 1.0, 0.0]/np.sqrt(2))]
        self.assertAlmostEqualVectors(np.array([0.0, 1.0, 1.0])/np.sqrt(2), coords.generateTNBVecXYZ(TNB, 3*np.pi/4, np.pi/4), delta=1e-12)

    def test_pickRandomTNBDirectionInAngRangeXYZ(self):
        for _ in range(0, self.maxRndTests):
            TNB = [ np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]) ]
            point = coords.pickRandomTNBDirectionInAngRangeXYZ(TNB, np.pi/2, 3.0*np.pi/4.0, -np.pi/3, np.pi/3)
            pointRTF = coords.XYZ2SphericalPolar(point)
            self.assertLessEqual(pointRTF[1], 3.0*np.pi/4.0 - np.pi/2, "beta out of range - too great")
            self.assertGreaterEqual(pointRTF[1], np.pi/2.0 - np.pi/2, "beta out of range - too low")
            self.assertLessEqual(pointRTF[2],  np.pi/3.0 + np.pi/2, "alpha out of range - too great")
            self.assertGreaterEqual(pointRTF[2],  -np.pi/3.0  + np.pi/2, "alpha out of range - too small")

    def test_pickRandomPointOnUnitSphereInAngRange(self):
        theta1 = -np.pi/3
        theta2 = np.pi/3
        phi1 = np.pi/4
        phi2 = np.pi/4
        for _ in range(0, self.maxRndTests):
            theta, phi = coords.pickRandomPointOnUnitSphereInAngRange(theta1, theta2, phi1, phi2)
            self.assertLessEqual(theta, theta2, "theta out of range - too great")
            self.assertGreaterEqual(theta, theta1, "theta out of range - too low")
            self.assertLessEqual(phi, phi2, "phi out of range - too great")
            self.assertGreaterEqual(phi, phi1, "phi out of range - too low")

    def test_pickRandomPointOnUnitSphere(self):
        theta1 = -np.pi/2
        theta2 = np.pi/2
        phi1 = -np.pi
        phi2 = np.pi
        for _ in range(0, self.maxRndTests):
            theta, phi = coords.pickRandomPointOnUnitSphere()
            self.assertLessEqual(theta, theta2, "theta out of range - too great")
            self.assertGreaterEqual(theta, theta1, "theta out of range - too low")
            self.assertLessEqual(phi, phi2, "phi out of range - too great")
            self.assertGreaterEqual(phi, phi1, "phi out of range - too low")

    def test_pickRandomPointUniformlyOnEllipsoid(self):
        # checks that the point picked is indeed on the ellipsoid specified. (satisfies the ellipsoid eqn)
        a=3
        b=4
        c=5
        for _ in range(0, self.maxRndTests):
            point = coords.pickRandomPointUniformlyOnEllipsoid(a, b, c)
            pointXYZ = coords.sphericalPolar2XYZ(point)
            d = np.power(pointXYZ[0],2) / (a * a) + np.power(pointXYZ[1],2) / (b * b) + np.power(pointXYZ[2],2) / (c * c) 
            self.assertAlmostEqual(1, d, delta=1e-12)


    def test_pickRandomPointUniformlyOnEllipsoidInAngRange(self):
        # checks the picked point is within angular range and satisfies ellipsoid equation
        a=3
        b=4
        c=5
        theta1 = -np.pi/3
        theta2 = np.pi/3
        phi1 = -np.pi/4
        phi2 = np.pi/4
        for _ in range(0, self.maxRndTests):
            point = coords.pickRandomPointUniformlyOnEllipsoidInAngRange(a, b, c, theta1, theta2, phi1, phi2)
            self.assertLessEqual(point[1], theta2, "theta out of range - too great")
            self.assertGreaterEqual(point[1], theta1, "theta out of range - too low")
            self.assertLessEqual(point[2], phi2, "phi out of range - too great")
            self.assertGreaterEqual(point[2], phi1, "phi out of range - too low")
            pointXYZ = coords.sphericalPolar2XYZ(point)
            d = np.power(pointXYZ[0],2) / (a * a) + np.power(pointXYZ[1],2) / (b * b) + np.power(pointXYZ[2],2) / (c * c) 
            self.assertAlmostEqual(1, d, delta=1e-12)

    def test_pickRandomPointInUnitCircleInAngRange(self):
        phi1 = -np.pi/4
        phi2 =  np.pi/4
        for _ in range(0, self.maxRndTests):
            phi, r = coords.pickRandomPointInUnitCircleInAngRange(phi1, phi2)
            self.assertLessEqual(phi, phi2, "phi out of range - too great")
            self.assertGreaterEqual(phi, phi1, "phi out of range - too low")
            self.assertLessEqual(r, 1, "r out of range - too great")
            self.assertGreaterEqual(r, 0.0, "r out of range - too low")

    def test_pickRandomPointOnUnitCircleInAngRange(self):
        phi1 = -np.pi/4
        phi2 =  np.pi/4
        for _ in range(0, self.maxRndTests):
            phi = coords.pickRandomPointOnUnitCircleInAngRange(phi1, phi2)
            self.assertLessEqual(phi, phi2, "phi out of range - too great")
            self.assertGreaterEqual(phi, phi1, "phi out of range - too low")

    def test_pickRandomPointInAnnulus(self):
        phi1 = -np.pi/4
        phi2 =  np.pi/4
        r1 = 5.0
        r2 = 10.0
        for _ in range(0, self.maxRndTests):
            r, phi = coords.pickRandomPointInAnnulus(r1, r2, phi1, phi2)
            self.assertLessEqual(phi, phi2, "phi out of range - too great")
            self.assertGreaterEqual(phi, phi1, "phi out of range - too low")
            self.assertLessEqual(r, r2, "r out of range - too great")
            self.assertGreaterEqual(r, r1, "r out of range - too low")


    def test_pickRandomPointInUnitCircle(self):
        for _ in range(0, self.maxRndTests):
            phi, r = coords.pickRandomPointInUnitCircle()
            self.assertLessEqual(phi, np.pi, "phi out of range - too great")
            self.assertGreaterEqual(phi, -np.pi, "phi out of range - too low")
            self.assertLessEqual(r, 1.0, "r out of range - too great")
            self.assertGreaterEqual(r, 0.0, "r out of range - too low")

    def test_pickRandomPointInEllipsoid(self):
        a = 3
        b = 4
        c = 5
        for _ in range(0, self.maxRndTests):
            point = coords.pickRandomPointInEllipsoid(a, b, c)
            self.assertLessEqual(point[0], a, "x out of range - too great")
            self.assertGreaterEqual(point[0], -a, "x out of range - too low")        
            self.assertLessEqual(point[1], b, "y out of range - too great")
            self.assertGreaterEqual(point[1], -b, "y out of range - too low")        
            self.assertLessEqual(point[0], a, "z out of range - too great")
            self.assertGreaterEqual(point[0], -a, "z out of range - too low")        
            d = np.power(point[0],2) / (a * a) + np.power(point[1],2) / (b * b) + np.power(point[2],2) / (c * c) 
            self.assertLessEqual(d, 1.0, "r out side of specified ellipse - too great")


    def test_pickRandomPointInSphericalRange(self):
        theta1 = -np.pi/2
        theta2 =  np.pi/2
        phi1 = -np.pi/4
        phi2 =  np.pi/4
        r1 = 5.0
        r2 = 10.0
        for _ in range(0, self.maxRndTests):
            point = coords.pickRandomPointInSphericalRange(r1, r2, theta1, theta2, phi1, phi2)
            pointPol = coords.XYZ2SphericalPolar(point)
            self.assertLessEqual(pointPol[1], theta2, "theta out of range - too great")
            self.assertGreaterEqual(pointPol[1], theta1, "theta out of range - too low")        
            self.assertLessEqual(pointPol[2], phi2, "phi out of range - too great")
            self.assertGreaterEqual(pointPol[2], phi1, "phi out of range - too low")        
            self.assertLessEqual(pointPol[0], r2, "r out of range - too great")
            self.assertGreaterEqual(pointPol[0], r1, "r out of range - too low")        

    def test_pickRandomPointInEllipsoidRange(self):
        a = 3
        b = 4
        c = 5
        theta1 = -np.pi/2
        theta2 =  np.pi/2
        phi1 = -np.pi/4
        phi2 =  np.pi/4

        for _ in range(0, self.maxRndTests):
            point = coords.pickRandomPointInEllipsoidRange(a, b, c, theta1, theta2, phi1, phi2)
            pointPol = coords.XYZ2SphericalPolar(point)
            self.assertLessEqual(pointPol[1], theta2, "theta out of range - too great")
            self.assertGreaterEqual(pointPol[1], theta1, "theta out of range - too low")        
            self.assertLessEqual(pointPol[2], phi2, "phi out of range - too great")
            self.assertGreaterEqual(pointPol[2], phi1, "phi out of range - too low")        
            self.assertLessEqual(point[0], a, "x out of range - too great")
            self.assertGreaterEqual(point[0], -a, "x out of range - too low")        
            self.assertLessEqual(point[1], b, "y out of range - too great")
            self.assertGreaterEqual(point[1], -b, "y out of range - too low")        
            self.assertLessEqual(point[0], a, "z out of range - too great")
            self.assertGreaterEqual(point[0], -a, "z out of range - too low")        
            d = np.power(point[0],2) / (a * a) + np.power(point[1],2) / (b * b) + np.power(point[2],2) / (c * c) 
            self.assertLessEqual(d, 1.0, "r outside of specified ellipse - too great")

    def test_pickRandomPointOnUnitCircle(self):
        for _ in range(0, self.maxRndTests):
            phi = coords.pickRandomPointOnUnitCircle()
            self.assertLessEqual(phi, np.pi, "phi out of range - too great")
            self.assertGreaterEqual(phi, -np.pi, "phi out of range - too low")        
        
    def test_pickRandomPointOnUnitCylinderInAngRange(self):
        phi1 = -np.pi/4
        phi2 =  np.pi/4
        for _ in range(0, self.maxRndTests):
            phi, z = coords.pickRandomPointOnUnitCylinderInAngRange(phi1, phi2)
            self.assertLessEqual(phi, phi2, "phi out of range - too great")
            self.assertGreaterEqual(phi, phi1, "phi out of range - too low")        
            self.assertLessEqual(z, 1.0, "z out of range - too great")
            self.assertGreaterEqual(z, 0.0, "z out of range - too low")        

    def test_pickRandomPointOnUnitCylinder(self):
        for _ in range(0, self.maxRndTests):
            phi, z = coords.pickRandomPointOnUnitCylinder()
            self.assertLessEqual(z, 1.0, "z out of range - too great")
            self.assertGreaterEqual(z, 0.0, "z out of range - too low")        
            self.assertLessEqual(phi, np.pi, "phi out of range - too great")
            self.assertGreaterEqual(phi, -np.pi, "phi out of range - too low")        

    def test_pickRandomPointInUnitCylinder(self):
        for _ in range(0, self.maxRndTests):
            pos = coords.pickRandomPointInUnitCylinder()
            self.assertLessEqual(pos[2], 1.0, "z out of range - too great")
            self.assertGreaterEqual(pos[2], 0.0, "z out of range - too low")        
            self.assertLessEqual(pos[1], np.pi, "phi out of range - too great")
            self.assertGreaterEqual(pos[1], -np.pi, "phi out of range - too low")
            self.assertLessEqual(pos[0], 1.0, "r out of range - too great")
            self.assertGreaterEqual(pos[0], 0.0, "r out of range - too low")
        

    def test_pickRandomPointOnUnitSquare(self):
        for _ in range(0, self.maxRndTests):
            point = coords.pickRandomPointOnUnitSquare()
            self.assertGreaterEqual(point[0], 0.0, "x out of range - too low")        
            self.assertLessEqual(point[0], 1.0, "x out of range - too great")
            self.assertGreaterEqual(point[1], 0.0, "y out of range - too low")        
            self.assertLessEqual(point[1], 1.0, "y out of range - too great")


class Test_Cartesian(additionalAsserts):

    def test_getCOM(self):  
        testVecs = [ np.array([-1.0,  0.0,  0.0]),
                     np.array([ 0.0, -1.0,  0.0]),
                     np.array([ 0.0,  0.0, -1.0]),
                     np.array([ 0.0,  0.0,  0.0]),
                     np.array([ 1.0,  0.0,  0.0]),
                     np.array([ 0.0,  1.0,  0.0]),
                     np.array([ 0.0,  0.0,  1.0])]

        self.assertAlmostEqualVectors(cart.getCentreOfMass(testVecs), np.array([0.0, 0.0, 0.0]), delta=1e-12)

 
    def test_alignTwoBlocksOfVectors(self):  
        Block1Vecs = [ np.array([ 0.0,  0.0,  0.0]),
                       np.array([ 1.0,  0.0,  0.0]),
                       np.array([ 1.0,  1.0,  0.0])]
        
        Block2Vecs = [ np.array([ 0.0,  0.0,  0.0]),
                       np.array([ 1.0,  0.0,  0.0]),
                       np.array([ 1.0,  1.0,  0.0])]

        alignedVecs = [ np.array([ 1.0,  0.0,  0.0]),
                        np.array([ 1.0,  1.0,  0.0]),
                        np.array([ 0.0,  1.0,  0.0])]                        
                        
        self.assertAlmostEqualListOfVectors(alignedVecs, cart.AlignTwoBlocksOfVectors(Block1Vecs, Block2Vecs), delta=1e-12)


    def test_translateBlock(self):
        basePosition = np.array([ 2.0, -1.0, -1.0])
        listOfVecs = [ np.array([ 3.0,  2.0,  1.0]),
                       np.array([ 1.0,  0.0,  0.0]),
                       np.array([ 1.0,  1.0,  0.0]),
                       np.array([ 1.0,  1.0,  1.0])]

        translateVecs = [ np.array([ 2.0, -1.0, -1.0]),
                          np.array([ 0.0, -3.0, -2.0]),
                          np.array([ 0.0, -2.0, -2.0]),
                          np.array([ 0.0, -2.0, -1.0])]
        
        self.assertAlmostEqualListOfVectors(translateVecs, cart.translateBlock(basePosition, listOfVecs), delta=1e-12)

    
    def test_translateBlockRelRef(self):
        basePosition = np.array([ 2.0, -1.0, -1.0])
        refVec = np.array([ 1.0,  1.0,  1.0])
        listOfVecs = [ np.array([ 3.0,  2.0,  1.0]),
                       np.array([ 1.0,  0.0,  0.0]),
                       np.array([ 1.0,  1.0,  0.0]),
                       np.array([ 1.0,  1.0,  1.0])]

        translateVecs = [ np.array([ 4.0, -0.0, -1.0]),
                          np.array([ 2.0, -2.0, -2.0]),
                          np.array([ 2.0, -1.0, -2.0]),
                          np.array([ 2.0, -1.0, -1.0])]
        
        self.assertAlmostEqualListOfVectors(translateVecs, cart.translateBlockRelRef(basePosition, listOfVecs, refVec), delta=1e-12)
    
    
    def test_alignBlockWithDirector(self):
        director = np.array([ 1.0,  0.0,  0.0])
        testVec =  np.array([ 0.0,  0.0,  1.0])
        listOfVecs = [ np.array([ 1.0,  0.0,  0.0]),
                       np.array([ 0.0,  1.0,  0.0]),
                       np.array([ 0.0,  0.0,  1.0])]

        alignedVecs = [ np.array([ 0.0,  0.0, -1.0]),
                          np.array([ 0.0,  1.0,  0.0]),
                          np.array([ 1.0,  0.0,  0.0])]
        
        self.assertAlmostEqualListOfVectors(alignedVecs, cart.alignBlockWithDirector(director, listOfVecs, testVec), delta=1e-12)
    
    
    def test_getPrincipalAxis(self):
        listOfVecs = [ np.array([ 2.0,  0.0,  0.0]),
                       np.array([ 0.0,  1.0,  0.0]),
                       np.array([ 0.0,  0.0,  1.0])]
        
        self.assertAlmostEqualVectors(np.array([1.0, 0.0, 0.0]), cart.getPrincipalAxis(listOfVecs), delta=1e-12)
    
    def test_rotPAboutAxisAtPoint(self):
        p = np.array([1.0, 2.0, 1.0])
        p0 = np.array([0.0, 2.0, 0.0])
        n  = np.array([0.0, 0.0, 1.0])
        angle = np.pi/2.0
        
        retVal = np.array([0.0, 3.0, 1.0])
        self.assertAlmostEqualVectors(retVal, cart.rotPAboutAxisAtPoint(p, p0, n, angle), delta=1e-12)
    
    def test_rotPAboutAxis(self):
        p = np.array([1.0, 0.0, 1.0])
        n  = np.array([0.0, 0.0, 1.0])
        angle = np.pi/2.0
        
        retVal = np.array([0.0, 1.0, 1.0])
        self.assertAlmostEqualVectors(retVal, cart.rotPAboutAxis(p, n, angle), delta=1e-12)
    
    
    def test_vectorBetweenTwoPoints(self):
        p1 = np.array([2.0, 1.0, 0.0])
        p2 = np.array([-2.0, 1.0, 0.0])
        
        self.assertAlmostEqualVectors(np.array([ 1.0, 0.0, 0.0]), cart.vectorBetweenTwoPoints(p1, p2, 'TwoToOne'), delta=1e-12)
        self.assertAlmostEqualVectors(np.array([-1.0, 0.0, 0.0]), cart.vectorBetweenTwoPoints(p1, p2, 'OneToTwo'), delta=1e-12)
        

    def test_clamp(self):
        self.assertAlmostEqual(1.5, cart.clamp(2.0, -1.5, 1.5), delta=1e-12)
        self.assertAlmostEqual(-1.5, cart.clamp(-2.0, -1.5, 1.5), delta=1e-12)

    def test_closestApproachTwoLineSegmentsSquared(self):
        p1 = np.array([-1.0, 0.0, -1.0])
        q1 = np.array([ 1.0, 0.0, -1.0])
        p2 = np.array([ 0.0, -1.0, 1.0])
        q2 = np.array([ 0.0,  1.0, 1.0])
        
        self.assertAlmostEqual(4.0, cart.closestApproachTwoLineSegmentsSquared( p1, q1, p2, q2, returnVec=False), delta=1e-12)
        
        retVec = cart.closestApproachTwoLineSegmentsSquared( p1, q1, p2, q2, returnVec=True)
        self.assertAlmostEqual(4.0, retVec[0])
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, -2.0]), retVec[1])
        
    def test_closestApproachPointToLineSegmentSquared(self):
        p = np.array([ -0.0, 0.0, -1.0])
        p2 = np.array([ 0.0, -1.0, 1.0])
        q2 = np.array([ 0.0,  1.0, 1.0])
        
        self.assertAlmostEqual(4.0, cart.closestApproachPointToLineSegmentSquared(p, p2, q2, returnVec=False), delta=1e-12)
        
        retVec = cart.closestApproachPointToLineSegmentSquared( p, p2, q2, returnVec=True)

        self.assertAlmostEqual(4.0, retVec[0])
        self.assertAlmostEqualVectors(np.array([0.0, 0.0, -2.0]), retVec[1])
    
    def test_distanceOfLineToAPlane(self):
        linePoint = np.array([-1.0, 0.0, -1.0])
        lineDirector = np.array([1.0, 0.0, 0.0])
        planePoint = np.array([1.0, 0.0, 1.0])
        planeDirector= np.array([1.0, 0.0, 0.0])
        
        self.assertAlmostEqual(2.0, cart.distanceOfLineToAPlane(linePoint, lineDirector, planePoint, planeDirector), delta=1e-12)
    
      
    
if __name__ == '__main__':
    unittest.main()
