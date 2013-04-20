#!/usr/bin/env python

# This python script runs a number of single iteration regression tests and diffs output with the expected output...
# Script should be run from the folder in which the cactus executable exists.
# Note: if no differences exist, the regression test output file is deleted for convenience.

import os
import shutil
import sys
import string
import filecmp
import difflib

# Note: To debug, call the python debugger module, pdb, as a script with this script as the argument (python -m pdb script.py)
# or insert: import pdb; pdb.set_trace() above, and call from terminal as usual

RegDir=os.path.dirname(os.path.realpath(__file__))
print 'Running runreg.py in ' + RegDir
print ''


# Run regression test 1
print 'Running regression test 1'
IFN=RegDir + '/RegTest1.in'
CCommand='cactus ' + IFN 
os.system(CCommand)

# clean up standard output files which are meaningless for this calculation
os.remove('RegTest1_Param.csv')
os.remove('RegTest1_RevData.csv')
os.remove('RegTest1_TimeData.csv')

# diff output
FN1='RegTest1_RegData.out'
FN2=RegDir + '/RegTest1_RegData_Ex.out'
if not filecmp.cmp(FN1,FN2):
    print 'Summary of differences between ' + FN1 + ' and ' + FN2 + ':'
    
    f=open(FN1,'r');
    F1=f.readlines()
    f.close()
    
    f=open(FN2,'r');
    F2=f.readlines()
    f.close()
    
    DOut=difflib.unified_diff(F1,F2,fromfile=FN1,tofile=FN2)
    for line in DOut:
        sys.stdout.write(line)
    #endfor
else:
    print 'No differences'
    os.remove(FN1)
#endif
print ''


# Run regression test 2
print 'Running regression test 2'
IFN=RegDir + '/RegTest2.in'
CCommand='cactus ' + IFN 
os.system(CCommand)

# clean up standard output files which are meaningless for this calculation
os.remove('RegTest2_Param.csv')
os.remove('RegTest2_RevData.csv')
os.remove('RegTest2_TimeData.csv')

# diff output
FN1='RegTest2_RegData.out'
FN2=RegDir + '/RegTest2_RegData_Ex.out'
if not filecmp.cmp(FN1,FN2):
    print 'Summary of differences between ' + FN1 + ' and ' + FN2 + ':'
    
    f=open(FN1,'r');
    F1=f.readlines()
    f.close()
    
    f=open(FN2,'r');
    F2=f.readlines()
    f.close()
    
    DOut=difflib.unified_diff(F1,F2,fromfile=FN1,tofile=FN2)
    for line in DOut:
        sys.stdout.write(line)
    #endfor
else:
    print 'No differences'
    os.remove(FN1)
#endif
print ''


# Run regression test 3
print 'Running regression test 3'
IFN=RegDir + '/RegTest3.in'
CCommand='cactus ' + IFN 
os.system(CCommand)

# clean up standard output files which are meaningless for this calculation
os.remove('RegTest3_Param.csv')
os.remove('RegTest3_RevData.csv')
os.remove('RegTest3_TimeData.csv')

# diff output
FN1='RegTest3_RegData.out'
FN2=RegDir + '/RegTest3_RegData_Ex.out'
if not filecmp.cmp(FN1,FN2):
    print 'Summary of differences between ' + FN1 + ' and ' + FN2 + ':'
    
    f=open(FN1,'r');
    F1=f.readlines()
    f.close()
    
    f=open(FN2,'r');
    F2=f.readlines()
    f.close()
    
    DOut=difflib.unified_diff(F1,F2,fromfile=FN1,tofile=FN2)
    for line in DOut:
        sys.stdout.write(line)
    #endfor
else:
    print 'No differences'
    os.remove(FN1)
#endif
print ''
