#!/usr/bin/env python

# DAKOTA analysis driver for CACTUS
# Reads DAKOTA input, writes input for CACTUS, executes, and compiles a DAKOTA return file.
# Analysis component inputs:
#   1. Path to CACTUS executable to use
#   2. CACTUS nominal input file to modify on each call
#   3. If the geometry file is to be modified, should be set to 'Y', otherwise 'N'
#   4. Path to nominal geometry generation script (if AC 3 set to 'Y')
#      Note: The name of the output geometry file should be assigned to a variable named 'FN' in the script.
#   Other optional analysis components in no particular order (KeepFiles, etc...). See those recognized below...

import os
import shutil
import sys
import string
import math

# Note: To debug, call the python debugger module, pdb, as a script with this script as the argument (python -m pdb SCRIPT)
# or insert: import pdb; pdb.set_trace() above, and call from terminal as usual

# Define DAKOTA data classes
class VarClass:
    def __init__(self,Tag,Val):
        self.Tag=Tag
        self.Val=Val
    #enddef
#endclass
    
class FuncClass:
    def __init__(self,Tag,Code,FnVal):
        self.Tag=Tag
        self.Code=Code
        self.FnVal=FnVal    
    #enddef
#endclass
    
class DVClass:
    def __init__(self,Tag,Code):
        self.Tag=Tag
        self.Code=Code 
    #enddef
#endclass   
    
class ACClass:
    def __init__(self,Tag,Code):
        self.Tag=Tag
        self.Code=Code    
    #enddef
#endclass
    
# Define functions
def TestNumString(TestStr):
    try:
        Num=float(TestStr)
    except ValueError:
        IsNumString=False
    else:
        IsNumString=True
    #endtry

    return IsNumString
#enddef    

def tail(fname, window):
    """Read last N lines from file fname."""
    try:
        f = open(fname, 'r')
    except IOError, err:
        if err.errno == errno.ENOENT:
            return []
        else:
            raise
        #endif
    else:
        BUFSIZ = 1024
        f.seek(0, 2)
        fsize = f.tell()
        block = -1
        data = ""
        exit = False
        while not exit:
            step = (block * BUFSIZ)
            if abs(step) >= fsize:
                f.seek(0)
                exit = True
            else:
                f.seek(step, 2)
            #endif
            data = f.read().strip()
            if data.count('\n') >= window:
                break
            else:
                block -= 1
            #endif
        #endwhile
            
        return data.splitlines()[-window:]
    #endtry
#enddef

def ReplaceScalar(FileLines,FileTag,VarVal,VarDes):
    # Replace scalar value in in template file...
    DoneLines=False
    LInd=0
    while not DoneLines:
        Line=FileLines[LInd]
        LL=Line.split('=')
        
        Cond=(LL[0].strip().lower() == FileTag.lower())  # Condition for identifying the correct line...
        if Cond:
            VS=LL[1]

            # remove comment if present (MATLAB)
            ComInd=VS.find('%')
            if ComInd>0:
                VS=VS[0:ComInd]
            #endif

            # strip ';' if present (MATLAB)
            HasSC=False
            if VS.find(';')>0:
                VS=VS.strip('; ')
                HasSC=True
            #endif

            # reset value
            if VarDes == 'F':
                LL[1]=str(float(VS)*VarVal) 
            elif VarDes == 'D':
                LL[1]=str(float(VS)+VarVal)   
            else:
                LL[1]=str(VarVal)        
            #endif

            if HasSC:
                LL[1]=LL[1] + ';'           
            #endif

            # rejoin LL and replace line in TempFileLines
            FileLines[LInd]='='.join(LL)

            DoneLines=True
        elif LInd == len(FileLines)-1:
            DoneLines=True
        #endif
        
        LInd+=1
    #endwhile    
#enddef

def ReplaceString(FileLines,FileTag,ValString):
    # Replace string value in in template file...
    DoneLines=False
    LInd=0
    while not DoneLines:
        Line=FileLines[LInd]
        LL=Line.split('=')
        
        Cond=(LL[0].strip().lower() == FileTag.lower())  # Condition for identifying the correct line...
        if Cond:
            VS=LL[1]
            
            # determine if ';' is present
            HasSC=False
            if VS.find(';')>0:
                HasSC=True
            #endif
            
            # reset value
            LL[1]="'" + ValString + "'"
            
            if HasSC:
                LL[1]=LL[1] + ';'           
            #endif
            
            # rejoin LL and replace line in TempFileLines
            FileLines[LInd]='='.join(LL)
            
            DoneLines=True
        elif LInd == len(FileLines)-1:
            DoneLines=True
        #endif
        
        LInd+=1
    #endwhile    
#enddef

def GetScalar(FileLines,FileTag):
    # Get default scalar value from template 
    DoneLines=False
    LInd=0
    Val=0.0
    while not DoneLines:
        Line=FileLines[LInd]
        LL=Line.split('=')
        
        Cond=(LL[0].strip().lower() == FileTag.lower())  # Condition for identifying the correct line...
        if Cond:
            VS=LL[1]
            
            # remove comment if present (MATLAB)
            ComInd=VS.find('%')
            if ComInd>0:
                VS=VS[0:ComInd]
            #endif
            
            # strip ';' if present (MATLAB)
            VS=VS.strip('; ')
            Val=float(VS)      
            DoneLines=True
        elif LInd == len(FileLines)-1:
            DoneLines=True
        #endif
        
        LInd+=1
    #endwhile 
    
    return Val   
#enddef

#########################

# Read input and output filenames from the command line
if len(sys.argv) is 3:
    DAKInFN=sys.argv[1]   
    DAKOutFN=sys.argv[2] 
else:
    print 'No files given...'
    sys.exit(1)
#endif

# Get DAKOTA filebase and run index (if exists)
DAKOutFB,DAKOutExt=DAKOutFN.split('.',1)

SplitList=DAKInFN.split('.')
DAKInd=''
DAKIndLast=None
for SL in SplitList:
    if TestNumString(SL):
        DAKInd=DAKInd + '_' + SL
        DAKIndLast=int(SL)
    #endif
#endfor

# Read DAKOTA input file into data classes
fDAKInput=open(DAKInFN,'r')

NumVars,Tag=string.split(fDAKInput.readline())
NumVars=int(NumVars)
Vars=[VarClass(None,None) for i in range(NumVars)]
for i in range(NumVars):
    Vars[i].Val,Vars[i].Tag=string.split(fDAKInput.readline())
#endfor

NumFuncs,Tag=string.split(fDAKInput.readline())
NumFuncs=int(NumFuncs)
Funcs=[FuncClass(None,None,None) for i in range(NumFuncs)]
for i in range(NumFuncs):
    Funcs[i].Code,Funcs[i].Tag=string.split(fDAKInput.readline())
#endfor

NumDVs,Tag=string.split(fDAKInput.readline())
NumDVs=int(NumDVs)
DVs=[DVClass(None,None) for i in range(NumDVs)]
for i in range(NumDVs):
    DVs[i].Code,DVs[i].Tag=string.split(fDAKInput.readline())
#endfor
    
NumACs,Tag=string.split(fDAKInput.readline())
NumACs=int(NumACs)
ACs=[ACClass(None,None) for i in range(NumACs)]
for i in range(NumACs):
    ACs[i].Code,ACs[i].Tag=string.split(fDAKInput.readline())
#endfor    
    
fDAKInput.close()  

#######################
# Check for additional analysis components
KeepFiles=False
RunCase=True
for i in range(NumACs):
    CodeStr=ACs[i].Code
    
    if CodeStr.find('KeepFiles') >= 0:
        # Check for keep files flag. If 'KeepFiles', keep files for all iterations.
        # If 'KeepFiles:1,2,5-8,...', files will be kept for the list of iteration numbers provided. Note the
        # iteration number is the last (lowest level) Dakota index attached to the analysis driver input file.
        
        # Check for specific iteration numbers
        if CodeStr.find(':') >= 0:
            CTag,IterList=CodeStr.split(':')
            
            # fill out - fields
            ILS=IterList.split(',')
            for j,ILSI in enumerate(ILS):
                Ind=ILSI.find('-')
                if Ind>0:
                    LBUB=ILSI.split('-')
                    ILS[j]=','.join(map(str,range(LBUB[0],LBUB[1]+1)))
                #endif
            #endfor
            IterList=','.join(ILS)

            IL=map(int,IterList.split(','))
            # Try to find current iteration in list
            try:
                IL.index(DAKIndLast)
            except ValueError:
                KeepFiles=False
            else:
                KeepFiles=True
            #endtry
        else:
            KeepFiles=True
        #endif
        
    elif CodeStr.find('SkipCase:') >= 0:    
        # Check for skip case flag. 
        # If 'SkipCase:1,2,5-8,...', the indicated cases will not be run. Note the
        # iteration number is the last (lowest level) Dakota index attached to the analysis driver input file.
        
        CTag,IterList=CodeStr.split(':')
        
        # fill out - fields
        ILS=IterList.split(',')
        for j,ILSI in enumerate(ILS):
            Ind=ILSI.find('-')
            if Ind>0:
                LBUB=ILSI.split('-')
                ILS[j]=','.join(map(str,range(int(LBUB[0]),int(LBUB[1])+1)))
            #endif
        #endfor
        IterList=','.join(ILS)

        IL=map(int,IterList.split(','))
        # Try to find current iteration in list
        try:
            IL.index(DAKIndLast)
        except ValueError:
            RunCase=True
        else:
            RunCase=False
        #endtry
        
    elif CodeStr.find('RunCase:') >= 0:    
        # Check for run case flag. 
        # If 'RunCase:1,2,5-8,...', the indicated cases will be run. RunCase overrides SkipCase if both present. Note the
        # iteration number is the last (lowest level) Dakota index attached to the analysis driver input file.
        
        CTag,IterList=CodeStr.split(':')
        
        # fill out - fields
        ILS=IterList.split(',')
        for j,ILSI in enumerate(ILS):
            Ind=ILSI.find('-')
            if Ind>0:
                LBUB=ILSI.split('-')
                ILS[j]=','.join(map(str,range(int(LBUB[0]),int(LBUB[1])+1)))
            #endif
        #endfor
        IterList=','.join(ILS)

        IL=map(int,IterList.split(','))
        # Try to find current iteration in list
        try:
            IL.index(DAKIndLast)
        except ValueError:
            RunCase=False
        else:
            RunCase=True
        #endtry            
        
    #endif 
    
#endfor

# End here if skipping this case
if not RunCase:
    fDAKOut = open(DAKOutFN, 'w')
    for Func in Funcs:
        OutLine = '-999 ' + Func.Tag + '\n'
        fDAKOut.write(OutLine)
    #endfor

    sys.exit(0)
#endif
######################

############## Add analysis specific pre-processing of top level inputs here #############
# Create dictionary from input variables...
VarDict={}
for Var in Vars:
    VarDict[Var.Tag]=float(Var.Val)
#endfor

##########################################################################################

# Get path to CACTUS executable
CACTUSExe=ACs[0].Code

# Get the CACTUS template file
CACTUSTempFN=ACs[1].Code
CACTUSFB,CACTUSExt=CACTUSTempFN.rsplit('.',1)
if CACTUSFB.find('/') >= 0:
    FF,CACTUSFB=CACTUSFB.rsplit('/',1) # keep temp files in local directory
#endif

CACTUSFN=CACTUSFB + '_Mod'
if DAKInd:
    CACTUSFN=CACTUSFN + DAKInd
#endif
CACTUSProbFN=CACTUSFN + '.' + CACTUSExt
CACTUSRevDataFN=CACTUSFN + '_RevData.csv'
CACTUSTimeDataFN=CACTUSFN + '_TimeData.csv'
CACTUSParamFN=CACTUSFN + '_Param.csv'

fTemp=open(CACTUSTempFN,'r')
TempFile=fTemp.read()  # read the file into a string
fTemp.close()

# Split template into lines
TempFileLines=TempFile.splitlines()


# Get geom file (if third component is 'Y')
CodeStr=ACs[2].Code
UseGeom=False
if CodeStr.strip().lower() == 'y':
    UseGeom=True
    GeomTempFN=ACs[3].Code
    GeomSFB,GeomSExt=GeomTempFN.rsplit('.',1)
    if GeomSFB.find('/') >= 0:
        FF,GeomSFB=GeomSFB.rsplit('/',1) # keep temp files in local directory
    #endif
    
    GeomSFN=GeomSFB + '_Mod'
    GeomFN='TurbGeom'
    if DAKInd:
        GeomSFN=GeomSFN + DAKInd
        GeomFN=GeomFN + DAKInd
    #endif
    GeomSFN=GeomSFN + '.' + GeomSExt
    GeomFN=GeomFN + '.geom'
    
    fGeom=open(GeomTempFN,'r')
    GeomFile=fGeom.read()  # read the file into a string
    fGeom.close()
    
    # Split template into lines
    GeomFileLines=GeomFile.splitlines()
    
    # Replace output filename
    ReplaceString(GeomFileLines,'FN',GeomFN) 
    
    # Turn off geom plotting
    ReplaceScalar(GeomFileLines,'PlotTurbine',0,None) 
    
    # Replace geom filename in CACTUS template
    ReplaceString(TempFileLines,'GeomFilePath',GeomFN) 

#endif



##### Apply modifications to template specific to each DAKOTA input tag. Input tags should generally have the form TAG:DES
##### where the TAG identifies the variable. The designator (DES) can be either "F" "D" or "V", indicating the value should
##### be applied as a factor on, delta over, or directly as the value of the identified variable...

for Var in Vars:
    VarTag=Var.Tag
    Val=float(Var.Val)
    
    if VarTag.find(':') < 0:
        Tag=VarTag
        Des=None
    else:
        Tag,Des=VarTag.split(':')
    #endif
    
    Tag=Tag.strip()
    # Recognize DAKOTA tags. If not a special tag defined below, Tag is assumed to be a scalar variable in the input file.
    if Tag == 'FnRRPM':
        # Set RPM to achieve a given FnR with a particular tip speed ratio
        FnR=Val
        
        # Check vars list for tip speed ratio (must be in there as 'Ut').
        # If not present, get default from template file
        if 'Ut' in VarDict:
            TSR=VarDict['Ut']
        else:
            TSR=GetScalar(TempFileLines,'Ut')    
        #endif
        
        # Check vars list for radius (must be in there as 'R').
        # If not present, get default from template file
        if 'R' in VarDict:
            R=VarDict['R']
        else:
            R=GetScalar(TempFileLines,'R')    
        #endif
        
        # Calc RPM
        g=32.174 # ft/s^2
        RPM=math.sqrt(FnR**2*TSR**2*g/R)*60.0/(2.0*math.pi)
        
        # Set new value
        ReplaceScalar(TempFileLines,'RPM',RPM,None)    
        
    else:
        # Assume Tag matches a variable name in either the CACTUS input file or the geometry script, and we are replacing a scalar...
        ReplaceScalar(TempFileLines,Tag,Val,Des) 
        if UseGeom:
            ReplaceScalar(GeomFileLines,Tag,Val,Des) 
        #endif
                      
    #endif
    
#endfor

############################################################################################

# Rejoin TempFileLines into file string (with carriage returns)
TempFileNew='\n'.join(TempFileLines) + '\n'

# Write the CACTUS input file
fCACTUSInput = open(CACTUSProbFN, 'w')
fCACTUSInput.write(TempFileNew)
fCACTUSInput.close()

if UseGeom:
    # Rejoin GeomFileLines into file string (with carriage returns)
    GeomFileNew='\n'.join(GeomFileLines) + '\n'

    # Write the CACTUS input file
    fGeomInput = open(GeomSFN, 'w')
    fGeomInput.write(GeomFileNew)
    fGeomInput.close()

    # Run geometry creation script
    GCLOutFN= 'Geom_CL_Output.' + DAKOutExt
    GCommand='octave ' + GeomSFN + ' ' + '&> ' + GCLOutFN
    os.system(GCommand)
#endif

# Run CACTUS and send output to file
CLOutFN= 'CACTUS_CL_Output.' + DAKOutExt
CCommand=CACTUSExe + ' ' + CACTUSProbFN + ' ' + '&> ' + CLOutFN
os.system(CCommand)

# Get CACTUS rev data
SIMFailFlag=False
fCACTUSDOut=open(CACTUSRevDataFN,'r')

RevHeader=fCACTUSDOut.readline().split(',')
RevData=[]
for line in fCACTUSDOut:
    DL=line.split(',')
    if DL: # not empty...
        RevData.append(map(float,DL))
    #endif
#endfor

fCACTUSDOut.close()

############## Add analysis specific post-processing here #############

# Check for non convergence of the NL iteration
NLConv=1
# Read last 4 lines
Lines=tail(CLOutFN,4)
for TL in Lines:
    if len(TL)>0:
        # Find error string
        if TL.find('NON-LINEAR ITERATION') >= 0:
            NLConv=0
            SIMFailFlag=True
        #endif
    else:
        DoneFile=True
    #endif
#endfor

# Check for potential presence of temporal odd-even non-physical mode
OECheck=0
# Read last 4 lines
Lines=tail(CACTUSTimeDataFN,4)
# Check if end-1 Cp value is consistent with end and end-2, to within tol...
if len(Lines) == 4:
    Tol=.01
    CPInd=3
    TData=[]
    for i in range(3):
        DL=Lines[i+1].split(',')
        if DL: # not empty...
            TData.append(float(DL[CPInd]))
        #endif
    #endfor
    if abs(TData[1]-(TData[2]+TData[0])/2.0) < Tol:
        OECheck=1
    #endif
#endif

#######################################################################

# Fill DAKOTA function output class with CACTUS data
for i,Func in enumerate(Funcs):
    OutTag=Func.Tag
    
    # Assign output based on DAKOTA tag names. Tags have the form 'ASV_n:TAG'
    A,Tag=OutTag.split(':')
    Tag=Tag.strip()
    Funcs[i].Tag=Tag
    if Tag == 'Cp':
        # Last Cp
        Ind=1
        if len(RevData)>0:
            Funcs[i].FnVal=RevData[len(RevData)-1][Ind]
        else:
            Funcs[i].FnVal=-999
        #endif   
        
    elif Tag == 'Ctr':
        # Last torque coeff
        Ind=3
        if len(RevData)>0:
            Funcs[i].FnVal=RevData[len(RevData)-1][Ind]
        else:
            Funcs[i].FnVal=-999
        #endif    
            
    elif Tag == 'Kp':
        # Last Kp
        Ind=2
        if len(RevData)>0:
            Funcs[i].FnVal=RevData[len(RevData)-1][Ind]
        else:
            Funcs[i].FnVal=-999
        #endif 
    
    elif Tag == 'Power':
        # Last Power
        Ind=3
        if len(RevData)>0:
            Funcs[i].FnVal=RevData[len(RevData)-1][Ind]
        else:
            Funcs[i].FnVal=-999
        #endif 
         
    elif Tag == 'Torque':
        # Last Torque
        Ind=4
        if len(RevData)>0:
            Funcs[i].FnVal=RevData[len(RevData)-1][Ind]
        else:
            Funcs[i].FnVal=-999
        #endif 
        
    elif Tag == 'CostFunc_MaxCP':
        # Cost function for optimizing for max CP with a method that can only do minimization.
        # Negative of last Cp...
        Ind=1
        if len(RevData)>0:
            Funcs[i].FnVal=-RevData[len(RevData)-1][Ind]
        else:
            Funcs[i].FnVal=-999
        #endif 
        
    elif Tag == 'Cp_3Rev':
        # Average Cp over last 3 revs
        Ind=1
        if len(RevData)>3:
            SumCp=0.0
            for RDInd in range(len(RevData)-3,len(RevData)):
                SumCp=SumCp+RevData[RDInd][Ind]
            #endfor
            Funcs[i].FnVal=SumCp/3
        else:
            SumCp=0.0
            for RD in RevData:
                SumCp=SumCp+RD[Ind]
            #endfor
            Funcs[i].FnVal=SumCp/len(RevData)
        #endif         
        
    elif Tag == 'NLConv':
        # Flag is 1 if NL iteration converged, 0 if the iteration failed to converge
        Funcs[i].FnVal=NLConv
        
    elif Tag == 'OECheck':
        # Flag is 1 if no temporal odd-even mode is present, 0 if odd-even mode is potentially present
        Funcs[i].FnVal=OECheck  
                     
    #endif
    
#endfor
        
# Write DAKOTA return file
fDAKOut = open(DAKOutFN, 'w')
for Func in Funcs:
    OutLine = str(Func.FnVal) + ' ' + Func.Tag + '\n'
    fDAKOut.write(OutLine)
#endfor

fDAKOut.close()

# Delete temporary files. If sim failed or keep files, keep temp files and save DAKOTA files with a .test ext...
if not SIMFailFlag and not KeepFiles:
    os.remove(CLOutFN)
    os.remove(CACTUSProbFN)
    os.remove(CACTUSParamFN)
    os.remove(CACTUSRevDataFN)
    os.remove(CACTUSTimeDataFN)
    if UseGeom:
        os.remove(GCLOutFN)
        os.remove(GeomSFN)
        os.remove(GeomFN)
    #endif
else:
    shutil.copy(DAKInFN,DAKInFN + '.test')  
    shutil.copy(DAKOutFN,DAKOutFN + '.test')
#endif
    

    
    
    
