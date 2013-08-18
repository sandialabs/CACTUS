#!/usr/bin/env python

# DAKOTA sub-level analysis driver. (Optimization) 
# Runs an initial optimization method and optionally refines the results with a second optimization method. Returns
# Typical application of this script is a coarse/fine hybrid approach in which a global method (like DAKOTA's DIRECT method) is used to find likely 
# global minimizer candidates, and the results are passed to a fast local gradient search method (like DAKOTA's DOTBFGS) for refinement.

# Analysis component inputs:
#   1. Input 'Y' if second refinement optimization is to be performed, otherwise 'N'
#   2. Template DAKOTA input file for the initial optimization analysis
#   3. Template DAKOTA input file for the refined optimization analysis
#   Other optional analysis components in no particular order (KeepFiles, etc...). See those recognized below...

# The template files should contain tags to be replaced with DAKOTA variables. Tags should be the top-level DAKOTA variable 
# descriptor name enclosed in curly braces (ex. {RPM}). In order to use the output of the first optimization as the initial conditions for the second,
# the tags for the optimization parameters in the refined optimization template should be the DAKOTA variable descriptor used for the corresponding parameter 
# in the initial optimization template (enclosed in curly braces).


import os
import shutil
import sys
import string
import math

# Note: To debug, call the python debugger module, pdb, as a script with this script as the argument (python -m pdb TAOS_Driver.py)
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
        
############################        
    
# Read input and output filenames from the command line
if len(sys.argv) is 3:
    DAKInFN=sys.argv[1]   
    DAKOutFN=sys.argv[2] 
else:
    print 'No files given...'
    sys.exit(1)
#endif

# Get DAKOTA filebase and run index (if exists)
DAKInFB,DAKInExt=DAKInFN.split('.',1)
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


########################################################################################
   
# If first AC is 'Y', the refinement optimization will be run, using the output of the initial
# optimization as an initial condition. The second and third AC must be template files for the initial
# and refinement optimization DAKOTA analyses.
RunFine=False
if ACs[0].Code == 'Y':
    RunFine=True
#endif
   
# Get the sub-level analysis template file for initial opt, and write a modified file for this call
SLATempFN=ACs[1].Code
TempFN=SLATempFN.rsplit('/',1) # split directory and filename
SLAFB,SLAExt=TempFN[len(TempFN)-1].split('.',1)

SLAFN=SLAFB + '_Mod'
if DAKInd:
    SLAFN=SLAFN + DAKInd
#endif
SLAInFN1=SLAFN + '.' + SLAExt
SLADataFN1=SLAFN + '_Data.dat'
SLADriverInFB=DAKInFB + '_INIT_SLA.' + DAKInExt
SLADriverOutFB=DAKOutFB + '_INIT_SLA.' + DAKOutExt

fTemp=open(SLATempFN,'r')
TempFile=fTemp.read()  # read the file into a string
fTemp.close()

# Replace tags with top-level DAKOTA values
for Var in Vars:
    TAGStr='{' + Var.Tag + '}'
    TempFile=TempFile.replace(TAGStr,Var.Val)
#endfor

# Write the analysis driver input and output filenames for the sub-level analysis...
# (these tags must exist in the sub-level analysis DAKOTA template)
TempFile=TempFile.replace('{SLADriverInFB}',SLADriverInFB) 
TempFile=TempFile.replace('{SLADriverOutFB}',SLADriverOutFB)   

# Write the sub-level analysis DAKOTA input file
fSLAInput = open(SLAInFN1, 'w')
fSLAInput.write(TempFile)
fSLAInput.close()

# Run DAKOTA sub-level analysis and save command line output to file
CLOutFN1= 'INIT_SLA_CL_Output.' + DAKOutExt
DCommand='dakota ' + SLAInFN1 + ' ' + '&> ' + CLOutFN1
os.system(DCommand)

# Get initial results data from the sub-level analysis (command line output)
fSLADOut=open(CLOutFN1,'r')

# Get variable tags and data from optimization output format, and insert into dictionary
OptOutDict={}
DoneFile=False
while not DoneFile:
    TL=fSLADOut.readline()
    if len(TL)>0:
        # Match header string
        if TL.strip() == '<<<<< Best parameters          =':
            # Collect available parameter tags, along with data for each
            DoneTags=False
            while not DoneTags:
                DL=string.split(fSLADOut.readline())
                if DL[0].strip() == '<<<<<': # done with parameters
                    DoneTags=True
                else:
                    OutTag=DL[1].strip()
                    OptOutDict[OutTag]=float(DL[0])
                #endif
            #endwhile
            
            # get objective function value
            DL=string.split(fSLADOut.readline())
            OptOutDict['Obj_Func']=float(DL[0])
        #endif
    else:
        DoneFile=True
    #endif
#endwhile

fSLADOut.close()


# Run refined optimization if requested
if RunFine:
    # Get the sub-level analysis template file, and write a modified file for this call
    SLATempFN=ACs[2].Code
    SLAFB,SLAExt=SLATempFN.split('.',1)
    
    SLAFN=SLAFB + '_Mod'
    if DAKInd:
        SLAFN=SLAFN + DAKInd
    #endif
    SLAInFN2=SLAFN + '.' + SLAExt
    SLADataFN2=SLAFN + '_Data.dat'
    SLADriverInFB=DAKInFB + '_FINE_SLA.' + DAKInExt
    SLADriverOutFB=DAKOutFB + '_FINE_SLA.' + DAKOutExt
    
    fTemp=open(SLATempFN,'r')
    TempFile=fTemp.read()  # read the file into a string
    fTemp.close()
    
    # Replace tags with top-level DAKOTA values
    for Var in Vars:
        TAGStr='{' + Var.Tag + '}'
        TempFile=TempFile.replace(TAGStr,Var.Val)
    #endfor
    
    # Replace initial conditions tags with outputs from initial opt. The initial condition tags
    # in the fine optimization template file should be the parameter tag name enclosed 
    # in curly braces (ex. {Rel_KCAS}) so that they can be matched with their corresponding 
    # opt. output dictionary entries...
    for DKey, DValue in OptOutDict.iteritems():
        TAGStr='{' + DKey + '}'
        TempFile=TempFile.replace(TAGStr,str(DValue))
    #endfor
    
    # Write the analysis driver input and output filenames for the sub-level analysis...
    # (these tags must exist in the sub-level analysis DAKOTA template)
    TempFile=TempFile.replace('{SLADriverInFB}',SLADriverInFB) 
    TempFile=TempFile.replace('{SLADriverOutFB}',SLADriverOutFB)   
    
    # Write the sub-level analysis DAKOTA input file
    fSLAInput = open(SLAInFN2, 'w')
    fSLAInput.write(TempFile)
    fSLAInput.close()
    
    # Run DAKOTA sub-level analysis and save command line output to file
    CLOutFN2= 'FINE_SLA_CL_Output.' + DAKOutExt
    DCommand='dakota ' + SLAInFN2 + ' ' + '&> ' + CLOutFN2
    os.system(DCommand)
    
    # Get fine results data from the sub-level analysis (command line output)
    fSLADOut=open(CLOutFN2,'r')
    
    # Get variable tags and data from optimization output format, and insert into dictionary
    DoneFile=False
    while not DoneFile:
        TL=fSLADOut.readline()
        if len(TL)>0:
            # Match header string
            if TL.strip() == '<<<<< Best parameters          =':
                # Collect available parameter tags, along with data for each
                DoneTags=False
                while not DoneTags:
                    DL=string.split(fSLADOut.readline())
                    if DL[0].strip() == '<<<<<': # done with parameters
                        DoneTags=True
                    else:
                        OutTag=DL[1].strip()
                        OptOutDict[OutTag]=float(DL[0])
                    #endif
                #endwhile
                
                # get objective function value
                DL=string.split(fSLADOut.readline())
                OptOutDict['Obj_Func']=float(DL[0])
            #endif
        else:
            DoneFile=True
        #endif
    #endwhile
    
    fSLADOut.close()
    
#endif

############## Add analysis specific post-processing here #############

# Check for sim failure
SIMFailFlag=False

############################################################################################

# Fill top-level DAKOTA function output class with values derived from sub-level data output
for i,Func in enumerate(Funcs):
    OutTag=Func.Tag
    
    # Assign output based on top-level and sub-level DAKOTA tag names.
    # The top-level output descriptor tag is first checked against the variable tag descriptor and assigned
    # the optimized value of that variable. Note that the objective function descriptor tag name in the DAKOTA optimization 
    # template file doesn't propagate to the optimizer output, so to request the objective function value, use the top-level output 
    # descriptor tag 'Obj_Func'. Any tags not recognized as sublevel analysis parameters are interpreted below...
    A,Tag=OutTag.split(':')
    Tag=Tag.strip()
    Funcs[i].Tag=Tag
    
    # Get sub-level output, if found
    if Tag in OptOutDict:
        Funcs[i].FnVal=OptOutDict[Tag]   
    #endif

#endfor
  
# Write DAKOTA return file
fDAKOut = open(DAKOutFN, 'w')
for Func in Funcs:
    OutLine = str(Func.FnVal) + ' ' + Func.Tag + '\n'
    fDAKOut.write(OutLine)
#endfor

fDAKOut.close()

# Delete temporary files. If sim failed, keep temp files and save DAKOTA input and output files with a .test ext...
if not SIMFailFlag and not KeepFiles:
    os.remove(CLOutFN1)
    os.remove(SLAInFN1)
    if RunFine:
        os.remove(CLOutFN2)
        os.remove(SLAInFN2)
    #endif
else:
    shutil.copy(DAKInFN,DAKInFN + '.test') 
    shutil.copy(DAKOutFN,DAKOutFN + '.test')  
#endif
    
    

    
    
    
