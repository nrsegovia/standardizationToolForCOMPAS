import numpy as np
import h5py as h5
from os import path as pth
import pandas as pd
import sys
import os

endConditionDict = {11 : 2,
                    12 : 9,
                    13 : 4, # how to distinguish between one and two massless remnants? this one could be 5, too
                    9 : 4,
                    8 : 9, # COMPAS says that "stars are touching"
                    3 : 1,
                    14 : 3}

compasToLucaRuggeroSNDIct = {1 : 2,
    2 : 3,
    4 : 4,
    8 : 5,
    16 : 7,#USSN : 16,
    32 : 8,#AIC : 32,
    64  : 1,
    128 : 9}#HeSD : 128}
# Solve conflict: some SNe types are not expected in the standardized output

def decodeType(hurleyType):
    hurleyType = int(hurleyType)
    if (hurleyType < 10) | (hurleyType == 16):
        out = 100
        if hurleyType < 7:
            out = out + 20 + hurleyType
            if hurleyType == 0:
                out += 1
            elif hurleyType > 4:
                out *= 10
                if hurleyType == 5:
                    out += 1
                else:
                    out -= 8
        elif hurleyType < 10:
            out = out + 30 + hurleyType - 6
        else:
            out += 50
            out /= 10
    elif (hurleyType > 9) & (hurleyType < 13):
        out = 20 + hurleyType - 9
    elif (hurleyType > 12) & (hurleyType < 16):
        out = (hurleyType - 10) if hurleyType != 15 else -1
    else:
        out = 9
    return(int(out))

def decodeStatus(cond, helper, others = []):
    out = cond * 10 + helper
    if cond == 2:
        out = out * 10 + others[0]
    return out

def whosRLOF(primary, companion):
    primary = primary.astype(int)
    companion = companion.astype(int) * 2
    return primary + companion

local = pth.abspath(sys.argv[1])
# first argument is COMPAS' output directory

toOpen = [x for x in os.listdir(local) if x.endswith(".h5")]
if len(toOpen) != 1:
    print("Either no h5 file or more than one. Check the directory,")
    exit()

configs = {}
etc = []
with open(pth.join(local, "Run_Details"), "r") as details:
    for line in details.readlines():
        crt = line.strip() # Current line
        if " = " in crt:
            stuff = crt.split(" = ")
            configs[stuff[0]] = stuff[1]
        else:
            etc.append(crt)

DataSum = h5.File(pth.join(local, toOpen[0]), 'r')

# Stuff from Sys params
sysPars = DataSum["BSE_System_Parameters"]
sysID = sysPars["ID"][()] # should be used as UIDs
sysTime = sysPars["Time"][()]
sysSMA = sysPars["SemiMajorAxis"][()]
sysEC = sysPars["Eccentricity"][()]
sysST1 = np.fromiter((decodeType(x) for x in sysPars["Stellar_Type(1)"][()]), int)
sysM1 = sysPars["Mass(1)"][()]
sysR1 = sysPars["Radius(1)"][()]
sysT1 = sysPars["Teff(1)"][()]
sysMHe1 = sysPars["Mass_He_Core(1)"][()]
sysST2 = np.fromiter((decodeType(x) for x in sysPars["Stellar_Type(2)"][()]), int)
sysM2 = sysPars["Mass(2)"][()]
sysR2 = sysPars["Radius(2)"][()]
sysT2 = sysPars["Teff(2)"][()]
sysMHe2 = sysPars["Mass_He_Core(2)"][()]
sysTerm = np.fromiter((endConditionDict[x] for x in sysPars["Evolution_Status"][()]), int)

# CEE
cee = DataSum["BSE_Common_Envelopes"]
ceID = cee["ID"][()]
ceTime = cee["Time"][()]
ceSMA = cee["SemiMajorAxis"][()]
ceEC = cee["Eccentricity"][()]
ceST1 = np.fromiter((decodeType(x) for x in cee["Stellar_Type(1)"][()]), int)
ceM1 = cee["Mass(1)"][()]
ceR1 = cee["Radius(1)"][()]
ceT1 = cee["Teff(1)"][()]
ceMHe1 = cee["Mass_He_Core(1)"][()]
ceST2 = np.fromiter((decodeType(x) for x in cee["Stellar_Type(2)"][()]), int)
ceM2 = cee["Mass(2)"][()]
ceR2 = cee["Radius(2)"][()]
ceT2 = cee["Teff(2)"][()]
ceMHe2 = cee["Mass_He_Core(2)"][()]

ceIs1 = cee["RLOF(1)"][()]
ceIs2 = cee["RLOF(2)"][()]
ceWho = whosRLOF(ceIs1, ceIs2)


# RLOF
rlof = DataSum["BSE_RLOF"]
rlID = rlof["ID"][()]
rlTimePre = rlof["Time<MT"][()]
rlTimePost = rlof["Time>MT"][()]
rlSMAPre = rlof["SemiMajorAxis<MT"][()]
rlSMAPost = rlof["SemiMajorAxis>MT"][()]
rlECPre = rlof["Eccentricity<MT"][()]
rlECPost = rlof["Eccentricity>MT"][()]
rlST1Pre = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(1)<MT"][()]), int)
rlST1Post = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(1)>MT"][()]), int)
rlM1Pre = rlof["Mass(1)<MT"][()]
rlM1Post = rlof["Mass(1)>MT"][()]
rlR1Pre = rlof["Radius(1)<MT"][()]
rlR1Post = rlof["Radius(1)>MT"][()]
rlT1 = rlof["Teff(1)"][()]
rlMHe1 = rlof["Mass_He_Core(1)"][()]
rlST2Pre = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(2)<MT"][()]), int)
rlST2Post = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(2)>MT"][()]), int)
rlM2Pre = rlof["Mass(2)<MT"][()]
rlM2Post = rlof["Mass(2)>MT"][()]
rlR2Pre = rlof["Radius(2)<MT"][()]
rlR2Post = rlof["Radius(2)>MT"][()]
rlT2 = rlof["Teff(2)"][()]
rlMHe2 = rlof["Mass_He_Core(2)"][()]

rlIs1 = rlof["RLOF(1)"][()]
rlIs2 = rlof["RLOF(2)"][()]
rlWho = whosRLOF(rlIs1, rlIs2)

# SNe
sne = DataSum["BSE_Supernovae"]
snID = sne["ID"][()]
snTime = sne["Time"][()]
snSMA = sne["SemiMajorAxis"][()]
snEC = sne["Eccentricity"][()]
snST1 = np.fromiter((decodeType(x) for x in sne["Stellar_Type(1)"][()]), int)
snM1 = sne["Mass(1)"][()]
snR1 = sne["Radius(1)"][()]
snT1 = sne["Teff(1)"][()]
snMHe1 = sne["Mass_He_Core(1)"][()]
snST2 = np.fromiter((decodeType(x) for x in sne["Stellar_Type(2)"][()]), int)
snM2 = sne["Mass(2)"][()]
snR2 = sne["Radius(2)"][()]
snT2 = sne["Teff(2)"][()]
snMHe2 = sne["Mass_He_Core(2)"][()]

snWho = sne["Supernova_State"][()]
snType1 = sne["SN_Type(1)"][()]
snType2 = sne["SN_Type(2)"][()]
snWhichType = {1 : snType1,
               2 : snType2}

# Switch
switchLog = DataSum["BSE_Switch_Log"]
swID = switchLog["ID"][()]
swTime = switchLog["Time"][()]
swSMA = switchLog["SemiMajorAxis"][()]
swEC = switchLog["Eccentricity"][()]
swST1 = np.fromiter((decodeType(x) for x in switchLog["Stellar_Type(1)"][()]), int)
swM1 = switchLog["Mass(1)"][()]
swR1 = switchLog["Radius(1)"][()]
swT1 = switchLog["Teff(1)"][()]
swMHe1 = switchLog["Mass_He_Core(1)"][()]
swST2 = np.fromiter((decodeType(x) for x in switchLog["Stellar_Type(2)"][()]), int)
swM2 = switchLog["Mass(2)"][()]
swR2 = switchLog["Radius(2)"][()]
swT2 = switchLog["Teff(2)"][()]
swMHe2 = switchLog["Mass_He_Core(2)"][()]

swWho = switchLog["Star_Switching"][()]

DataSum.close()

# Physical properties needed for L0 standardized output. Problems: Temperature is the same pre and post RLOF
dictAll = {"Sys" : [sysID,sysTime, 8,sysSMA,sysEC,sysST1,sysM1,sysR1,sysT1,sysMHe1,sysST2,sysM2,sysR2,sysT2,sysMHe2],
           "CEE" : [ceID,ceTime, 51,ceSMA,ceEC,ceST1,ceM1,ceR1,ceT1,ceMHe1,ceST2,ceM2,ceR2,ceT2,ceMHe2],
           "RLOF" : [rlID,[rlTimePre, rlTimePost], [3, 4],[rlSMAPre,rlSMAPost],[rlECPre, rlECPost],[rlST1Pre, rlST1Post],[rlM1Pre, rlM1Post],[rlR1Pre, rlR1Post],[rlT1, rlT1],[rlMHe1, rlMHe1],[rlST2Pre, rlST2Post],[rlM2Pre, rlM2Post],[rlR2Pre, rlR2Post],[rlT2, rlT2],[rlMHe2, rlMHe2]],
           "SNe" : [snID,snTime, 2,snSMA,snEC,snST1,snM1,snR1,snT1,snMHe1,snST2,snM2,snR2,snT2,snMHe2],
           "Switch" : [swID,swTime, 1,swSMA,swEC,swST1,swM1,swR1,swT1,swMHe1,swST2,swM2,swR2,swT2,swMHe2]}

# Additional info: which star is undergoing the process, termination condition.
helperDict = {"Sys" : sysTerm,
           "CEE" : ceWho,
           "RLOF" : [rlWho, rlWho],
           "SNe" : snWho,
           "Switch" : swWho}

allLines = []
totSys = len(sysID)
for x in range(totSys):
    current = sysID[x]
    tempLine = f"{x},{current},"
    times = []
    eventPriority = []
    allTemps = []
    # dummy = np.argwhere(current).flatten()
    for log in dictAll.keys():
        crt = dictAll[log]
        mask = crt[0] == current
        args = np.argwhere(mask)
        if any(mask):
            if log == "RLOF":
                for r in range(2):
                    for arg in args:
                        thisTime = crt[1][r][arg][0]
                        others = ",".join(["{0:.4E}".format(crt[x][r][arg][0]) if isinstance(crt[x][r][arg][0], float) else str(crt[x][r][arg][0]) for x in range(3,15)])
                        event = crt[2][r]
                        eventPriority.append(event)
                        ending = decodeStatus(event, helperDict[log][r][arg][0])
                        thisLine = tempLine+f"{thisTime:.4E},{ending},{others}\n"
                        allTemps.append(thisLine)
                        times.append(thisTime)
            else:
                for arg in args:
                    thisTime = crt[1][arg][0]
                    others = ",".join(["{0:.4E}".format(crt[x][arg][0]) if isinstance(crt[x][arg][0], float) else str(crt[x][arg][0]) for x in range(3,15)])
                    extras = [compasToLucaRuggeroSNDIct[snWhichType[snWho[arg][0]][arg][0]]] if log == "SNe" else []
                    event = crt[2]
                    eventPriority.append(event)
                    ending = decodeStatus(event, helperDict[log][arg][0], extras)
                    thisLine = tempLine+f"{thisTime:.4E},{ending},{others}\n"
                    allTemps.append(thisLine)
                    times.append(thisTime)
    sortedByTime = np.lexsort((eventPriority, times))
    for s in sortedByTime:
        allLines.append(allTemps[s])




# Needed: Header with physical constants, COMPAS version and stuff required by the standard output.

cofVer = 0.2 # Common output format
level = "L0" # Cof Level as shown in the PDF document
ext = "" # Extension name, if used. Otherwise empty string
bps = "COMPAS" # Code name
ver = etc[1].split()[1][1:]
contact = "n.rsegovia@adfa.edu.au" # Contact email
nsys = totSys
nlines = len(allLines) + 3 # Should this be at the end. somehow? I will suppose that it counts the header as well.
metal = configs["metallicity"].split(", ")[0]

header = ["cofVer,cofLevel,cofExtension,bpsName,bpsVer,contact,NSYS,NLINES,Z\n",
          f"{cofVer},{level},{ext},{bps},{ver},{contact},{nsys},{nlines},{metal}\n",
          "ID,UID,time,event,semiMajor,eccentricity,type1,mass1,radius1,Teff1,massHecore1,type2,mass2,radius2,Teff2,massHecore2\n"]

# Output is saved in the COMPAS directory, under the name "standardOutputV0.csv"
output = open(pth.join(local, "standardOutputV0.csv"), "w")
output.writelines(header+allLines)
output.close()
