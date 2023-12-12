# %matplotlib inline

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

import os, sys
import numpy as np
import h5py as h5
import pandas as pd

# +
# ## COMPAS of course has many different output file types, 
# ## and does not print a chronological set of events.
# ## To retrieve this for the UCB standardized outputs, 
# ## it is easiest to process each entry in 'BSE_RLOF', 
# ## 'BSE_Supernovae', and 'BSE_Switch_Log' separately, 
# ## and then combine and reorder these afterwards.  Notably, 
# ## ordering needs to be done on seeds first, and time second.
#
# ## Star IDs should be assigned at the end, COMPAS IDs should not be used
# ## This is because COMPAS IDs are not guarunteed to be unique (if multiple hdf5 runs are combied together), and there may be some missing if, e.g., the sampled initial conditions of a system were invalid. 
#
# ## Events require:
# ## ID UID time event semiMajor eccentricity type1 mass1 radius1 Teff1 massHecore1 type2 mass2 radius2 Teff2 massHeCore2
# -

#

# +
# TODO: check the whole thing
# TODO: what about events that double up? Say MT and SN in the same timestep?
# -



# +
def main():
    filepath = './COMPAS_Output/COMPAS_Output.h5'
    return load_COMPAS_data(filepath)
    
def load_COMPAS_data(filepath):
    ucb_events_obj = UCB_Events(filepath)
    getUCBEventForSupernova(ucb_events_obj)
    getUCBEventForMassTransfer(ucb_events_obj)
    #getUCBEventForStellarTypeChanges(ucb_events_obj) # TODO
    #getUCBEventForEndCondition(ucb_events_obj) # TODO
    return ucb_events_obj.getEvents()

class UCB_Events(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.Data = h5.File(filepath, 'r')
        self.all_UCB_events = None
        self.initialiaze_header()
        
    def initialiaze_header(self):
        self.header = {
            "cofVer" : 1.0, 
            "cofLevel": "L0",
            "cofExtension": "None", 
            "bpsName": "COMPAS",
            "bpsVer": self.Data['Run_Details']['COMPAS-Version'][()][0].decode('UTF-8'),
            "contact": "reinhold.willcox@gmail.com", 
            "NSYS": 0, 
            "NLINES": 0,
            #"Z": metallicity
        }
                       
    def update_header(self):
        df = self.all_UCB_events
        ids = df.loc[:,'ID']
        self.header.update({
            "NSYS": len(np.unique(ids)),
            "NLINES": len(ids)
        })

    def addEvents(self, uid=None, time=None, event=None, semiMajor=None, eccentricity=None, 
                        stellarType2=None, mass2=None, radius2=None, teff2=None, massHeCore2=None,
                        stellarType1=None, mass1=None, radius1=None, teff1=None, massHeCore1=None, 
                        scrapSeeds=None):
        columns   = [  "UID", "time", "event", "semiMajor", "eccentricity", 
                       "stellarType1", "mass1", "radius1", "teff1", "massHeCore1", 
                       "stellarType2", "mass2", "radius2", "teff2", "massHeCore2", 
                       "scrapSeeds" ] 
        data_list = [   uid,   time,   event,   semiMajor,   eccentricity, 
                        stellarType1,   mass1,   radius1,   teff1,   massHeCore1, 
                        stellarType2,   mass2,   radius2,   teff2,   massHeCore2, 
                        scrapSeeds   ] 
    
        # Want to enter data using name keywords, but all of them are required
        if np.any([ ii is None for ii in data_list ]):
            raise Exception("Can't skip any of the required input values. Currently missing {}".format(columns[ii]))
        new_events = pd.DataFrame(np.vstack(data_list).T, columns=columns)    
        if self.all_UCB_events is None:
            self.all_UCB_events = new_events
        else:
            self.all_UCB_events = pd.concat([self.all_UCB_events, new_events])

    def getEvents(self):
        df = self.all_UCB_events                                        # Convert to df for convenience
        
        # Clean up events - remove bad seeds
        allSeeds = df.loc[:,"UID"]                                      # Identify all seeds (or UIDs)
        scrapSeedMask = df.loc[:,"scrapSeeds"] == 1                     # Identify and create mask from the scrapSeeds column
        badSeedMask = np.in1d(allSeeds, allSeeds[scrapSeedMask])        # Create mask for all seeds to be scrapped, including rows with scrappable seeds that were not masked for it
        df = df[~badSeedMask]                                           # Return df without scrapped seeds
        df = df.drop(columns='scrapSeeds')                              # Remove scrap seeds column

        # Reorder the cells, add ID column
        df = df.sort_values(["UID", 'time'])                            # Reorder by uid (seed) first, then time second        
        uid_arr = df.loc[:,"UID"]                                       # Get list of UIDs
        uniq_uid = np.unique(uid_arr)                                   # Get the unique sorted list of UIDs
        uniq_id = np.arange(len(uniq_uid)) + 1                          # Start IDs counter at 1
        dict_uid_id = dict(zip(uniq_uid, uniq_id))                      # Map the UIDs to the IDs
        id_arr = np.vectorize(dict_uid_id.__getitem__)(uid_arr)         # Apply the map to the list of UIDs (with repeats)
        df.insert(0, "ID", id_arr)                                      # Insert the full IDs list at the front of the df
        
        self.all_UCB_events = df                                        # Convert back from df
        self.update_header()                                            # Update the header with new information
        return self.all_UCB_events



# -

main()


def verifyAndConvertCompasDataToUcbUsingDict(compasData, conversionDict):
    """
    General convenience function to verify and convert compas data arrays to their 
    equivalent values in UCB format, using the dictionaries defined variously below
    """
    try: 
        valid_types = np.array(list(conversionDict.keys()))
        assert np.all(np.in1d(compasData, valid_types))
    except:
        raise Exception('Invalid input array')
    return np.vectorize(conversionDict.get)(compasData)  # Quickly process entire input vector through converter dict




# ## Convert COMPAS output to UCB events 

#

# +
# Stellar type conversion dictionary

compasStellarTypeToUCBdict = {
    # COMPAS : UCB
    0: 121,
    1: 121,
    2: 122,
    3: 123,
    4: 124,
    5: 1251,
    6: 1252,
    7: 131,
    8: 132,
    9: 133,
    10: 21,
    11: 22,
    12: 23,
    13: 3,
    14: 4,
    15: -1,
    16: 9,  # CHE star - doesn't exist in UCB notation
}
# -

#


# +
## Supernova output processing

# Supernova conversion dictionary
compasSupernovaToUCBdict = {
    1:   2, # CCSN
    2:   3, # ECSN
    4:   4, # PISN
    8:   5, # PPISN
    16:  7, # USSN 
    32:  8, # AIC 
    64:  1, # Type Ia
    128: 9, # HeSD
}
UCB_SN_TYPE_DIRECT_COLLAPSE = 6 # Need to separately treat failed SNe (i.e direct collapse)


# -

def getUCBEventForSupernova(ucb_events_obj):

    Data = ucb_events_obj.Data
    SN = Data["BSE_Supernovae"]
    
    # Direct output
    uid = SN["SEED"][()]
    time = SN["Time"][()]
    semiMajorAxis = SN["SemiMajorAxis"][()]
    eccentricity = SN["Eccentricity"][()]
    mass1 = SN["Mass(1)"][()]
    mass2 = SN["Mass(2)"][()]
    radius1 = SN["Radius(1)"][()]
    radius2 = SN["Radius(2)"][()]
    teff1 = SN["Teff(1)"][()]
    teff2 = SN["Teff(2)"][()]
    massHeCore1 = SN["Mass_He_Core(1)"][()]
    massHeCore2 = SN["Mass_He_Core(2)"][()]
    stellarType1 = verifyAndConvertCompasDataToUcbUsingDict(SN["Stellar_Type(1)"][()], compasStellarTypeToUCBdict)
    stellarType2 = verifyAndConvertCompasDataToUcbUsingDict(SN["Stellar_Type(2)"][()], compasStellarTypeToUCBdict)
    
    # Processed output
    whichStar = SN["Supernova_State"][()] 
    assert np.all(np.in1d(whichStar, np.array([1, 2, 3]))) # TODO: need to address State 3 systems somehow.

    snType = verifyAndConvertCompasDataToUcbUsingDict(SN["SN_Type(SN)"][()], compasSupernovaToUCBdict)
    fb = SN['Fallback_Fraction(SN)'][()]
    snType[fb == 1] == UCB_SN_TYPE_DIRECT_COLLAPSE
    
    event = 2*100 + whichStar*10 + snType
    
    scrapSeeds = whichStar == 3 # need to remove these seeds at the end
    
    ucb_events_obj.addEvents(  uid=uid, time=time, event=event, semiMajor=semiMajorAxis, eccentricity=eccentricity, 
                               stellarType1=stellarType1, mass1=mass1, radius1=radius1, teff1=teff1, massHeCore1=massHeCore1, 
                               stellarType2=stellarType2, mass2=mass2, radius2=radius2, teff2=teff2, massHeCore2=massHeCore2,
                               scrapSeeds=scrapSeeds)


# +
## BSE_RLOF output processing

def getUCBEventForMassTransfer(ucb_events_obj):

    Data = ucb_events_obj.Data
    MT = Data["BSE_RLOF"]
    # Need to distinguish:
    # 1. Start of RLOF
    # 2. End of RLOF
    # 3. CEE events
    # 4. Mergers
    # 5. Contact phase (do we do this?)
    
    # Direct output
    uid = MT["SEED"][()]
    time = MT["Time>MT"][()]
    semiMajorAxis = MT["SemiMajorAxis>MT"][()]
    eccentricity = MT["Eccentricity>MT"][()]
    mass1 = MT["Mass(1)>MT"][()]
    mass2 = MT["Mass(2)>MT"][()]
    radius1 = MT["Radius(1)>MT"][()]
    radius2 = MT["Radius(2)>MT"][()]
    teff1 = MT["Teff(1)"][()]
    teff2 = MT["Teff(2)"][()]
    massHeCore1 = MT["Mass_He_Core(1)"][()]
    massHeCore2 = MT["Mass_He_Core(2)"][()]
    stellarType1 = verifyAndConvertCompasDataToUcbUsingDict(MT["Stellar_Type(1)>MT"][()], compasStellarTypeToUCBdict)
    stellarType2 = verifyAndConvertCompasDataToUcbUsingDict(MT["Stellar_Type(2)>MT"][()], compasStellarTypeToUCBdict)
    
    # Indirect output
    isRlof1 = MT["RLOF(1)>MT"][()] == 1
    isRlof2 = MT["RLOF(2)>MT"][()] == 1
    wasRlof1 = MT["RLOF(1)<MT"][()] == 1
    wasRlof2 = MT["RLOF(2)<MT"][()] == 1
    isCEE = MT["CEE>MT"][()] == 1
    isMerger = MT["Merger"][()] == 1
    scrapSeeds = np.zeros_like(uid).astype(bool) # TODO Scrap seeds if start of RLOF for both in the same timestep - is there any way to work with these??

    # Every mask in allmasks corresponds to an event in allevents
    allmasks = []
    allevents = []
    
    # Could make an events array of Nones, and then fill as they come up
    # The advantage of this is that for timesteps that qualify as 2 different events, you overwrite the wrong one...
    # Maybe I should just include the flags explicitly, that's probably more careful
    # So instead of doing a bunch of final MT timesteps and overwriting with any CEEs, I just include ~CEE in the condition.
    
    # 1. Start of RLOF.
    maskStartOfRlof1 = isRlof1 & ~wasRlof1 & ~isCEE
    maskStartOfRlof2 = isRlof2 & ~wasRlof2 & ~isCEE

    for ii in range(2):
        whichStar = ii+1 # either star 1 or 2
        allmasks.append([ maskStartOfRlof1, maskStartOfRlof2 ][ii])
        allevents.append( 3*10 + whichStar )

    # 2. End of RLOF
    maskFirstMtInParade1 = isRlof1 & ~wasRlof1
    maskFirstMtInParade2 = isRlof2 & ~wasRlof2
    for ii in range(2):
        whichStar = ii+1 # either star 1 or 2
        maskFirstMtInParade = [ maskFirstMtInParade1, maskFirstMtInParade2][ii]
        idxLastMtInParade = maskFirstMtInParade.nonzero()[0] - 1
        maskLastMtInParade = np.zeros_like(uid).astype(bool)
        maskLastMtInParade[idxLastMtInParade] = True
        allmasks.append(maskLastMtInParade & ~isCEE)
        allevents.append(4*10 + whichStar)

    # 3. CEE events - Process each CEE donor separately, plus double CEE for both
    maskAnyCEE = isCEE & ~isMerger
    whichStar = 1
    maskCEE1 = isRlof1 & ~isRlof2 & maskAnyCEE
    allmasks.append(maskCEE1)
    allevents.append(510 + whichStar)
    whichStar = 2
    maskCEE2 = isRlof2 & ~isRlof1 & maskAnyCEE
    allmasks.append(maskCEE2)
    allevents.append(510 + whichStar)
    whichStar = 3
    maskCEE3 = isRlof2 & isRlof1 & maskAnyCEE
    allmasks.append(maskCEE3)
    allevents.append(510 + whichStar)

    # 4. Mergers
    allmasks.append(isMerger)
    allevents.append(52)
    
    # 5. Contact phase (do we do this?)
    # TBD
    
    
    # Use masks to add all the events back into the array
    for mask, event in zip(allmasks, allevents):

        ucb_events_obj.addEvents(  uid=uid[mask], time=time[mask], event=event*np.ones_like(uid)[mask], semiMajor=semiMajorAxis[mask], eccentricity=eccentricity[mask], 
                                   stellarType1=stellarType1[mask], mass1=mass1[mask], radius1=radius1[mask], teff1=teff1[mask], massHeCore1=massHeCore1[mask], 
                                   stellarType2=stellarType2[mask], mass2=mass2[mask], radius2=radius2[mask], teff2=teff2[mask], massHeCore2=massHeCore2[mask],
                                   scrapSeeds=scrapSeeds[mask])


# +
## BSE_RLOF output processing

def getUCBEventForStellarTypeChanges(ucb_events_obj):
    
    Data = ucb_events_obj.Data
    SL = Data["BSE_Switch_Log"]
    
    # TODO
    
    for mask, event in zip(allmasks, allevents):

        ucb_events_obj.addEvents(  uid=uid[mask], time=time[mask], event=event*np.ones_like(uid)[mask], semiMajor=semiMajorAxis[mask], eccentricity=eccentricity[mask], 
                                   stellarType1=stellarType1[mask], mass1=mass1[mask], radius1=radius1[mask], teff1=teff1[mask], massHeCore1=massHeCore1[mask], 
                                   stellarType2=stellarType2[mask], mass2=mass2[mask], radius2=radius2[mask], teff2=teff2[mask], massHeCore2=massHeCore2[mask],
                                   scrapSeeds=scrapSeeds[mask])
# +
## BSE_RLOF output processing

def getUCBEventForEndCondition(ucb_events_obj):
    
    Data = ucb_events_obj.Data
    SP = Data["BSE_System_Parameters"]
    
    # TODO
    
    for mask, event in zip(allmasks, allevents):

        ucb_events_obj.addEvents(  uid=uid[mask], time=time[mask], event=event*np.ones_like(uid)[mask], semiMajor=semiMajorAxis[mask], eccentricity=eccentricity[mask], 
                                   stellarType1=stellarType1[mask], mass1=mass1[mask], radius1=radius1[mask], teff1=teff1[mask], massHeCore1=massHeCore1[mask], 
                                   stellarType2=stellarType2[mask], mass2=mass2[mask], radius2=radius2[mask], teff2=teff2[mask], massHeCore2=massHeCore2[mask],
                                   scrapSeeds=scrapSeeds[mask])


# -


if __name__ == "__main__":
    # for testing
    Data = h5.File('./COMPAS_Output/COMPAS_Output.h5', 'r')
    print(Data.keys())
    MT = Data['BSE_RLOF']
    SN = Data['BSE_Supernovae']
    SP = Data['BSE_System_Parameters']
    SL = Data['BSE_Switch_Log']



























# +
#endConditionDict = {11 : 2,
#                    12 : 9,
#                    13 : 4, # how to distinguish between one and two massless remnants? this one could be 5, too
#                    9 : 4,
#                    8 : 9, # COMPAS says that "stars are touching"
#                    3 : 1,
#                    14 : 3}
#
#
#
#
#
#def decodeStatus(cond, helper, others = []):
#    out = cond * 10 + helper
#    if cond == 2:
#        out = out * 10 + others[0]
#    return out
#
#def whosRLOF(primary, companion):
#    primary = primary.astype(int)
#    companion = companion.astype(int) * 2
#    return primary + companion
#
#local = os.path.abspath(sys.argv[1])
## first argument is COMPAS' output directory
#
#toOpen = [x for x in os.listdir(local) if x.endswith(".h5")]
#if len(toOpen) != 1:
#    print("Either no h5 file or more than one. Check the directory,")
#    exit()
#
#configs = {}
#etc = []
#with open(os.path.join(local, "Run_Details"), "r") as details:
#    for line in details.readlines():
#        crt = line.strip() # Current line
#        if " = " in crt:
#            stuff = crt.split(" = ")
#            configs[stuff[0]] = stuff[1]
#        else:
#            etc.append(crt)
#
#Data = h5.File(os.path.join(local, toOpen[0]), 'r')
#
## Stuff from Sys params
#sysPars = Data["BSE_System_Parameters"]
#sysID = sysPars["ID"][()] # should be used as UIDs
#sysTime = sysPars["Time"][()]
#sysSMA = sysPars["SemiMajorAxis"][()]
#sysEC = sysPars["Eccentricity"][()]
#sysST1 = np.fromiter((decodeType(x) for x in sysPars["Stellar_Type(1)"][()]), int)
#sysM1 = sysPars["Mass(1)"][()]
#sysR1 = sysPars["Radius(1)"][()]
#sysT1 = sysPars["Teff(1)"][()]
#sysMHe1 = sysPars["Mass_He_Core(1)"][()]
#sysST2 = np.fromiter((decodeType(x) for x in sysPars["Stellar_Type(2)"][()]), int)
#sysM2 = sysPars["Mass(2)"][()]
#sysR2 = sysPars["Radius(2)"][()]
#sysT2 = sysPars["Teff(2)"][()]
#sysMHe2 = sysPars["Mass_He_Core(2)"][()]
#sysTerm = np.fromiter((endConditionDict[x] for x in sysPars["Evolution_Status"][()]), int)
#
## CEE
#cee = Data["BSE_Common_Envelopes"]
#ceID = cee["ID"][()]
#ceTime = cee["Time"][()]
#ceSMA = cee["SemiMajorAxis"][()]
#ceEC = cee["Eccentricity"][()]
#ceST1 = np.fromiter((decodeType(x) for x in cee["Stellar_Type(1)"][()]), int)
#ceM1 = cee["Mass(1)"][()]
#ceR1 = cee["Radius(1)"][()]
#ceT1 = cee["Teff(1)"][()]
#ceMHe1 = cee["Mass_He_Core(1)"][()]
#ceST2 = np.fromiter((decodeType(x) for x in cee["Stellar_Type(2)"][()]), int)
#ceM2 = cee["Mass(2)"][()]
#ceR2 = cee["Radius(2)"][()]
#ceT2 = cee["Teff(2)"][()]
#ceMHe2 = cee["Mass_He_Core(2)"][()]
#
#ceIs1 = cee["RLOF(1)"][()]
#ceIs2 = cee["RLOF(2)"][()]
#ceWho = whosRLOF(ceIs1, ceIs2)
#
#
## RLOF
#rlof = Data["BSE_RLOF"]
#rlID = rlof["ID"][()]
#rlTimePre = rlof["Time<MT"][()]
#rlTimePost = rlof["Time>MT"][()]
#rlSMAPre = rlof["SemiMajorAxis<MT"][()]
#rlSMAPost = rlof["SemiMajorAxis>MT"][()]
#rlECPre = rlof["Eccentricity<MT"][()]
#rlECPost = rlof["Eccentricity>MT"][()]
#rlST1Pre = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(1)<MT"][()]), int)
#rlST1Post = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(1)>MT"][()]), int)
#rlM1Pre = rlof["Mass(1)<MT"][()]
#rlM1Post = rlof["Mass(1)>MT"][()]
#rlR1Pre = rlof["Radius(1)<MT"][()]
#rlR1Post = rlof["Radius(1)>MT"][()]
#rlT1 = rlof["Teff(1)"][()]
#rlMHe1 = rlof["Mass_He_Core(1)"][()]
#rlST2Pre = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(2)<MT"][()]), int)
#rlST2Post = np.fromiter((decodeType(x) for x in rlof["Stellar_Type(2)>MT"][()]), int)
#rlM2Pre = rlof["Mass(2)<MT"][()]
#rlM2Post = rlof["Mass(2)>MT"][()]
#rlR2Pre = rlof["Radius(2)<MT"][()]
#rlR2Post = rlof["Radius(2)>MT"][()]
#rlT2 = rlof["Teff(2)"][()]
#rlMHe2 = rlof["Mass_He_Core(2)"][()]
#
#rlIs1 = rlof["RLOF(1)"][()]
#rlIs2 = rlof["RLOF(2)"][()]
#rlWho = whosRLOF(rlIs1, rlIs2)
#
## SNe
#
## Switch
#switchLog = Data["BSE_Switch_Log"]
#swID = switchLog["ID"][()]
#swTime = switchLog["Time"][()]
#swSMA = switchLog["SemiMajorAxis"][()]
#swEC = switchLog["Eccentricity"][()]
#swST1 = np.fromiter((decodeType(x) for x in switchLog["Stellar_Type(1)"][()]), int)
#swM1 = switchLog["Mass(1)"][()]
#swR1 = switchLog["Radius(1)"][()]
#swT1 = switchLog["Teff(1)"][()]
#swMHe1 = switchLog["Mass_He_Core(1)"][()]
#swST2 = np.fromiter((decodeType(x) for x in switchLog["Stellar_Type(2)"][()]), int)
#swM2 = switchLog["Mass(2)"][()]
#swR2 = switchLog["Radius(2)"][()]
#swT2 = switchLog["Teff(2)"][()]
#swMHe2 = switchLog["Mass_He_Core(2)"][()]
#
#swWho = switchLog["Star_Switching"][()]
#
#Data.close()
#
## Physical properties needed for L0 standardized output. Problems: Temperature is the same pre and post RLOF
#dictAll = {"Sys" : [sysID,sysTime, 8,sysSMA,sysEC,sysST1,sysM1,sysR1,sysT1,sysMHe1,sysST2,sysM2,sysR2,sysT2,sysMHe2],
#           "CEE" : [ceID,ceTime, 51,ceSMA,ceEC,ceST1,ceM1,ceR1,ceT1,ceMHe1,ceST2,ceM2,ceR2,ceT2,ceMHe2],
#           "RLOF" : [rlID,[rlTimePre, rlTimePost], [3, 4],[rlSMAPre,rlSMAPost],[rlECPre, rlECPost],[rlST1Pre, rlST1Post],[rlM1Pre, rlM1Post],[rlR1Pre, rlR1Post],[rlT1, rlT1],[rlMHe1, rlMHe1],[rlST2Pre, rlST2Post],[rlM2Pre, rlM2Post],[rlR2Pre, rlR2Post],[rlT2, rlT2],[rlMHe2, rlMHe2]],
#           "SNe" : [snID,snTime, 2,snSMA,snEC,snST1,snM1,snR1,snT1,snMHe1,snST2,snM2,snR2,snT2,snMHe2],
#           "Switch" : [swID,swTime, 1,swSMA,swEC,swST1,swM1,swR1,swT1,swMHe1,swST2,swM2,swR2,swT2,swMHe2]}
#
## Additional info: which star is undergoing the process, termination condition.
#helperDict = {"Sys" : sysTerm,
#           "CEE" : ceWho,
#           "RLOF" : [rlWho, rlWho],
#           "SNe" : snWho,
#           "Switch" : swWho}
#
#allLines = []
#totSys = len(sysID)
#for x in range(totSys):
#    current = sysID[x]
#    tempLine = f"{x},{current},"
#    times = []
#    eventPriority = []
#    allTemps = []
#    # dummy = np.argwhere(current).flatten()
#    for log in dictAll.keys():
#        crt = dictAll[log]
#        mask = crt[0] == current
#        args = np.argwhere(mask)
#        if any(mask):
#            if log == "RLOF":
#                for r in range(2):
#                    for arg in args:
#                        thisTime = crt[1][r][arg][0]
#                        others = ",".join(["{0:.4E}".format(crt[x][r][arg][0]) if isinstance(crt[x][r][arg][0], float) else str(crt[x][r][arg][0]) for x in range(3,15)])
#                        event = crt[2][r]
#                        eventPriority.append(event)
#                        ending = decodeStatus(event, helperDict[log][r][arg][0])
#                        thisLine = tempLine+f"{thisTime:.4E},{ending},{others}\n"
#                        allTemps.append(thisLine)
#                        times.append(thisTime)
#            else:
#                for arg in args:
#                    thisTime = crt[1][arg][0]
#                    others = ",".join(["{0:.4E}".format(crt[x][arg][0]) if isinstance(crt[x][arg][0], float) else str(crt[x][arg][0]) for x in range(3,15)])
#                    if log == "SNe":
#                        whoIsGoing = snWho[arg][0]
#                        extras = [compasToLucaRuggeroSNDIct[snWhichType[whoIsGoing][arg][0]]]
#                        # The following lines do nothing for now, but we should use them in case we have to do something different when both stars go SN.
#                        # if whoIsGoing != 3:
#                        #     extras = [compasToLucaRuggeroSNDIct[snWhichType[whoIsGoing][arg][0]]]
#                        # else:
#                        #     pass 
#                    else:
#                        extras = []              
#                    event = crt[2]
#                    eventPriority.append(event)
#                    ending = decodeStatus(event, helperDict[log][arg][0], extras)
#                    thisLine = tempLine+f"{thisTime:.4E},{ending},{others}\n"
#                    allTemps.append(thisLine)
#                    times.append(thisTime)
#    sortedByTime = np.lexsort((eventPriority, times))
#    for s in sortedByTime:
#        allLines.append(allTemps[s])
#
#
#
#
## Needed: Header with physical constants, COMPAS version and stuff required by the standard output.
#
#cofVer = 0.2 # Common output format
#level = "L0" # Cof Level as shown in the PDF document
#ext = "" # Extensiona name, if used. Otherwise empty string
#bps = "COMPAS" # Code name
#ver = etc[1].split()[1][1:]
#contact = "n.rsegovia@adfa.edu.au" # Contact email
#nsys = totSys
#nlines = len(allLines) + 3 # Should this be at the end. somehow? I will suppose that it counts the header as well.
#metal = configs["metallicity"].split(", ")[0]
#
#header = ["cofVer,cofLevel,cofExtension,bpsName,bpsVer,contact,NSYS,NLINES,Z\n",
#          f"{cofVer},{level},{ext},{bps},{ver},{contact},{nsys},{nlines},{metal}\n",
#          "ID,UID,time,event,semiMajor,eccentricity,type1,mass1,radius1,Teff1,massHecore1,type2,mass2,radius2,Teff2,massHecore2\n"]
#
## Output is saved in the COMPAS directory, under the name "standardOutputV0.csv"
#output = open(os.path.join(local, "standardOutputV0.csv"), "w")
#output.writelines(header+allLines)
#output.close()
