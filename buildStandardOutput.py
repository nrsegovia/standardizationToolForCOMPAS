import os, sys
import numpy as np
import h5py as h5
import pandas as pd


# +
# TODO: check the whole thing
# TODO: what about events that double up? Say MT and SN in the same timestep?
# -



# +
def load_COMPAS_data(filepath, testing=False):
    ucb_events_obj = UCB_Events(filepath, testing)
    return ucb_events_obj.getEvents()

class UCB_Events(object):

    """
    COMPAS has many different output file types, and does not print a chronological set of events.
    To retrieve this for the UCB standardized outputs, it is easiest to process each entry in the 
    output files 'BSE_Switch_Log', 'BSE_RLOF', 'BSE_Supernovae', and 'BSE_System_Parameters' separately, 
    and then combine and reorder these afterwards. Notably, ordering needs to be done on seeds first, 
    and time second. This is done in getEvents().
    
    Star IDs are assigned at the end, COMPAS IDs (i.e SEED) is not used. This is because COMPAS IDs are 
    not guarunteed to be unique (if multiple hdf5 runs are combied together), and there may be some missing 
    if, e.g., the sampled initial conditions of a system were invalid. 
    
    Events require:
    ID UID time event semiMajor eccentricity type1 mass1 radius1 Teff1 massHecore1 type2 mass2 radius2 Teff2 massHeCore2
    """

    def __init__(self, filepath, testing=False):
        self.filepath = filepath
        self.testing = testing
        self.Data = h5.File(filepath, 'r')
        self.all_UCB_events = None
        self.initialiaze_header()
        # Calculate the events as prescribed for the UCBs
        self.getUCBEventForStellarTypeChanges()
        self.getUCBEventForSupernova()
        self.getUCBEventForMassTransfer()
        self.getUCBEventForEndCondition()
        
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

    def verifyAndConvertCompasDataToUcbUsingDict(self, compasData, conversionDict):
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

    
    ########################################################################################
    ### 
    ### Convert COMPAS output to UCB events 
    ### 
    ########################################################################################
    
    ############################################
    ### 
    ### BSE_Switch_Log output processing
    ### 
    ############################################
    
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
    
    def getUCBEventForStellarTypeChanges(self):
        
        SL = self.Data["BSE_Switch_Log"]
        
        # Direct output
        uid = SL["SEED"][()]
        time = SL["Time"][()]
        semiMajorAxis = SL["SemiMajorAxis"][()]
        eccentricity = SL["Eccentricity"][()]
        mass1 = SL["Mass(1)"][()]
        mass2 = SL["Mass(2)"][()]
        radius1 = SL["Radius(1)"][()]
        radius2 = SL["Radius(2)"][()]
        teff1 = SL["Teff(1)"][()]
        teff2 = SL["Teff(2)"][()]
        massHeCore1 = SL["Mass_He_Core(1)"][()]
        massHeCore2 = SL["Mass_He_Core(2)"][()]
        stellarType1 = self.verifyAndConvertCompasDataToUcbUsingDict(SL["Stellar_Type(1)"][()], self.compasStellarTypeToUCBdict)
        stellarType2 = self.verifyAndConvertCompasDataToUcbUsingDict(SL["Stellar_Type(2)"][()], self.compasStellarTypeToUCBdict)
        
        # Indirect output
        whichStar = SL['Star_Switching'][()]
        event = 10 + whichStar
        scrapSeeds = np.zeros_like(uid).astype(bool) # TODO Scrap seeds if start of RLOF for both in the same timestep - is there any way to work with these??
    
        self.addEvents(  uid=uid, time=time, event=event, semiMajor=semiMajorAxis, eccentricity=eccentricity, 
                         stellarType1=stellarType1, mass1=mass1, radius1=radius1, teff1=teff1, massHeCore1=massHeCore1, 
                         stellarType2=stellarType2, mass2=mass2, radius2=radius2, teff2=teff2, massHeCore2=massHeCore2,
                         scrapSeeds=scrapSeeds)
    
    ############################################
    ### 
    ### BSE_RLOF output processing
    ### 
    ############################################
    
    def getUCBEventForMassTransfer(self):
    
        MT = self.Data["BSE_RLOF"]
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
        stellarType1 = self.verifyAndConvertCompasDataToUcbUsingDict(MT["Stellar_Type(1)>MT"][()], self.compasStellarTypeToUCBdict)
        stellarType2 = self.verifyAndConvertCompasDataToUcbUsingDict(MT["Stellar_Type(2)>MT"][()], self.compasStellarTypeToUCBdict)
        
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
    
            self.addEvents(  uid=uid[mask], time=time[mask], event=event*np.ones_like(uid)[mask], semiMajor=semiMajorAxis[mask], eccentricity=eccentricity[mask], 
                             stellarType1=stellarType1[mask], mass1=mass1[mask], radius1=radius1[mask], teff1=teff1[mask], massHeCore1=massHeCore1[mask], 
                             stellarType2=stellarType2[mask], mass2=mass2[mask], radius2=radius2[mask], teff2=teff2[mask], massHeCore2=massHeCore2[mask],
                             scrapSeeds=scrapSeeds[mask])
    
    ############################################
    ### 
    ### BSE_Supernovae output processing
    ### 
    ############################################
    
    # Supernova conversion dictionary
    compasSupernovaToUCBdict = {
        # COMPAS : UCB
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
    
    def getUCBEventForSupernova(self):
    
        SN = self.Data["BSE_Supernovae"]
        
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
        stellarType1 = self.verifyAndConvertCompasDataToUcbUsingDict(SN["Stellar_Type(1)"][()], self.compasStellarTypeToUCBdict)
        stellarType2 = self.verifyAndConvertCompasDataToUcbUsingDict(SN["Stellar_Type(2)"][()], self.compasStellarTypeToUCBdict)
        
        # Indirect output
        whichStar = SN["Supernova_State"][()] 
        assert np.all(np.in1d(whichStar, np.array([1, 2, 3]))) # TODO: need to address State 3 systems somehow.
        scrapSeeds = whichStar == 3 # need to remove these seeds at the end
    
        snType = self.verifyAndConvertCompasDataToUcbUsingDict(SN["SN_Type(SN)"][()], self.compasSupernovaToUCBdict)
        fb = SN['Fallback_Fraction(SN)'][()]
        snType[fb == 1] == self.UCB_SN_TYPE_DIRECT_COLLAPSE
        event = 2*100 + whichStar*10 + snType    
        
        self.addEvents(  uid=uid, time=time, event=event, semiMajor=semiMajorAxis, eccentricity=eccentricity, 
                         stellarType1=stellarType1, mass1=mass1, radius1=radius1, teff1=teff1, massHeCore1=massHeCore1, 
                         stellarType2=stellarType2, mass2=mass2, radius2=radius2, teff2=teff2, massHeCore2=massHeCore2,
                         scrapSeeds=scrapSeeds)
    
    ############################################
    ### 
    ### BSE_System_Parameters output processing
    ### 
    ############################################
    
    # End condition conversion dictionary
    compasOutcomeToUCBdict = {
        # COMPAS : UCB
        1:  -1, # simulation completed?
        2:  9, # error
        3:  1, # max time reached
        4:  1, # max timesteps reached, kind of the same as above
        5:  9, # error
        6:  9, # error
        7:  -1, # time exceeded dco merger time?
        8:  -1, # stars touching ?
        9:  4, # merger
        10: 4, # merger
        11: 2, # dco formed
        12: 2, # dwd formed
        13: 4, # massless remnant
        14: 3, # unbound    
    }
    
    # -1 applies if the compas output description is unclear, just toss these seeds for now.
    # Q: how does 85 happen and not 84? Wouldn't simulations stop at 84?
        
    # UCBs
    # 81 - max time reached
    # 82 - both components are compact remnants - RTW: including WDs?
    # 83 - the binary system is dissociated
    # 84 - only one object is left (e.g. due to a merger or because the companion has been disrupted)
    # 85 - nothing left (both components are massless remnants)
    # 89 - other: a terminating condition different from any previous one
    
    # COMPAS
    # Simulation completed = 1     
    # Evolution stopped because an error occurred = 2
    # Allowed time exceeded = 3
    # Allowed timesteps exceeded = 4
    # SSE error for one of the constituent stars = 5
    # Error evolving binary = 6
    # Time exceeded DCO merger time = 7
    # Stars touching = 8
    # Stars merged = 9
    # Stars merged at birth = 10
    # DCO formed = 11
    # Double White Dwarf formed = 12
    # Massless Remnant formed = 13
    # Unbound binary = 14
    
    def getUCBEventForEndCondition(self):
        
        SP = self.Data["BSE_System_Parameters"]
        
        # Direct output
        uid = SP["SEED"][()]
        time = SP["Time"][()]
        semiMajorAxis = SP["SemiMajorAxis"][()]
        eccentricity = SP["Eccentricity"][()]
        mass1 = SP["Mass(1)"][()]
        mass2 = SP["Mass(2)"][()]
        radius1 = SP["Radius(1)"][()]
        radius2 = SP["Radius(2)"][()]
        teff1 = SP["Teff(1)"][()]
        teff2 = SP["Teff(2)"][()]
        massHeCore1 = SP["Mass_He_Core(1)"][()]
        massHeCore2 = SP["Mass_He_Core(2)"][()]
        stellarType1 = self.verifyAndConvertCompasDataToUcbUsingDict(SP["Stellar_Type(1)"][()], self.compasStellarTypeToUCBdict)
        stellarType2 = self.verifyAndConvertCompasDataToUcbUsingDict(SP["Stellar_Type(2)"][()], self.compasStellarTypeToUCBdict)
        
        # Indirect output
        evolStatus = self.verifyAndConvertCompasDataToUcbUsingDict(SP["Evolution_Status"][()], self.compasOutcomeToUCBdict)
        assert np.all(np.in1d(evolStatus, np.array([1, 2, 3, 4, 5, 9, -1]))) 
        scrapSeeds = evolStatus == -1 # -1 means I don't understand the compas outcome
        if self.testing:
            if np.any(scrapSeeds):
                print("There were {} strange evolutionary outcomes".format(np.sum(scrapSeeds)))
                print("Their seeds were:")
                print("[" + ", ".join(uid[scrapSeeds].astype(str)) + "]")
        event = 8*10 + evolStatus
        
        self.addEvents(  uid=uid, time=time, event=event, semiMajor=semiMajorAxis, eccentricity=eccentricity, 
                         stellarType1=stellarType1, mass1=mass1, radius1=radius1, teff1=teff1, massHeCore1=massHeCore1, 
                         stellarType2=stellarType2, mass2=mass2, radius2=radius2, teff2=teff2, massHeCore2=massHeCore2,
                         scrapSeeds=scrapSeeds)


# -

