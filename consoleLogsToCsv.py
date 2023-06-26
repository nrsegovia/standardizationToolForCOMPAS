# This file needs the COMPAS output to be saved locally, i.e. $COMPAS/src/COMPAS --myParams &> storeConsoleLogs.txt,
# and then use it as first argument. By default it will create a new csv file named console.csv in the same directory
# as this file. Note that you have to specify --evolve-unbound-systems FALSE for unbound systems to show up in the logs
# as well as terminate the evolution at that point.

from os import path as pth
import sys

local = pth.split(pth.abspath(__file__))[0]
cont = 0
converted = open(pth.join(local, "console.csv"), "w")
converted.write("Index,LineNumber,Termination,EndStage1,EndStage2\n")
# allCauses = []
endDict = {'Allowed time exceeded' : "Time",
           'An error occurred' : "Error",
           'Double White Dwarf' : "DWD",
           'Evolution stopped' : "DCO",
           'Massless Remnant formed' : "MasslessRemnant",
           'Stars merged' : "Merger",
           'Unbound Binary' : "Unbound"}

idx = 0
with open(pth.abspath(sys.argv[1])) as file:
    for line in file:
        if line[0].isdigit():
            impList = line.strip().split(": ")
            cause = endDict[impList[1]]
            end = impList[2].split(" + ")
            end1 = end[0].split(" -> ")[1].strip(")")
            end2 = end[1].split(" -> ")[1].strip(")")
            if cause == "Error":
                converted.write(f"-99,{impList[0]},{cause},{end1},{end2}\n")
            else:
                converted.write(f"{idx},{impList[0]},{cause},{end1},{end2}\n")
                idx += 1

converted.close()
