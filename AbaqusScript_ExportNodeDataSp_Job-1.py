# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:49:20 2019

@author: AK115190
"""
#OutputSpoiler Info

from odbAccess import *
from abaqusConstants import *
import numpy as np
import csv
import sys
import datetime

jobName = 'Job-1'
odb = openOdb(path=jobName+'.odb', readOnly=True)
#lastFrame = odb.steps['Step-1'].frames[-1]

exportSets = []
exportSets.append({'Name':'Part-1', 'Step':'Step-1', 'Instance':'PART-1-1', 'NodeSet':'SET-1', 'sectionPoint': 'Top'})

def updateProgressbar(i): # i in range(101)=[0 - 100]
    sys.__stdout__.write("\r%3d%%" % i)
    barStr = '#'*int(i/2)
    barStr = barStr + ' '*(50-len(barStr))
    sys.__stdout__.write('|'+barStr+'|')
    sys.__stdout__.flush()

for exportSet in exportSets:
    
    print >> sys.__stdout__, 'Export data for %s (%s, %s, %s)' % (exportSet['Name'], exportSet['Step'], exportSet['Instance'], exportSet['NodeSet'])
    
    for frame in odb.steps[exportSet['Step']].frames:
        
        print >> sys.__stdout__, 'Frame: %d' % (frame.frameId)
        allFields = frame.fieldOutputs
        
        ##############################
        ####Export Field Outputs at NODES
        ##############################
        
        #Define nodeSet
        nodeSet = odb.rootAssembly.instances[exportSet['Instance']].nodeSets[exportSet['NodeSet']]
        #Define elementSet
        elemSet = odb.rootAssembly.instances[exportSet['Instance']].elementSets[exportSet['NodeSet']]
        
        if exportSet['sectionPoint']=='Bottom':
            spNr = 0
        elif exportSet['sectionPoint']=='Top':
            spNr = -1
        else:
            spNr = None
        
        #get section point of top surface (= last section point) --> assume that all elements have same sectionCategory and sectionPoints
        sP = elemSet.elements[0].sectionCategory.sectionPoints[spNr]
        
        # Get Stresses at nodes (position = ELEMENT_NODAL) for each node
        # Note: for each node the extrapolated stress of neighhor elements is in this stressSet, therefor the results for each node must be averaged (see below)
        stressSet = allFields['S'].getSubset(region=nodeSet,position=ELEMENT_NODAL, sectionPoint=sP)
        # --> Attention: for shell elements top and bottom results are exported to stressSet
        # important to use sectionPoint option
        nStresses = len(stressSet.values)
        allStresses = np.empty([nStresses,10])
        for i,v in enumerate(stressSet.values):
            allStresses[i,0] = v.nodeLabel 
            allStresses[i,1] = v.data[0]  #S11
            allStresses[i,2] = v.data[1]  #S22
            allStresses[i,3] = v.data[2]  #S33
            allStresses[i,4] = v.data[3]  #S12
            #allStresses[i,1:5] = v.data #alternative to export S11, S22, S33 and S12 in one line
            allStresses[i,5] = v.mises
            allStresses[i,6] = v.maxPrincipal
            allStresses[i,7] = v.minPrincipal
            allStresses[i,8] = v.maxInPlanePrincipal
            allStresses[i,9] = v.minInPlanePrincipal
        
        strainSet = allFields['E'].getSubset(region=nodeSet,position=ELEMENT_NODAL, sectionPoint=sP)
        # --> Attention: for shell elements top and bottom results are exported to strainSet
        # important to use sectionPoint option
        nStrains = len(strainSet.values)
        allStrains = np.empty([nStrains,10])
        for i,v in enumerate(strainSet.values):
            allStrains[i,0] = v.nodeLabel 
            allStrains[i,1] = v.data[0]  #E11
            allStrains[i,2] = v.data[1]  #E22
            allStrains[i,3] = v.data[2]  #E33
            allStrains[i,4] = v.data[3]  #E12
            allStrains[i,5] = v.mises
            allStrains[i,6] = v.maxPrincipal  
            allStrains[i,7] = v.minPrincipal
            allStrains[i,8] = v.maxInPlanePrincipal
            allStrains[i,9] = v.minInPlanePrincipal 
            
        dispSet = allFields['U'].getSubset(region=nodeSet)
        nDisps = len(dispSet.values)
        allDisps = np.empty([nDisps,5])
        for i,v in enumerate(dispSet.values):
            allDisps[i,0] = v.nodeLabel 
            allDisps[i,1] = v.data[0]  #U1
            allDisps[i,2] = v.data[1]  #U2
            allDisps[i,3] = v.data[2]  #U3
            allDisps[i,4] = v.magnitude  #U Magnitude
        
        with open('%s_NodalFO_%s_frameId%02d_sp%s.csv' %(jobName,exportSet['Name'],frame.frameId,exportSet['sectionPoint']), 'wb') as csvfile:
            csvWriter = csv.writer(csvfile, delimiter=';')
            csvWriter.writerow(['nodeLabel', 'x', 'y', 'z',
                                'U_1', 'U_2', 'U_3', 'U_mag',
                                'S_11', 'S_22', 'S_33', 'S_12', 'S_mises', 'S_maxPrincipal', 'S_minPrincipal', 'S_maxInPlanePrincipal', 'S_minInPlanePrincipal', 
                                'E_11', 'E_22', 'E_33', 'E_12', 'E_mises', 'E_maxPrincipal', 'E_minPrincipal', 'E_maxInPlanePrincipal', 'E_minInPlanePrincipal'])
            for i,n in enumerate(nodeSet.nodes):
                indicesStress = np.where(allStresses[:,0] == float(n.label)) #at which indices are stress results for the current node
                stressesAtIndices = allStresses[indicesStress,:] #get this results and average for each stress invariant (mises, maxPrincipal)
                S_11 = np.mean(stressesAtIndices[0][:,1])
                S_22 = np.mean(stressesAtIndices[0][:,2])
                S_33 = np.mean(stressesAtIndices[0][:,3])
                S_12 = np.mean(stressesAtIndices[0][:,4])
                S_mises = np.mean(stressesAtIndices[0][:,5])
                S_maxPrincipal = np.mean(stressesAtIndices[0][:,6])
                S_minPrincipal = np.mean(stressesAtIndices[0][:,7])
                S_maxInPlanePrincipal = np.mean(stressesAtIndices[0][:,8])
                S_minInPlanePrincipal = np.mean(stressesAtIndices[0][:,9])
                
                indicesStrain = np.where(allStrains[:,0] == float(n.label)) #at which indices are strain results for the current node
                strainsAtIndices = allStrains[indicesStrain,:] #get this results and average to get single value for the current node
                E_11 = np.mean(strainsAtIndices[0][:,1])
                E_22 = np.mean(strainsAtIndices[0][:,2])
                E_33 = np.mean(strainsAtIndices[0][:,3])
                E_12 = np.mean(strainsAtIndices[0][:,4])
                E_mises = np.mean(strainsAtIndices[0][:,5])
                E_maxPrincipal = np.mean(strainsAtIndices[0][:,6])
                E_minPrincipal = np.mean(strainsAtIndices[0][:,7])
                E_maxInPlanePrincipal = np.mean(strainsAtIndices[0][:,8])
                E_minInPlanePrincipal = np.mean(strainsAtIndices[0][:,9])
                
                indicesDisp = np.where(allDisps[:,0] == float(n.label)) #at which indices are strain results for the current node
                dispsAtIndices = allDisps[indicesDisp,:] # averaging is accually not neccessary for displacements (node results!)
                U_1 = np.mean(dispsAtIndices[0][:,1])
                U_2 = np.mean(dispsAtIndices[0][:,2])
                U_3 = np.mean(dispsAtIndices[0][:,3])
                U_mag = np.mean(dispsAtIndices[0][:,4])
                
                csvWriter.writerow([n.label, n.coordinates[0], n.coordinates[1], n.coordinates[2],
                                    U_1, U_2, U_3, U_mag,
                                    S_11, S_22, S_33, S_12, S_mises, S_maxPrincipal, S_minPrincipal, S_maxInPlanePrincipal, S_minInPlanePrincipal,
                                    E_11, E_22, E_33, E_12, E_mises, E_maxPrincipal, E_minPrincipal, E_maxInPlanePrincipal, E_minInPlanePrincipal])
                # Updating progress bar
                updateProgressbar(int((i+1)*100.0/len(nodeSet.nodes)))
        print >> sys.__stdout__, ''

odb.close()