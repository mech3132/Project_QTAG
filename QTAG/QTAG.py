#!bin/bash python3

#==========================================
# QUANTIFICATION OF TURNOVER ALONG GRADIENTS
# By Melissa Chen
# Version 2022.03.22
# Built under python 3.8.5
# Created April 10th, 2017
#==========================================

# REQUIREMENTS:
	# Imports:
		# numpy
		# math
		# scipy stats
		# argparse
		# os
		# subprocess
		# time sleep
		# sys
		# progressbar
		# matplotlib pyplot
		# statsmodels stats.multitest
import math 
import numpy # For standard deviation and quantiles
	# Required for BINNING and BINNING SALINITY sections
from scipy import stats # For welch's t-test
	# Required TYPETAXA sections
# from graphics import *
import argparse # for user import
import os
import subprocess
from time import sleep
import sys
import copy # to copy lists and dictionaries
from progressbar import ProgressBar
import matplotlib.pyplot as plt # For plotting
# from scipy.stats.distributions import chi2 # for LLR test
import statsmodels.stats.multitest as multitest
#import QTAG_func.py 
#===================== ARG PARSE ========================#

parser = argparse.ArgumentParser(
	description="Bins and classifies OTUs according to gradient specialization\n NOTE: WRITE MORE LATER \n\n\n")
parser.add_argument(
	'-t',
	'--taxaTable',
	help = "Path to taxa table file, text format. Header for OTUs must be '#OTU ID'. Preceding lines beginning with # are ignored. ",
	required = True)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'Path to metadata file. First column are site names that match taxaTable sites. Metadata must have at least one other column that indicates gradient values (see metadata_name). Header required.',
	required = True)
parser.add_argument(
	'-M',
	'--metadata_name',
	help = 'Header name of column in metadata file that specifies the gradient variable to analyze',
	required = True)
parser.add_argument(
	'-o',
	'--output_dir',
	help = 'Output directory [default: GradientProcessing]',
	required = False,
	default = 'GradientProcessing')
parser.add_argument(
	'--minX',
	help = "Minimum value for boundary X. IMPORTANT: default rounds the minimum gradient value down to nearest whole number. There may be situations where this value should be set manually. [default is minVal, the minimum gradient value]",
	required = False,
	default = 'Check')
parser.add_argument(
	'--maxY',
	help = 'Maximum value for boundary Y IMPORTANT: default rounds the maximum gradient value up to nearest whole number. There may be situations where this value should be set manually. [default is maxVal,the maximum gradient value]',
	required = False,
	default = 'Check')
parser.add_argument(
	'--XYdiff',
	help = 'Minimum difference between X and Y when iterating through X/Y combinations. Note: XYdiff should always be larger than your unit_Size. [default 2]',
	required = False,
	type = float,
	default = 2)
parser.add_argument(
	'--unitSize',
	help = 'Unit size of gradient to iterate through when fitting model. [default 1]',
	required = False,
	default = 1)
parser.add_argument(
	'--minCountBin',
	help = 'Minimum number of taxa in each bin (A, B, C) during model selection in order for that bin combination to be considered. (Equal or greater than minCountBin) [Default: 3]',
	required = False,
	default = 3)
parser.add_argument(
	'--minCountTable',
	help = 'Minimum number of reads found in entire taxa table for an taxa to be kept. [Default: 50]',
	required = False,
	default = 50)
parser.add_argument(
	'--minCountOTUinSample',
	help = 'Minimum number of reads of any taxa in a sample in order for that taxa to be considered "non-zero". If the number of reads of that taxa in a sample is LESS than this value, then it will be changed to have an abundance of zero. [Default: 5]',
	required = False,
	default = 5)
parser.add_argument(
	'--minSamplePres',
	help = 'Threshold for number of samples an taxa must appear in. [default 3]',
	required = False,
	default = 3)
parser.add_argument(
	'--pvalthresh',
	help = 'Significance threshold to use when comparing differences in bins. P-values are FDR-corrected during multiple comparisons. [Default: 0.05]',
	required = False,
	default = 0.05)
parser.add_argument(
	'--qquant',
	help = 'Quantile to use during quantile regression when finding best-fit model for taxa distribution across salinity. [Default: 0.95]',
	required = False,
	default = 0.95)	
parser.add_argument(
	'--bandwidth_percent',
	help = 'The bandwidth of samples to consider when creating quantile regression model; unit is percent of total gradient breadth (e.g. if the gradient extends from 0 to 100, a bandwidth percent of 0.1 would use a bandwidth of 10) [Default: 0.2]',
	required = False,
	default = 0.2)	
parser.add_argument(
	'--qregression',
	help = 'When calculating MSE, should you weight errors by rho (quantile weights) (True) or consider +/- errors using the same weight (False)? [Default: True]',
	required = False,
	default = True)
parser.add_argument(
	'--tolThresh',
	help = 'When calculating tolerance ranges, the total summed values of both sides of a PDF when calcuating the HDI (Highest density interval). A value of 0.05 means each tail is 0.025. [Default: 0.05]',
	required = False,
	default = 0.05)		
parser.add_argument(
	'--tolThreshAbund',
	help = 'When calculating tolerance ranges, individual observations must be less than this percent of maximum relative abundance observation to be excluded from the HDI value. [Default: 0.05]',
	required = False,
	default = 0.05)	
	
args = parser.parse_args()

taxaTablePWD = args.taxaTable
metadataPWD = args.metadata
metadata_name = args.metadata_name
minX = args.minX
maxY = args.maxY
XYdiff = float(args.XYdiff)
unitSize = float(args.unitSize)
output_dir = args.output_dir

minCountTable = int(args.minCountTable)
minCountBin = int(args.minCountBin)
minCountOTUinSample = int(args.minCountOTUinSample)
minSamplePres = int(args.minSamplePres)
pvalthresh = float(args.pvalthresh)
qquant = float(args.qquant)
bandwidth_percent = float(args.bandwidth_percent)
qregression = float(args.qregression)
tolThresh = float(args.tolThresh)
tolThreshAbund = float(args.tolThreshAbund)


# ======================================== FUNCTIONS ==============================++++++++= #
#!bin/bash python3

class QTAGCommunity:
    name = ''
    taxaTablePWD = ''
    metadataPWD = ''
    taxaTable = ''
    taxaIDs = ''
    
    def __init__(self, taxaTablePWD, metadataPWD, output_dir):
        self.taxaTablePWD = taxaTablePWD
        self.metadataPWD = metadataPWD
        self.output_dir = output_dir
        self.func_mkdir(output_dir)
        
    ## Load table
    def func_makeTaxaSummaries(self, metadata_name, minCountOTUinSample, minCountTable, minSamplePres):
        # Input is file path for OTU table and metadata
        # Output is 
        # (a) Dictionary where keys are OTU IDs and content is two lists. (taxaTable)
        # 		The first list contains the relative abundance of that OTU across salinity
        # 		The second list contains the corresponding salinities in the first list
        # (b) Dictionary where keys are OTU IDs and content is taxonomy (taxaIDs)
        # 		Taxonomies are taken from the observation metadata from the OTU table
        print("Making and loading taxa table and metadata table...")
        taxaOpen = open(self.taxaTablePWD, 'r')
        taxaOpenTemp = []
        for i in taxaOpen:	# Read each line and concatenate into single file
            taxaOpenTemp += [i] # A single string; OTU table should be Unix(LF) format (\n at ends)
        taxaOpen.close()
        tempList =[] # to split lines
        for j in taxaOpenTemp: # Each line split by "\n"
            tempLine = j.strip()
            tempList += [tempLine.split('\t')] # Make a list of lists; each smaller list is abundance data for each OTU
        while '#OTU ID' not in tempList[0]:
            del tempList[0] # Deletes first row until the first row is #OTU ID
        # Sort information in relavent piles
        taxaIDs = {}
        taxaCountsTemp = {}
        first = True
        for y in tempList: # tempList is a list of lists; smaller list is OTUID,+ abundance data
            if y[0] == '#OTU ID': # The first line should be #OTU ID
                sites = y[1:len(y)-1] # So we will take the site names from the first row
            else: # Every other row is abundance data
                taxaIDs[y[0]] = y[len(y)-1] # Make dictionary of taxaIDs for later
                for x in range(len(y)): # Make file of 'total' counts for each OTU
                    if (x != 0) and (x != (len(y)-1)): # If it's not the first element (which is the OTU ID) or last element (which is the taxonomy), then add to list
                        if first: # First is made 'equal', but everything else is just appended.
                            taxaCountsTemp[str(x)] = float(y[x])
                        else:
                            taxaCountsTemp[str(x)] += float(y[x])
                first = False
            # Output of this loop is taxaCountsTemp (a dictionary of abundance data), taxaIDs, and sites
        headers = sites[:]
        taxaTable = {}
        for i in range(len(tempList)):
            taxaTable[tempList[i][0]] = [[]] # Make empty list of lists for each OTU ID
            if tempList[i][0] == '#OTU ID': # 
                pass
            else:
                for j in range(1,len(tempList[i])-1):
                    # Sum of all abundances to make relative abundances
                    value = float(tempList[i][j])
                    if value < float(minCountOTUinSample):
                        value = 0.0	
                    taxaTable[tempList[i][0]][0].append(value)
                    # Save values as relative abundances instead of absolute ones

        # Now get rid of low abundance in table
#         taxaTableFilt = self.deleteLowAbund(taxaTable,minCountTable,minSamplePres)
        taxaTableFilt = self.func_deleteLowAbund(taxaTable,minCountTable,minSamplePres)
        # Get total counts
        totalCounts = [ 0 for i in range(len(sites))]
        for i in range(len(sites)):
            for taxa in taxaTableFilt:
                totalCounts[i] += taxaTableFilt[taxa][0][i]
        # Convert to relative abundance
        taxaTableFinal =copy.deepcopy(taxaTableFilt)
        for OTU in taxaTableFilt.keys():
            tempOTUlist = [0 for i in range(len(sites))]
            for i in range(len(sites)):
                tempOTUlist[i] = float(taxaTableFilt[OTU][0][i])/float(totalCounts[i])
            taxaTableFinal[OTU][0] = tempOTUlist
        metadataOpen = open(self.metadataPWD, 'r')
        metadataOpenTemp = []
        for i in metadataOpen:
            metadataOpenTemp += [i]
        metadataOpen.close()
        tempMeta =[]
        for j in metadataOpenTemp:
            tempLine = j.strip()
            tempMeta += [tempLine.split('\t')]
        positionSal = tempMeta[0].index(metadata_name)
        metadata = []
        for line in tempMeta:
            metadata.append([line[0],line[positionSal]])
        # Now, change key names so they're not longer sites; they're gradient values
        for site in metadata:
            sites = [site[1] if x == site[0] else x for x in sites]
        # Make proper format, but with site names instead of numbers
        for taxa in taxaTableFinal:
            taxaTableFinal[taxa].append(sites)
        # Make abundance values integer as well
        for x in taxaTableFinal:
            taxaTableFinal[x] = [[float(strings) for strings in keys] for keys in taxaTableFinal[x]]
        self.func_printOTUTable(taxaTableFinal,headers,taxaIDs,output_dir)
        self.taxaTableRA = taxaTableFinal
        self.taxaTableCounts = taxaTableFilt
        self.taxaIDs = taxaIDs
    
    def QTAG(self, qtag_params, overwrite=False):
        
        # Check if overwrite, run fit models
        if (hasattr(self, 'allbestMSE')) & (not overwrite):
            print("Models already fit; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_fit_quantile_model(qtag_params)

        if (hasattr(self, 'classifications') & (not overwrite)):
            print("Classifications already determined; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_get_specialist_classifications(qtag_params)
        if (hasattr(self, 'toleranceRanges')) & (not overwrite):
            print("Tolerance already calculated; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_get_tolerance_range(qtag_params.tolThresh,qtag_params.tolThreshAbund)
        if (hasattr(self, 'classSummary') & (not overwrite)):
            print("Class Summary already calculated; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_get_class_summary()
        if ( os.path.isfile(self.output_dir+'/QTAG_output.txt') & (not overwrite) ):
            print("QTAG output exists; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_print_output()
        if ( os.path.isfile(self.output_dir+'/QTAG_quantile_curves.txt') & (not overwrite) ):
            print("QTAG quantile curves output exists; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_print_qCurves()
        if ( os.path.isfile(self.output_dir+'/QTAG_class_summary.txt') & (not overwrite) ):
            print("QTAG class summary output exists; not overwriting. Set overwrite=True if you wish to overwrite previous.")
        else:
            self.QTAG_print_classSummary()
            
    def QTAG_plot_all(self, overwrite=False):
        self.QTAG_plot_fit(overwrite=overwrite)
        self.QTAG_plot_class_summary(overwrite=overwrite)

    def QTAG_inherit(self, class_object, list_attributes):
        # To inherit previously calculated fits; to save time
        for att in list_attributes:
            setattr(self, att, getattr(class_object, att))
    
    def QTAG_fit_quantile_model(self, qtag_params):
        #==========================================
        # Calculate smoothed quantile values
        self.allMSE=dict()
        self.qCurves=dict()
        self.allListAbund=dict()
        self.allbestMSE=dict()
        
        # Extracting param values
        bandwidth_percent = qtag_params.bandwidth_percent
        qquant = qtag_params.qquant
        minX = qtag_params.minX
        maxY = qtag_params.maxY
        XYdiff = qtag_params.XYdiff
        unitSize = qtag_params.unitSize
        minCountBin = qtag_params.minCountBin
        
        print("Iterating through all combinations...")
        # Go through each taxa and iterate through all combinations of boundaries X and Y
        sleep(0.2)
        pbar = ProgressBar()
        for taxa in pbar(list(self.taxaTableRA.keys())):
            sleep(0.2)
            listAbund=self.taxaTableRA[taxa]
            qcurves=self.func_qcurve(minX, maxY, unitSize, bandwidth_percent,qquant,listAbund)
            MSE = numpy.array([0,0,0,0,0,0])
            for X in self.func_sequenceGenerator(minX,(maxY-XYdiff),unitSize):
                for Y in self.func_sequenceGenerator((X+XYdiff),maxY,unitSize):
                    # Original verison; raw data
                    binValues = self.func_sortBins(X,Y,listAbund)
                    binValuesQ = self.func_sortBinsQuantile(X,Y,qcurves)
                    if len(binValues['lista']) <= minCountBin or len(binValues['listb']) <= minCountBin or len(binValues['listc']) <= minCountBin: # if the bins have nothing in them, then don't use that bin combination
                        pass
                    else:
                        combined = self.func_makeAModelQ(X,Y,binValuesQ, binValuesQ)
                        modVal=combined['modelabbr']
                        MSE = numpy.vstack((MSE, numpy.array([X,Y,self.func_calcError(combined,qquant,qregression=qregression),modVal[0],modVal[1],modVal[2]])))
            # Delete first entry, which is 0,0,0
            MSET=numpy.delete(MSE,0,0).T

            self.allMSE[taxa] = MSET
            self.qCurves[taxa] = qcurves
            self.allListAbund[taxa] = listAbund
            # Order is: X, Y, Amodel,Bmodel,Cmodel
            self.allbestMSE[taxa] = MSET.T[numpy.where(MSET[2]==min(MSET[2]))[0].tolist()][0].tolist()
    
    def QTAG_get_specialist_classifications(self, qtag_params):
        print("Classifying taxa by specialist type")
        sleep(0.2)
        pbar = ProgressBar()
        classifications = dict()
        pvalthresh=qtag_params.pvalthresh
        for taxa in pbar(list(self.qCurves.keys())):
            sleep(0.2)
            currqcurves = self.qCurves[taxa]
            currmodel = self.allbestMSE[taxa]
            currBins = self.func_sortBinsQuantile(currmodel[0],currmodel[1],currqcurves)
            currListAbund = self.allListAbund[taxa]
            ## AB, BC, AC quantiles
            lista = currBins['lista']
            listb = currBins['listb']
            listc = currBins['listc']
            AQuant = self.func_average(currBins['lista'])
            BQuant = self.func_average(currBins['listb'])
            CQuant = self.func_average(currBins['listc'])
            # AB, BC, AC ranked test
            sig_ab = stats.mannwhitneyu(lista, listb)[1]
            sig_bc = stats.mannwhitneyu(listb, listc)[1]
            sig_ac = stats.mannwhitneyu(lista, listc)[1]
            # Correction using fdr for multiple comparisons
            sig_abadj, sig_bcadj, sig_acadj = list(multitest.fdrcorrection([sig_ab, sig_bc, sig_ac])[1])

            spclass = ''
            if (BQuant>AQuant) & (BQuant>CQuant) & (sig_abadj<pvalthresh) & (sig_bcadj<pvalthresh):
                spclass='intermediate'
            elif (AQuant>CQuant) & (sig_acadj<pvalthresh):
                spclass='low'
            elif (CQuant>AQuant) & (sig_acadj<pvalthresh):
                spclass='high'
            else:
                spclass='unclassified'
            classifications[taxa] = {'spClass':spclass,'fitModel':[AQuant, BQuant, CQuant],'pvals_ab_bc_ac':[sig_abadj, sig_bcadj, sig_acadj]}
        self.classifications = classifications
    
    def QTAG_get_tolerance_range(self, thresh, threshAbund):
        toleranceRanges=dict()
        for taxa in list(self.allListAbund.keys()):
            tempListAbund=self.allListAbund[taxa]
            # Sum total abundance
            sumRA=sum(tempListAbund[0])
            abundLim=max(tempListAbund[0])*threshAbund
            lwrRA=''
            lwrOccur=''
            uprRA=''
            uprOccur=''
            currRA=0
            for i in set(sorted(tempListAbund[1])):
                currAbund=[ tempListAbund[0][x[0]] for x in enumerate(tempListAbund[1]) if x[1]==i]
                # Need to fix instance where either list is empty
                prevAbund=[ tempListAbund[0][x[0]] for x in enumerate(tempListAbund[1]) if x[1]<i]+[0] # add 0 in case list is empty
                postAbund=[ tempListAbund[0][x[0]] for x in enumerate(tempListAbund[1]) if x[1]>i]+[0] # add 0 in case list is empty
                prevAbundMax = max(prevAbund)
                futureAbundMax = max(postAbund)
                currRA+=sum(currAbund)
                if lwrRA=='':
                    if ((currRA/sumRA) >= thresh/2) | (prevAbundMax>=abundLim):
                        lwrRA=i
                if uprRA=='':
                    if (currRA/sumRA >= 1-thresh/2) & (futureAbundMax<abundLim):
                        uprRA=i
                if lwrOccur=='':
                    if (currRA>0):
                        lwrOccur=i
                if uprOccur=='': 
                    if (currRA >= sumRA*0.99999): # needs to be not quite equal due to math errors in python
                        uprOccur=i
#                     elif (i==max(tempListAbund[1])):
#                         uprOccur=i
            toleranceRanges[taxa]=[lwrOccur, lwrRA, uprRA, uprOccur]
        self.toleranceRanges= toleranceRanges
    
    def QTAG_get_class_summary(self):
        lowSum = []
        interSum = []
        highSum = []
        unclassifiedSum=[]
        sals = []
        for taxa in list(self.taxaTableRA.keys()):
            #initilize
            if len(sals)==0:
                sals = self.taxaTableRA[taxa][1]
                sals_set = list(set(sals))
                sals_order = [sals_set.index(x) for x in sorted(sals_set)]
                lowSum = [0]*len(sals_set)
                interSum = [0]*len(sals_set)
                highSum = [0]*len(sals_set)
                unclassifiedSum=[0]*len(sals_set)
            sumRA=numpy.bincount(sals, weights=self.taxaTableRA[taxa][0])
            sumN=numpy.bincount(sals)
            sumTemp=[x[0]/x[1] if x[1]>0 else 0 for x in zip(sumRA, sumN)]
            # sum each taxa
            if self.classifications[taxa]['spClass']=='low':
                lowSum = [sum(x) for x in zip(lowSum, sumTemp)]
            elif self.classifications[taxa]['spClass']=='intermediate':
                interSum = [sum(x) for x in zip(interSum, sumTemp)]
            elif self.classifications[taxa]['spClass']=='high':
                highSum = [sum(x) for x in zip(highSum, sumTemp)]
            else:
                unclassifiedSum = [sum(x) for x in zip(unclassifiedSum, sumTemp)]

        sals_ordered=[sals_set[i] for i in sals_order]
        lowSum_ordered=[lowSum[i] for i in sals_order ]
        interSum_ordered=[interSum[i] for i in sals_order]
        highSum_ordered=[highSum[i] for i in sals_order]
        unclassifiedSum_ordered=[unclassifiedSum[i]for i in sals_order]
        classSummary=numpy.array([sals_ordered,lowSum_ordered, interSum_ordered, highSum_ordered,unclassifiedSum_ordered])
        self.classSummary=classSummary
        
    def QTAG_plot_fit(self, overwrite=False):
        if (os.path.isdir(self.output_dir + '/plots')) & (not overwrite):
            print("Plot directory already exists. To overwrite, include overwrite=True.")
        else:
            print("Plotting all fits")
            sleep(0.2)
            pbar = ProgressBar()
            self.func_mkdir(self.output_dir+'/plots')
            if hasattr(self, 'classifications') & hasattr(self, 'toleranceRanges'):
                for taxa in pbar(list(self.allListAbund.keys())):
                    sleep(0.2)
                    listAbund=self.allListAbund[taxa]
                    XYbestMSE=self.allbestMSE[taxa]
                    spclass=self.classifications[taxa]
                    qcurves=self.qCurves[taxa]
                    toleranceRanges=self.toleranceRanges[taxa]
                    underdist=-1*max(listAbund[0])*0.05
                    if (spclass['spClass']=='low'):
                        col='blue'
                    elif (spclass['spClass']=='intermediate'):
                        col='purple'
                    elif (spclass['spClass']=='high'):
                        col='red'
                    else:
                        col='black'
                    # Start plot
                    plt.scatter(listAbund[1],listAbund[0])
                    plt.ylim(top=max(listAbund[0])*1.1, bottom=-max(listAbund[0])*0.1)
                    plt.axhline(y=0, color="black", alpha=0.1)
                    plt.plot(qcurves[0], qcurves[3])
                    plt.plot([qtag_params.minX,XYbestMSE[0]], [spclass['fitModel'][0],spclass['fitModel'][0]], color="blue") # h line segment, low
                    plt.plot([XYbestMSE[0],XYbestMSE[1]], [spclass['fitModel'][1],spclass['fitModel'][1]], color="purple")# h line segment, mid
                    plt.plot([XYbestMSE[1],qtag_params.maxY], [spclass['fitModel'][2],spclass['fitModel'][2]], color="red")# h line segment, high
                    plt.plot([toleranceRanges[0],toleranceRanges[3]],[underdist,underdist], color=col, alpha=0.5)
                    plt.plot([toleranceRanges[1],toleranceRanges[2]],[underdist,underdist], color=col, alpha=1)
                    plt.axvline(x=XYbestMSE[0], color="black", linestyle="dashed") #X
                    plt.axvline(x=XYbestMSE[1], color="black", linestyle="dashed") #Y
                    plt.xlabel('Gradient value')
                    plt.ylabel('Relative abundance')
                    plt.title(taxa+": "+spclass['spClass'])
                    plt.subplots_adjust(left=0.2)
                    plt.savefig(self.output_dir + '/plots/'+taxa+'.png', transparent=False)
                    plt.clf()

            else:
                print("Error: QTAG not run yet, cannot plot best fit models.")
    
    def QTAG_plot_class_summary(self, overwrite=False):
        if (os.path.isfile(self.output_dir + '/class_summary_plot.png')) & (not overwrite):
            print("Class summary plot already exists.To overwrite, include overwrite=True.")
        else:
            if hasattr(self, 'classSummary'):
                sals=self.classSummary[0]
                low=self.classSummary[1]
                inter=self.classSummary[2]
                high=self.classSummary[3]
                noclass=self.classSummary[4]

                fig, ax = plt.subplots()
                ax.bar(sals, height=low, label='low', color='blue',bottom=0)
                ax.bar(sals, height=inter, label='inter', color='purple',bottom=low)
                ax.bar(sals, height=high, label='high', color='red',bottom=low+inter)
                ax.bar(sals, height=noclass, label='unclassified', color='grey',bottom=low+inter+high)
                ax.set_ylabel('Relative Abundance')
                ax.set_title('Class composition across gradient')
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.set_xlabel("Gradient Value")
                plt.subplots_adjust(right=0.75)
                plt.savefig(self.output_dir + '/class_summary_plot.png', transparent=False)
                plt.clf()
            else:
                print("Error: Class summary attribute not created yet.")
                
    
    def QTAG_print_output(self):
        modelsToPrint='TaxaID\tclass\tX\tY\tMSE\tmodelA\tmodelB\tmodelC\tpval_ab\tpval_bc\tpval_ac\tlwrTol\tlwr95\tupr95\tuprTol\n'
        for taxa in list(self.taxaTableRA.keys()):
            modelsToPrint+= taxa + '\t' + self.classifications[taxa]['spClass']
            for i in self.allbestMSE[taxa]:
                modelsToPrint+='\t'+str(i)
            for j in self.classifications[taxa]['pvals_ab_bc_ac']:
                modelsToPrint+='\t'+str(j)
            for k in self.toleranceRanges[taxa]:
                modelsToPrint+='\t'+str(k)
            modelsToPrint+='\n'
        modelsToPrint=modelsToPrint.strip()
        QtagOutput = open(self.output_dir + '/QTAG_output.txt', 'w')
        QtagOutput.write(modelsToPrint)
        QtagOutput.close()
    
    def QTAG_print_qCurves(self):
        sals=[]
        toPrintQCurves='Gradient'
        for taxa in list(self.qCurves.keys()):
            # If first loop, save sals
            if len(sals)==0:
                sals = self.qCurves[taxa][0]
                qcurveMat=numpy.array(sals)
            # Print all taxa headers
            toPrintQCurves+='\t'+taxa
            qcurveMat=numpy.vstack([qcurveMat,self.qCurves[taxa][3]])
        toPrintQCurves+='\n' #Add end to header line
        for r in range(qcurveMat.shape[1]): # Loop through array and print table
            for c in list(qcurveMat[:,r]):
                toPrintQCurves+=str(c)+'\t'
            toPrintQCurves=toPrintQCurves.strip()
            toPrintQCurves+='\n'
        toPrintQCurves=toPrintQCurves.strip()
        qCurvePrint = open(self.output_dir + '/QTAG_quantile_curves.txt','w')
        qCurvePrint.write(toPrintQCurves)
        qCurvePrint.close()
    
    def QTAG_print_classSummary(self):
        toPrintSummary='Gradient\tlow\tinter\thigh\tunclassified\n'
        for r in range(self.classSummary.shape[1]):
            for c in list(self.classSummary[:,r]):
                toPrintSummary+= str(c)+'\t'
            toPrintSummary = toPrintSummary.strip() + '\n'
        classSummaryToPrint = open(self.output_dir + '/QTAG_class_summary.txt','w')
        classSummaryToPrint.write(toPrintSummary)
        classSummaryToPrint.close()
        
        
    #### INTERNALLY CALLED FUNCITONS ####
    ### Delete low abund taxa
    def func_deleteLowAbund(self, taxaTable, minCountTable,minSamplePres): # Change to absolute threshold, or get rid of entirely. Get rid of show ups in < 3 samples
        taxaTableFilt = {}
        for taxa in taxaTable:
            nonzeroCount = 0
            presCount = 0
            for i in taxaTable[taxa][0]:
                nonzeroCount += i
                if i > 0:
                    presCount += 1
            if int(nonzeroCount) >= int(minCountTable) and int(presCount) >= int(minSamplePres):
                taxaTableFilt[taxa] = taxaTable[taxa]
        return(taxaTableFilt)

    ### Save OTU table
    def func_printOTUTable(self,taxaTableFinal,headers,taxaIDs,output_dir):
        OTUtabletoPrint = open(output_dir+'/OTUTableText.txt','w')
        toPrint = '#OTUID'
        for head in headers:
            toPrint += '\t' + head
        toPrint += '\t'+'taxonomy'+'\n'
        for taxa in taxaTableFinal:
            toPrint += taxa
            for val in taxaTableFinal[taxa][0]:
                toPrint += '\t' + str(val)
            toPrint += '\t' + taxaIDs[taxa]
            toPrint += '\n'
        OTUtabletoPrint.write(toPrint)
        OTUtabletoPrint.close()
        print("Printed OTU table")
    
    # Fit a quantile curve
    def func_qcurve(self, minX, maxY, unitSize, bandwidth_percent, qquant, listAbund):
        # Get % band width to nearest unit size
        bandwidth = round(((float(maxY)-float(minX))*float(bandwidth_percent))/float(unitSize))*float(unitSize)
        smooth=numpy.array([0,0,0,0])
        for U in self.func_sequenceGenerator(minX,maxY,unitSize):
            minU=max(U-bandwidth/2,0)
            maxU=min(U+bandwidth/2,maxY)
            sals=numpy.array(listAbund[1])
            keep=numpy.where((sals>=minU) & (sals<=maxU))[0].tolist()
            if len(keep)>0:
                q=numpy.quantile(a=[listAbund[0][i] for i in keep],q=qquant,interpolation='linear')
                smooth=numpy.vstack((smooth,numpy.array([U,minU,maxU,q])))
        smoothT=numpy.delete(smooth, 0,0).T
        return(smoothT)
    
    def func_sequenceGenerator(self, min,max,binSize): # Make sequence from max value, min value, and binsize. List includes min and max bin values. Integers.
        if max < min:
            print("ERROR, max is less than min when trying to generate sequence in sequenceGenerator")
            return
        breadth = max-min
        nbins = math.ceil(float(breadth)/float(binSize))
        nbins = int(nbins) # Make into integer from float
        seq = []
        for i in range(nbins+1):
            current = float(min) + float(binSize)*i
            seq.append(current)
        return seq

    def func_sortBins(self, X, Y, listAbund): # Sorts list of abundances into bins according to X and Y and their associated salinities
        lista = []
        listb = []
        listc = []
        for i in range(len(listAbund[1])):
            if listAbund[1][i] < X:  # if X = Y, value will be sorted into X
                lista.append(listAbund[0][i])
            if listAbund[1][i] >= X and listAbund[1][i] <= Y:
                listb.append(listAbund[0][i])
            if listAbund[1][i] > Y:
                listc.append(listAbund[0][i])
        binValues = {'lista':lista, 'listb':listb, 'listc':listc}
        return binValues # output is a dictionary with abundance observations sorted into lists according to what salinity they were found in

    def func_sortBinsQuantile(self, X, Y, arr): # Sorts list of abundances into bins according to X and Y and their associated salinities
        lista = [arr[3][i] for i in numpy.where((arr[0] <= X))[0].tolist()]
        listb = [arr[3][i] for i in numpy.where((arr[0] > X) & (arr[0] <Y))[0].tolist()]
        listc = [arr[3][i] for i in numpy.where((arr[0] >= Y))[0].tolist()]
        binValues = {'lista':lista, 'listb':listb, 'listc':listc}
        return binValues # output is a dictionary with abundance observations sorted into lists according to what salinity they were found in

    def func_makeAModelQ(self, X,Y, binValues_data, binValues_model_fit): # output is a dictionary with the abundances at each salinity, and then the 'model', which is just means
        aModel = self.func_average(binValues_model_fit['lista'])
        bModel = self.func_average(binValues_model_fit['listb'])
        cModel = self.func_average(binValues_model_fit['listc'])
        model = [aModel]*len(binValues_data['lista']) + [bModel]*len(binValues_data['listb']) + [cModel]*len(binValues_data['listc'])
        # Make a combined 'piecewise' list that has each value and its corresponding model value
        abundance = binValues_data['lista'] + binValues_data['listb'] + binValues_data['listc']
        combined = {'abundance': abundance, 'model': model, 'modelabbr': [aModel,bModel,cModel]}
        return combined

    def func_calcError(self, combined,qquant,qregression): # Calculated the MES of data vs model
        errorSum = 0
        n = range(len(combined['abundance']))
        for i in n:
            error = (combined['model'][i] - combined['abundance'][i])
            if qregression!=True :
                indivError = error**2
            else:
                rhoi = float(qquant)*max(error,0) + (1-float(qquant))*max(-error,0)
                indivError = abs(error)*rhoi
            errorSum += indivError
            err = float(errorSum)/len(combined['abundance'])
        return err

    def func_mkdir(self, newdir):
    #     directory = os.path.dirname(file_path)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
    
    def func_average(self,listValues): # shortcut for 'average'
        # finds average for list of numbers
        if type(listValues) in [int,float]:
            listValues = [listValues]
        if len(listValues) == 0:
            return None
        else:
            average = float(sum(listValues))/len(listValues)
            return average

            
class QTAG_parameters:
    unitSize = ''
    XYdiff = ''
    bandwidth_percent = ''
    qquant = ''
    minCountBin = ''
    pvalthresh = ''
    tolThresh = ''
    tolThreshAbund = ''
    minX = ''
    maxY = ''
    
    ##### Set up ######
    def __init__(self, unitSize, XYdiff, bandwidth_percent, qquant, minCountBin, pvalthresh, tolThresh, tolThreshAbund):
        self.unitSize = unitSize
        self.XYdiff = XYdiff
        self.bandwidth_percent = bandwidth_percent
        self.qquant = qquant
        self.minCountBin = minCountBin
        self.pvalthresh = pvalthresh
        self.tolThresh = tolThresh
        self.tolThreshAbund = tolThreshAbund

    def func_set_minmax(self, taxaTableRA, minX, maxY):
        # Find min/max gradient values
        maxVal = max(taxaTableRA[list(taxaTableRA.keys())[1]][1])
        minVal = min(taxaTableRA[list(taxaTableRA.keys())[1]][1])

        # Assign min/max gradient values
        if str(minX) == 'Check':
            minX = math.floor(minVal)
        else:
            minX = minX

        if str(maxY) == 'Check':
            maxY = math.ceil(maxVal)
        else:
            maxY = maxY
        self.minX = minX
        self.maxY = maxY
    
    def print_params(self, output_dir):
        toPrint='parameter\tvalue\n'
        for a in ['unitSize', 'XYdiff', 'bandwidth_percent', 'qquant', 'minCountBin', 'pvalthresh', 'tolThresh', 'tolThreshAbund', 'minX', 'maxY' ]:
            toPrint+= a + '\t' + str(getattr(self, a)) + '\n'
        toPrint=toPrint.strip()
        paramPrint = open(output_dir+'/QTAG_parameters.txt', 'w')
        paramPrint.write(toPrint)
        paramPrint.close()




# ========================= RUN ===========================#


# LOAD #

community = QTAGCommunity(taxaTablePWD, metadataPWD, output_dir)
# Make community file
community.func_makeTaxaSummaries(metadata_name, minCountOTUinSample, minCountTable, minSamplePres)
# Set parameters, include new file
qtag_params = QTAG_parameters(unitSize, XYdiff, bandwidth_percent, qquant, minCountBin, pvalthresh, tolThresh, tolThreshAbund)
qtag_params.func_set_minmax(community.taxaTableRA, minX, maxY)
# Run QTAG
community.QTAG(qtag_params)
# Plot QTAG
community.QTAG_plot_all(qtag_params)
# Print parameters
qtag_params.print_params(output_dir)



