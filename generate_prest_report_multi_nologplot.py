import sys
import csv
import os
from cycler import cycler
import matplotlib.pyplot as plt
from mpltools import color
import numpy as np

def main(argv):
  if len(argv) >= 2 and argv[0].upper() in ['A', 'B', 'AB']:
    protFiles = argv[1:]
    vial = argv[0]
  else:
    sys.exit("\nusage: python generate_prest_report_multi.py <vial> <results_file1> (<results_file2> ...)\n  where vial is one of 'A', 'B' or 'AB'.")
  
  color_cycle = color.colors_from_cmap(len(protFiles), 'tab10')
  plt.rc('axes', prop_cycle= cycler('color',color_cycle.tolist()))
  #plt.rcParams['axes.color_cycle'] = 

  p = plt.figure(figsize=(16, 12))
  i = 0
  for protFile in protFiles:
    plotPrestReport(protFile, vial.upper(), p)
    i += 1
  handles, labels = p.get_axes()[-1].get_legend_handles_labels()
  p.legend(handles, labels,bbox_to_anchor=(0.5, 0., 0.45, 0.45), shadow=True, fontsize='x-large',  loc='upper left')
  #plt.show()
  plt.savefig("result.svg")
    
def plotPrestReport(protFile, vial, p):
  vials = ["A","B","AB"]
  
  poolAFastaFile = "prest_pool_a.fasta"
  poolBFastaFile = "prest_pool_b.fasta"
  poolRandomFastaFile = "prest_1000_random.fasta"
  
  poolAProteins = list(getProteinIds(open(poolAFastaFile,'r')))
  poolBProteins = list(getProteinIds(open(poolBFastaFile,'r')))
  poolRandomProteins = list(getProteinIds(open(poolRandomFastaFile,'r')))
  totalDbProteins = len(poolAProteins) + len(poolBProteins) + len(poolRandomProteins)
  
  if vial == vials[0]:
    presentProteins = poolAProteins
  elif vial == vials[1]:
    presentProteins = poolBProteins
  elif vial == vials[2]:
    presentProteins = poolAProteins + poolBProteins

  scorename = "Method"
  if not os.path.isfile(protFile) or os.path.getsize(protFile) <= 0:
    print("Warning: " + protFile + " empty or does not exist. Skipping..")
    for i in range(1,4):
      plt.subplot(2, 2, i)
      plt.plot([],[])
    return
   
  with open(protFile) as f:
    scorename = f.readline().split('\t')[1].rstrip('\n').strip('\"')
  f.close()
  reportedQvals, entrapmentFDRs, tpfp = getQvalues(protFile, totalDbProteins, presentProteins)
  numSignificant = sum(1 if qval < 0.05 else 0 for qval in reportedQvals)
  print(scorename + ":")
  print(("#significant protein groups Reported 5% FDR:"), (numSignificant))
  print(("#significant protein groups Observed 5% FDR:"), (sum(1 if qval < 0.05 else 0 for qval in entrapmentFDRs)))
  print(("Observed FDR at reported 5% FDR:"), (entrapmentFDRs[numSignificant-1]))
  numRows, numCols = 2, 2

  csv = open("reportedEntrapmentFDRs_"+scorename+".csv", "w")
  for rep, trap, tf in zip(reportedQvals,entrapmentFDRs,tpfp):
    csv.write(str(rep) + "," + str(trap) + "," + str(tf[0]) + "," + str(tf[1]) + "\n")

  csv.close()
  plotQvalues(reportedQvals, entrapmentFDRs, tpfp, color, '*', numRows, numCols, scorename)
  
  upperMargin = 1.5
  lowerMargin = 1.0/upperMargin
  
  xlabel = "Decoy FDR"
  ylabel = "Entrapment FDR"
  labelFontSize = 24
  titleFontSize = 28
  axisFontSize = 20
  
  x = np.linspace(1e-20, 1, num=1000)
  
  figIdx = 1
  
  plt.subplot(numRows, numCols, figIdx)
  plt.title("Calibration (all)", fontsize = titleFontSize, fontweight = 'bold', y = 1.02)
  plt.axis([0, 1, 0, 1])
  plt.plot(x, x, 'k-')
  plt.plot(x, [a*upperMargin for a in x], 'k--')
  plt.plot(x, [a*lowerMargin for a in x], 'k--')
  plt.xlabel(xlabel, fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  setAxisFontSize(axisFontSize)
  
  figIdx += 1
  
  plt.subplot(numRows, numCols, figIdx)
  plt.title("Calibration (zoom)", fontsize = titleFontSize, fontweight = 'bold', y = 1.02)
  plt.plot(x, x, 'k-')
  plt.plot(x, [a*upperMargin for a in x], 'k--')
  plt.plot(x, [a*lowerMargin for a in x], 'k--')
  plt.axis([0, 0.1, 0, 0.1])
  plt.xlabel(xlabel, fontsize = labelFontSize)
  plt.ylabel(ylabel, fontsize = labelFontSize)
  setAxisFontSize(axisFontSize)
  
  figIdx += 1
  
  maxY = (np.floor(len(presentProteins) / 100) + 1) * 100
  plt.subplot(numRows, numCols, figIdx)
  plt.title("Performance", fontsize = titleFontSize, fontweight = 'bold', y = 1.02)
  plt.plot([0.05, 0.05], [0, maxY], 'k', linestyle = 'dotted')
  plt.plot([0, 0.1], [len(presentProteins), len(presentProteins)], 'k', linestyle = 'dashed')
  plt.xlim([0, 0.1])
  plt.ylim([0, maxY])
  plt.xlabel(ylabel, fontsize = labelFontSize)
  plt.ylabel("Number of protein groups", fontsize = labelFontSize)
  setAxisFontSize(axisFontSize)
  
  plt.subplots_adjust(top=0.94, left=0.08, right=0.96, bottom=0.08, hspace=0.3, wspace=0.3)
  #handles, labels = ax.get_legend_handles_labels()

def setAxisFontSize(size):
  for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(size)
  for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(size)
    
def getQvalues(fileName, totalProteins, presentProteins):
  # The input is a file starting with a header line, followed by lines: protein1 <tab> q-value1
  #   protein groups can be specified as "protein1a,protein1b" (without quotes)
  reader = csv.reader(open(fileName, 'r'), delimiter='\t')
  next(reader) # skip header
  x = np.linspace(0, 1, num=1000)
  fp = 0
  tp = 0

  reportedQvals, entrapmentFDRs, tpfp = [], [], []

  fpSeen = False
  pSeen = False
  
  numAbsentProteins = totalProteins - len(presentProteins)
  piA = 1.0 - float(len(presentProteins)) / totalProteins

  auc = 0
  thresh = 0.05
  maxarea = 0
  oldQval = -1
  for row in reader:
    reportedQvalue = float(row[1])
    
    proteinName = row[0].strip('\"') # Protein name (unquoted)
    
    for proteinId in proteinName.split(","):
      proteinAbsent = proteinId not in presentProteins
      if proteinAbsent:
        fpSeen = True
      else:
        pSeen = True
    
    if pSeen:
      tp = tp + 1
    elif fpSeen:
      fp = fp + 1
      if fp/(tp+fp) <= thresh:
        auc += tp
    if fp/(tp+fp) <= thresh:
      maxarea = (len(presentProteins)**2 * thresh)

    fpSeen = False
    pSeen = False

    if False:
      entrapmentFdr = float(fp+1) / (tp + fp + 1)
    else:
      entrapmentFdr = float(fp) / (tp + fp)

    tpfp.append([tp,fp])
    reportedQvals.append(piA * reportedQvalue)
    entrapmentFDRs.append(entrapmentFdr)
  
  entrapmentFDRs = fdrsToQvals(entrapmentFDRs)
  if maxarea == 0:
    rocx = 0
  else:
    rocx = float(auc)/float(maxarea)

  print (("AUC at 0.05:"),(rocx))

  return reportedQvals, entrapmentFDRs, tpfp

def fdrsToQvals(fdrs):
  qvals = [0] * len(fdrs)
  qvals[len(fdrs)-1] = fdrs[-1]
  for i in range(len(fdrs)-2, -1, -1):
    qvals[i] = min(qvals[i+1], fdrs[i])
  return qvals
  
def plotQvalues(reportedQvals, entrapmentFDRs, tpfp, color, marker, numRows, numCols, label = ""):    
  figIdx = 1
  
  plt.subplot(numRows, numCols, figIdx)
  plt.plot(reportedQvals, entrapmentFDRs, linewidth = 2,marker = marker, markevery = 0.05, markeredgecolor = 'none', label = label)
  
  figIdx += 1
  
  plt.subplot(numRows, numCols, figIdx)
  plt.plot(reportedQvals, entrapmentFDRs, linewidth = 2,marker = marker, markevery = 0.05, markeredgecolor = 'none', label = label)
  
  figIdx += 1
  
  plt.subplot(numRows, numCols, figIdx)
  plt.plot(entrapmentFDRs, [x[0]+x[1] for x in tpfp], linewidth = 2, marker = marker, markevery = 0.05, markeredgecolor = 'none', label = label)

def getProteinIds(fp):
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield name
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name: yield name
     
if __name__ == "__main__":
  main(sys.argv[1:])

