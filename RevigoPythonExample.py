#!/usr/bin/env python
# -*- coding: utf-8 -*-

# An example of using a REVIGO core library for your own projects.
# To Run this example you need RevigoCore library and a 
# set of database files available at: http://revigo.irb.hr/RevigoDatabases.zip
# 
# Authors:
#		Rajko Horvat (rhorvat at irb.hr)
#	
# License:
# 	MIT License
#		Copyright(c) 2011-2023 Ruðer Boškoviæ Institute
#		
#		Permission is hereby granted, free of charge, to any person obtaining a copy
#		of this software and associated documentation files (the "Software"), to deal
#		in the Software without restriction, including without limitation the rights
#		to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#		copies of the Software, and to permit persons to whom the Software is
#		furnished to do so, subject to the following conditions:
#
#		The above copyright notice and this permission notice shall be included in all
#		copies or substantial portions of the Software.
#
#		THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#		IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#		FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#		AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#		LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#		OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#		SOFTWARE.

import time
import math

from clr_loader import get_coreclr
from pythonnet import set_runtime

set_runtime(get_coreclr(runtime_config="PythonRuntimeConfig.json"))

import clr

clr.AddReference("RevigoCore")
from IRB.Revigo.Core.Worker import RevigoWorker, ValueTypeEnum, RequestSourceEnum
from IRB.Revigo.Core import SemanticSimilarityTypeEnum, RevigoTerm, RevigoTermCollection, Utilities
from IRB.Revigo.Core.Databases import GeneOntology, SpeciesAnnotationList

clr.AddReference("mscorlib")
from System.IO import StreamWriter

def main():
	dCutoff = 0.7
	eValueType = ValueTypeEnum.PValue
	iSpeciesTaxon = 0
	eMeasure = SemanticSimilarityTypeEnum.SIMREL
	bRemoveObsolete = True
	print("Loading Ontology")
	oOntology = GeneOntology.Deserialize("C:\\Revigo\\Databases\\Current\\GeneOntology.xml.gz")
	print("Loading Species Annotations")
	oAnnotations = SpeciesAnnotationList.Deserialize("C:\\Revigo\\Databases\\Current\\SpeciesAnnotations.xml.gz")
	sExample1 = None
	sExample2 = None
	sExample3 = None

	with open('Example1.csv', 'r') as file:
		sExample1 = file.read()

	with open('Example2.csv', 'r') as file:
		sExample2 = file.read()

	with open('Example3.csv', 'r') as file:
		sExample3 = file.read()

	# Create worker 1
	oWorker1 = RevigoWorker(
		# JobID
		1, 
		# Ontology
		oOntology, 
		# Annotations for a given dataset
		oAnnotations.GetByID(iSpeciesTaxon), 
		# Timeout in minutes
		20, 
		# Job source
		RequestSourceEnum.JubSubmitting,
		# Dataset
		sExample1, 
		# Job parameters
		dCutoff, eValueType, eMeasure, bRemoveObsolete)

	# Create worker 2
	oWorker2 = RevigoWorker(2, oOntology, oAnnotations.GetByID(9606), 20, RequestSourceEnum.JubSubmitting,
		sExample2, 0.9, eValueType, SemanticSimilarityTypeEnum.LIN, bRemoveObsolete);

	# Create worker 3
	oWorker3 = RevigoWorker(3, oOntology, oAnnotations.GetByID(iSpeciesTaxon), 20, RequestSourceEnum.JubSubmitting,
		sExample3, 0.4, eValueType, eMeasure, bRemoveObsolete);

	# Workers will notify when the are finished processing the data
	oWorker1.OnFinish += OWorker_OnFinish
	oWorker2.OnFinish += OWorker_OnFinish
	oWorker3.OnFinish += OWorker_OnFinish

	# Start Workers and wait for their completion
	# They will automatically be assigned to different CPU core, if available
	print("Starting Workers...")
	oWorker1.Start()
	oWorker2.Start()
	oWorker3.Start()

	while oWorker1.IsFinished != True or oWorker2.IsFinished != True or oWorker3.IsFinished != True:
		time.sleep(0.1)

	print("All Workers have finished processing.")

	# export our results
	print("Exporting data.")
	ExportTable(oOntology, oWorker1, oWorker1.BPVisualizer, "..\\Example1_BPTable.tsv")
	ExportScatterplot(oOntology, oWorker2, oWorker2.CCVisualizer, "..\\Example2_CCScatterplot.tsv")
	ExportTreeMap(oOntology, oWorker3, oWorker3.MFVisualizer, "..\\Example3_MFTreeMap.tsv")
	ExportCytoscapeXGMML(oWorker1.BPVisualizer, "..\\Example1_BPCytoscape.xgmml")
	ExportSimMat(oWorker1.BPVisualizer, "..\\Example1_BPSimilarityMatrix.tsv")
	ExportWordClouds(oWorker1, "..\\Example1_WordClouds.json")

	# We are finished

def ExportTable(ontology, worker, visualizer, fileName):
	if visualizer.IsEmpty != True:
		with open(fileName, 'w') as oWriter:
			oTerms = visualizer.Terms.FindClustersAndSortByThem(ontology, worker.CutOff)

			oWriter.write("TermID\tName\tValue\t")
			c = 1
			while c < worker.MinNumColsPerGoTerm:
				oWriter.write("UserValue_{0}\t", c - 1)
				c += 1
			
			oWriter.write("LogSize\tFrequency\tUniqueness\tDispensability\tRepresentative\n")

			# print the data
			i = 0
			while i < oTerms.Count:
				term = oTerms[i]

				oWriter.write("\"{}\"\t".format(term.GOTerm.FormattedID))
				oWriter.write("\"{}\"\t".format(term.GOTerm.Name))
				oWriter.write("{}\t".format(term.Value))

				c = 1
				while c < worker.MinNumColsPerGoTerm:
					oWriter.write("{}\t".format(term.UserValues[c - 1]))
					c += 1

				oWriter.write("{}\t".format(term.LogAnnotationSize))
				oWriter.write("{}\t".format(term.AnnotationFrequency * 100.0))
				oWriter.write("{}\t".format(term.Uniqueness))
				oWriter.write("{}\t".format(term.Dispensability))

				if term.RepresentativeID > 0:
					oWriter.write("{}".format(term.RepresentativeID));
				else:
					oWriter.write("null");

				oWriter.write("\n")
				i += 1

def ExportScatterplot(ontology, worker, visualizer, fileName):
	if visualizer.IsEmpty != True:
		with open(fileName, 'w') as oWriter:
			oTerms = visualizer.Terms.FindClustersAndSortByThem(ontology, worker.CutOff)

			oWriter.write("TermID\tName\tValue\tLogSize\tFrequency\tUniqueness\tDispensability\tPC_0\tPC_1\tRepresentative\n")

			# print the data
			i = 0
			while i < oTerms.Count:
				term = oTerms[i]

				oWriter.write("\"{}\"\t".format(term.GOTerm.FormattedID))
				oWriter.write("\"{}\"\t".format(term.GOTerm.Name))
				oWriter.write("{}\t".format(term.Value))

				oWriter.write("{}\t".format(term.LogAnnotationSize))
				oWriter.write("{}\t".format(term.AnnotationFrequency * 100.0))
				oWriter.write("{}\t".format(term.Uniqueness))
				oWriter.write("{}\t".format(term.Dispensability))

				# 2D
				oWriter.write("{}\t".format(term.PC[0] if (term.PC.Count > 0) else "null"))
				oWriter.write("{}\t".format(term.PC[1] if (term.PC.Count > 1) else "null"))

				oWriter.write("{}".format(term.RepresentativeID if (term.RepresentativeID > 0) else "null"))

				oWriter.write("\n")
				i += 1

def ExportTreeMap(ontology, worker, visualizer, fileName):
	if visualizer.IsEmpty != True:
		with open(fileName, 'w') as oWriter:
			terms = visualizer.Terms.FindClustersAndSortByThem(ontology, 0.1)

			oWriter.write("# WARNING - This exported Revigo data is only useful for the specific purpose of constructing a TreeMap visualization.\n")
			oWriter.write("# Do not use this table as a general list of non-redundant GO categories, as it sets an extremely permissive \n")
			oWriter.write("# threshold to detect redundancies (c=0.10) and fill the 'representative' column, while normally c>=0.4 is recommended.\n")
			oWriter.write("# To export a reduced-redundancy set of GO terms, go to the Scatterplot or Table tab, and export from there.\n")

			oWriter.write("TermID\tName\tFrequency\tValue\t")
			c = 1
			while c < worker.MinNumColsPerGoTerm:
				oWriter.write("UserValue_{}\t".format(c - 1))
				c += 1
			oWriter.write("Uniqueness\tDispensability\tRepresentative\n")

			# print the data
			i = 0
			while i < terms.Count:
				term = terms[i]
				
				isTermEliminated = term.Dispensability > worker.CutOff
				if (isTermEliminated):
					i += 1
					continue # will not output terms below the dispensability threshold at all

				oWriter.write("\"{}\"\t".format(term.GOTerm.FormattedID))
				oWriter.write("\"{}\"\t".format(term.GOTerm.Name))
				oWriter.write("{}\t".format(term.AnnotationFrequency * 100.0))

				oWriter.write("{}\t".format(term.Value))

				c = 1
				while c < worker.MinNumColsPerGoTerm:
					oWriter.write("{}\t".format(term.UserValues[c - 1]))
					c += 1

				oWriter.write("{}\t".format(term.Uniqueness))
				oWriter.write("{}\t".format(term.Dispensability))
				if term.RepresentativeID > 0:
					oWriter.write("\"{}\"".format(ontology.GetValueByKey(term.RepresentativeID).Name))
				else:
					oWriter.write("null")

				oWriter.write("\n")
				i += 1

def ExportCytoscapeXGMML(visualizer, fileName):
	if visualizer.IsEmpty != True:
		oWriter = StreamWriter(fileName)
		visualizer.SimpleOntologram.GraphToXGMML(oWriter)
		oWriter.Close()

def ExportSimMat(visualizer, fileName):
	if visualizer.IsEmpty != True:
		with open(fileName, 'w') as oWriter:
			i = 0
			while i < visualizer.Terms.Length:
				oWriter.write("\t{}".format(visualizer.Terms[i].GOTerm.FormattedID))
				i += 1

			oWriter.write("\n")
			i = 0
			while i < visualizer.Terms.Length:
				oWriter.write(visualizer.Terms[i].GOTerm.FormattedID)
				j = 0
				while j < visualizer.Terms.Length:
					oWriter.write("\t{}".format(visualizer.Matrix.GetSimilarity(i, j)))
					j += 1

				oWriter.write("\n")
				i += 1

def ExportWordClouds(worker, fileName):
	with open(fileName, 'w') as oWriter:
		oWriter.write("{")
		if worker.Enrichments.Count > 0:
			oWriter.write("\"Enrichments\":[")

			MIN_UNIT_SIZE = 1.0
			MAX_UNIT_SIZE = 9.0
			RANGE_UNIT_SIZE = MAX_UNIT_SIZE - MIN_UNIT_SIZE
			minFreq = 999999.0
			maxFreq = 0.0

			i = 0
			while i < worker.Enrichments.Count:
				dFrequency = math.sqrt(worker.Enrichments[i].Value)
				if dFrequency > 0.0:
					minFreq = min(minFreq, dFrequency)
					maxFreq = max(maxFreq, dFrequency)
				i += 1
			

			if minFreq > maxFreq:
				dTemp = minFreq
				minFreq = maxFreq
				maxFreq = dTemp
			
			if minFreq == maxFreq:
				maxFreq += 1
			
			range = maxFreq - minFreq
			bFirst = True

			i = 0
			while i < worker.Enrichments.Count:
				sWord = worker.Enrichments[i].Key.replace("'", "")
				dFrequency = math.sqrt(worker.Enrichments[i].Value)

				if dFrequency > 0.0:
					if bFirst == False:
						oWriter.write(",")

					size = math.ceil(MIN_UNIT_SIZE + round(((dFrequency - minFreq) * RANGE_UNIT_SIZE) / range))
					oWriter.write("{{\"Word\":\"{0}\",\"Size\":{1}}}".format(Utilities.StringToJSON(sWord), size))
					bFirst = False
				i += 1
			
			oWriter.write("]")

		if worker.Correlations.Count > 0:
			if worker.Enrichments.Count > 0:
				oWriter.write(",")

			oWriter.write("\"Correlations\":[")

			MIN_UNIT_SIZE = 1.0
			MAX_UNIT_SIZE = 9.0
			RANGE_UNIT_SIZE = MAX_UNIT_SIZE - MIN_UNIT_SIZE
			minFreq = 999999.0
			maxFreq = 0.0
			i = 0
			while i < worker.Correlations.Count:
				dFrequency = worker.Correlations[i].Value;
				if dFrequency > 0.0:
					minFreq = min(minFreq, dFrequency)
					maxFreq = max(maxFreq, dFrequency)
				i += 1

			if minFreq > maxFreq:
				dTemp = minFreq
				minFreq = maxFreq
				maxFreq = dTemp
			
			if minFreq == maxFreq:
				maxFreq += 1
			
			range = maxFreq - minFreq
			bFirst = True

			i = 0
			while i < worker.Correlations.Count:
				sWord = worker.Correlations[i].Key.replace("'", "")
				dFrequency = worker.Correlations[i].Value

				if dFrequency > 0.0:
					if bFirst == False:
						oWriter.write(",")

					size = math.ceil(MIN_UNIT_SIZE + round(((dFrequency - minFreq) * RANGE_UNIT_SIZE) / range))
					oWriter.write("{{\"Word\":\"{0}\",\"Size\":{1}}}".format(Utilities.StringToJSON(sWord), size))
					bFirst = False
				i += 1
			
			oWriter.write("]");
		
		oWriter.write("}");

def OWorker_OnFinish(sender, e):
	if sender != None:
		print("Worker {} has finished processing the data in {} seconds.".format(sender.JobID, sender.ExecutingTime.TotalSeconds))

if __name__ == '__main__':
	main()
