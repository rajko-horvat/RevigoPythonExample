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
import clr

clr.AddReference("RevigoCore")
from IRB.Revigo.Worker import RevigoWorker, ProgressEventArgs, ValueTypeEnum, RequestSourceEnum
from IRB.Revigo.Core import SemanticSimilarityScoreEnum, GOTermList
from IRB.Revigo.Databases import SpeciesAnnotationsList, GeneOntology

def main():
	dCutoff = 0.7
	eValueType = ValueTypeEnum.PValue
	iSpeciesTaxon = 0
	eMeasure = SemanticSimilarityScoreEnum.SIMREL
	bRemoveObsolete = True
	print("Loading Ontology")
	oOntology = GeneOntology("go.obo")
	print("Loading Species Annotations")
	oAnnotations = SpeciesAnnotationsList.Deserialize("SpeciesAnnotations.xml")
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
		sExample2, 0.9, eValueType, SemanticSimilarityScoreEnum.LIN, bRemoveObsolete);

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
	ExportTable(oOntology, oWorker1, oWorker1.BPVisualizer, "Example1_BPTable.tsv")
	#ExportScatterplot(oOntology, oWorker2, oWorker2.CCVisualizer, "Example2_CCScatterplot.tsv")
	#ExportTreeMap(oOntology, oWorker3, oWorker3.MFVisualizer, "Example3_MFTreeMap.tsv")
	#ExportCytoscapeXGMML(oWorker1.BPVisualizer, "Example1_BPCytoscape.xgmml")
	#ExportSimMat(oWorker1.BPVisualizer, "Example1_BPSimilarityMatrix.tsv")
	#ExportWordClouds(oWorker1, "Example1_WordClouds.json")

	# We are finished

def ExportTable(ontology, worker, visualizer, fileName):
	with open(fileName, 'w') as oWriter:
		oTerms = GOTermList(visualizer.Terms)
		oTerms.FindClustersAndSortByThem(ontology, worker.AllProperties, worker.CutOff)

		oWriter.write("TermID\tName\tValue\t")
		c = 1
		while c < worker.MinNumColsPerGoTerm:
			oWriter.write("UserValue_{0}\t", c - 1)
			c += 1
			
		oWriter.write("LogSize\tFrequency\tUniqueness\tDispensability\tRepresentative\n")

		# print the data
		i = 0
		while i < oTerms.Count:
			oTerm = oTerms[i]
			oProperties = worker.AllProperties.GetValueByKey(oTerm.ID)

			oWriter.write("\"{}\"\t".format(oTerm.FormattedID))
			oWriter.write("\"{}\"\t".format(oTerm.Name))
			oWriter.write("{}\t".format(oProperties.Value))

			c = 1
			while c < worker.MinNumColsPerGoTerm:
				oWriter.write("{}\t".format(oProperties.UserValues[c - 1]))
				c += 1

			oWriter.write("{}\t".format(oProperties.LogAnnotationSize))
			oWriter.write("{}\t".format(oProperties.AnnotationFrequency * 100.0))
			oWriter.write("{}\t".format(oProperties.Uniqueness))
			oWriter.write("{}\t".format(oProperties.Dispensability))

			if oProperties.Representative > 0:
				oWriter.write("{}".format(oProperties.Representative));
			else:
				oWriter.write("null");

			oWriter.write("\n")
			i += 1

def OWorker_OnFinish(sender, e):
	print("Worker {} has finished processing the data in {} seconds.".format(sender.JobID, sender.ExecutingTime.TotalSeconds))

if __name__ == '__main__':
	main()
