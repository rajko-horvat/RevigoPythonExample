## Repository description
<p>This is the example Python project that illustrates how to use a REVIGO core library for your own projects.</p>

## How to run this example
<p>To compile from command line:
<ul>
	<li>git clone https://github.com/rajko-horvat/RevigoPythonExample</li>
	<li>Install pythonnet package (pip install pythonnet)</li>
	<li>Compile <a href="https://github.com/rajko-horvat/RevigoCore">RevigoCore library</a> (copy generated binaries to RevigoPythonExample directory)</li>
	<li>Download a set of precompiled databases:
	<a href="http://revigo.irb.hr/Databases/GeneOntology.xml.gz" target="_blank">Gene Ontology</a> and 
	<a href="http://revigo.irb.hr/Databases/SpeciesAnnotations.xml.gz" target="_blank">Species annotations</a>, 
	or build your own databases with <a href="https://github.com/rajko-horvat/RevigoGenerateDatabases">RevigoGenerateDatabases</a> command line utility</li>
	<li>Adjust the path to databases in python code</li>
	<li>Run the RevigoPythonExample.py script.</li>
</ul></p>

## About REVIGO (REduce + VIsualize Gene Ontology) project
<p>Outcomes of high-throughput biological experiments are typically interpreted by statistical testing
for enriched gene functional categories defined by the Gene Ontology (GO). The resulting lists of GO terms 
may be large and highly redundant, and thus difficult to interpret.<p>
<p>REVIGO is a successful project to summarize long, unintelligible lists of Gene Ontology terms by finding a representative subset 
of the terms using a simple clustering algorithm that relies on semantic similarity measures.</p>
<p>For any further information about REVIGO project please see  
<a href="https://dx.doi.org/10.1371/journal.pone.0021800" target="_blank">published paper</a> and  
<a href="http://revigo.irb.hr/FAQ.aspx" target="_blank">Frequently Asked Questions page</a></p>
