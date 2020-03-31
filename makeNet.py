#importing libraries

import argparse
import time
from os import system
import sys 
try:
	import networkx as nx
except:
	flag == False
	while flag == False:
		option = input("Would you like to install networkx in your system? Y/N")
		if option.lower() == "n":
			exit()
		else:
			if option.lower() == "y":
				try:
					system("pip3 --user install networkx")
					import networkx as nx
				except:
					print("Error installing networkx. Please install it manually. Exiting")
					exit()
	
	#print("You do not have the library networkx in you system. Please install it before use makeNet. Exiting")
	#exit()
	
#defining functions
def done():
	print ("Done (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+")")
	exit()

#to save graph
def save(G, Format, outFile):
	if Format == "edge_list":
		nx.write_edgelist(G, outFile)
	elif Format == "GEXF":
		nx.write_gexf(G, outFile)
	elif Format == "GML":
		nx.write_gml(G, outFile)
	elif Format == "Pickle":
		nx.write_gpickle(G, outFile)
	elif Format == "GraphML":
		nx.write_graphml(G, outFile)
	elif Format == "YAML":
		nx.write_yaml(G, outFile)
	elif Format == "Pajek":
		nx.write_pajek(G, outFile)
	else:
		print("Your format selection can not be processed. Exiting")
		done()

	return
	
if __name__ == "__main__":


	print("\n\n##################################")
	print("##                              ##")
	print("##     makeNet version 2.0.     ##")
	print("##  Data taken from regNetwork  ##")
	print("##       on March 31, 2020      ##")
	print("##                              ##")
	print("##          Created by          ##")
	print("## Sebastian Contreras-Riquelme ##")
	print("## Search me on github as Cold7 ##")
	print("##                              ##")
	print("##################################\n\n")
	
	#checking python version
	if (sys.version_info < (3, 0)):
		print("You are using a python version lower than 3.0. Please use python3 to execute this script. Exiting")
		exit()

	#input
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--organism", help="Organism to use. Default: human", default="human", choices=['human','mouse'])
	parser.add_argument("-d", "--database", help="Database to use. Default: All", default="All", choices=["All", "biogrid", "ensembl", "hprd", "intact", "kegg", "microT", "miRanda", "miRBase", "miRecords", "miRTarBase", "PicTar", "string", "Tarbase", "TargetScan", "transmir", "tred", "ucsc"])
	parser.add_argument("-e", "--evidence", help="Evidence to use. Default: All", default="All", choices=["All", "Experimental", "Predicted"])
	parser.add_argument("-c", "--confidence", help="Confidence to use. Default: All", default="All", choices=["All", "High", "Medium", "Low"])
	parser.add_argument("-t", "--tfList", help="File with TF list.", required=True)
	parser.add_argument("-O","--output", help="Output file to save results", required=True)
	
	args = parser.parse_args()


	print("Starting makeNet (Date: "+str(time.strftime("%y/%m/%d"))+" Time: "+str(time.strftime("%H:%M:%S"))+") using the following parameters:")
	print("  Organism: "+args.organism)
	print("  Database: "+args.database)
	print("  Evidence: "+args.evidence)
	print("  Confidence: "+args.confidence)
	print("  Expression files: "+args.tfList)
	print("  Output File: "+args.output)
	print("\n")
	
	#getting net name
	netName = ""
	if args.organism.lower() == "human":
		netName = "./dat/human_all_regnetwork.csv"
	elif args.organism.lower() == "mouse":
		netName = "./dat/mouse_all_regnetwork.csv"
	else:
		print("Your organism currently is not supported. Exiting")
		done()
	#######################
	## making source graph
	#######################
	G =nx.MultiDiGraph()
	f = open(netName,"r")
	for line in f:
		line = line[:-1].split("\",")
		line2 = []
		for item in line:
			line2.append(item.replace("\"",""))
		line = line2
		node1 = line[0].split("::")
		node2 = line[2].split("::")
		for node in node1:
			if node not in G.nodes():
				G.add_node(node)
		for node in node2:
			if node not in G.nodes():
				G.add_node(node)
		
		for n1 in node1: #some nodes from regnetwork are in XXX::YYY format (coregulating a gene) so i had splitted it and now i will loop over list
			for n2 in node2:
				dbase = line[4].split(",")
				ev = line[5]
				confid = line[6]
				for db in dbase:
					if (args.database == "All" or args.database == db) and (args.evidence == "All" or args.evidence == ev) and (args.confidence == "All" or args.confidence == confid):
						G.add_edge(n1,n2,database=db, evidence=ev, confidence=confid)
	f.close()

	
	##########################
	## filtering by gene list
	##########################
	f = open(args.tfList,"r")
	geneList = []
	for gene in f:
		gene = gene[:-1]
		if gene not in geneList and gene != "":
			geneList.append(gene)
	f.close()

	#getting TFs in G using outdegree
	tfList = []
	for node in G.nodes():
		if G.out_degree(node)>0:
			tfList.append(node)
	
	edgesToDel = []
	for tf in tfList:
		flag = False
		for gene in geneList:
			if gene == tf:
				flag = True
		
		if flag == False: #tf not in the sample
			neigh = G.neighbors(tf)
			for n in neigh:
				edgesToDel.append([tf, n])
	
	for edges in edgesToDel:
		G.remove_edge(edges[0],edges[1])

	##########################
	## saving graph
	##########################
	#creating digraph to remove multiple edges between two nodes in order to save the final graph
	H = nx.DiGraph()
	for edges in G.edges():
		H.add_edge(edges[0],edges[1])

	output = open(args.output,"w")
	for edges in H.edges():
		output.write(edges[0]+"\t"+edges[1]+"\t1\n")
	output.close()
		
	done()
	

