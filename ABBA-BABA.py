#load packages
import json
import ipyrad.analysis as ipa
import ipyparallel as ipp
import toytree
import toyplot
#read the text
f = open("taxa_assignment.txt","r")#assign P1,P2,P3,P4
data = f.read().splitlines()
result = []
for i in range(len(data)):
    temp = data[i].strip(",")
    user_dict = json.loads(temp)
    result.append(user_dict)
#muti cores
#In a terminal on your computer you must launch the ipcluster instance by hand, like this:
#ipcluster start -n 10 --daemonize
ipyclient=ipp.Client()
## ipyrad and raxml output files
locifile="sequence.loci" #sequence file
newick="quartet.tre"#species tree file
## create a baba object linked to a data file and newick tree
aa=ipa.baba(data=locifile,newick=newick)
# do not need copy,close ipyclient
aa.tests = result
aa.run(ipyclient)
sorted_results=aa.results_table.sort_values(by="Z", ascending=False)
sorted_taxa=aa.taxon_table.iloc[sorted_results.index]
sorted_taxa.to_csv("taxa.abba-baba.csv")
aa.results_table.to_csv("result.abba-baba.csv")

ipyclient.shutdown(hub=True)
