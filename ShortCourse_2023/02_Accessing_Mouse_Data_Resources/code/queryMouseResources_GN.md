
**This notebook/tutorial shows how one can query MouseMine and GeneNetwork to retrieve data programmatically**  

# Using GeneNetwork to retrieve a phenotype

***import requests*** : imports the Python library that contains the necessary functions to query a web resource.

***resp = requests.get("url")*** : Runs the query and gets results back. Here, we will make use of GeneNetwork.org's programmable API to retrieve trait data for a single trait collected in the BXD RI lines.

***print()*** : The results are contained in the resp variable and stored in the "text" attribute.

***with open("filename", 'w') as fl**: Writes the contents to a file. Here the file is named "queryGeneNetwork_SingleTrait.tsv"

***fl.writelines(resp.text)*** : The function that write the contents of resp.text to the file.


```python
import requests
resp = requests.get("http://robot.genenetwork.org/webqtl/main.py?cmd=trait&probeset=18495&db=BXDPublish&format=col") 
print (resp.text)

with open("queryGeneNetwork_SingleTrait.tsv", 'w') as fl:    
    fl.writelines(resp.text)

```

# Using GeneNetwork to retrieve a phenotype

***import requests***  : imports the Python library that contains the necessary functions to query a web resource.

***resp = requests.get("url")*** : Runs the query and gets results back. Here, we will make use of GeneNetwork.org's programmable API to retrieve gene expression for Syn1 (gene=Syn1) in the hippocampus (tissue=hip) collected in the BXD RI lines.

***print()***  : The results are contained in the resp variable and stored in the "text" attribute.

***with open("filename", 'w') as fl**: Writes the contents to a file. Here the file is named "queryGeneNetwork_SingleGene_Tissue.tsv"

***fl.writelines(resp.text)*** : The function that write the contents of resp.text to the file.


```python
import requests
geneData = requests.get("http://www.genenetwork.org/webqtl/main.py?cmd=sch&gene=Syn1&tissue=hip&format=text") 
print (geneData.text)

with open("queryGeneNetwork_SingleGene_Tissue.tsv", 'w') as fl:    
    fl.writelines(geneData.text)

```
