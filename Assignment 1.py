#!/usr/bin/env python
# coding: utf-8

# On this assignment I worked with Shyan Polman

# In[61]:


import numpy as np

def FastQC (fastqc):
    
    scores =[]
    file = open(fastqc, "r")

    line_num = 0 
    
    #retrieve all scores from file
    for line in file:
        line_num += 1
        if line_num%4 == 0:
            scores.append(line.strip())
            
    diction = {}

    for score in scores: 
        key = 0
        
        for letter in score:
            Q = ord(letter) - 33

            if key in diction: 
                value = diction[key]
                value.append(Q)
                diction[key] = value
                
            else:
                diction[key]= [Q]
            key += 1

    print("Position\tMean")
    for position, scores in diction.items():
        mean = np.mean(scores)
        print(f"{position + 1}\t\t{mean:.2f}")
        
    #write out file for distribution scores
    with open('qscores.txt', 'w') as file:
        file.write("Position\tMean"+'\n')
        for position, scores in diction.items():
            mean = np.mean(scores)
            file.write(f"{position + 1}\t\t{mean:.2f}"+'\n')

    file.close()
    return 


# In[62]:


FastQC('assignment1.fastq')


# In[45]:


from Bio import SeqIO
def deBruijn_kmer (fastqc, k):
    
    reads =[]
    graph = {}
    with open(fastqc, 'r') as f:
        for record in SeqIO.parse(f, 'fastq'):
            reads.append(str(record.seq))
            
            
    for read in reads:
        for i in range(len(read) - k + 1):
             #creating kmers = nodes
            node = read[i:i+k]
            
            #find edges through kmers using suffix and prefix
            prefix = node[:-1]
            suffix = node[1:]
            
            if prefix in graph:
                graph[prefix].add(suffix)
                
            else:
                graph[prefix] =set()
                graph[prefix].add(suffix)
    
    #write out file for nodes
    with open('nodes.txt', 'w') as file:
        for node in graph.keys():
            file.write(str(node) +'\n')
            
    #write out file for edges
    with open('edges.txt', 'w') as file:
        for edge in graph.values():
            file.write(str(edge) +'\n')    
  
            
    return graph


# In[46]:


deBruijn_kmer('assignment1.fastq', 8)

