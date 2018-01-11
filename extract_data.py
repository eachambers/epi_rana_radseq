
# coding: utf-8

# In[2]:


import pandas as pd
import numpy
import itertools
import csv
import glob


# In[3]:


# Make a list of the file names I want to iterate through using glob
epifilenames = glob.glob('epi2brad/*outfiles/*stats.txt')
print epifilenames


# In[4]:


# The *stats.txt files produced by iPyrad have two blocks of data we want to retrieve: the distribution of SNPs
# (vars and pis) per locus, and the final sample stats summary.
# Start with a list that will contain the first set of data, distribution of SNPs
snp_dist = []

# To just see line numbers within the file:
#with open("epi2brad/clust_84_outfiles/clust_84_stats.txt", "rt") as infile:
#   for line in infile:
#        lines.append(line)
#print(list(enumerate(lines)))

#for name in epifilenames:
#    open(name,'rt') # iterate through filenames generated from glob

# for now, just leave hard-coded for epifilenames[0]; happens to be clust_threshold=91

infile=open(epifilenames[0], 'rt')
counter=0

for line in infile:
    counter+=1
    if 'var  sum_var' in line:
        varline=counter
        print varline # not necessary but want to just check
        print line
    elif '## Final' in line:
        endline=counter
        print endline
        print line
    else:
        continue

infile.close() # because it's iterated through the whole file already so we need to close and reopen

#infile=open('epi2brad/clust_84_outfiles/clust_84_stats.txt','rt')
infile=open(epifilenames[0], 'rt')

#assert itertools.islice(infile, endline-2, endline-1) # throw an error if endline-2 isn't a blank line

for lines in itertools.islice(infile, varline, endline-3):
    lines2 = lines.strip().split()
    lines2.append('91') # eventually change to whatever file it's on
    snp_dist.append(lines2)
        
print snp_dist

infile.close()

snpdist_labels=['number','variable','sum_var','pis', 'sum_pis', 'clust_threshold']
df_snpdist = pd.DataFrame.from_records(snp_dist, columns=snpdist_labels)

print df_snpdist

df_snpdist.to_csv("./epi2brad_clust_91_snpdist.csv")


# In[5]:


print epifilenames


# In[75]:


# Another list containing the second set, the summary stats
sum_stats = []

infile=open(epifilenames[0], 'rt') # again, hard-coded for clustering threshold 91
counter=0
for line in infile:
    counter+=1
    if 'state  reads_raw' in line:
        stateline=counter
        print varline # not necessary but want to just check
        print line
    else:
        continue
     
infile.close() # because it's iterated through the whole file already so we need to close and reopen

infile=open(epifilenames[0], 'rt')

for lines in itertools.islice(infile, stateline, stateline+12):
    lines2 = lines.strip().split()
    lines2.append('91') # eventually change to whatever file it's on
    sum_stats.append(lines2)
        
print sum_stats

infile.close()

sumstats_labels=['sample', 'state', 'reads_raw', 'reads_passed', 'clust_total', 'clust_hidepth','hetero_est','error_est','reads_consens','loci_assembly','clust_threshold']
df_sumstats = pd.DataFrame.from_records(sum_stats, columns=sumstats_labels)

print df_sumstats

df_sumstats.to_csv("./epi2brad_clust_91_sumstats.csv")

