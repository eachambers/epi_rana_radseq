
# coding: utf-8

# In[2]:


import pandas as pd
import numpy
import itertools
import csv
import glob


# In[3]:


# Make a list of the file names I want to iterate through using glob
ranafilenames = glob.glob('../rana2brad/*outfiles/*stats.txt')
print ranafilenames


# In[110]:


def get_file_list(directory):
    return glob.glob(directory)

def get_line_numbers(filename):
    counter=0
    infile=open(filename, 'rt')
    for line in infile:
        counter+=1
        if 'var  sum_var' in line:
            varline=counter
#            print name,varline # not necessary but want to just check
#            print line
        elif '## Final' in line:
            endline=counter
#            print name,endline
#            print line
        else:
            continue
    infile.close()
    return [varline,endline]

def slice_data(start,end,filename):
    snp_dist=[]
    infile=open(filename,'rU')
    for lines in itertools.islice(infile, start, end-3):
        lines2 = lines.strip().split()
#        print lines2
        lines2.append(filename) # eventually change to whatever file it's on
        snp_dist.append(lines2)
    infile.close()
    return snp_dist
        

def pd_conversion(filename):
    nums = get_line_numbers(filename)
    snp_dist = slice_data(nums[0],nums[1],filename)
    snpdist_labels=['number','variable','sum_var','pis', 'sum_pis', 'clust_threshold']
    df_snpdist = pd.DataFrame.from_records(snp_dist, columns=snpdist_labels)
    return df_snpdist  

def main():
    directory='../epiddrad/*outfiles/*stats.txt'
    #filename='../epi2brad/clust_92_outfiles/clust_92_stats.txt'
    #directory=filename
    file_list = get_file_list(directory)
    dfs=[]
    for filename in file_list:
        print filename
        pd_df = pd_conversion(filename)
        dfs.append(pd_df)
        pd_df.to_csv(filename+"_snpdist.csv")
    return dfs

print main()


# In[12]:


# MESSY CODE; INCLUDES NOTES I TOOK WHILE WRITING IT; IGNORE AND DON'T RUN

# The *stats.txt files produced by iPyrad have two blocks of data we want to retrieve: the distribution of SNPs
# (vars and pis) per locus, and the final sample stats summary.
# Start with a list that will contain the first set of data, distribution of SNPs
snp_dist = []

# To just see line numbers within the file:
#with open("epi2brad/clust_84_outfiles/clust_84_stats.txt", "rt") as infile:
#   for line in infile:
#        lines.append(line)
#print(list(enumerate(lines)))

# for now, just leave hard-coded for epifilenames[0]; happens to be clust_threshold=91

#infile=open(ranafilenames[0], 'rt')

def get_line_numbers(filename):
    counter=0
    infile=open(filename, 'rt')
    for line in infile:
        counter+=1
        if 'var  sum_var' in line:
            varline=counter
            print name,varline # not necessary but want to just check
            print line
        elif '## Final' in line:
            endline=counter
            print name,endline
            print line
        else:
            continue
    return varline,endline

#get_line_numbers('../rana2brad/clust_95_outfiles/clust_95_stats.txt')
                
def slice_data(start,end,filename)
    for lines in itertools.islice(filename, start, end-3):
        lines2 = lines.strip().split()
        print lines2
        lines2.append(filename) # eventually change to whatever file it's on
        snp_dist.append(lines2)


#print snp_dist
#for name in ranafilenames:
 #   infile=open(name, 'rt')
  #  for lines in itertools.islice(infile, varline, endline-3):
   #     lines2 = lines.strip().split()
    #    lines2.append(name) # eventually change to whatever file it's on
     #   snp_dist.append(lines2)

    #infile.close() # because it's iterated through the whole file already so we need to close and reopen

#snpdist_labels=['number','variable','sum_var','pis', 'sum_pis', 'clust_threshold']
df_snpdist = pd.DataFrame.from_records(snp_dist, columns=snpdist_labels)

print df_snpdist

#df_snpdist.to_csv("./rana2brad_snpdist.csv")


# In[111]:


# The SECOND SET of data we're interested in, the summary stats

def get_file_list(directory):
    return glob.glob(directory)

def get_line_numbers(filename):
    counter=0
    infile=open(filename, 'rt')
    for line in infile:
        counter+=1
        if 'state  reads_raw' in line:
            varline=counter
#            print name,varline # not necessary but want to just check
#        print line
        else:
            continue
    infile.close()
    return [varline,varline]

def slice_data(start,end,filename):
    sum_stats=[]
    infile=open(filename,'rU')
    for lines in itertools.islice(infile, start, end+12):
        lines2 = lines.strip().split()
        print lines2
        lines2.append(filename) # eventually change to whatever file it's on
        sum_stats.append(lines2)
    infile.close()
    return sum_stats
        

def pd_conversion(filename):
    nums = get_line_numbers(filename)
    sum_stats = slice_data(nums[0],nums[1],filename)
    sumstats_labels=['sample', 'state', 'reads_raw', 'reads_passed', 'clust_total', 'clust_hidepth','hetero_est','error_est','reads_consens','loci_assembly','clust_threshold']
    df_sumstats = pd.DataFrame.from_records(sum_stats, columns=sumstats_labels)
    return df_sumstats  

def main():
    directory='../epiddrad/*outfiles/*stats.txt'
    #filename='../rana2brad/clust_95_outfiles/clust_95_stats.txt'
    #directory=filename
    file_list = get_file_list(directory)
    dfs=[]
    for filename in file_list:
        print filename
        pd_df = pd_conversion(filename)
        dfs.append(pd_df)
        pd_df.to_csv("./" +filename+ "_sumstats.csv")
    return dfs

print main()


# In[114]:


# The LAST SET that we're interested in is the LOCUS COVERAGE; this will tell us about missing data and shared loci.

def get_file_list(directory):
    return glob.glob(directory)

def get_line_numbers(filename):
    counter=0
    infile=open(filename, 'rt')
    for line in infile:
        counter+=1
        if 'locus_coverage' in line:
            varline=counter
#            print name,varline # not necessary but want to just check
#        print line
        elif '## The distribution' in line:
            endline=counter
            print name,endline
        else:
            continue
    infile.close()
    return [varline,endline]

def slice_data(start,end,filename):
    coverage=[]
    infile=open(filename,'rU')
    for lines in itertools.islice(infile, start, end-3):
        lines2 = lines.strip().split()
        print lines2
        lines2.append(filename) # eventually change to whatever file it's on
        coverage.append(lines2)
    infile.close()
    return coverage

def pd_conversion(filename):
    nums = get_line_numbers(filename)
    coverage = slice_data(nums[0],nums[1],filename)
    cov_labels=['number','locus_coverage', 'sum_coverage','clust_threshold']
    df_cov = pd.DataFrame.from_records(coverage, columns=cov_labels)
    return df_cov  

def main():
    directory='../rana2brad/*outfiles/*stats.txt'
    #filename='../rana2brad/clust_95_outfiles/clust_95_stats.txt'
    #directory=filename
    file_list = get_file_list(directory)
    dfs=[]
    for filename in file_list:
        print filename
        pd_df = pd_conversion(filename)
        dfs.append(pd_df)
        pd_df.to_csv("./" +filename+ "_sumstats.csv")
    return dfs

print main()


# In[109]:


# Now, let's work on merging these dataframes together from each parameter run.
# They're contained within their _outfile folder, so we need to pull them out and merge them somehow
# To make life easier, I just did this in terminal using the following:
# $ mv *outfiles/*.csv .

# Now they're all present within our directory of choice (../rana2brad/*.csv)

df1 = pd.read_csv("../epi2brad/clust_80_stats.txt_sumstats.csv")
df2 = pd.read_csv("../epi2brad/clust_81_stats.txt_sumstats.csv")
df3 = pd.read_csv("../epi2brad/clust_82_stats.txt_sumstats.csv")
df4 = pd.read_csv("../epi2brad/clust_83_stats.txt_sumstats.csv")
df5 = pd.read_csv("../epi2brad/clust_84_stats.txt_sumstats.csv")
df6 = pd.read_csv("../epi2brad/clust_85_stats.txt_sumstats.csv")
df7 = pd.read_csv("../epi2brad/clust_86_stats.txt_sumstats.csv")
df8 = pd.read_csv("../epi2brad/clust_87_stats.txt_sumstats.csv")
df9 = pd.read_csv("../epi2brad/clust_88_stats.txt_sumstats.csv")
df10 = pd.read_csv("../epi2brad/clust_89_stats.txt_sumstats.csv")
df11 = pd.read_csv("../epi2brad/clust_90_stats.txt_sumstats.csv")
df12 = pd.read_csv("../epi2brad/clust_91_stats.txt_sumstats.csv")
df13 = pd.read_csv("../epi2brad/clust_92_stats.txt_sumstats.csv")
df14 = pd.read_csv("../epi2brad/clust_93_stats.txt_sumstats.csv")
df15 = pd.read_csv("../epi2brad/clust_94_stats.txt_sumstats.csv")
df16 = pd.read_csv("../epi2brad/clust_95_stats.txt_sumstats.csv")

sumstats_concat = df1.append([df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16])
sumstats_concat.to_csv("epi2brad_sumstats.csv")

