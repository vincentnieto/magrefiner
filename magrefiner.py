#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Interactive MAG analysis program by Vincent Nieto

print('Interactive MAG Analysis Program by Vincent Nieto, 2021')
print('Made specifically for Compiled Pasolli Data')


# In[3]:


#Importing tools for data manipulation (pandas) and mathematics (numpy).
print('Importing pandas and numpy...')
import pandas as pd
import numpy as np
print('Import complete!')


# In[5]:


#Making data frames using pandas.
#Replace paths here!

print('Please set the paths to your tab-deliminated files. Please change the file type if necessary (e.g., .txt, .csv).')
df1 = pd.read_csv ('/Users/vin176/Desktop/Scripts/pasolli_analysis/pasolli_eliminases_header.txt', sep='\t')
df2 = pd.read_csv ('/Users/vin176/Desktop/Scripts/pasolli_analysis/mmc3.txt', sep='\t')
df3 = pd.read_csv ('/Users/vin176/Desktop/Scripts/pasolli_analysis/mmc4.txt', sep='\t')
df4 = pd.read_csv ('/Users/vin176/Desktop/Scripts/pasolli_analysis/metadata.txt', sep='\t')
print('Construction of data frames complete!')


# In[6]:


#Merging (joining) data frames using pandas.

df5 = pd.merge(df1, df2, how='left')
print('df5 Done!')
print('Working...')
df6 = pd.merge(df5, df3, how='left')
print('df6 Done!')
print('Working...')
df7 = pd.merge(df6, df4, how='left')
print('df7 Done!')
print("Merging of data frames complete!")


# In[7]:


#Checks the rows and column counts.
print('These are the dimensions (rows, columns) of your merged and unrefined dataframe:')
df7.shape


# In[8]:


#Makes csv file of merged dataframes. 

df7.to_csv('unrefined_dataframe.csv', index=False)
print('Your merged dataframes have been saved as \"unrefined_dataframe.csv\" to your data drive.')


# In[11]:


#Refining your merged dataframe by gene of interest.

print('Let\'s refine your unrefined dataframe.')
df7['QSEQID'] = df7['QSEQID'].str.lower()
gene_name = input("Please type your gene-of-interest: ").lower()
df_gene_filt = df7[df7['QSEQID'].str.contains(gene_name)]
df_gene_filt_len = len(df_gene_filt)
while df_gene_filt_len <= 0:
    gene_name = input("Input invalid, please try again: ")
    df_gene_filt = df7[df7['QSEQID'].str.contains(gene_name)]
    df_gene_filt_len = len(df_gene_filt)
print('These are the dimensions (rows, columns) of your gene-refined dataframe:')
df_gene_filt.shape


# In[12]:


#Refining your dataframe by mismatches.

print('Let\'s now set a cutoff for how many mismatches you would like to keep.')
mismatch_cutoff = input("Keep mismatches greater-than or equal-to: ")
df_mismatch_filt = df_gene_filt[df_gene_filt.MISMATCH >= int(mismatch_cutoff)]
print("These are the dimensions (rows, columns) of your mismatch-refined dataframe:")
df_mismatch_filt.shape


# In[14]:


#Refining your dataframe by percent ID.

print('Let\'s set a range for your data using percent ID.')
piden_low = input("Low range cutoff for percent ID: ")
df_piden_low_filt = df_mismatch_filt[df_mismatch_filt.PIDEN >= int(piden_low)]
piden_high = input("High range cutoff for percent ID: ")
df_piden_high_filt = df_piden_low_filt[df_piden_low_filt.PIDEN <= int(piden_high)]
print("These are the dimensions (rows, columns) of your perent ID-refined dataframe:")
df_piden_high_filt.shape


# In[15]:


#Drop duplicates of sample sequences.

df_sseq_dropped = df_piden_high_filt.drop_duplicates(subset=['SSEQ'])
print("Automatically removing sample sequence duplicates from your dataframe...")
print("Done!")
print("These are the dimensions (rows, columns) of your sample sequence-refined dataframe:")
df_sseq_dropped.shape


# In[16]:


#Refining by length (coverage).

print("Let's refine the dataframe by length (coverage).")
print("We will refine the dataframe one baceteria at a time.")
print("Please calculate the percent cutoff based on the gene size from each bacteria you want to process.")


b_one = input("Enter a bacterial species (e.g., ddesulfuricans): ").lower()
b_one_filt = df_sseq_dropped[df_sseq_dropped.QSEQID.str.contains(b_one)]
b_one_len = input("Lower cutoff for length (coverage): ")
b_one_len_filt = b_one_filt[b_one_filt.LENGTH >= int(b_one_len)]
print("Done!")
print(b_one_len_filt.shape)

b_two = input("Enter a bacterial species (e.g., ecoli): ").lower()
b_two_filt = df_sseq_dropped[df_sseq_dropped.QSEQID.str.contains(b_two)]
b_two_len = input("Lower cutoff for length (coverage): ")
b_two_len_filt = b_two_filt[b_two_filt.LENGTH >= int(b_two_len)]
print("Done!")
print(b_two_len_filt.shape)

b_three = input("Enter a bacterial species (e.g., dalaskensis): ").lower()
b_three_filt = df_sseq_dropped[df_sseq_dropped.QSEQID.str.contains(b_three)]
b_three_len = input("Lower cutoff for length (coverage): ")
b_three_len_filt = b_three_filt[b_three_filt.LENGTH >= int(b_three_len)]
print("Done!")
print(b_three_len_filt.shape)

#To include more bacteria, please remove the hash tags, below.

#b_four = input("Enter a bacterial species (e.g., difficile): ").lower()
#b_four_filt = df_sseq_dropped[df_sseq_dropped.QSEQID.str.contains(b_four)]
#b_four_len = input("Lower cutoff for length (coverage): ")
#b_four_len_filt = b_four_filt[b_four_filt.LENGTH >= int(b_four_len)]
#print("Done!")

#b_five = input("Enter a bacterial species (e.g., bwadsworthia): ").lower()
#b_five_filt = df_sseq_dropped[df_sseq_dropped.QSEQID.str.contains(b_five)]
#b_five_len = input("Lower cutoff for length (coverage): ")
#b_five_len_filt = b_five_filt[b_four_filt.LENGTH >= int(b_five_len)]
#print("Done!")

df_concat1 = pd.concat([b_one_len_filt, b_two_len_filt], ignore_index = True)
print("Working...")
df_concat2 = pd.concat([df_concat1, b_three_len_filt], ignore_index = True)
print("Working...")
#df_concat3 = pd.concat([df_concat2, b_four_len_filt], ignore_index = True)
#print("Working...")
#df_concat4 = pd.concat([df_concat3, b_five_len_filt], ignore_index = True)
#print("Working...")
print("Concatenation done!")
print("These are the dimensions (rows, columns) or your length-refined dataframe:")

#Change number, below, depending on how many dataframes you will concatenate.
df_concat2.shape


# In[17]:


#Makes CSV file of your refined dataframe.
#Change the number, below, depending on how many dataframes you concatenated.

df_concat2.to_csv('refined_dataframe.csv', index=False)
print('Your refined dataframe have been saved as \"refined_dataframe.csv\" to your data drive.')


# In[ ]:




