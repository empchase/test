import sys
import argparse
import re 
import os
import csv
import pandas as pd
import numpy as np

# duckdb_dir = '~/bin/duckdb'
# # Append the DuckDB directory to sys.path
# sys.path.append(duckdb_dir)

import duckdb
conn = duckdb.connect('/global/scratch/users/empchase/A10_sequencing/v2/analysis2.db') # connect to database with TBB map


# define functions for data processing

def find_designed(des):
    """ takes a file where every line is a designed tile and creates a lookup dictionary of all designed tiles """
    dt = [] # initialize a list of designed tiles from design file
    with open(des, 'r') as f_des:
        # open the design file for reading
        for line in f_des:
            if "ArrayDNA" in line: #skip the "arrayDNA" line if applicable
                pass
            else:
                dt.append(line.strip()) #add tile to list without whitespace
    
    dt_dict = {} # initialize dictionary to lookup tiles in
    for i in dt:
        dt_dict[i] = 1 # add tiles into a diction ary
        
    print('Number of designed tiles:', len(dt_dict)) 

    return dt_dict

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', 'X':'X'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)
def reverse_complement(s):
        return complement(s[::-1])
    
    
def getmid(seq, pre, post, bclen):
    # seq = the sequence to parse
    # pre = substring that precedes piece of interest
    # post = substring that follows piece of interest
    # returns piece of interest

    re_key = pre + "(.*)"+ post 
    poi_search = re.search(re_key, seq)
    if poi_search is None:
        #the barcode will be called X
        poi = "X"
        
        #then we search for which restriction site had the error
        #regex for the bc we want to ignore
        w = "(.{"+str(bclen)+"})" 
        pre_re = pre + w + "(.{7})"
        pre_search = re.search(pre_re, seq)
        post_re = "(\w{7})" + w + post
        post_search = re.search(post_re, seq)
        
        if pre_search is None and post_search is None:
            a = 'X'
            z = 'X'
        elif pre_search is None:
            poi = post_search.group(2)
            a = post_search.group(1)
            z = post
        elif post_search is None:
            poi = pre_search.group(1)
            z = pre_search.group(2)
            a = pre
        else:
            a = 'Z'
            z = 'Z'            
    else:
        poi = poi_search.group(1)
        a = pre
        z = post
    
    return poi, a, z

#putative consensus sequences ***reverse complement of snapgene***
adp = 'CGGGCCC'#7 bp ; beforeAD barcode in read1
adf = 'GGCGCGC' #7bp ; after AD barcode in read1

tp = "TAGTCA" # 6 bp ; before AD Tile in read1
tf = "GCTAGC"

# function that just looks for Tile/AD bc
def adbc_mapper(readfile, dtiles, tile_pre = tp, tile_post = tf, 
                  adBC_pre = adp, adBC_post = adf):
    #readfile = fastq (ADBC fastqs are typically paired)
    # dtiles = dictionary of designed tiles, output of find_designed
    # *_pre or *_post = the consensus sequences before or after each feature, defaults defined above

    # make lists of reads
    readlist = []    
    with open(readfile, 'r') as fin:
        for line in fin:
            if line.startswith('@'):
                #look at next line to get read sequence, add to list
                seq = next(fin)
                seq = seq.strip()
                readlist.append(seq)
                
    #make lists of tiles/BCs from list of reads
    tile_list = []
    tile_lengths= []
    
    des_query = [] # tells us if tile matches design or not
    
    adBC_list = []
    adBC_lengths = []
    
    tile_adbc = []
   
    
    for read in readlist:
        tile, tpre, tpost = getmid(read, tile_pre, tile_post, 120) #use consensus seq to find tile
        tile = reverse_complement(tile)
        tile_list.append(tile) #add tile to list
        tile_len = len(tile) #find length of tile
        tile_lengths.append(tile_len) #add length to list
            
        if tile in dtiles:
            des_query.append(1)
        else:
            des_query.append(0)
#         print(tile)
        
        adBC, prea, posta = getmid(read, adBC_pre, adBC_post, 11)
        adBC = reverse_complement(adBC)
        adBC_list.append(adBC)
        adBC_len = len(adBC)
        adBC_lengths.append(adBC_len)
        
        tile_adbc.append(tile+'-'+adBC)

            
    # make the df
    tileBC_dict = {"Tiles":tile_list, "T Len" : tile_lengths, "Designed": des_query, 
                  "AD BCs":adBC_list, "A Len": adBC_lengths, 'Tile-AD': tile_adbc}
    tileBC_df = pd.DataFrame.from_dict(tileBC_dict)
    
    libname = '_'.join(readfile.split('/')[-1].split('_')[0:6])
    tileBC_df['Library'] = libname
    
    return tileBC_df


def SQLanalyze_tiles_adbcs (df, bc_len=11, tile_len=120):
    # df = barcode containing df, parsed from fastq
    # bc_len = int, expected barcode length
    # tbb_dictkey = str, either 'ADbc' or 'RPTRbc'
    
    print(df.loc[0,'Library'])

    tr = df.shape[0] #print total reads
    print(f'Total Reads {tr}')
    cls = df[(df['A Len']== bc_len) & (df['T Len']==tile_len)] #cl = correct length
    print('Reads w BC and Tile of correct length')
    clcount = cls.shape[0]
    print(clcount) #print total reads w correct len bc/ad
    
    print('% Reads w correct length Tiles/BCs')
    clpct = cls.shape[0]/df.shape[0]
    print (clpct)
    
    
    #filter for only designed tiles:
    cls = cls[cls['Designed']==1]
    descount = cls.shape[0]
    print(f'# designed {descount}') #print total reads w correct len bc/ad
    
    
    #df of BC coverage
    cl_covdf = cls['Tile-AD'].value_counts().to_frame().reset_index() 
    print('# Unique Tile-ADbc pairs')
    uniqbccount = cl_covdf.shape[0]
    print(uniqbccount)
    print('SUM Unique BCs')
    print(cl_covdf.sum(numeric_only=True)['Tile-AD'])

    #copy down unique BCs into a list
    bcs = cl_covdf['index'].tolist()
    matchlist = [] #list of TBBs that match the BCs
    for x in bcs:
        try:
            both = x.split('-')
            ak = both[1]
            tk = both[0]
            myquery = conn.sql("SELECT RPTR_BC FROM A10_2_T_NODBLMAP_20240422 WHERE AD_BC='{}' AND AD='{}'".format(ak, tk)).to_df() #try ("SELECT RPTRBC FROM TBB WHERE ADBC={adfill} AND TILES={tilefill}".format(adfill=ak, tilefill=tk))
            pr = myquery.loc[0,'RPTR_BC'] #pr as in putative reporter
        except KeyError:
            matchlist.append(0)
        else:
            matchlist.append(pr)
    
            
    cl_covdf['PutativeRPTR'] = matchlist
    
    matchesonly = cl_covdf.replace(0, np.nan)
    matchesonly = matchesonly.dropna()
    totmatches = matchesonly.sum(numeric_only=True)['Tile-AD']
    print('# BC matches to A10 deep seq map')
    mapmatches = matchesonly.shape[0]
    print(mapmatches)
    
    print('TOT BC matches to A10 deep seq map')
    print(totmatches)    
    
    print()
    
    #label the df with library name
    libname = df.loc[0,'Library']
    matchesonly['Library'] = libname
    
    statslist = [libname, tr, clcount, clpct, descount, uniqbccount, mapmatches, totmatches]
    print(statslist)
    
    return matchesonly
    


# set up main code to be run
def main(input_file, output_file, design_file):
    # get list of designed tiles
    designed_tile_dictionary = find_designed(design_file)

    # process input file into  df
    rawmap = adbc_mapper(input_file, designed_tile_dictionary)

    # Filter raw df 
    countsdf = SQLanalyze_tiles_adbcs(rawmap)


    # write the output file to csv
    countsdf.to_csv(output_file)





if __name__ == "__main__":
    # Check if the script is being run as the main program

    # Create ArgumentParser object with a description
    parser = argparse.ArgumentParser(description='Build Tile-BC map from step1 fastq file')

    # Add arguments to the parser
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file path')
    # Add an argument for the input file path. -i or --input specifies the argument name.
    # It's required (required=True), of type string (type=str), and has a help message.

    parser.add_argument('-o', '--output', type=str, required=True, help='Raw data output file path')
    # Add an argument for the output file path. -o or --output specifies the argument name.
    # It's required (required=True), of type string (type=str), and has a help message.

    # parser.add_argument('-c', '--clean_output', type=str, required=True, help='Clean data output file path')
    # # Add an argument for the output file path. -o or --output specifies the argument name.
    # # It's required (required=True), of type string (type=str), and has a help message.

    parser.add_argument('-d', '--design', type=str, required=True, help='Design file (each line in file is a designed tile) path')
    # Add an argument for the output file path. -o or --output specifies the argument name.
    # It's required (required=True), of type string (type=str), and has a help message.

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args.input, args.output, args.design)
    conn.close()