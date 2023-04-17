##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
## If using Espalier or this code, please cite:
##
##      Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.
##
############################

import subprocess
import sys
import glob
import re
import pandas as pd
import dendropy
from Bio import SeqIO

"""
    Wrappers and utility functions for working with ARGweaver
"""

def run_arg_sampler(seq_file,out_file,Ne,recombination_rate,mut_rate,verbose=False):
    
    """
        Run ARGweaver MCMC sampler to reconstruct trees        
        -ntimes controls # of discretized time steps
        -c is a compression param
        -n is the number of mcmc samples, default is to sample every 10 mcmc iterations.
        --no-compress-output can be used if don’t want output gzipped
    """
    
    # Compute expected tree length in order to set reasonable max_time
    sample_size = len(list(SeqIO.parse(seq_file, "fasta")))
    tree_height = 0.0
    for k in range(2,sample_size+1):
        tree_height += Ne / (k*(k-1)/2)
    max_time = tree_height * 100.0 # multiply by 100 to make sure MRCA is below max_time
    
    flag = False
    cmd = 'arg-sample -f ' + seq_file + ' -o ' + out_file + ' -n 100 -N ' + str(Ne) + ' -r ' + str(recombination_rate) + ' -m ' + str(mut_rate) + ' --ntimes 40 --maxtime ' + str(max_time) + ' --no-compress-output'
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        if verbose:
            sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        flag = True #sys.exit(1)

    return flag


def parse_smc(smc_files):
    
    """
        Parse trees from .smc files created by ARGWeaver
    """
    
    # Set up lists to hold data
    samples = []
    start_positions = []
    end_positions = []
    newicks = []
    for smc in smc_files:
        
        sample_num = int(smc.split('.')[-2])
        print(sample_num)
        f = open(smc)
        line = f.readline()
        while line:
            if 'NAMES' in line:
                "Parse taxon naming scheme"
                fields = re.split(r'\t+', line.strip('\n'))
                fields.pop(0) # remove first
                taxa_nums = list(range(len(fields)))
                taxa_map = {str(key):value for key,value in zip(taxa_nums,fields)} # dict maps taxa numbers back to original taxa labels
            if 'TREE' in line: #if a tree
                fields = re.split(r'\t+', line)
                samples.append(sample_num)
                start_positions.append(int(fields[1]))
                end_positions.append(int(fields[2]))
                newicks.append(fields[3])
            line = f.readline()
        f.close()
    
    # Convert to pandas dataframe
    d = {'Sample':samples,'StartPosition':start_positions,'EndPosition':end_positions,'Newick':newicks}
    df = pd.DataFrame(d)
    
    return df, taxa_map

def local_consensus_tree(df,site,root=False):

    """
    Get a local consensus tree based on all sampled ARG local trees for a particular site
    """
    
    # NOTE: Sites are indexed from 1 here
    newicks = df['Newick'][(df['StartPosition'] <= site) & (df['EndPosition'] >= site)]
    
    # TODO: May want to write this to a file rather than joining as one big string"
    trees_as_str = "".join(newicks)

    """
    The following block generates inputs and runs sumtrees.py on smc trees
    """

#    file1 = open('myfile.txt', 'w')

#    file1.writelines(trees_as_str)

#    file1.close()

#    cladefreq = '0.95'

#    burnin = '20'	# only 11 trees for each site

#    outfile = 'ConsensusTree' + str(site) + '.tre'

#    print(outfile)

#    subprocess.run(['sumtrees.py', '-f', cladefreq, '-b', burnin, '-o', outfile, 'myfile.txt', '--no-analysis-metainformation', '-Fnewick'])

#    os.remove(outfile)

    """
    end of sumtrees.py block
    """


    """
    dendropy consensus tree building code block
    """

    trees = dendropy.TreeList()
    trees.read(data=trees_as_str, schema="newick",rooting="default-rooted")
#    print(trees)
    
    # Get MCMC consensus
    consensus = trees.maximum_product_of_split_support_tree()
#    print("\nTree {} maximizes the product of split support (log product = {}): {}".format(trees.index(mcct), mcct.log_product_of_split_support, mcct))
    
    # Should not need to root trees since ARGs should be rooted
    #if root:
        #try:
        #    consensus.reroot_at_midpoint()
        #except:
        #    print()
    
    return consensus

def concat_alignments(seq_files,file_out):

    for idx, file in enumerate(seq_files):
        seq_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta")) # one line alternative
        if idx == 0:
            concat_seqs = seq_dict
        else:
            for key in concat_seqs:
                concat_seqs[key] += seq_dict[key]            

    # Convert to list and write seq records to fasta file
    concat_records = [concat_seqs[key] for key in concat_seqs]
    SeqIO.write(concat_records, file_out, "fasta")


if  __name__ == '__main__':
    
    print("Running...")
    smc_file_path = '/Users/mfarman/SMC_trees/' # dir where smc files live    
    print(smc_file_path)
    smc_list = glob.glob(smc_file_path + '*.smc') # get all smc files in dir as list
    print(smc_list)
    df,taxa_map = parse_smc(smc_list) # parse smc files into pandas df
    df.sort_values(['Sample', 'StartPosition'], inplace=True)
    treenum = 0
    for site in (500000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 7800000):
        treenum = treenum + 1
        consensusTree = local_consensus_tree(df,site,root=False)
        outfile = 'MCCT' + str(treenum) + '.tre'
        file1 = open(outfile, 'w')
        file1.writelines(str(consensusTree) + ";")
        file1.close()        
        #print(consensusTree, ";", sep="")
    """
        Assign mrca times to trees using dendropy
    """

    df['tMRCA'] = 0.0
    for index, row in df.iterrows():
        tree = dendropy.Tree.get(data=row['Newick'], schema="newick", rooting="default-rooted") #, taxon_namespace=tree.taxon_namespace)
        tree.calc_node_ages()
        df.at[index,'tMRCA'] = tree.seed_node.age
        
    df.to_csv('sampled_arg_tMRCAs.csv',index=False) 
