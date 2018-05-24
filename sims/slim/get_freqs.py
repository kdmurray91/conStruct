#!/usr/bin/env python3
description = '''
Get, and dump to text, the allele frequencies of some individuals from the .trees output from a SLiM simulation.
'''

import sys, os
import gzip
import glob
import re
import argparse
import struct
import numpy as np
import random

import msprime

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--tree_file", "-t", type=str, nargs="*", dest="tree_file", 
                    help="name of file to load tree sequences from [default: .trees files in basedir.")
parser.add_argument("--basedir", "-o", type=str, dest="basedir", 
                    help="name of directory to save output files to.")
parser.add_argument("--popfile", "-p", type=str, nargs="*", dest="popfile", 
                    help="name of population output files [default: as trees but with .pop.tsv]")
parser.add_argument("--indivfile", "-i", type=str, nargs="*", dest="indivfile", 
                    help="name of individual output files [default: as trees but with .indiv.tsv]")
parser.add_argument("--freqfile", "-f", type=str, nargs="*", dest="freqfile", 
                    help="name of allele frequency output files [default: as trees but with .freq.tsv]")
parser.add_argument("--mut_rate", "-u", type=float, dest="mut_rate", 
                    help="mutation rate", default=1e-9)
parser.add_argument("--num_pops", "-n", type=int, dest="num_pops", 
                    help="Number of 'populations' to sample.")
parser.add_argument("--popsize", "-s", type=int, dest="popsize", 
                    help="Number of individuals in each 'population'.")
parser.add_argument("--logfile", "-g", type=str, dest="logfile", 
                    help="name of log file")
parser.add_argument("--seed", "-d", dest="seed", type=int, 
                    help="random seed", default=random.randrange(1,1000))

args = parser.parse_args()
argdict = vars(args)

if args.basedir is None and (args.indivfile is None or args.popfile):
    print(description)
    raise ValueError("Must specify at least basedir and indivfile (run with -h for help).")

if args.num_pops is None or args.popsize is None:
    print(description)
    raise ValueError("Must specify num_pops and popsize (run with -h for help).")

if args.tree_file is None or len(args.tree_file) == 0:
    args.tree_file = glob.glob(os.path.join(args.basedir, "*.trees"))

if args.popfile is None or len(args.popfile) == 0:
    args.popfile = [os.path.join(args.basedir, re.sub("[.]trees$", "", os.path.basename(x))) 
                        + ".pop.tsv" for x in args.tree_file]

if args.indivfile is None or len(args.indivfile) == 0:
    args.indivfile = [os.path.join(args.basedir, re.sub("[.]trees$", "", os.path.basename(x))) 
                        + ".indiv.tsv" for x in args.tree_file]

if args.freqfile is None or len(args.freqfile) == 0:
    args.freqfile = [os.path.join(args.basedir, re.sub("[.]trees$", "", os.path.basename(x))) 
                        + ".freq.tsv" for x in args.tree_file]

if args.logfile is None:
    args.logfile = os.path.join(args.basedir, "add_muts.log")

assert len(args.popfile) == len(args.tree_file)
assert len(args.indivfile) == len(args.tree_file)

logfile = open(args.logfile, "w")

random.seed(args.seed)
seeds = [random.randrange(1,10000000) for _ in range(len(args.tree_file))]

# typedef struct __attribute__((__packed__)) {
#   slim_pedigreeid_t pedigree_id_; // 8 bytes (int64_t): the SLiM pedigree ID for this individual, assigned by pedigree rec
#   slim_age_t age_; // 4 bytes (int32_t)
#   slim_objectid_t subpopulation_id_; // 4 bytes (int32_t)
# } IndividualMetadataRec;

class slimIndividual(object):
    def __init__(self, table_row):
        ped_id, age, subpop = struct.unpack("qii", table_row.metadata)
        location = table_row.location[:3]
        self.ped_id = ped_id
        self.age = age
        self.subpop = subpop
        self.location = location


header = ["id", "age", "subpop", "location_x", "location_y", "location_z", "node_ma", "node_pa"]


def write_info(chrom):
    treefile = args.tree_file[chrom]
    indiv_out = open(args.indivfile[chrom], "w")
    pop_out = open(args.popfile[chrom], "w")
    freq_out = open(args.freqfile[chrom], "w")
    logfile.write("Reading trees from " + treefile + "\n")
    ts = msprime.load(treefile)
    # add mutations
    mut_rate = args.mut_rate
    seed = seeds[chrom]
    logfile.write("Simulating mutations on" + treefile + "\n")
    tables = ts.dump_tables()
    rng = msprime.RandomGenerator(seed)
    mutgen = msprime.MutationGenerator(rng, mut_rate)
    mutgen.generate(tables.nodes, tables.edges, tables.sites, tables.mutations)
    ts = msprime.load_tables(**tables.asdict())
    # get individuals
    node_inds = [np.where(ts.tables.nodes.individual == u)[0] 
                 for u in range(ts.tables.individuals.num_rows)]
    individuals = [slimIndividual(ts.tables.individuals[k]) 
                   for k in range(ts.tables.individuals.num_rows)]
    ## get populations
    ind_x = np.array([ind.location[0] for ind in individuals])
    ind_y = np.array([ind.location[1] for ind in individuals])
    pop_lims = [(min(ind_x), max(ind_x)),
                (min(ind_y), max(ind_y))]
    pops = []
    allpops = []
    pop_locs = []
    j = 0
    while j < args.num_pops:
        pop_loc = (pop_lims[0][0] + np.random.random(1)[0] * (pop_lims[0][1] - pop_lims[0][0]),
                    pop_lims[0][0] + np.random.random(1)[0] * (pop_lims[0][1] - pop_lims[0][0]))
        dists = np.sqrt((ind_x - pop_loc[0])**2 + (ind_y - pop_loc[1])**2)
        pop = list(set(np.argpartition(dists, 2*args.popsize)[:2*args.popsize]) - set(allpops))[:args.popsize]
        if len(pop) > 0:
            pops.append(pop)
            pop_locs.append(pop_loc)
            allpops.extend(pop)
            j += 1
    # nodes corresponding to those individuals
    pop_nodes = [[u for j in pop for u in node_inds[j]] for pop in pops]
    allnodes = [u for x in pop_nodes for u in x]
    all_inds = [individuals[j] for j in allpops]
    # write out pop info
    logfile.write("Saving population info to" + args.popfile[chrom] + "\n")
    pop_out.write("\t".join(["pop", "x", "y", "num_genomes"]) + "\n")
    for k, (x, y) in enumerate(pop_locs):
        pop_out.write("\t".join(map(str,[k, x, y, len(pops[k])])) + "\n")
    pop_out.close()
    # write out indivdual info
    logfile.write("Saving individual info to" + args.indivfile[chrom] + "\n")
    indiv_out.write("\t".join(header) + "\n");
    for k, ind in enumerate(all_inds):
        data = [ind.ped_id, ind.age, ind.subpop] + list(ind.location) + list(node_inds[k])
        indiv_out.write("\t".join(map(str, data)) + "\n")
    indiv_out.close()
    # write out allele frequency info
    ## hack to remove individual info
    tables = ts.dump_tables()
    tables.nodes.set_columns(flags=tables.nodes.flags,
            time=tables.nodes.time, population=tables.nodes.population,
            metadata=tables.nodes.metadata, metadata_offset=tables.nodes.metadata_offset)
    sts = msprime.load_tables(**tables.asdict())
    # recover indices after simplification
    n = 0
    sub_pop_nodes = [[-1 for a in pop] for pop in pop_nodes]
    for k in range(len(sub_pop_nodes)):
        for j in range(len(sub_pop_nodes[k])):
            sub_pop_nodes[k][j] = n
            n += 1
    ## just so we can simplify
    sts = sts.simplify(allnodes)
    freq_out.write("\t".join(['pop_' + str(j) for _ in pops]) + "\n")
    for var in sts.variants():
        freqs = [sum(var.genotypes[nodes])/len(nodes) for nodes in sub_pop_nodes]
        freq_out.write("\t".join(map(str, freqs)) + "\n")
    freq_out.close()
    return True

for k in range(len(args.tree_file)):
    write_info(k)

logfile.write("Done!\n")
logfile.close()



