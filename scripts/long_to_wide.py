#!/usr/bin/env python3

## sys.argv:
#   1. AMR_aggregated_output.csv
#   2. nextflow_output folder

import sys
import glob
import numpy as np

taxa_level = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5, 'S': 6}
taxa_level_names = {1: 'Domain', 2: 'Phylum', 3: 'Class', 4: 'Order',
                    5: 'Family', 6: 'Genus', 7: 'Species', 8: 'Unclassified'}
amr_level_names = {0: 'Class', 1: 'Mechanism', 2: 'Group'}


def dict_to_matrix(D):
    ncol = len(D.keys())
    unique_nodes = []
    samples = []
    for sample, tdict in D.items():
        for taxon in tdict.keys():
            if taxon not in unique_nodes:
                unique_nodes.append(taxon)
    nrow = len(unique_nodes)
    ret = np.zeros((nrow, ncol), dtype=np.float)
    for j, (sample, tdict) in enumerate(D.items()):
        samples.append(sample)
        for i, taxon in enumerate(unique_nodes):
            if taxon in tdict:
                ret[i, j] = np.float(tdict[taxon])
    return ret, unique_nodes, samples


def amr_load_long_data(file):
    samples = {}
    labels = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')[1:]
        for entry in data:
            if not entry:
                continue
            entry = entry.split(',')
            try:
                samples[entry[0]][entry[1]][entry[2]] = float(entry[3])
            except KeyError:
                try:
                    samples[entry[0]][entry[1]].setdefault(entry[2], float(entry[3]))
                except KeyError:
                    try:
                        samples[entry[0]].setdefault(entry[1], {entry[2]: float(entry[3])})
                    except KeyError:
                        samples.setdefault(entry[0], {entry[1]: {entry[2]: float(entry[3])}})
            try:
                if entry[2] not in labels[entry[0]]:
                    labels[entry[0]] += (entry[2],)
            except KeyError:
                labels.setdefault(entry[0], (entry[2],))
    return samples, labels


def output_amr_analytic_data(outdir, S, L):
    with open(outdir, 'w') as amr:
        for flevel, sdict in S.items():
            local_sample_names = []
            for sample in sdict.keys():
                local_sample_names.append(sample)
            if flevel == 'Gene':
                amr.write(','.join(local_sample_names) + '\n')
            for label in L[flevel]:
                local_counts = []
                if flevel == 'Gene':
                    amr.write(label + ',')
                for sample in local_sample_names:
                    if label in sdict[sample]:
                        local_counts.append(str(sdict[sample][label]))
                    else:
                        local_counts.append(str(0))
                if flevel == 'Gene':
                    amr.write(','.join(local_counts) + '\n')


if __name__ == '__main__':
    S, L = amr_load_long_data(sys.argv[1])
    output_amr_analytic_data(sys.argv[2], S, L)

