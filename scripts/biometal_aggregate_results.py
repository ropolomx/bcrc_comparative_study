#!/usr/bin/env python3

import sys

biomet_level_names = {0: 'Gene', 2: 'Class', 3: 'Mechanism', 1: 'Group'}


def load_megabio_annotations(file):
    annot = {}
    levels = {}
    with open(file, 'r') as f:
        data = f.read().split('\n')[1:]
        for entry in data:
            if not entry:
                continue
            if entry[0] == "\"":
                temp = entry.lstrip("\"").split('\"')[0]
                repl = temp.replace(',', '|')
                entry = ''.join([repl] + entry.lstrip("\"").split('\"')[1:])
            entry = entry.split(',')
            annot.setdefault(entry[0], entry[1:])
            for e, level in enumerate(entry):
                levels.setdefault(level, e)
    return annot, levels


def biomet_aggregate_results(file, A, L, N):
    long_results = {}
    sample_id = file.split('/')[-1].replace('_coverage_sampler_biomet.tab', '')
    with open(file, 'r') as f:
        data = f.read().split('\n')[1:]
        for result in data:
            if not result:
                continue
            result = result.split('\t')
            if result[2].count(',') > 0:
                result[2] = result[2].replace(',', '|')
            if result[2] not in A:
                continue
            for l in A[result[2]]:
                if l not in long_results:
                    long_results[l] = float(result[4])
                else:
                    long_results[l] += float(result[4])
            sys.stdout.write('Gene,{},{},{}\n'.format(
                sample_id,
                result[2],
                result[4]
            ))
    for name, count in long_results.items():
        sys.stdout.write('{},{},{},{}\n'.format(
            N[L[name]],
            sample_id,
            name,
            count
        ))

if __name__ == '__main__':
    A, L = load_megabio_annotations(sys.argv[1])
    biomet_aggregate_results(sys.argv[2], A, L, biomet_level_names)
