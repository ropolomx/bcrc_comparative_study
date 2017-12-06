#!/usr/bin/env python3

import sys

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


if __name__ == '__main__':
    A, L = load_megabio_annotations(sys.argv[1])
    for k, v in A.items():
        sys.stdout.write('{},{},{},{}\n'.format(k, v[1], v[2], v[0]))

