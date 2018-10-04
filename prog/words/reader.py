#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io, os, sys, glob
from math import *

blksz = 1000


def tokenize(fn):
    s = io.open(fn, encoding="utf-8").read()
    s = s.lower()
    punctuations = u"\"',;.:!?()[]{}-_—`“”‘’…£°˙/\\&\n\t\r"
    for c in punctuations:
        s = s.replace(c, " ")
    arr = s.split()
    for i in range(len(arr)):
        w = arr[i]
        if not w.isalnum():
            print w.encode("utf-8"), arr[i-1], arr[i+1], arr[i+2], arr[i+23], arr[i+24]
            raw_input()
    return arr


def num_tokens(arr, fnseq):
    tab = {}
    n = len(arr)
    for i in range(n):
        w = arr[i]
        if not w in tab:
            tab[w] = 0
        tab[w] += 1

    wmap = {}
    i = 0
    ent = 0
    for w in tab:
        freq = tab[w]
        x = 1.0 * freq / n
        ent -= x * log(x)
        wmap[w] = i
        i += 1
    print "%s words mapped, entropy %s" % (i, ent)

    return wmap, ent


def traj(arr, blksz, fnout):
    n = len(arr)
    nb = n / blksz
    fp = open(fnout, "w")
    for ib in range(nb):
        ent, tab = get_entropy(arr[0:(ib+1)*blksz])
        print (ib+1)*blksz, ent, "\r",
        fp.write("%s %s\n" % ((ib+1)*blksz, ent))
    fp.close()


def analyze(folder, blksz):
    # get the name of the folder
    if folder.endswith("/"):
        folder = folder[:-1]
    name = os.path.split(folder)[-1]

    # set the output files
    fnseq = os.path.join(folder, name + ".seq")
    fnout = name + ".log"

    # get a list of files to analyze
    ls = glob.glob(os.path.join(folder, "*.txt"))
    warr = []
    merged = []
    for fn in ls:
        wseq = tokenize(fn)
        warr += [wseq,]
        merged += wseq

    # convert words to numbers
    wmap, ent = num_tokens(merged, fnseq)

    narr = []
    for i in range(len(warr)):
        wseq = warr[i]
        nseq = [wmap[w] for w in wseq]
        narr += [nseq,]

    nwords = len(wmap)
    print name, ent, len(warr), len(merged), nwords, "words"

    slog = "# %s %s\n" % (nwords, len(warr))
    for i in range(len(narr)):
        nseq = narr[i]
        fn = ls[i]
        slog += "# %s %s %s\n" % (i, len(nseq), fn)
        for j in range(len(nseq)):
            slog += "%s\n" % nseq[j]
        slog += "\n"
    print "saving sequence to %s" % fnseq
    open(fnseq, "w").write(slog)


if __name__ == "__main__":
    folder = "../../data/words/shakespeare"
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    analyze(folder, blksz)

