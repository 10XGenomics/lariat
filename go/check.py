#!/usr/bin/env python

import tenkit.bam
import subprocess
import os
import sys
import numpy as np
import itertools

rfa = os.path.join(os.environ["GOPATH"], "bin", "RFA")
test_fq = os.path.join(os.environ["GOPATH"], "src/test/inputs", "segdups.fastq.gz")

hg19 = os.path.join(os.environ['TENX_REFDATA_PATH'], "fasta/hg19/hg19.fa")

if len(sys.argv) < 2:
    print "Running RFA"
    args = [rfa, "-reads", test_fq, "-genome", hg19]
    print args
    subprocess.check_call(args)
    bam_fn = "default.bam"
else:
    print "Skipping RFA. Analyzing %s" % sys.argv[1]
    bam_fn = sys.argv[1]
    try:
        subprocess.check_call("go tool pprof -dot  bin/RFA cpu.pprof  > prof.dot", shell=True)
        subprocess.check_call("dot -Tpng cpu.dot  > dot.png", shell=True)
    except Exception:
        print "failed making profiling graph. continuing"


#subprocess.check_call(['samtools','sort','default.bam','default_sorted'])
#subprocess.check_call(['samtools','index','default_sorted.bam'])
bam = tenkit.bam.create_bam_infile(bam_fn)
reads = list(bam)

logf = open("check_log.txt", "a")
date = subprocess.check_output("date")
git = subprocess.check_output("git describe --tags --dirty", shell=True)
log_lines = ['', '----------', date, git]

def show_frac(label, r, f):
    n_match = len([x for x in r if f(x)])
    frac = float(n_match) / len(r)
    s = "{0:15}: {1:3f}".format(label, frac)
    log_lines.append(s)

show_frac("Unmapped", reads, lambda x: x.is_unmapped)
show_frac("Proper pair", reads, lambda x: x.is_proper_pair)
show_frac("mapq = 0", reads, lambda x: x.mapq == 0)
show_frac("mapq < 30", reads, lambda x: x.mapq < 30)
show_frac("mapq >= 30", reads, lambda x: x.mapq >= 30)



def correct_alignment(bam,r):
    # qname = mol:blah:chr3:number:number:80771341:80771523
    parts = r.qname.split(":")
    al_pos = int(parts[5])
    return parts[2] == bam.references[r.tid] and abs(r.pos - al_pos) < 200


def round_closest(val, opts):
    diff = np.abs(val - opts)
    return opts[np.argmin(diff)]

def report_mapq(res):
    mapqs = np.array([x[1] for x in res])
    correct = np.array([int(x[2]) for x in res])
    return {'med_map': np.median(mapqs), 'emp_mapq': -10.0*np.log10(1.0-correct.mean()), 'n': len(res)}


def analyze_mapqs(bam):
    opts = np.array([5, 15, 30, 45])
    obs = []
    bam.reset()
    for r in bam:
        rqv = round_closest(r.mapq, opts)
        obs.append((rqv, r.mapq, correct_alignment(bam, r)))

    obs.sort()

    print obs[:5]
    print obs[-5:]

    # Split into 5 groups and report mapqs
    results = []
    for (k, vals) in itertools.groupby(obs, lambda x: x[0]):
        stats = report_mapq(list(vals))
        results.append(stats)

    return results



rname_parts = reads[0].qname.split(":")

if len(rname_parts) == 7 and rname_parts[0] == 'mol':
    r = analyze_mapqs(bam)
    print r
    for l in r:
        log_lines.append("{}".format(l))


logf.writelines(log_lines)
for l in log_lines:
    print l
