#!/usr/bin/env python3

__prog__ = "fasta2reads"
__version__ = "0.0.1"
__date__ = "2016/03/25"
__author__ = "Kevin Drees"
__email__ = "kevin.drees@unh.edu"

'''

Creates artificial reads from a fasta sequence.
No attempt is made to simulate actual sequence output.
Quality at each locus is set to Q40. 

version history:

0.0.1 2016/03/25 New

'''

import logging
import os


def _parse_args():
    import argparse
    from fasta2reads import __version__ as fasta2reads_version

    parser = argparse.ArgumentParser(description="fasta2reads, version %s" % fasta2reads_version)
    parser.add_argument("fasta", 
                        default="", 
                        help="Path to the input fasta") 
    parser.add_argument("-r",
                        default="100",
                        type=int,
                        help="Simulated read length")
    parser.add_argument("-d",
                        default="20",
                        type=int,
                        help="Depth of read coverage")
    return parser.parse_args()


def _parse_filename(fasta):
    import re
    match_object = re.search(r"(.*)\.f[ast|n]?a", fasta)
    if match_object:
        prefix = match_object.group(1)
    else:
        print("Cannot parse input filename %s" % ( fasta ))
        raise SystemExit
    return prefix


def _create_read(seq, header):
    fake_qual = "I"
    seq = seq.replace("\n",'')
    length = len(seq)
    newseq = ''
    qual = ''
    for i in range(1, length):
        newseq+= seq[i]
        qual += fake_qual
    read = "%s\n%s\n+\n%s\n" % (header, newseq, qual)
    return read  


def _app2gz(outfile, string):
    import gzip

    bstring = bytearray(string, 'utf-8')
    with gzip.open(outfile, 'a') as stream:
        stream.write(bstring)

def _write_read(seq, header, outfile, read_id):
    header += "_%i" % (read_id)
    read_id = read_id + 1
    read = _create_read(seq, header)
    _app2gz(outfile, read) 
    return read_id
    

def main():
    import re

    cli_args = _parse_args()
    fasta = cli_args.fasta
    fasta_path = os.path.abspath(fasta)    
    readlength = cli_args.r
    depth = cli_args.d
    slide = int(readlength / depth)
    prefix = _parse_filename(fasta)
    outfile = "%s.fastq.gz" % prefix
    try:
        os.remove(outfile)
    except(Exception):
        pass    
    fastq_header = "@%s_simulated_read" % prefix 
    read_id = 1
    with open(fasta_path, 'r') as fasta_stream:
        seq = ''
        while True:
            curr_seq_len = len(seq)
            if curr_seq_len < readlength:
               line = fasta_stream.readline()
               line = line.replace("\n","")
               if line:        
                   if re.match('>', line) is None:
                       seq += line
                   elif curr_seq_len > 0:
                       read_id = _write_read(seq, fastq_header, outfile, read_id)
                       seq = ''
               elif curr_seq_len > 0:
                   read_id = _write_read(seq, fastq_header, outfile, read_id)
                   break
               else:
                   break
            else:
                start = 0
                end = readlength
                trimmings = seq[start:end]
                read_id = _write_read(trimmings, fastq_header, outfile, read_id)
                start = start + slide  
                seq = seq[start:]
 
'''
old version
with open(fasta_path, 'r') as fasta_stream:
        seq = ''
        while True:
           curr_seq_len = len(seq)
           if curr_seq_len < readlength:
               line = fasta_stream.readline()
               if line:        
                   if re.match('>', line) is None:
                        seq += line 
               elif curr_seq_len > 0:
                   read_id = _write_read(seq, fastq_header, outfile, read_id)
                   break
               else:
                   break
           else:
               trimmings, seq = seq[:readlength], seq[readlength:]
               read_id = _write_read(trimmings, fastq_header, outfile, read_id)
'''

if __name__ == "__main__":
    main()
