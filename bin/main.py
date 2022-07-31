#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-07-30
Purpose: To run an in-silico PCR using bacterial genomes.
"""

# argument parsing
import argparse
from typing import NamedTuple

# required modules
from Bio import Entrez, SeqIO
import pathlib
import time
import subprocess
import shutil

class Args(NamedTuple):
    """ Command-line arguments """
    init: bool
    run: bool
    name: str
    overwrite: bool

    db_fetch_num: int

    f_seq: str
    r_seq: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        usage='%(prog)s [--init/--run] [name] [options]',
        description='Fetch sequences from NCBI Entrez with a list of search queries',
        formatter_class=argparse.MetavarTypeHelpFormatter)

    group_main = parser.add_mutually_exclusive_group(required=True)
    group_main.add_argument(
        '--init',
        action='store_true',
        help="initialize 16S database")
    group_main.add_argument(
        '--run',
        action='store_true',
        help="run PCR with a specified database")

    parser.add_argument(
        'name',
        metavar='name',
        type=str,
        default=None,
        help="name of database to either initialize (--init) or run (--run)")


    group_database = parser.add_argument_group(
        title='database options',
        description='use with --init')
    group_database.add_argument(
        '--overwrite',
        action='store_true',
        help="overwrite an existing database")
    group_database.add_argument(
        '-n',
        dest='db_fetch_num',
        metavar='<int>',
        type=int,
        default=10,
        help="number of sequences to fetch (default=1000)")

    group_ipcr = parser.add_argument_group('PCR options')
    group_ipcr.add_argument(
        '-f',
        dest='f_seq',
        metavar='FWD',
        type=str,
        default=None,
        help="How deep sohuld this database go")
    group_ipcr.add_argument(
        '-r',
        dest='r_seq',
        metavar='REV',
        type=str,
        default=None,
        help="How deep sohuld this database go")

    args = parser.parse_args()

    if args.init: print_runtime(f'Running prokaryotic-ipcr with the --init setting.')
    if args.run: print_runtime(f'Running prokaryotic-ipcr with the --run setting.')

    # if initializing the database and it already exists, make sure overwrite flag is active
    if all((args.init, pathlib.Path.cwd().joinpath('db').joinpath(args.name).exists(), not args.overwrite)):
        parser.error('Database already exists. Use --overwrite to overwrite this database.')

    # if running the database, make sure the forward and reverse primer sequences are specified
    if args.run and not (args.f_seq and args.r_seq):
        parser.error('Specify the forward and reverse primer sequences with -f and -r, respectively. ')

    # if primers are specified, make sure that they're actually
    allowed_chars = 'ACGTUWSMKRYBDHVN'
    if (args.f_seq and args.r_seq):
        for nucleotide in (args.f_seq.upper() + args.r_seq.upper()):
            if nucleotide not in allowed_chars:
                parser.error(f'Invalid character {nucleotide} in primer sequences.')

    return Args(args.init, args.run, args.name, args.overwrite, args.db_fetch_num, args.f_seq, args.r_seq)
# --------------------------------------------------
def main() -> None:
    args = get_args()
    home=pathlib.Path.cwd()

    # create_database directory
    db_output = home.joinpath('db')
    db_output.mkdir(parents=True, exist_ok=True)
    
    if args.init:
        named_db=db_output.joinpath(args.name)
        create_db(named_db, args.db_fetch_num)

    elif args.run:
        named_db=db_output.joinpath(args.name)
        total_count = total_bacterial_count = total_archaeal_count = 0
        iter_represented_accession_count = {'Bacteria': 0, 'Archaea': 0}
        iter_total_represented_accession_count = {'Bacteria': 0, 'Archaea': 0}
        actual_represented_accession_count = {'Bacteria': 0, 'Archaea': 0}
        actual_total_accession_count = {'Bacteria': 0, 'Archaea': 0}
        
        # realistically 10?
        for mismatch_iter in range(0, 7):
            print_runtime(f'Iteratively allowing more mismatches in the primer sequence... (Iteration {mismatch_iter+1})')
            iPCR_results = iPCR(named_db, (args.f_seq, args.r_seq), mismatch_iter)

            total_count += iPCR_results[0]
            total_bacterial_count += iPCR_results[1]
            total_archaeal_count += iPCR_results[2]

            for taxon in iter_represented_accession_count:
                iter_represented_accession_count[taxon] += len(iPCR_results[3][taxon])
                iter_total_represented_accession_count[taxon] += iPCR_results[4][taxon]
                actual_total_accession_count[taxon] = iPCR_results[4][taxon]
                #print(iter_represented_accession_count[taxon], iter_total_represented_accession_count[taxon])
                try:
                    if len(iPCR_results[3][taxon]) > actual_represented_accession_count[taxon]:
                        actual_represented_accession_count[taxon] = len(iPCR_results[3][taxon])
                except: pass
        print_runtime("Completed job. \n")
        print_report("Population", f"Total: {actual_total_accession_count['Bacteria']+actual_total_accession_count['Archaea']}\tBacteria: {actual_total_accession_count['Bacteria']}\tArchaea: {actual_total_accession_count['Archaea']}")
        print_report("Total 'reads'", f'Total: {total_count}\tBacteria: {total_bacterial_count}\tArchaea: {total_archaeal_count}')
        for taxon in iter_represented_accession_count:
            approx_efficiency = iter_represented_accession_count[taxon]/iter_total_represented_accession_count[taxon]
            percent_represented = actual_represented_accession_count[taxon]/actual_total_accession_count[taxon]
            print_report(f"{taxon}", f'Approximate efficiency: {round(approx_efficiency*100, 2)}\tPercent representation of accessions: {round(percent_represented*100, 2)} % ({actual_represented_accession_count[taxon]}/{actual_total_accession_count[taxon]})')
# --------------------------------------------------'
def create_db(db_name_arg: pathlib.Path, db_fetch_num_arg: int) -> None:
    """
    Create a named database and fetch a number of sequences.

    Parameters:
        db_name_arg (path): the path to the named database
        db_fetch_num_arg (int): the number of sequences to fetch'

    Returns:
        None
    """
    #query='(16s[All Fields] AND rRNA[All Fields] NOT "partial sequence"[All Fields] NOT "operon"[All Fields]) AND 00000000001[SLEN] : 00000005000[SLEN] AND ("Bacteria"[Primary Organism] OR "Archaea"[Primary Organism])'
    #query='(16s[All Fields] AND rRNA[All Fields] NOT "partial sequence"[All Fields] NOT "operon"[All Fields]) AND 00000000001[SLEN] : 00000005000[SLEN] AND ("Bacteria"[Primary Organism] OR "Archaea"[Primary Organism]) AND ("cow"[All Fields] OR "rumen"[All Fields] OR "ruminant"[All Fields])'
    #query='(16s[Title] AND rRNA[Title] NOT "partial sequence"[All Fields] NOT "operon"[All Fields]) AND 00000000001[SLEN] : 00000005000[SLEN] AND ("Bacteria"[Primary Organism] OR "Archaea"[Primary Organism]) AND "ruminant"[All Fields] NOT "whole genome"[All Fields]'
    query='(biomol_genomic[PROP] AND 00000001000[SLEN] : 00000005000[SLEN]) AND ("bacteria"[Primary Organism] OR "archaea"[Primary Organism]) AND (“16S ribosomal RNA”[All Fields] or “SSU”[All Fields]) NOT ("whole genome"[All Fields] OR "assembly"[All Fields] OR "partial"[All Fields])'
    #query='("16S"[title] AND "rRNA"[title] NOT "partial sequence"[All Fields] NOT "operon"[All Fields]) AND 00000001000[SLEN] : 00000005000[SLEN] AND (((("Bacteria"[Primary Organism] OR "Bacteria Latreille et al. 1825"[Primary Organism]) OR "Bacteria Latreille et al. 1825"[Primary Organism]) OR "Bacteria Latreille et al. 1825"[Primary Organism]) OR "Archaea"[Primary Organism]) AND biomol_genomic[PROP]'

    # create directory for the named database
    if db_name_arg.exists(): 
        shutil.rmtree(db_name_arg)
        print_runtime(f"Overwrote database db/{db_name_arg.name}.")
    db_name_arg.mkdir(parents=True, exist_ok=True)
    db_name_arg.joinpath('fasta').mkdir(parents=True, exist_ok=True)
    print_runtime(f"Created database db/{db_name_arg.name}. Now attempting to retrieve {db_fetch_num_arg} sequences.")

    # Entrez wants your email
    Entrez.email='erick.samera@kpu.ca'

    # perfom Entrez esearch
    handle = Entrez.esearch(db="nucleotide", term=query, idtype="acc", retmax=db_fetch_num_arg, sort='relevance')
    record = Entrez.read(handle)

    running_counts = {'Bacteria': 0, 'Archaea': 0, 'Methanogens': 0}

    # variables for the progress bar
    bar_width = 15
    step = (1/bar_width)
    current_value = 0
    max_value = record["RetMax"]
    current_progress_step = 0 + step

    # try to retrieve sequences via accession number
    with open(db_name_arg.joinpath('taxonomy_map.txt'), 'w') as taxonomy_file:
        print_runtime(f"Retrieving sequences: [{''.join(['-']*bar_width)}] {round(current_value/int(max_value)*100)}%")
        for accession_num in record['IdList']:
            handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="gb", retmode="text")

            # progress bar
            current_value += 1
            #print(current_value, current_value/int(max_value), current_progress_step)
            if  (current_value/int(max_value)) >= current_progress_step:
                progress_bar = ''.join(['*' for i in range(int(current_progress_step*bar_width))])
                unprogress_bar = ''.join(['-' for i in range(int(bar_width-len(progress_bar)))])
                print_runtime(f'Retrieving sequences: [{progress_bar+unprogress_bar}] {round(current_value/int(max_value)*100)}%')
                current_progress_step += step

            try:
                parsed_handle = SeqIO.read(handle, "gb")

                if parsed_handle.annotations['taxonomy']:
                    try:
                        # running counts of bacterial and archaeal sequences
                        running_counts[parsed_handle.annotations['taxonomy'][0]] += 1
                        # if archaea, check whether it's a methanogen 
                        if parsed_handle.annotations['taxonomy'][0] == 'Archaea' and any(['meth' in taxon.lower() for taxon in parsed_handle.annotations['taxonomy']]): 
                            running_counts['Methanogens'] += 1

                        try: taxonomy_file.write(f"{parsed_handle.annotations['accessions'][0]}.{parsed_handle.annotations['sequence_version']}\t{parsed_handle.annotations['taxonomy']+[parsed_handle.annotations['organism']]}\n")
                        except: print_runtime(f'Failed to write to assession {accession_num} to taxonomy_map.txt !')
                        try: SeqIO.write(parsed_handle, db_name_arg.joinpath('fasta').joinpath(f'{accession_num}.fasta'), 'fasta')
                        except: print_runtime(f'Failed to write assession {accession_num} FASTA !')

                    except: print_runtime(f'Failed to process taxonomy for accession {accession_num} !')
            except: print_runtime(f'Failed to retrieve accession {accession_num} !')
    #print(running_counts)
    print_runtime(f"Created db/{db_name_arg.name}/taxonomy_map.txt")
    print_runtime(f"Retrieved {running_counts['Bacteria']+running_counts['Archaea']} total sequences: {running_counts['Bacteria']} bacterial; {running_counts['Archaea']} archaeal; specifically {running_counts['Methanogens']} methanogens.")
def iPCR(db_name_arg, primers_arg: tuple, mismatch_arg: int = 0):
    """
    Function will perform in-silico PCR with ipcress and output products.

    Parameters:
        db_name_arg (Path?): the path to the named database to run in-silico PCR
        primers_arg (tuple): the forward and reverse sequences of the primers to use for in-silico PCR
        mismatch_arg (int): the number of mismatches to allow for PCR amplification

    Returns:
        lots
    """

    # ambiguous nucleotides are not considered matches? 
    #print(sum([sum([primer.count(ambiguous_nucleotide) for ambiguous_nucleotide in 'UWSMKRYBDHVN']) for primer in primers_arg]))

    # variables to declare for something?
    represented_accessions = {'Bacteria': [], 'Archaea': []}
    total_accession_count = {'Bacteria': 0, 'Archaea': 0}
    count_bacteria = count_archaea = 0
    size_filter = (0, 600)

    with open(db_name_arg.joinpath('taxonomy_map.txt'), 'r') as taxonomy_map_file:

        # process the file again
        taxonomy_data = {}  
        taxonomy_map = [line for line in taxonomy_map_file.readlines()]
        for result in taxonomy_map:
            accession = result.split('\t')[0]
            taxonomy_data[accession] = [taxon.strip()[1:-1] for taxon in result.split('\t')[1].strip()[1:-1].split(',')]
            total_accession_count[taxonomy_data[accession][0]] += 1

        with open(db_name_arg.parent.joinpath('ipcr_file.ipcress'), 'w') as ipcr_file:
            ipcr_file.write(f'ID0001 {primers_arg[0]} {primers_arg[1]} 0 600')

        args = [
            'ipcress',
            '--input', 'db/ipcr_file.ipcress',
            '--sequence', f"{db_name_arg.joinpath('fasta')}",
            '--mismatch', str(mismatch_arg),
            '--pretty', 'False',
            '--products', 'True',
            '--seed', '15s',
            '--memory', '24000'
        ]

        iPCR_results = subprocess.run(args, check=False, capture_output=True)
        PCR_products_list = str(iPCR_results.stdout).split('>')[1:]
        for PCR_result in [i.split('\\n') for i in PCR_products_list]:
            PCR_assession = PCR_result[0].split(' ')[2].split(':')[0]
            PCR_sequence = ''.join(PCR_result[1:-2])
            if (size_filter[0] < len(PCR_sequence) < size_filter[1]):
                if taxonomy_data[PCR_assession][0] == 'Bacteria': count_bacteria += 1
                if taxonomy_data[PCR_assession][0] == 'Archaea': count_archaea += 1
                represented_accessions[taxonomy_data[PCR_assession][0]].append(PCR_assession)
    for taxon in represented_accessions:
        represented_accessions[taxon] = sorted(set(represented_accessions[taxon]))

    return count_archaea+count_bacteria, count_bacteria, count_archaea, represented_accessions, total_accession_count
def print_report(heading, action) -> None:
    print(f'{heading}:\t\t{action}')
def print_runtime(action) -> None:
    """Prints an action with the time it was performed."""
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')
# --------------------------------------------------
if __name__ == '__main__':
    main()