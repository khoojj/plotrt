#!/usr/bin/env python
# coding: utf-8
#Version 1.0

from Bio import SeqIO #call biopython to import seq info from fastq
import pandas as pd #call pandas for managing dataframe to allow us to sort based on variables
import matplotlib.pyplot as plt
import datetime, sys, os, argparse, glob, re

def import_seq_id_desc(input_file):
    seqid = []
    desc = []
    #import sequence id and description from fastq file into separate lists
    for seq_record in SeqIO.parse(input_file, "fastq"):
        seqid.append(seq_record.id)
        desc.append(seq_record.description)
    #create df with seqids and descriptions
    df_list = zip(seqid, desc) #combine both lists
    column_names = ["seqid", "description"] #set column names
    #insert each list into a column and add corresponding column name
    df = pd.DataFrame(df_list, columns = column_names)
    return df

def extract_and_sort_date_time(df):
    df["date_time"] = df.description.str.replace("T", " ").str.extract(r"(\d{4}\-\d{2}\-\d{2}\s\d{2}\:\d{2}\:\d{2})").replace("T", " ")
    df["date_time"] = pd.to_datetime(df.date_time) #convert from string to datetime data type
    #sort df based on datetime column in ascending order
    df_sorted = df.sort_values(by=['date_time']).reset_index(drop=True)
    return df_sorted

#create columns of values for x(time in hrs) and y(num of reads) axis
def create_x_y(df_sorted):
    time_list = [0]
    time_hrs = []
    for i in range(len(df_sorted)-1):
        time_list.append((df_sorted.loc[i+1, "date_time"] - df_sorted.loc[i, "date_time"]).total_seconds() + time_list[i])

    for t in time_list:
        time_hrs.append(t / 3600)

    read_count = list(i + 1 for i in list(range(len(df_sorted))))
    return time_hrs, read_count

#plot x and y values on scatter plot
def plot_scatter(time_reads, out_path, output_file_name):
    output = out_path + "/" + output_file_name + ".png"
    plt.figure(figsize=(12,8))
    plt.plot(time_reads[0], time_reads[1])
    plt.xlabel("hours")
    plt.ylabel("num of reads", labelpad =12)
    plt.savefig(output)

#function to produce one plot for each fastq file
def plot_read_vs_time(input_file, out_path, output_file_name):
    df = import_seq_id_desc(input_file)
    df_sorted = extract_and_sort_date_time(df)
    time_reads = create_x_y(df_sorted)
    plot_scatter(time_reads, out_path, output_file_name)

#create a directory to hold all output files
def create_out_dir(out_dir):
    parent_dir = os.getcwd()
    directory = out_dir
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    out_path = parent_dir + "/" + out_dir
    return out_path

#function to loop through all input files provided as fastqs in current directory
#or path to directory containing fastqs to be processed
def repeat(input_files, out_dir):
    start_time = datetime.datetime.now()
    out_path = create_out_dir(out_dir)
    for file in input_files:
        if ".fastq" in file:
            output_name = file.replace(".fastq", "")
            plot_read_vs_time(file, out_path, output_name)
        else:
            path = input_files[0]
            files = [f for f in glob.glob(path + "**/*.fastq", recursive=True)]
            for f in files:
                output_name = re.findall("\/.+\/(.+)\.fastq", f)
                plot_read_vs_time(f, out_path, output_name[0])
    end_time = datetime.datetime.now()
    print('Process complete. Total processing time: {}'.format(end_time - start_time))

#use argparser to set inputs from command line
parser = argparse.ArgumentParser(
description='''Sort reads in fastq produced from Guppy basecalling based on timestamp and
produce scatter plot for num of reads produced vs time in hours.
Plots will be saved in a new directory created in the current directory.
''')

parser.add_argument('--input_files', '-i',
                       type=str,
                       action='store',
                       nargs='+',
                       help='''input fastq file in current directory, can be single file name
                               or multiple files given as *.fastq, or path to directory containing
                               all the input fastq files (do not include the file names in path),
                               example: 'file1.fastq file2.fastq' or '*.fastq'
                                         or '/home/username/directory' ''')

parser.add_argument('--out_dir', '-o',
                       type=str,
                       action='store',
                       help='''name for directory created to hold all output files,
                               do not include the entire path,
                               example: 'out_folder' ''')

args = parser.parse_args()

if args.input_files:
    repeat(args.input_files, args.out_dir)
