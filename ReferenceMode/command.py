import os
import multiprocessing
from subprocess import Popen, PIPE, STDOUT


def runCommand(command):
    #p = Popen(command + " > log 2>&1", stdin=PIPE, stdout=PIPE, shell=True).communicate()
    p = Popen(command, stdin=PIPE, stdout=PIPE, shell=True).communicate()

def run_bowtie2(repeat_contig_path, reference_path, thread_num, tools_dir):
    (dir, filename) = os.path.split(repeat_contig_path)
    (name, extension) = os.path.splitext(filename)
    sam_path = dir + '/' + name + '.sorted.sam'
    (ref_dir, ref_filename) = os.path.split(reference_path)
    (ref_name, ref_extension) = os.path.splitext(ref_filename)
    index_path = ref_dir + '/' + ref_name
    bowtie2_index_command = tools_dir + '/bowtie2-build ' + reference_path + ' ' + index_path
    print(bowtie2_index_command)
    runCommand(bowtie2_index_command)
    bowtie2_align_command = tools_dir + '/bowtie2 -p ' + str(thread_num) + ' -f --very-fast -a -x ' + index_path + ' -U ' + repeat_contig_path + ' -S ' + sam_path
    print(bowtie2_align_command)
    runCommand(bowtie2_align_command)
    return sam_path

def run_bwa(repeat_contig_path, reference_path, thread_num, tools_dir):
    repeat_contig_dir = os.path.abspath(os.path.dirname(repeat_contig_path))
    (dir, filename) = os.path.split(repeat_contig_path)
    (name, extension) = os.path.splitext(filename)
    sam_path = repeat_contig_dir + '/' + name + '.sam'

    if not os.path.exists(reference_path+'.bwt'):
        bwa_index_command = tools_dir + '/bwa index -a bwtsw ' + reference_path
        print(bwa_index_command)
        runCommand(bwa_index_command)
    bwa_align_command = tools_dir + '/bwa mem -t ' + str(thread_num) + ' -x intractg ' + reference_path + ' -a ' + repeat_contig_path + ' > ' + sam_path
#    '/bwa mem -t 48 -x intractg /public/home/hpc194701009/Ref/dmel-all-chromosome-r5.43.fasta -a /public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/no_simple_repeats/dmel_curated.fasta > /public/home/hpc194701009/KmerRepFinder_test/library/curated_lib/no_simple_repeats/dmel_curated.sam'
    print(bwa_align_command)
    runCommand(bwa_align_command)
    return sam_path

def run_minimap2(repeat_contig_path, reference_path, thread_num, tools_dir):
    (dir, filename) = os.path.split(reference_path)
    (name, extension) = os.path.splitext(filename)

    index_file = dir + '/' + name + '.mmi'
    minimap2_index_command = tools_dir + '/minimap2 -d ' + index_file + ' ' + reference_path
    print(minimap2_index_command)
    runCommand(minimap2_index_command)

    (dir, filename) = os.path.split(repeat_contig_path)
    (name, extension) = os.path.splitext(filename)
    sam_path = dir + '/' + name + '.sam'
    minimap2_align_command = tools_dir + '/minimap2 -a -t ' + str(thread_num) + ' ' + index_file + ' ' + repeat_contig_path + ' > ' + sam_path
    print(minimap2_align_command)
    runCommand(minimap2_align_command)
    return sam_path


def run_minimap2_self_align(repeat_contig_path, tools_dir):
    (dir, filename) = os.path.split(repeat_contig_path)
    (name, extension) = os.path.splitext(filename)
    sam_path = dir + '/' + name + '.self_align.sam'
    minimap2_align_command = tools_dir + '/minimap2 -ax map-pb --eqx -t ' + str(multiprocessing.cpu_count()) + ' ' + repeat_contig_path + ' ' + repeat_contig_path + ' > ' + sam_path
    print(minimap2_align_command)
    runCommand(minimap2_align_command)
    return sam_path