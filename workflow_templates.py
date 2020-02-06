# Author: adapted from Skov/Lucotte by CMK.
# common for X and Y chromosomes.
# The templates are ordered by use in the protocol.

import os



# 0: initialize

def initialize(title, description): # lav den her funktion i python i stedet for.
    inputs = []
    outputs = [title + "/batch_description.txt"] # Jeg forstår ikke hvornår den her fil bliver lavet?
    options = {'cores': 1, 'memory': '1g', 'walltime': '00:01:00', 'account': "simons"}
    spec = '''
        mkdir -p {title}
        cd {title}
        echo $(date)'\n{description}\n' >> batch_description.txt
        mj > mj.txt'''.format(
            title = title,
            description = description)
    return inputs, outputs, options, spec    


# 1: Index the ref. genome
def index_genome(title, refgenome):
    refgenome_stem = os.path.splitext(refgenome)[0]
    
    inputs = [refgenome]
    outputs = [title+'/'+refgenome_stem+extension for extension in ['.amb', '.ann', '.pac', '.sa', '.bwt', '.fa']]
    #options = {'cores': 1, 'memory': 2000, 'walltime': '00:25:00', 'account': "simons"}
    options = {'cores': 1, 'memory': '10g', 'walltime': '10:00:00', 'account': "simons"}
    
    #options = {}
    spec =  """sleep 1; cd {title}; cp ../{refgenome} .

#source /com/extra/bwa/0.7.5a/load.sh
bwa index -p {refgenome_stem} -a bwtsw {refgenome}""".format(title=title, refgenome_stem=refgenome_stem, refgenome=refgenome)
    #spec =  """mkdir {title}; cd {title}; touch test; """.format(title=title)

    return inputs, outputs, options, spec


#2
def bwa_map_pe(title, refgenome, read1, read2, subject):
    '''
    Maps the hap. genomes to the reference

    1st sambamba view -F : Filters for the paired reads
    2nd sambamba view -F : specify the reference genome for writing the output
    sambamba sort : sorts the reads per coordinate (default)
    sambamba markdup : removes the duplicate reads in bam file
    sambamba flagstat : gives a flag stat to BAM files
    sambamba view -F : filter the bam files
    '''
    refgenome_stem = os.path.splitext(refgenome)[0]
    inputs = [read1,
              read2,
              refgenome]

    for extension in ['.amb', '.ann', '.pac']: inputs.append(title+'/'+refgenome_stem+extension) # kan godt flattenes pænere
    # ooutputs = [
    #   title+'/'+subject+'_sorted.bam',
    #   title+'/'+subject+'_sorted.bam.bai',
    #   title+'/'+subject+'_unsorted.bam',
    #   title+'/'+subject+'_sort_dedup.bam',
    #   title+'/'+subject+'_sort_dedup.bam.bai',
    #   title+'/'+subject+'_sort_dedup.bam.flagstat',
    #   title+'/'+subject+'.COMPLETED']

    outputs = [title+'/'+subject+extension for extension in [
        '_sorted.bam',
        '_sorted.bam.bai',
        '_unsorted.bam',
        '_sort_dedup.bam',
        '_sort_dedup.bam.bai',
        '_sort_dedup.bam.flagstat']]

    #options = {'cores': 16, 'memory': 16000, 'walltime': '24:00:00', 'account': "simons"} # changed memory to 16, # back to 12, and with double time (24h)
    options = {'cores': 16, 'memory': 64000, 'walltime': '120:00:00', 'account': "simons"} # changed memory to 16, # back to 12, and with double time (24h)
    spec = '''
        #source /com/extra/bwa/0.7.5a/load.sh
        #source /com/extra/sambamba/0.5.1/load.sh
        cd {title}
        echo hvade
        bwa mem -M -t 16 -a {refgenome_stem} {R1} {R2} \
        | sambamba view -f bam -F "proper_pair" -S -t 8 /dev/stdin \
        | sambamba view -f bam -T {refgenome} /dev/stdin > {ind}_unsorted.bam
          
        sambamba sort -t 16 -m 16GB --out={ind}_sorted.bam --tmpdir=/scratch/$GWF_JOBID/ {ind}_unsorted.bam
        sleep 10

        #rm -f {ind}_unsorted.bam
        # I'm trying to use a different sambamba for the mardup stage, as the conda 0.7.0 install gives a segmentation fault.
        /home/cmkobel/software/sambamba71/sambamba-0.7.1-linux-static markdup -t 16 {ind}_sorted.bam --tmpdir=/scratch/$GWF_JOBID/ {ind}_sort_dedup.bam
        sleep 10

        #rm -f {ind}_sorted.bam
        #rm -f {ind}_sorted.bam.bai
        sambamba flagstat -t 16 {ind}_sort_dedup.bam > {ind}_sort_dedup.bam.flagstat
        sleep 10'''.format(
            title = title, 
            refgenome = refgenome,
            refgenome_stem = refgenome_stem,
            R1 = read1,
            R2 = read2,
            ind = subject)

    return inputs, outputs, options, spec



# 3
def merge_bams_new(title, subject, infiles, outfile, input):
    inputs = [title+'/'+i for i in infiles]
    outputs = [title+'/'+subject+'/'+outfile] # must be a list
    if len(inputs) > 1:
        options = {'cores': 4, 'memory': 8000, 'walltime': '02:45:00', 'account': "simons"}
        spec = '''
            #source /com/extra/sambamba/0.5.1/load.sh
            cd {title}

            mkdir -p {subject}

            sambamba merge -t 4 {subject}/{outfile} {inbams}


            sleep 30
            #rm -f {inbams}
            echo "All done at "$(date) > {subject}.COMPLETED'''.format(
                title = title,
                outfile  = outfile,
                inbams = ' '.join(infiles),
                subject = subject)
    else:
        options = {'cores': 1, 'memory': 1000, 'walltime': '00:10:00', 'account': "simons"}
        spec = '''
            #source /com/extra/sambamba/0.5.1/load.sh
            cd {title}

            mkdir -p {subject}
            touch {subject}_was_here.txt
            cp {inbams} {subject}/{outfile}

            sleep 30
            #rm -f {inbams}
            echo "All done at "$(date) > {subject}.COMPLETED'''.format(
                title = title,
                outfile  = outfile,
                inbams = ' '.join(infiles),
                subject = subject)


    #print(spec)
    return inputs, outputs, options, spec


# 4
def filter_bam_file(title, individual):
    '''
    Filter bam files

    This functions filters the BAM file so we only get the reads that map uniquely and are paired
    '''
    inputs = [title+'/'+individual+'/'+individual+'_merged.bam'] # er det nødvendigt at kalde den individual to gange?
    outputs = [title+'/'+individual+'/'+individual+'_filtered.bam',
               title+'/'+individual+'/'+individual+'_filtered.bam.bai']
    options = {'cores': 1, 'memory': 8000, 'walltime': '10:00:00', 'account': "simons"} # 2 timer er nok til 80% af chimp38
    spec = '''
        cd {title}
        #source /com/extra/sambamba/0.5.1/load.sh

        # original
        sambamba view -F "not (duplicate or secondary_alignment or unmapped) and mapping_quality >= 50 and cigar =~ /100M/ and [NM] < 3" -f bam {ind}/{ind}_merged.bam > {ind}/{ind}_filtered.bam

	# different tries
        #sambamba view -F "mapping_quality >= 50 and cigar =~ /100M/ and [NM] < 3" -f bam {ind}/{ind}_merged.bam > {ind}/{ind}_filtered.bam #0
        #sambamba view -F "cigar =~ /100M/ and [NM] < 3" -f bam {ind}/{ind}_merged.bam > {ind}/{ind}_filtered.bam #0
        
        #sambamba view -F "not (duplicate or secondary_alignment or unmapped)" -f bam {ind}/{ind}_merged.bam > {ind}/{ind}_filtered.bam # this gave 2GB coverage
        




        sambamba index ./{ind}/{ind}_filtered.bam 

        #rm -f {ind}_merged.bam'''.format(title=title, ind=individual)
    #print(spec)
    return inputs, outputs, options, spec


# 5
def get_coverage(title, individual):
    """
    saves the files *_cov.txt inside each individuals directory
    This file contains the number of reads touching each position in the ar.chromosome.
    """
    inputs = [title+'/'+individual+'/'+individual+'_filtered.bam']
    outputs = [title+'/'+individual+'/'+individual+'_cov.txt']
    options = {'cores': 4, 'memory': 4000, 'walltime': '04:00:00', 'account': "simons"}
    spec = '''
        cd {title}
        #source /com/extra/sambamba/0.5.1/load.sh
        #source /com/extra/samtools/1.3/load.sh
        cd {ind}

        samtools depth {ind}_filtered.bam > {ind}_cov.txt
        mail -s "chimpx done {ind}" kobel@pm.me <<< "$(du -h {ind}_cov.txt)"
        '''.format(title=title, ind=individual)
        

    return inputs, outputs, options, spec


# 6:
def get_cnv(title, individual, chrom):
    
    
    inputs = [title+'/'+individual+'/'+individual+'_cov.txt']
    #ooutputs = [inputs[0]+'_pd_median.txt']
    outputs = [title+'/cn_medians/'+chrom+'_'+individual+'_cn_median.csv']
    options = {'cores': 4, 'memory': 64000, 'walltime': '04:00:00', 'account': "simons"}
    spec = '''
    if [ ! -d {dir} ]; then mkdir {dir}; fi
    ~/miniconda3/bin/python compute_cov_median.py {input} {output}
    '''.format(dir = title+'/cn_medians/', input = inputs[0], output = outputs[0])
    return inputs, outputs, options, spec
    #fejl ved manglende fil? fordi den ligge en aterproccessing.
