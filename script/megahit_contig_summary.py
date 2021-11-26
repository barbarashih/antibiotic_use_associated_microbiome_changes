import os
import numpy as np

## Setting
in_dir="analysis/megahit"
in_megares_dir = "analysis/abricate/MegaresDB_megahit/"
out_dir="analysis/megahit_contigsummary/"
out_dir_contig800_2000="analysis/megahit_contigsummary/contig_800_2000/"
out_dir_contig800="analysis/megahit_contigsummary/contig_800/"
out_fp_histo = "analysis/megahit_contigsummary/histogram.csv"
out_fp_binning = "analysis/megahit_contigsummary/contig_binning.csv"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
if not os.path.exists(out_dir_contig800_2000):
    os.mkdir(out_dir_contig800_2000)

# Read in the contig lengths
all_seq_len=np.zeros(500000, dtype=int)
binned_contig_summary = {}
all_samples = os.listdir(in_dir)

for sample in all_samples:
    megares_contigs = {} # contigs with a match with megares db
    binned_contig_names = {}
    count_binned = 0
    count_notbinned = 0
    count_above2000 = 0
    count_between800and2000_withMegares = 0
    count_above2000_withMegares = 0
    contig_name_between_800_2000 = {}
    contig_name_above_800 = {}
    current_seq = ""
    current_seqname = ""
    current_contig_name = ""
    next_contig_name = ""

    # Read in the list of contig_ids that have an entry in megares
    in_megares_fp =  in_megares_dir+ "/" + sample + ".csv"
    with open(in_megares_fp, 'r') as f_in:
        for line in f_in:
            line=line.strip()
            vals = line.split(",")
            contig_id = vals[1]
            megares_contigs[contig_id] = 1
            
    # read in binned contigs names in metabat
    binned_fp = os.listdir("analysis/metabat/megahit/%s/bin" %(sample))
    for current_bin in binned_fp:
        current_bin_fp = "analysis/metabat/megahit/%s/bin/%s" %(sample, current_bin)
        with open(current_bin_fp, 'r') as f_in:
            for line in f_in:
                line=line.strip()
                if line.startswith(">") :
                    binned_contig_names[line.replace(">", "")] = 1
    # read in contigs from megahit
    in_fp = in_dir+ "/" + sample + "/final.contigs.fa"
    with open(in_fp, 'r') as f_in:
        for line in f_in:
            line=line.strip()
            if line.startswith(">") :
                next_contig_name = line.replace(">", "")
                next_contig_name = next_contig_name.split(" ")[0]
                if next_contig_name in binned_contig_names:
                    count_binned +=1
                else:
                    count_notbinned +=1
                # record the number of samples with each length in a numpy array
                if len(current_seq) > 0:
                    current_len = len(current_seq)
                    all_seq_len[current_len] = all_seq_len[current_len]  + 1
                    if 800 <= current_len:
						contig_name_above_800[current_contig_name]=1
                    if 800 <= current_len < 2000:
                        contig_name_between_800_2000[current_contig_name]=1
                        if current_contig_name in megares_contigs:
                            count_between800and2000_withMegares +=1
                    elif current_len >= 2000:
                        count_above2000 += 1
                        if current_contig_name in megares_contigs:
                            count_above2000_withMegares += 1
                    current_seq = ""
                    current_contig_name = next_contig_name
                    
            else:
                current_seq = current_seq + line
    # add the last record
    if len(current_seq) > 0:
        current_len = len(current_seq)
        all_seq_len[current_len] = all_seq_len[current_len]  + 1
        if 800 < current_len < 2000:
            contig_name_between_800_2000[current_contig_name]=1
    # save the list of contig names that is between 800 and 2000 bp
    out_fp= "%s/%s.csv" %(out_dir_contig800_2000 , sample)
    f_out = open(out_fp, 'w')
    for contig in contig_name_between_800_2000:
        f_out.write(contig + "\n")
    f_out.close()
    # save the list of contig names that is between 800 and 2000 bp
    out_fp= "%s/%s.csv" %(out_dir_contig800 , sample)
    f_out = open(out_fp, 'w')
    for contig in contig_name_above_800:
        f_out.write(contig + "\n")
    f_out.close()
    # keep the summary for binned/unbinned contigs
    binned_contig_summary[sample] = {}
    binned_contig_summary[sample]['binned'] = count_binned
    binned_contig_summary[sample]['notbinned'] = count_notbinned
    binned_contig_summary[sample]['count_above2000'] = count_above2000
    binned_contig_summary[sample]['count_between800and2000_withMegares'] = count_between800and2000_withMegares
    binned_contig_summary[sample]['count_above2000_withMegares'] = count_above2000_withMegares

# write the output
f_out = open(out_fp_histo, 'w')
f_out.write("length,count\n")
for idx in range(len(all_seq_len)):
    current_count=all_seq_len[idx]
    if current_count>0:
        f_out.write("%d,%d\n" %(idx, current_count))
f_out.close()


#### Work out the number of binned and not-binned contigs
f_out = open(out_fp_binning, 'w')
f_out.write("sample,binned,notbinned,above2000,between800and2000withMegares,eqabove2000withMegares\n")
for sample in binned_contig_summary:
    f_out.write("%s,%d,%d,%d,%d,%d\n"%(sample,  binned_contig_summary[sample]['binned'] ,  binned_contig_summary[sample]['notbinned'], binned_contig_summary[sample]['count_above2000'] , binned_contig_summary[sample]['count_between800and2000_withMegares'], binned_contig_summary[sample]['count_above2000_withMegares']))
f_out.close()

