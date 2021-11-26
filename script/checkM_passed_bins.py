# this script go through checkM result and print out a list of bins that pass the checkM completeness
# and prints out an annotation file for each bin for plotting the phylophlan tree using graphlan
import os, sys

# takes input
# 1] checkm_batch directory
# 2] checkm output prefix
# 3] number of batches in checkm
# 4] file path for printing out a list of passed bins
# 5] file path for printing out summarised data for all checkm
# 6] completeness_threshold
# 7] contamination threshold

in_dir=sys.argv[1]
in_file_prefix=sys.argv[2]
in_range=int(sys.argv[3]) # the total number of batches in checkm
out_fp_passed_bins = sys.argv[4]
out_fp_all_bins = sys.argv[5]
completeness_threshold = float(sys.argv[6]) # only keeping bins with a completeness over this threshold
contamination_threshold = float(sys.argv[7]) # only keeping bins with a completeness over this threshold
passed_bins = []
bin_annotation = {}
passed=0
out_fp="file_list/bins_annotation.txt" # this records the annotation for the bins using checkm info
f_out = open(out_fp, 'w', newline='')
# Go through each checkM output
for i in range(in_range +1 ):
    current_fp=os.path.join(in_dir, in_file_prefix + str(i) + "_plus.txt")
    if os.path.exists(current_fp):
        with open(current_fp, 'r') as f_in:
            for line_num, line in enumerate(f_in):
                if line_num>0:
                    line=line.strip()
                    vals=line.split("\t")
                    current_bin=vals[0]
                    marker_linage = vals[1]
                    marker_linage_uid = vals[2]
                    completeness = vals[12]
                    contamination = vals[13]
                    current_sample = current_bin.split("_bin")[0]
                    # check if the sample is pig or human
                    if current_sample.endswith("P"):
                        species='pig'
                        species_col="#ef8a62"
                    else:
                        species='human'
                        species_col="#999999"
                    #f_out.write("%s\t%s\t%s\n"%(current_bin, "clade_marker_color", species_col)) # write the annotation file for graphlan
                    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_shape", 1, 'v')) # write the annotation file for graphlan
                    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 2, species_col)) # write the annotation file for graphlan
                    # record the marker linage
                    marker_linage_out =marker_linage + " " + marker_linage_uid
                    # put the bins, and its corresponding sample name, into a separate folder
                    if float(completeness) > completeness_threshold and float(contamination) <= contamination_threshold:
                        passed_bins.append(current_bin)
                        passed +=1
                        last_passed = current_bin # record what was the last bin that passed the completeness so I can double check that is correct
                    bin_annotation[current_bin]={}
                    bin_annotation[current_bin]['completeness']=completeness
                    bin_annotation[current_bin]['contamination']=contamination
                
f_out.close()
# print out passed bins
f_out = open(out_fp_passed_bins, 'w', newline="")
for file in passed_bins:
    f_out.write(file + "\n")
f_out.close()

# print out checkm annotation
f_out = open(out_fp_all_bins, 'w', newline="")
out_line = "%s\t%s\t%s\n"%("bin",'completeness', 'contamination')
f_out.write(out_line)
for current_bin in bin_annotation:
    out_line = "%s\t%s\t%s\n"%(current_bin, bin_annotation[current_bin]['completeness'], bin_annotation[current_bin]['contamination'])
    f_out.write(out_line)
f_out.close()
