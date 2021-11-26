# encoding=utf8

# this script go through checkM result and print out a list of bins that pass the checkM completeness
# and prints out an annotation file for each bin for plotting the phylophlan tree using graphlan
import os, sys, re

# takes input
# 1] file list with bins that have passed checkm completeness
# 2] diamond_plus 
# 3] abricate plasmidfinder results directory
# 4] abricate megaresdb results directory
# 5] megaresdb annotation (the csv file for converting gene name to amr class)
# 6] this indicate if the script should ignore coassembly for pig/human
# 7] checkm result
# 8] in_fp_phylophlan_input
# 9] output file prefix

########################
##### settings - variables
in_fp_checkm_passed=sys.argv[1]
in_fp_diamond_plus=sys.argv[2]
in_dir_plasmidfinder=sys.argv[3]
in_dir_megaresdb=sys.argv[4]
in_fp_megaresdb_annotation=sys.argv[5]
in_coassembly=sys.argv[6] # this value is yes/no - yes to include co-assembly, no to exclude
in_checkm_results=sys.argv[7]
in_fp_phylophlan_input=sys.argv[8]
out_fp=sys.argv[9] # this is the prefix of the output files
in_fp_annotation="sample_annotation/metadata_20210122.csv"

amr_coverage_threshold=0
plasmidfinder_coverage_threshold=0

passed_bins={}
unclassified_bins={}
all_phylum={}

# the phylum details for some of the bins were not successfully converted in MAGpy. These are manually converted by matching their taxid
manually_added_phylum={
'33039':['Bacteria','Terrabacteria group','Firmicutes','Clostridia','Eubacteriales','Lachnospiraceae','Mediterraneibacter'],
'39491':['Bacteria','Terrabacteria group','Firmicutes','Clostridia','Eubacteriales','Lachnospiraceae','unclassified Lachnospiraceae'],
'59620':['Bacteria','Terrabacteria group','Firmicutes','Clostridia','Eubacteriales','Clostridiaceae','Clostridium'],
'165186':['Bacteria','Terrabacteria group','Firmicutes','Clostridia','Eubacteriales','Oscillospiraceae','Ruminococcus'],
'446043':['Bacteria','Terrabacteria group','Firmicutes','Clostridia','Eubacteriales','Lachnospiraceae','Lachnospira'],
'1076179':[],
'1318617':['Bacteria','Terrabacteria group','Tenericutes','Mollicutes','Mycoplasmatales','Mycoplasmataceae','Mycoplasma'],
'1643353':['Bacteria','Bacteria incertae sedis','Bacteria candidate phyla','Candidatus Absconditabacteria','unclassified Candidatus Absconditabacteria','',''],
'1768112':[],
'1898205':['Bacteria','Terrabacteria group','Firmicutes','Clostridia','Eubacteriales','Oscillospiraceae','unclassified Oscillospiraceae'],
'1916229':['Bacteria','Terrabacteria group','Cyanobacteria/Melainabacteria group','Candidatus Melainabacteria','Candidatus Gastranaerophilales','',''] ,
'2320115':[]}


phylum_hex_col={"Not classified":'#ffffff'}
phylum_hex_col_list=[
        '#1f78b4',
        '#b2df8a',
        '#33a02c',
        '#fb9a99',
        '#e31a1c',
        '#fdbf6f',
        '#ff7f00',
        '#cab2d6',
        '#6a3d9a',
        '#ffff99',
        '#b15928',
        '#301403',
        '#696969',
        '#7FFFD4',
        '#ccccff',
        '#ff007f',
        '#c8a2c8',
        '#1f78b4',
        '#b2df8a',
        '#33a02c',
        '#fb9a99',
        '#e31a1c',
        '#fdbf6f',
        '#ff7f00',
        '#cab2d6',
        '#6a3d9a',
        '#ffff99',
        '#702963']
sample_col_list=[
        '#4575b4',
        '#543005',
        '#bf812d',
        '#f6e8c3',
        '#003c30',
        '#35978f',
        '#c7eae5',
        '#f781bf',
        '#999999', ]
megares_hex_col_list=['#1E78B5',
        '#AEC7E7',
        '#F07E20',
        '#F9BA7A',
        '#2AA137',
        '#9DCB88',
        '#D62729',
        '#F29697',
        '#8F67A9',
        '#C5B0D5',
        '#8D574C',
        '#C39C94',
        '#D279AF',
        '#F5B6D2',
        '#7F7F7F',
        '#C8C7C7',
        '#BCBD20',
        '#DCDC8D',
        '#28B8CB',
        '#9FD6E2']
origin_hex_col={
        'Pig':"#0000ff",
        'Human':"#ff0000"}
sample_type_hex_col={
        'Dry Sows':'#D83C50',
        'Farrowing Sows':'#F2EB91',
        'Pooled Sows':'#228B22',
        'Piglets':'#3488BF',
#        'Sow_faeces_pooled':'#6a3d9a',
#        'Soil':'#b15928',
#        'Weaner_faeces_6_wks':'#ece2f0',
#        'Weaner_faeces_8_wks':'#a6bddb',
#        'Finisher_faeces_24_wks':'#1c9099',
#        'Sow_faeces':'#cab2d6',
#        'Sow_faeces_pooled':'#6a3d9a',
#        'coassembly':'#ffffff'
        }
sample_timepoint_hex_col={
        'T1':'#edf8fb',
        'T2':'#b2e2e2',
        'T3':'#66c2a4',
        'T4':'#238b45',
        }
sample_amrstatus_hex_col={
        'NA':'#ffffff',
        'Post_tylosin':'#ff7f00',
        'Post_chlorotetracycline':'#1f78b4',
        'Pre_chlorotetracycline':'#a6cee3',
        'Pre_chlorotetracycline_2':'#b2df8a',
        'coassembly':'#ffffff',
        }

amr_class_grouping={
    'betalactams':1,
    'Aminoglycosides':9,
    'Trimethoprim':13,
    'Sulfonamides':11,
    'Fluoroquinolones':4,
    'Elfamycins':9,
    'MLS':8,
    'Rifampin':9,
    'Tetracyclines':12,
    'Multi-drug resistance':3,
    'Metronidazole':7,
    'Phenicol':10,
    'Bacitracin':9,
    'Glycopeptides':5,
    'Cationic antimicrobial peptides':2,
    'Aminocoumarins':9,
    'Lipopeptides':6,
    'Mycobacterium tuberculosis-specific Drug':10,
        }

amr_class_grouping_to_name={
#        1:'Aminoglycosides',
        1:'Betalactams',
        2:'CAP',
        3:'CCR',
        4:'Fluoroquinolones',
        5:'Glycopeptides',
        6:'Lipopeptides',
        7:'Metronidazole',
        8:'MLS',
        9:'Others',
        10:'Phenicol',
        11:'Sulfonamides',
        12:'Tetracyclines',
        13:'Trimethoprim',

#        2:'Aminoglycosides',
#        3:'TMPS',
#        5:'Elfamycins',
#        7:'Rifampin',
#        11:'Metronidazole',
#        12:'Glycopeptides',
#        13:'Phenicol',
#        14:'Bacitracin',
#        15:'Cationic antimicrobial peptides',
#        16:'Aminocoumarins',
#        17:'Lipopeptides',
#        18:'Mycobacterium tuberculosis-specific Drug',
 }

amr_frequency_col={
        '0':"#ffffff",
        '1':"#cccccc",
        '2':'#969696',
        '3-10':'#636363',
        '10+':'#252525',        
        }
amr_diversity_col={
        '0':"#ffffff",
        '1':"#bae4b3",
        '2':'#74c476',
        '3-10':'#31a354',
        '10+':'#006d2c',
        }
amr_hex_col={
        1:'#ff0000', 
        2:'#1a1aff',
        3:'#00b386',
        4:'#86592d',
        5:'#e60073',
        6:'#0099cc',
        7:'#009900',
        8:'#ff8c1a',
        9:'#9900cc',
        10:'#000000',
        11:'#ff0000',
        12:'#1a1aff',
        13:'#00b386',
        14:'#86592d',
        }
# these bins are in cluster 3 for the amr genes plot Adrian produced. 
plot_specific_samples = False # change this to false if it is not to be plotted
specific_samples=[
        'MBDBAG011P',
        'KLAKYN063P',
        'KLAKYN064P',
        'MBDMAD033P',
        'MBDBAG012P',
        'KLAKYN065P',
        'KLAMAK051P',
        'KLAKYN061H',
        'MBDMADCTR02-5',
        'KLANAN062H',
        'KLANANCTR02-3',
        'MBDMADCTR01−5',
        'KLANANCTR02−4',
        'MBDMADCTR01−2',
        ]
#############################################            
# read and construct a list of all the phylophlan input
phylophlan_input={}
with open(in_fp_phylophlan_input, 'r') as f_in:
    for line in f_in:
        if line.startswith(">"):
            line=line.strip()
            phylophlan_input[line.replace(">","")] =1

#### Base annotation 
#### Read in sample annotation

sample_annotation={}
annotation_header=[]
with open(in_fp_annotation, 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        line=line.split(",")
        if line_num==0:
            annotation_header=line

        else:
            sample_id=line[0]
            sample_annotation.update({sample_id:{}}) if sample_id not in sample_annotation else 0
            for idx, val in enumerate(annotation_header):                
                current_annotation_value = line[idx].strip()
                annotation_swap ={
                        'Farrowing_Sows':'Farrowing Sows',
                        'Dry_Sows':'Dry Sows',
                        'piglets':'Piglets',
                        'Pooled':'Pooled Sows',
                        }
                sample_annotation[sample_id][val]=current_annotation_value if current_annotation_value not in annotation_swap else annotation_swap[current_annotation_value]



            
#### checkm
# Go through checkm passed samples
all_samples_col={}
present_sample_type={}
with open(in_fp_checkm_passed, 'r') as f_in:
    for line in f_in:
        current_bin=line.strip()
        if current_bin in phylophlan_input:
            current_sample=current_bin.split("_bin")[0]
            # check if the sample is pig or human
            passed_bins[current_bin]={}
            passed_bins[current_bin]['sample_id']=current_sample
            passed_bins[current_bin]['sample_type']=sample_annotation[current_sample]['sampletype']
            passed_bins[current_bin]['timepoint']=sample_annotation[current_sample]['timepoint']
            present_sample_type[sample_annotation[current_sample]['sampletype']]=1
            all_samples_col[current_sample]=""

#### checkm result
with open(in_checkm_results, 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        if line_num > 0:
            vals=line.split("\t")
            current_bin=vals[0]
            if current_bin in passed_bins:
                passed_bins[current_bin]['checkm_completeness']=vals[1]
                passed_bins[current_bin]['checkm_contamination']=vals[2]

#### Go through specific samples and label bins from them; these are in Cluster 3 of ARM heatmap)
for current_bin in passed_bins:
    current_sample = current_bin.split("_bin")[0]
    if current_sample in specific_samples:
        passed_bins[current_bin]['specific_bins']='#000000'

#### Diamond 
# Go through diamond_bin_report_plus.tsv
not_classified_count={}
tax_classification_colidx={'class':12, 'order':13, 'family':14, 'genus':15}
with open(in_fp_diamond_plus, 'r') as f_in:
    for line in f_in:
        line=line.strip()
        vals = line.split("\t")
        current_bin=vals[0]
        if current_bin in passed_bins:
            if len(vals) > 11:
                passed_bins[current_bin]['phylum']=vals[11]
            else:
                # find the taxid and check it against the manually curated dictionary
                current_taxid=re.match('.*?OX=(.*?)$', vals[6])
                tax_id = current_taxid.group(1)
                tax_id=tax_id.strip()
                if len(manually_added_phylum[tax_id]) > 0:
                    passed_bins[current_bin]['phylum']=manually_added_phylum[tax_id][2]
                else:
                    passed_bins[current_bin]['phylum']="Not classified"
                    not_classified_count.update({tax_id:0}) if tax_id not in not_classified_count else 0
                    not_classified_count[tax_id] +=1
                    print("%s\t%s"%(current_bin, tax_id))
                
            for tax_key in tax_classification_colidx:
                if len(vals) > (tax_classification_colidx[tax_key]):
                    passed_bins[current_bin][tax_key]=vals[tax_classification_colidx[tax_key]]                    
            passed_bins[current_bin]['species'] = vals[6]
            all_phylum[passed_bins[current_bin]['phylum']]=1
print("not classified" )
print(not_classified_count)
current_hex_col_idx=0
for current_phylum in all_phylum:
    if current_phylum not in phylum_hex_col:
        phylum_hex_col[current_phylum]=phylum_hex_col_list[current_hex_col_idx]
        current_hex_col_idx+=1
    
# Print out the annotations for genomes with 
current_out_fp=out_fp+"_annotation_base.txt"
f_out = open(current_out_fp, 'w', newline='')

f_out.write("annotation_background_alpha\t0.1\n")
f_out.write("start_rotation\t270\n")
f_out.write("class_legend_font_size\t10\n")
f_out.write("annotation_legend_font_size\t10\n")
f_out.write("annotation_background_separation\t-0.03\n")
f_out.write("annotation_background_offset\t0.0\n")
f_out.write("annotation_background_width\t270\n")
f_out.write("start_rotation\t270\n")
# print the ring annontation settings

ring_labels=['phylum']
for ring_idx in range(len(ring_labels)):
    f_out.write("ring_label_font_size\t%d\t10\n"%(ring_idx +1))
    f_out.write("ring_internal_separator_thickness\t%d\t0.2\n"%(ring_idx +1))
    f_out.write("ring_separator_color\t%d\t#FFFFFF\n"%(ring_idx +1))
    f_out.write("ring_label\t%d\t%s\n"%(ring_idx +1 ,ring_labels[ring_idx]))
    f_out.write("ring_label_color\t%d\t#0000FF\n"%(ring_idx +1) )

for current_bin in passed_bins:
    if 'phylum' in passed_bins[current_bin]:
        current_phylum=passed_bins[current_bin]['phylum']
        current_phylum_col=phylum_hex_col[current_phylum]
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 1, current_phylum_col)) # write the annotation file for graphlan
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", 1, 0.5)) # write the annotation file for graphlan

all_phylum_keys=list(phylum_hex_col.keys())
all_phylum_keys.sort()
all_phylum_keys = [val for val in all_phylum_keys if val != "Not classified"]
#all_phylum_keys.append("Not classified")
for current_phylum in all_phylum_keys:
    current_phylum_col=phylum_hex_col[current_phylum]
    f_out.write("%s\t%s\t%s\n"%(current_phylum, "clade_marker_color", current_phylum_col)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_phylum, "clade_marker_size", 50)) # write the annotation file for graphlan
f_out.close()

####################################
#### AMR class + plasmid annotation
if True:
    all_plasmidfinder_fp=os.listdir(in_dir_plasmidfinder)
    for file in all_plasmidfinder_fp:
        current_fp= os.path.join(in_dir_plasmidfinder,file)
        with open(current_fp, 'r') as f_in:
            for line_num, line in enumerate(f_in):
                if line_num > 0:
                    line=line.strip()
                    line=line.split(",")
                    current_bin=os.path.basename(line[0]).replace(".fa", "")
                    current_coverage=float(line[8])
                    current_gene=line[4]
                    if current_coverage > plasmidfinder_coverage_threshold and current_bin in passed_bins:
                        passed_bins[current_bin].update({'plasmidfinder':[]})
                        passed_bins[current_bin]['plasmidfinder'].append(current_gene)
            
#### Megares + abricate - antibiotic resistance annotation
# Get the megaresdb annotation for annotating abricate data
megares_gene2class={}
with open(in_fp_megaresdb_annotation, 'r') as f_in:
    for line_num, line in enumerate(f_in):
        line=line.strip()
        line=line.split(",")
        if line_num==0:
            continue
        else:
            megares_gene2class[line[0]]=line[1]
# Go through all results in the abricate folder
all_genes={}
all_class={}
if True:
    all_abricate_fp = os.listdir(in_dir_megaresdb)
    for current_abricate_fp in all_abricate_fp:
        headerline=[]
        with open(os.path.join(in_dir_megaresdb, current_abricate_fp), 'r') as f_in:
            for line_num, line in enumerate(f_in):
                line=line.strip()
                line=line.split(",")
                if line_num==0:
                    headerline=line
                    gene_idx=[val_idx for val_idx,val in enumerate(headerline) if val == "GENE"][0]
                    coverage_idx=[val_idx for val_idx,val in enumerate(headerline) if val == "%COVERAGE"][0]
                else:
                    # only reccord it if the coverage is > the amr_coverage_threshold
                    current_coverage=float(line[coverage_idx])
                    if current_coverage > amr_coverage_threshold:
                        current_bin=line[0]
                        current_bin=current_bin.split("/")[-1]
                        current_bin=current_bin.replace(".fa","")
                        current_gene=line[gene_idx]
                        # get the class annotation by using the gene name. Record the amr classes present for each bin
                        # note that each bin can have multiple class
                        current_class=megares_gene2class[current_gene]
                        amr_class_grouping.update({current_class:10}) if current_class not in amr_class_grouping else 0
                        # this just records the total number of times a class has appeared in the full dataset
                        all_genes.update({current_gene:0}) if current_gene not in all_genes else 0
                        all_genes[current_gene]+=1
                        # this records the amr calss for each of the bin, the following ar recorded
                        # - the frequency (number of times a gene of a class has occurred) 
                        # - the diversity (number of different genes for the particular amr class) 
                        if current_bin in passed_bins:
                            passed_bins[current_bin].update({'amr_class':{}}) if 'amr_class' not in passed_bins[current_bin] else 0
                            passed_bins[current_bin]['amr_class'].update({current_class:{}}) if current_class not in passed_bins[current_bin]['amr_class'] else 0
                            passed_bins[current_bin]['amr_class'][current_class].update({current_gene:0}) if current_gene not in passed_bins[current_bin]['amr_class'][current_class] else 0
                            passed_bins[current_bin]['amr_class'][current_class][current_gene] += 1
                            all_class.update({current_class:0}) if current_class not in all_class else 0
                            all_class[current_class]+=1
# generate reverse keys for amr class
amr_class_grouping_rev={}
for key in amr_class_grouping:
    amr_class_grouping_rev.update({amr_class_grouping[key]:[]}) if amr_class_grouping[key] not in amr_class_grouping_rev else 0
    amr_class_grouping_rev[amr_class_grouping[key]].append(key)

#### Print annotation
current_out_fp=out_fp+"_annotation_amr.txt"
f_out = open(current_out_fp, 'w', newline='') 
# print out the base graph annotation
f_out.write("annotation_background_alpha\t0.1\n")
f_out.write("start_rotation\t270\n")
f_out.write("class_legend_font_size\t10\n")
f_out.write("annotation_legend_font_size\t10\n")
f_out.write("annotation_background_separation\t-0.03\n")
f_out.write("annotation_background_offset\t0.0\n")
f_out.write("annotation_background_width\t270\n")
f_out.write("start_rotation\t270\n")  
# print the setting for each of the amr class group ring. There are 8 groups, 1-10 ,in total. 
# depending on if any other annotations are printed, a number, ring_idx_modifier, was added to each group
ring_labels=['phylum','sample_type', 'timepoint', 'plasmid']
base_ring_num=len(ring_labels)
ring_idx_modifier=len(ring_labels)
amr_group2ringid={}
amr_ringid2group={}
# check which amr group are actually present
amr_ring_id=1 + len(ring_labels)
for group_idx in list(amr_class_grouping_to_name.keys()):
    amr_group_members=amr_class_grouping_rev[group_idx]
    current_amr_group_name=amr_class_grouping_to_name[group_idx]
    amr_recorded=False
    for current_amr_group in amr_group_members:
        if current_amr_group in all_class and len(current_amr_group) > 0:
            amr_recorded=True
    if amr_recorded:
        amr_group2ringid[current_amr_group_name]=amr_ring_id
        amr_ringid2group[amr_ring_id]=group_idx
        amr_ring_id+=1

for ring_id in range(len(ring_labels)+1, amr_ring_id):
    # get the classes present in each group and use it as a label
    current_ring_label_frequency = amr_class_grouping_to_name[amr_ringid2group[ring_id]] # the name of the black ring
    ring_labels.append(current_ring_label_frequency)

    #current_ring_label_diversity = amr_class_grouping_to_name[group_idx] + "_diversity" # the name of the but let me know 
    #ring_labels.append(current_ring_label_diversity)
if plot_specific_samples:
    ring_labels.append("specific_bin")
# print ring index
if True:
    for ring_idx in range(len(ring_labels)):
        f_out.write("ring_label_font_size\t%d\t5\n"%(ring_idx+1))
        f_out.write("ring_internal_separator_thickness\t%d\t0.2\n"%(ring_idx+1))
        f_out.write("ring_separator_color\t%d\t#D3D3D3\n"%(ring_idx+1))
        f_out.write("ring_label\t%d\t%s\n"%(ring_idx+1,ring_labels[ring_idx]))
        if (ring_idx >= base_ring_num) and (ring_idx < amr_ring_id):
            ring_label_col=amr_hex_col[(ring_idx-base_ring_num + 1)]
            f_out.write("ring_label_color\t%d\t%s\n"%(ring_idx+1, ring_label_col))
        else:
            f_out.write("ring_label_color\t%d\t#000000\n"%(ring_idx+1))
# go through each of the passed bins and add their annotation
for current_bin in passed_bins:
    current_sample_id=passed_bins[current_bin]['sample_id']
    current_sample_col=all_samples_col[current_sample_id]
    current_sample_type=passed_bins[current_bin]['sample_type']
    current_sample_timepoint=passed_bins[current_bin]['timepoint']
    current_sample_type_col=sample_type_hex_col[current_sample_type]
    current_sample_timepoint_col=sample_timepoint_hex_col[current_sample_timepoint]
    
    
    # check if there is any result from plasmidfinder
    current_plasmid_presence=True if ('plasmidfinder' in passed_bins[current_bin] and (len(passed_bins[current_bin]['plasmidfinder']) > 0)) else False
    if current_plasmid_presence:
        print(current_bin)
    current_specific_bins=True if ('specific_bins' in passed_bins[current_bin] and (len(passed_bins[current_bin]['specific_bins']) > 0)) else False
    # base colour
    if 'phylum' in passed_bins[current_bin]:
        current_phylum=passed_bins[current_bin]['phylum']
        current_phylum_col=phylum_hex_col[current_phylum]
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 1, current_phylum_col)) # write the annotation file for graphlan
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", 1, 0.5)) # write the annotation file for graphlan
        f_out.write("%s\t%s\t%d\n"%(current_bin, "clade_marker_size", 12)) # write the annotation file for graphlan
        f_out.write("%s\t%s\t%s\n"%(current_bin, "clade_marker_color", current_phylum_col)) # write the annotation file for graphlan
    else:
        unclassified_bins[current_bin]=1
    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 2, current_sample_type_col)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", 2, 0.5)) # write the annotation file for graphlan
#    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 3, current_sample_armstatus_col)) # write the annotation file for graphlan
#    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", 3, 0.5)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 3, current_sample_timepoint_col)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", 3, 0.5)) # write the annotation file for graphlan

    # annotation for plasmidfinder
    if current_plasmid_presence:
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", 4, "#fc9272")) # write the annotation file for graphlan
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", 4, 0.2)) # write the annotation file for graphlan        
    # annotation for specific bins
    if current_specific_bins and plot_specific_samples:
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", len(ring_labels), "#FF0000")) # write the annotation file for graphlan
        f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", len(ring_labels), 0.2)) # write the annotation file for graphlan        
    # annotation for amr class rings
    if 'amr_class' in passed_bins[current_bin]:
        current_amr=passed_bins[current_bin]['amr_class']
        for current_amr_class in current_amr:
            amr_number_frequency=sum(current_amr[current_amr_class].values())
            amr_number_diversity=len(current_amr[current_amr_class].keys())
            if current_amr_class in amr_class_grouping:
                current_amr_group=amr_class_grouping[current_amr_class]
                current_amr_grouping_name=amr_class_grouping_to_name[current_amr_group]
                current_ring_frequency_idx=amr_group2ringid[current_amr_grouping_name]

    #            current_ring_frequency_idx=((current_amr_group*2)-1)+ring_idx_modifier
    #            current_ring_diversity_idx=(current_amr_group*2)+ring_idx_modifier
                def count_to_bin(f_count):
                    if f_count < 1 or not f_count or f_count==None:
                        f_count='0'
                        f_alpha=0
                    elif f_count<3:
                        f_alpha=(f_count* 0.2) + 0.2
                        f_count=str(f_count)
                    elif 2 < f_count < 11:
                        f_count="3-10"
                        f_alpha=0.8
                    else:
                        f_count="10+"
                        f_alpha=1
                    return([f_count, f_alpha])
                current_fill_col= amr_hex_col[current_amr_group]
                frequency_bin=count_to_bin(amr_number_frequency) # the lists returned by the function are the count and alpha catogeries
                frequency_alpha=frequency_bin[1]
    #            diversity_bin=count_to_bin(amr_number_diversity)
    #            current_fill_col_frequency=amr_frequency_col[frequency_bin]
    #            current_fill_col_diversity=amr_diversity_col[diversity_bin]
                f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", current_ring_frequency_idx, current_fill_col)) # write the annotation file for graphlan
                f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", current_ring_frequency_idx, 0.2)) # write the annotation file for graphlan
                f_out.write("%s\t%s\t%d\t%.2f\n"%(current_bin, "ring_alpha", current_ring_frequency_idx, frequency_alpha)) # write the annotation file for graphlan
    #            f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_color", current_ring_diversity_idx, current_fill_col_diversity)) # write the annotation file for graphlan
    #            f_out.write("%s\t%s\t%d\t%s\n"%(current_bin, "ring_height", current_ring_diversity_idx, 0.1)) # write the annotation file for graphlan
            else:
                pass

for current_phylum in all_phylum_keys:
    current_phylum_col=phylum_hex_col[current_phylum]
    #f_out.write("%s\t%s\t%s\n"%(current_phylum, "annotation_font_size ", 3)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_phylum, "clade_marker_color", current_phylum_col)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_phylum, "clade_marker_size", 30)) # write the annotation file for graphlan
for current_key in present_sample_type:
    current_key_col=sample_type_hex_col[current_key]
    #f_out.write("%s\t%s\t%s\n"%(current_key, "annotation_font_size ", 3)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_key, "clade_marker_color", current_key_col)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_key, "clade_marker_size", 30)) # write the annotation file for graphlan
for current_key in sample_timepoint_hex_col:
    current_key_col=sample_timepoint_hex_col[current_key]
    #f_out.write("%s\t%s\t%s\n"%(current_key, "annotation_font_size ", 3)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_key, "clade_marker_color", current_key_col)) # write the annotation file for graphlan
    f_out.write("%s\t%s\t%s\n"%(current_key, "clade_marker_size", 30)) # write the annotation file for graphlan
    
f_out.close()  
# print out the colour annotation for each ring for each bin.
current_out_fp=out_fp+"_flat_annotation.txt"
out_parameters=['phylum','sample_id']
f_out=open(current_out_fp, 'w+')
f_out.write("sample\t" + "\t".join(out_parameters) + "\n")
for current_bin in passed_bins:
    out_line=[current_bin]
    for current_parameter in out_parameters:
        if current_parameter  in passed_bins[current_bin]:
            if current_parameter in ['plasmidfinder']:
                current_append = "||".join(passed_bins[current_bin][current_parameter])
            elif current_parameter in ['amr_class']:
                current_append = "||".join(list(passed_bins[current_bin][current_parameter].keys()))
            else:
                current_append = passed_bins[current_bin][current_parameter]
        elif current_parameter in ["household"]:
            current_basename = current_bin.split("_bin")[0]
            if "-" in current_basename:
                current_household = current_basename.split("-")[0]
            elif current_basename in ['human', 'pig']:
                current_household = current_basename
            else:
                current_household = current_basename[0:-2]
            current_append = current_household        
        else:
            current_append = ""
        out_line.append(current_append)
    f_out.write("\t".join(out_line) + "\n")
f_out.close()

# print out gene hit count
current_out_fp1=out_fp+"_flat_annotation_amr_pergene.csv"
current_out_fp2=out_fp+"_flat_annotation_amr_perclass.csv"
current_out_fp3=out_fp+"_flat_annotation_tax.csv"
f_out1=open(current_out_fp1, 'w+')
f_out1.write("%s\n"%(','.join(['sample_id', 'bin', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'gene', 'count', 'amr_class'])))
f_out2=open(current_out_fp2, 'w+')
f_out2.write("%s\n"%(','.join(['sample_id', 'bin', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'amr_class', 'count'])))
f_out3=open(current_out_fp3, 'w+')
f_out3.write("%s\n"%(','.join(['sample_id', 'bin', 'phylum', 'class', 'order', 'family', 'genus', 'species'])))
for current_bin in passed_bins:
    sample_id=current_bin.split("_bin")[0]
    out_tax=[]
    for tax_key in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
        if tax_key in passed_bins[current_bin]:
            out_tax.append(passed_bins[current_bin][tax_key])
        else:
            out_tax.append("")
    out_line=[sample_id, current_bin] + out_tax
    f_out3.write("%s\n" %(",".join(out_line)))
    if 'amr_class' in passed_bins[current_bin]:
        for current_amr_class in passed_bins[current_bin]['amr_class']:
            current_amr_count=str(len(passed_bins[current_bin]['amr_class'][current_amr_class]))
            out_line=[sample_id, current_bin] + out_tax + [current_amr_class, current_amr_count]
            f_out2.write("%s\n" %(",".join(out_line)))
            for current_gene in passed_bins[current_bin]['amr_class'][current_amr_class]:
                current_gene_count=str(passed_bins[current_bin]['amr_class'][current_amr_class][current_gene])
                out_line=[sample_id, current_bin] + out_tax + [current_gene, current_gene_count, current_amr_class]
                f_out1.write("%s\n" %(",".join(out_line)))
    else:
        pass
#        out_line=[sample_id, current_bin] + out_tax + ["no_amr", "0"]
#        f_out2.write("%s\n" %(",".join(out_line)))
#        out_line=[sample_id, current_bin] + out_tax + ["no_amr", "0", "no_amr"]
#        f_out1.write("%s\n" %(",".join(out_line)))
            
f_out1.close()
f_out2.close()
f_out3.close()
# print out the number of unclassified bin + total bin
current_out_fp=out_fp+"_classification_summary.txt"
f_out=open(current_out_fp, 'w+')
all_bin=len(passed_bins)
unclassified_bins=len(unclassified_bins)
f_out.write("Total number of bins: %d\n"%all_bin)
f_out.write("Total number of unclassified bins: %d\n"%unclassified_bins)
f_out.close()
