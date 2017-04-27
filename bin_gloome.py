import csv, sys
from collections import defaultdict

def bin_gloome(gloome_in, outfile):

    # Input should be only the potential gain nodes of interest (e.g. N1-N6)
    with open(gloome_in, 'r') as infile:
        dictreader = csv.DictReader(infile)
        # Put all the GLOOME data in a list
        gloome_dat = [(row['POS'],(row['Node'], row['Prob'])) for row in dictreader]

    # Consolidate list data into non-redundant dictionary {sRNA:[node, probabiliy]}
    gloome_nr = defaultdict(list)
    for pos, node_prob in gloome_dat:
        gloome_nr[pos].append(node_prob)

    srna_gain_nodes = []
    srna_ct = 0
    for srna, node_prob_list in gloome_nr.items():
        srna_ct += 1
        # Sort the nodes in descending order (i.e. youngest first)
        node_prob_sort = sorted(node_prob_list, key=lambda tup: tup[0], reverse=True)
        if any(float(prob) >= 0.70 for node, prob in node_prob_sort):
            for node, prob in node_prob_sort:
                if float(prob) >= 0.70:
                    gain_node = node
                    gn_prob = prob
                else:
                    break
        else:
            gain_node = 'leaf'
            gn_prob = 1
        srna_gain_nodes.append((srna, gain_node, gn_prob))

    age_bin_dict = {'N6':'nascent', 'N5':'nascent','N4':'nascent','N3':'nascent','N2':'established','N1':'established', 'leaf': 'nascent'}
    print(age_bin_dict)

    nascent_ct = 0
    intermed_ct = 0
    est_ct = 0
    with open(outfile, 'w') as outfile:
        fieldnames = ['sRNA_id', 'gain_node', 'prob', 'age_bin']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for srna, gn, pp in srna_gain_nodes:
            results_dict = {}
            results_dict['sRNA_id'] = srna
            results_dict['gain_node'] = gn
            results_dict['prob'] = pp
            results_dict['age_bin'] = age_bin_dict[gn]
            if results_dict['age_bin'] == 'nascent':
                nascent_ct += 1 
            if results_dict['age_bin'] == 'intermediate':
                intermed_ct += 1
            if results_dict['age_bin'] == 'established':
                est_ct += 1
            writer.writerow(results_dict)

    print('%i sRNAs analyzed. RESULTS: nascent: %i, intermediate: %i, established: %i' % (srna_ct, nascent_ct, intermed_ct, est_ct) )

if __name__ == '__main__':
    if len(sys.argv) == 3:
         bin_gloome(sys.argv[1], sys.argv[2])
    else:
         print("Usage: bin_gloome.py gloome_in.csv age_bin_out.csv")
         sys.exit(0)






