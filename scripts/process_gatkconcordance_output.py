import argparse
import pandas as pd

def parse_args():
    """
    Parsing command-line arguments
    """
    parser = argparse.ArgumentParser(description='This script processes output from GATK Genotype Concordance an returns '
                                                 'a file with genotype concordance and a file with allelic concordance for '
                                                 'each individual')

    parser.add_argument('--gatk_concordance', required=True,
                        help='REQUIRED. Input the path to the GATK Genotype Concordance output.')

    parser.add_argument('--num_indiv', required=True,
                        help='REQUIRED. Input the number of individuals you have.')

    parser.add_argument('--count_out', required=True,
                        help='Required. This is the output that lists the counts.')

    parser.add_argument('--genotype_concordance_out', required=True,
                        help='Required. This is the output for genotype concordance')

    parser.add_argument('--allelic_concordance_out', required=True,
                        help='Required. This is the output for genotype concordance')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    # Obtain the counts and output it as a file from GATK Genotype Concordance
    file = open(args.gatk_concordance, "r")
    lines = file.read().splitlines()

    count_outfile = open(args.count_out, "w")
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith("#:GATKTable:GenotypeConcordance_Counts:Per-sample concordance tables: comparison counts"):
            for sample in range(1, int(args.num_indiv) + 2):
                items = lines[i + sample].rstrip("\n").split()
                print("\t".join(items), file=count_outfile)
    count_outfile.close()

    # Calculate genotype and allelic concordance
    data = pd.read_table(args.count_out)
    genotype_out = open(args.genotype_concordance_out, "w")
    allelic_out = open(args.allelic_concordance_out, "w")

    for index, row in data.iterrows():
        match_total = row.HOM_REF_HOM_REF + row.HET_HET + row.HOM_VAR_HOM_VAR
        comparison_total = row.HOM_REF_HOM_REF + row.HOM_REF_HET + row.HOM_REF_HOM_VAR + row.HET_HOM_REF + row.HET_HET + row.HET_HOM_VAR + row.HOM_VAR_HOM_REF + row.HOM_VAR_HOM_VAR + row.HOM_VAR_HET
        print (str(match_total/comparison_total), file=genotype_out)
        allele_score = 2 * row.HOM_REF_HOM_REF + row.HOM_REF_HET + 0 * row.HOM_REF_HOM_VAR + row.HET_HOM_REF + 2 * row.HET_HET + row.HET_HOM_VAR + 0 * row.HOM_VAR_HOM_REF + 2 * row.HOM_VAR_HOM_VAR + row.HOM_VAR_HET
        print (str(allele_score / (comparison_total * 2)), file=allelic_out)

if __name__ == '__main__':
    main()