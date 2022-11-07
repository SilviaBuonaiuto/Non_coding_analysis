import re, sys, argparse
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="path to input FASTA file ",type=str,required=True)
    parser.add_argument("-o", help="path to output table",type=str, required= True)
    args = parser.parse_args()

    input_file = open(args.i, 'r')
    output_file = open(args.o,'w')
    output_file.write('Gene\tA\tC\tG\tT\tCG\tLength\tCG%\n')
    from Bio import SeqIO
    for cur_record in SeqIO.parse(input_file, "fasta"):
        gene_name = cur_record.name
        A_count = cur_record.seq.upper().count('A')
        C_count = cur_record.seq.upper().count('C')
        G_count = cur_record.seq.upper().count('G')
        T_count = cur_record.seq.upper().count('T')
        CG_count = cur_record.seq.upper().count('CG')
        length = len(cur_record.seq)
        cg_percentage = float(C_count + G_count) / length
        output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n' % (gene_name, A_count, C_count, G_count, T_count, CG_count, length, cg_percentage)
        output_file.write(output_line)
    output_file.close()
    input_file.close()

if __name__ == "__main__":
    main()