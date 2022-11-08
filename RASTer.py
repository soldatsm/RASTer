# TODO Переделать логику чтобы не не чиать фасту два раза

from Bio import SeqIO
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-gff_path', help='Path to RAST gff3 table')
parser.add_argument('-id_column_index',
                    help='Number of column in gff, where ids '
                         'is store, probably it is 8th column.'
                         'However', type=int)
parser.add_argument('-faa_path', help='Path to faa file from RAST')
parser.add_argument('-delim', help='Delimiter for table splitting 1 - tab, 2 - comma',
                    choices=[1, 2], type=int)
parser.add_argument('-out_tab', help='Path for output table')
parser.add_argument('-out_faa', help='Path for output FASTA file')
parser.add_argument('-rna', help='Delete rna sequences from table')

args = parser.parse_args()


def write_fasta(fasta_path: str,
                names_list: list[str],
                fasta_output: str):
    """
    :param fasta_path:
    :param names_list:
    :param fasta_output:
    :return: None

    Write faa file with sequence and proteins info from gff table
    :param names_list - names of proteins from gff table
    """

    for string in SeqIO.parse(open(fr'{fasta_path}', 'r'), 'fasta'):
        for names in names_list:
            if string.description == names[0].replace('ID=', ''):
                with open(fr'{fasta_output}', 'a') as write_file:
                    write_file.write('>' + str(string.description)
                                     + ' '
                                     + str(names[1].replace('Name=', '').replace('%', '; '))
                                     + '\n'
                                     + str(string.seq)
                                     + '\n'
                                     + '\n')


def _fasta_reread(path: str):
    protein_names = []
    sequences = []

    for string in SeqIO.parse(open(fr'{path}', 'r'),
                              'fasta'):
        protein_names.append(str(string.description).replace('%', '; '))
        sequences.append(str(string.seq))

    return protein_names, sequences


def make_table(out_faa_path: str, gffdata,
               table_output: str):
    names, seqs = _fasta_reread(out_faa_path)

    # print(len([i.split(' ')[1] for i in names]))
    # print(len(names))
    # print(len(gffdata.iloc[:, 2]))
    # print(len([len(i) for i in seqs]))
    # print(len(seqs))

    table = pd.DataFrame({
        'IDs': [i.split(' ')[0] for i in names],
        'Description': [' '.join(i.split(' ')[1::]) for i in names],
        'Type': gffdata.iloc[:, 2],
        'Lenght_AA': [len(i) for i in seqs],
        'AAs': seqs
                })

    table.to_csv(fr'{table_output}', header=True, index=False)


if __name__ == '__main__':
    if args.delim == 1:
        gff_data = pd.read_csv(fr'{args.gff_path}',
                               delimiter='\t')
    else:
        gff_data = pd.read_csv(fr'{args.gff_path}')

    if args.rna is not None:
        gff_data = gff_data[gff_data.iloc[:, args.id_column_index].str.contains('rna') == False]

    nested_annot_list = [i.split(';')
                         for i in
                         gff_data.iloc[:, args.id_column_index]]

    write_fasta(args.faa_path,
                nested_annot_list,
                args.out_faa)
    if args.out_tab is not None:
        make_table(args.out_faa,
                   gff_data,
                   args.out_tab)
    print('Done')
