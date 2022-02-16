import subprocess
from Bio import SeqIO
import plotly.express as px

trans_cdna = "Homo_sapiens.GRCh38.cdna.all.fa"
trans_ncrna = "Homo_sapiens.GRCh38.ncrna.fa"
genome_file = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# IF YOU ARE USING WINDOWS COMPLETE THESE STRINGS WITH SUITABLE PATHS
path_blastn = "C:/Users/andrz/PycharmProjects/pythonProject/blast-2.12.0+/bin/blastn"  # path to your blastn program
path_makeblastdb = "C:/Program Files/NCBI/magicblast-1.6.0/bin/makeblastdb"  # path to your makeblastdb program

'''We will look for DNA contamination from organism (ref_organism) in human/bacterial DNA. In order to do that, we will
check the first 10 000 sequences'''

ref_organism = []
records = SeqIO.parse("DRR008476.fastq", "fastq")

for record in records:
    if int(record.id.split('.')[1]) <= 10000:
        ref_organism.append(record)
    else:
        break

SeqIO.write(ref_organism, "ref_organism.fasta", "fasta")

# Here you put the downloaded transcriptome into fasta file "transcriptome", do the same with the genome
transcriptome = []

records = SeqIO.parse(trans_cdna, "fasta")
for record in records:
    transcriptome.append(record)
records = SeqIO.parse(trans_ncrna, "fasta")
for record in records:
    transcriptome.append(record)
SeqIO.write(transcriptome, "transcriptome.fasta", "fasta")

genome = []
records = SeqIO.parse(genome_file, "fasta")

for record in records:
    genome.append(record)
SeqIO.write(genome, "genome.fasta", "fasta")

human_genome = []
human_transcriptome = []
bacteria_genome = []
bacteria_transcriptome = []
# Here you can choose where you want to look for contamination
print("1 - human_genome \n2 - human_transcriptome \n3 - bacteria_genome \n4 - bacteria_transcriptomee")
x = int(input("Choose number of operation you want to perform: "))  

if x == 1:
    records = SeqIO.parse("genome.fasta", "fasta")

    for record in records:
        human_genome.append(record)
    SeqIO.write(human_genome, "human_genome.fasta", "fasta")

    subprocess.call(path_makeblastdb + " -in human_genome.fasta -dbtype nucl " +
                    " -out human_genome.fasta ")

    subprocess.call(path_blastn + " -task megablast " +
                    " -db human_genome.fasta -outfmt 6 -query ref_organism.fasta -max_target_seqs 1 " +
                    "-out human_genome_result.out ")

    hg_result = []
    z = open('human_genome_result.out')

    for line in z:
        if float(line.split('\t')[3]) >= 95 and float(line.split('\t')[2]) >= 95:
            hg_result.append(line)

    num_lines = sum(1 for line in open('human_genome_result.out'))
    uniqueperc = round(((len(hg_result)/num_lines)*100), 2)

    print("The number of unique readings is: " + str(len(hg_result)))
    print("It is " + str(uniqueperc) + "% of all readings")

    labels = ['Unique readings', 'Not unique readings']
    values = [uniqueperc, 100 - uniqueperc]
    fig = px.pie(values=values, names=labels, color_discrete_sequence=px.colors.sequential.RdBu)
    fig.show()


elif x == 2:
    records = SeqIO.parse("cdna.fasta", "fasta")
    for record in records:
        human_transcriptome.append(record)
    SeqIO.write(human_transcriptome, "human_transcriptome.fasta", "fasta")
    subprocess.call(
        path_makeblastdb + " -in human_transcriptome.fasta -dbtype nucl " +
        " -out human_transcriptome.fasta ")

    subprocess.call(path_blastn + " -task megablast " +
                    " -db human_transcriptome.fasta -outfmt 6 -query ref_organism.fasta -max_target_seqs 1 " +
                    "-out human_transcriptome_result.out ")

    ht_result = []
    z = open('human_transcriptome_result.out')

    for line in z:
        if float(line.split('\t')[3]) >= 95 and float(line.split('\t')[2]) >= 95:
            ht_result.append(line)

    num_lines = sum(1 for line in open('human_transcriptome_result.out'))
    uniqueperc = round(((len(ht_result) / num_lines) * 100), 2)

    print("The number of unique readings is: " + str(len(ht_result)))
    print("It is " + str(uniqueperc) + "% of all readings")

    labels = ['Unique readings', 'Not unique readings']
    values = [uniqueperc, 100 - uniqueperc]
    fig = px.pie(values=values, names=labels,
             color_discrete_sequence=px.colors.sequential.RdBu)
    fig.show()

elif x == 3:
    records = SeqIO.parse("genome.fasta", "fasta")
    for record in records:
        bacteria_genome.append(record)
    SeqIO.write(bacteria_genome, "bacteria_genome.fasta", "fasta")
    subprocess.call(path_makeblastdb + " -in bacteria_genome.fasta -dbtype nucl " +
                    "-out bacteria_genome.fasta")

    subprocess.call(path_blastn + " -task megablast " +
                    " -db bacteria_genome.fasta -outfmt 6 -query ref_organism.fasta -max_target_seqs 1" +
                    "-out bacteria_genome_result.out ")

    bg_result = []
    z = open('bacteria_genome_result.out')
    for line in z:
        if float(line.split('\t')[3]) >= 95 and float(line.split('\t')[2]) >= 95:
            bg_result.append(line)

    num_lines = sum(1 for line in open('bacteria_genome_result.out'))
    uniqueperc = round(((len(bg_result) / num_lines) * 100), 2)

    print("The number of unique readings is: " + str(len(bg_result)))
    print("It is " + str(uniqueperc) + "% of all readings")

    labels = ['Unique readings', 'Not unique readings']
    values = [uniqueperc, 100 - uniqueperc]
    fig = px.pie(values=values, names=labels,
                 color_discrete_sequence=px.colors.sequential.RdBu)
    fig.show()

elif x == 4:
    records = SeqIO.parse("cdna.fasta", "fasta")
    for record in records:
        bacteria_transcriptome.append(record)
    SeqIO.write(bacteria_transcriptome, "bacteria_transcriptome.fasta", "fasta")
    subprocess.call(path_makeblastdb + " -in bacteria_transcriptome.fasta " +
                    " -dbtype nucl -out bacteria_genome.fasta")

    subprocess.call(path_blastn + " -task megablast " +
                    " -db bacteria_transcriptome.fasta -outfmt 6 -query ref_organism.fasta -max_target_seqs 1" +
                    " -out bacteria_transcriptome_result.out ")

    bt_result = []
    z = open('bacteria_transcriptome_result.out')
    for line in z:
        if float(line.split('\t')[3]) >= 95 and float(line.split('\t')[2]) >= 95:
            bt_result.append(line)

    num_lines = sum(1 for line in open('bacteria_transcriptome_result.out'))
    uniqueperc = round(((len(bt_result) / num_lines) * 100), 2)

    print("The number of unique readings is: " + str(len(bt_result)))
    print("It is " + str(uniqueperc) + "% of all readings")

    labels = ['Unique readings', 'Not unique readings']
    values = [uniqueperc, 100 - uniqueperc]
    fig = px.pie(values=values, names=labels,
                 color_discrete_sequence=px.colors.sequential.RdBu)
    fig.show()

else:
    print('wrong number')
    pass
