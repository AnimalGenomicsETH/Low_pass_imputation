import fnmatch
import pandas as pd

sumcov_file = open(snakemake.output[0],"w+")
sumcov_file.write("bam_id\tITB_id\tautosome\tx\ty\tPAR\tMT\n")
info_file = pd.read_table(snakemake.input.info)
info_file.index = info_file['bam_id']
coverage_files=snakemake.input.file_list
for file in coverage_files:
    segments = file.split("/")
    seg_file = segments[10].split(".")
    ITB_id = (info_file.loc[info_file['bam_id'] == seg_file[0], 'Interbull-ID'].iloc[0])
    covy = pd.read_table(file, header=None , sep = ' ')
    AUTOS = round(sum((covy.iloc[0:28,3])/28),2)
    AUTOS = str(AUTOS)
    X = covy.iloc[29,3]
    X = str(X)
    PAR = covy.iloc[30,3]
    PAR = str(PAR)
    Y = covy.iloc[31,3]
    Y = str(Y)
    MT = covy.iloc[32,3]
    MT = str(MT)
    sumcov_file.write(seg_file[0] + ".bam\t" + ITB_id + "\t" + AUTOS + "\t" + X + "\t" + Y + "\t" + PAR + "\t" + MT + "\n")
sumcov_file.close()
