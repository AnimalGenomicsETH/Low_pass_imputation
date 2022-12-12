import gzip

myfile = snakemake.input[0]
vcf_file = gzip.open(myfile, "rt")
mypath = myfile.split("/")
mysamp = mypath[12].split(".")[0]
outfile = open(snakemake.output[0], "w+")

outfile.write("stage\tvar_cal\tsample\tvariants_chip\tmatch_variants\tdiscrepant_variants\tALT2HET\tHET2ALT\tALT2NOCALL\tHET2NOCALL\tHOMREF2HET\tHOMREF2HOMALT\tNOCALLTRUTH\n")

nvar = 0
match = 0
discr = 0
ALT2HET = 0
HET2ALT = 0
ALT2NOCALL = 0
HET2NOCALL = 0
HOMREF2HET = 0
HOMREF2HOMALT = 0
NOCALL = 0

for line in vcf_file:
    if line [0] != "#":
        nvar+=1
        spl=line.rstrip().split()
        truth=spl[9].split(":")
        query=spl[10].split(":")
        if truth[5] == query[5]:
            match+=1
        else:
            discr+=1
            if truth[5] == "homalt" and query[5] == "het":
                ALT2HET +=1
            elif truth[5] == "het" and query[5] == "homalt":
                HET2ALT +=1
            elif truth[5] == "homalt" and query[5] == "nocall":
                ALT2NOCALL +=1
            elif truth[5] == "het" and query[5] == "nocall":
                HET2NOCALL +=1
            elif truth[5] == "homref" and query[5] == "het":
                HOMREF2HET +=1
            elif truth[5] == "homref" and query[5] == "homalt":
                HOMREF2HOMALT +=1
            elif truth[5] == "nocall":
                NOCALL +=1

results = (mypath[10], mypath[11], mysamp, str(nvar), str(match), str(discr), str(ALT2HET), str(HET2ALT), str(ALT2NOCALL), str(HET2NOCALL), str(HOMREF2HET), str(HOMREF2HOMALT), str(NOCALL))
towrite="\t".join(results)
outfile.write (f"{towrite}\n")
