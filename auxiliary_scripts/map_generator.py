samples=[]
with open ("selected_samples.txt", "rt") as inf:
    next(inf, None) #skips the header
    for line in inf:
        spl=line.rstrip().split()
        samples.append(spl[0])

chromosomes= list (range (1,30)) + ["X", "Y"]
for mychr in chromosomes:
    with open ("maps/" + str(mychr) + ".map", "w" ) as out:
        for sample in samples:
            out.write (f"{sample}\t{sample}_{mychr}.g.vcf.gz\n")
    print (f"done for chr : {mychr}")
