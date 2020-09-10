import re
import os

# requires python 3.8
fname = "test_db"
organisms = []
regex = re.compile('[^a-zA-Z]')

if not os.path.exists(fname+"_split"):
    os.makedirs(fname+"_split")

with open(fname+".fasta", "r") as f:
    outname = ""
    while (r := f.readline()):
        if ">" in r:
            with open((outname := fname+"_split/"+r.strip().replace(">", "")+".fasta"), "a") as wf:
                if outname not in organisms:
                    organisms.append(outname)
                wf.write(r)
        else:
            with open(outname, "a") as wf:
                wf.write(r)

counter = 0
with open(fname+"_split/label_file.txt", "w") as wf:
    for o in organisms:
        wf.write(o.split("/")[1]+"\t"+str(counter)+"\n")
        counter += 1