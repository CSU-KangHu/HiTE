import os, subprocess, shutil,sys
from collections import OrderedDict


def buildCompleteTransposon(database,combined):
    ''' database = FASTA formatted reference sequence with LTRS and interiors separate
        combined = FASTA formatted reference sequence with LTRS and interiors combined into single sequence

        returns: dictionary mapping names of LTRS to sequences
    '''
   

    with open(database, 'r') as infile:

        lines = infile.readlines()

        interiors = {}
        ltrs = {}
        repeats = set()

        i = 0

        while i < len(lines)-1:

            curr = lines[i]

            if curr[0] == ">":

                label = curr.split()[0][1:]

                seq = ""
                next = ""


            while i < len(lines)-1:
                i += 1
                next = lines[i][0]

                if next == ">":
                    break
                seq += lines[i].rstrip()

            if label[-3:] == "LTR":

                label = label[:-3]
                label = label[:-1] if label[-1] == "_" else label

                if label not in ltrs:
                    repeats.add(label)

                ltrs[label] = seq

            elif label[-1] == "I":

                label = label[:-1]
                label = label[:-1] if label[-1] == "_" else label
                if label not in interiors:
                    repeats.add(label)
                interiors[label] = seq

        storage = []

        for r in repeats:
            if r in interiors and r in ltrs:
                inside = interiors[r]
                outside = ltrs[r]
                total = outside+inside+outside
                data = (r, total)
                storage.append(data)

        with open(combined, 'w') as outfile:
            print "building "+combined
            for s in storage:
                outfile.write(">"+s[0]+"\n")
                outfile.write(s[1]+"\n")
        
        return ltrs


def queryDatabase(database, query, name):
    ''' database = FASTA formatted chromosome
        query = FASTA formatted representative sequences
    '''
    dbDirectory = "dbTemp"+database[:5]

    if not os.path.isdir(dbDirectory):
        os.mkdir(dbDirectory)
    out = dbDirectory + "/"+name
    

    db_params = ['makeblastdb', '-in', database, '-dbtype', 'nucl', '-out', out]
    log =subprocess.check_output(db_params)
    identity = '0.70'
    format = "7 qseqid sseqid qstart qend sstart send length pident qcovhsp qcovs"
    blast_params = ['blastn', '-db', out, '-query', query,'-perc_identity', identity, '-outfmt', format]
    locAlign = subprocess.check_output(blast_params)


    output = name+".blast"

    with open(output, 'w') as outfile:
        outfile.write(locAlign)

    shutil.rmtree(dbDirectory)
    #os.remove(query)
    return output


def parseBlast(name):
    types = {}
    storage = []
    blastFile = name+".blast"
    with open(blastFile, 'r') as infile:

            lines = infile.readlines()
            print "LEN BLAST ="
            print len(lines)
            for l in lines:

                fields = l.split()

                if len(fields) > 0:

                    head = fields[0]

                    if head != "#":  # is valid entry line
                        # identity and query coverage above threshold
                        if float(fields[7]) > 70.0 and float(fields[8]) == 100.0:
                            data = ""
                            if int(fields[4]) <= int(fields[5]):
                                data = (fields[4], fields[5],fields[0])

                            else:
                                data = (fields[5], fields[4], fields[0])

                            storage.append(data)
    os.remove(blastFile)
    storage.sort(key=lambda tup: int(tup[0]))
    storage = OrderedDict((x, True) for x in storage).keys()
    bedFile = name+".bed"

    start = name.rfind('/')+1
    chromName = name[start:]

    with open(bedFile, "w") as outfile:
        for s in storage:
            #print chromName+"\t"+str(s[0])+"\t"+str(s[1])+"\n"
            outfile.write(chromName+"\t"+str(s[0])+"\t"+str(s[1])+"\n")
            types[s[0]] = s[2]

    return bedFile,types


def parseRepeatMasker(file, out):

    storage = []
    interiors = {}

    with open(file, 'r') as infile:

        lines = infile.readlines()

        for l in lines:

            fields = l.split()

            if len(fields) > 0:

                head = fields[0]

                if head.isdigit():

                    ltr = fields[9][-3:] == "LTR"
                    interior = fields[9][-5:] == "I-int" or  fields[9][-7:] == "int-int" or fields[9][-1] == "I"

                    begin = fields[11]
                    left = fields[13]

                    full = (begin == "(0)" and left == "1") or ( begin == "1" and left == "(0)")

                    if ltr and full:
                       

                        label = fields[9][:-3]
                        label = label[:-1] if label[-1] == "_" else label


                        info = (label, fields[5], fields[6])

                        storage.append(info)

                    elif interior:
                        
                        if fields[9][-1] == "I":
                            label = fields[9][:-1]
                        elif fields[9][-5:] == "I-int":
                            label = fields[9][:-5]
                        else:
                            label = fields[9][:-7]

                        label = label[:-1] if label[-1] == "_" else label

                        info = (label, fields[5], fields[6],fields[11],fields[12])


                        if label not in interiors:

                            interiors[label] = [info]
                        else:
                            interiors[label].append(info)

    repeats = {}

    for s in storage:
        name, start, end = s

        if name not in repeats:
            repeats[name] = [(name, start, end)]
        else:
            repeats[name].append((name, start, end))

    next = []

    for r in repeats:

        copies = repeats[r]

        i = 1

        while i < len(copies):

            prev = copies[i-1]
            curr = copies[i]

            if tupleCompare(prev, curr):
                info = (prev[0], int(prev[1]), int(curr[2]), int(prev[1]), int(prev[2]), int(curr[1]), int(curr[2]))
                next.append(info)

                i += 2
            else:
                i += 1

    next.sort()

    final = []

    for f in next:
        name = f[0]


        if name in interiors:

            interior = interiors[name]

        
            interior_len = f[5]-f[4]
            interior_space = [0]*interior_len

            for i in interior:

                if checkInterior(f, i):

                    x = i[3]

                    if x[0] == "(":

                        x = x[1:-1]
                    
                    y = i[4]

                    if y[0] == "(":

                        y = y[1:-1]                    
                    
                    interior_start = int(x)
                    interior_end = int(y)

                    if interior_start > interior_end:

                        temp = interior_start
                        start = end
                        end = interior_start

                    if interior_start > interior_len-1:
                        interior_start = interior_len-1
                    if interior_end > interior_len -1:
                        interior_end = interior_len-1


                    for q in range(interior_start,interior_end):
                        interior_space[q] = 1
            
            coverage = sum(interior_space)

            pct = 0.0

            if float(interior_len) !=0.0:
                pct = float(coverage)/float(interior_len)

            if pct >=0.5:

                final.append((f[1], f[2], f[3], f[4], f[5], f[6]))
                    

    start = out.rfind('/')+1
    chromName = out[start:]

    final.sort(cmp = lambda x,y:compare(x,y))

    out += "Masked.bed"

    with open(out, 'w') as outFile:
        for f in final:
            outFile.write(chromName+"\t"+str(f[0]) + "\t"+str(f[1])+"\t"+str(f[2])+"\t"+str(f[3]) + "\t"+str(f[4]) + "\t"+str(f[5]) + "\n")
    return out

def compare(x,y):

    first = int(x[0])
    second = int(y[0])

    diff = first- second

    val = -1 if diff <0 else 1
    if diff == 0:
        val =0
    return val



def tupleCompare(first, second):

    distance = int(second[2]) - int(first[1])

    return distance > 400 and distance < 22000


def checkInterior(candidate, interior):

    return candidate[1] < int(interior[1]) and candidate[2] > int(interior[2])


def runBedTools(truthFile, bedFile, overlap):
    ''' Runs bedtools intersect utility to calculate the entries in truthFile that overlap with entries in bedFile above some similarity. Returns number of overlaps
    '''
    similarity = "0.95"

    mode = "-u" if overlap else "-v" 
    params = ['bedtools', 'intersect','-sorted', '-a', truthFile, '-b',bedFile, mode, '-f', similarity, '-F', similarity]
    #mode = "-e" if overlap else "-n"
    #params = ['bedops', mode, "95%" ,truthFile,bedFile]

    result= subprocess.check_output(params)

    return result


def getSequence(file, start, end):
  
    with open(file, 'r') as infile:
        label = infile.readline()
        lines = [l.rstrip() for l in infile.readlines()]

        i = 0
        seq = ""
        begin = 0

        for j in range(len(lines)):
            l = lines[j]
            length = len(l)

            if i+length > start:
    
                remainder = length - abs(i+length - start)
                seq += l[remainder-1:]
                i += length

                begin = j+1
                break

            else:
                i += length

        for j in range(begin, len(lines)):
            l = lines[j]
            length = len(l)

            if i+length > end:
                remainder = length - abs(i+length - end)
                seq += l[:remainder]
                i += length
                begin = j
                break

            else:
                seq += l
                i += length
    return seq



def parseBedOutput(bedString):

    lines = bedString.split("\n")
    storage = []

    for l in lines:
        fields = l.split("\t")
        info = tuple([int(x) for x in fields[1:]])
        if len(info) >=1:
            storage.append(info)
    return storage


def merge(first,second,third):

    all = first+second+third
    all.sort()
    return all


def buildGroundTruth(repeats, repeatMaskerFile, chromFile,queryFile, name,outDir):
    print "Querying Database"
    queryDatabase(chromFile, queryFile, name)
    print("Intital Blast query passed")
    blast_output, types = parseBlast(name)
    rm_output = parseRepeatMasker(repeatMaskerFile, name)

    print("finding Intersections")
    masker = runBedTools(rm_output, blast_output, False)
    shared = runBedTools(rm_output, blast_output, True)
    blast = runBedTools(blast_output, rm_output, False)

    os.remove(blast_output)
    os.remove(rm_output)

   

    incompletes = blast.split("\n")
    blast_processed = []

    print "+++++++++++++++++++++++HOW MANY ARE UNIQUE TO BLAST?++++++++++++++++++++++++++++"
    print(len(incompletes))

    mask_count = len(masker.split("\n"))
    shared_count = len(shared.split("\n"))

    print "+++++++++++++++++++++++HOW MANY ARE UNIQUE TO RMASKER?++++++++++++++++++++++++++++"
    print mask_count

    print "+++++++++++++++++++++++HOW MANY ARE Shared?++++++++++++++++++++++++++++"
    print shared_count


    print "Entering Blast corrections"
    
    for i in incompletes:

        fields = i.split()

        if len(i) >= 1:

            start = int(fields[1])
            end = int(fields[2])

            subject_loc = fields[1]+"-"+fields[2]

            seq = getSequence(chromFile, start, end)

            origin = types[fields[1]]
           

            query = repeats[origin]

            queryFile = "temp.fast"

            with open(queryFile, "w") as tempFile:
                tempFile.write(">chrTemp"+"\n")
                tempFile.write(query)

            outfmt = "10 qstart qend sstart send"
            params = ['blastn', '-query', queryFile, '-subject', chromFile, '-subject_loc',
                      subject_loc, '-perc_identity', '70', '-max_hsps', '2', '-outfmt', outfmt]
            locAlign = subprocess.check_output(params).split("\n")[:2]

           
            coords = [y.split(',')[2:] for y in locAlign]


            flattened = [int(item) for l in coords for item in l]
            flattened.sort()
            bounds = [start, end]

            output = tuple(bounds+flattened)

            if len(output) ==6:
                blast_processed.append(output)

    masker_processed = parseBedOutput(masker)

    shared_processed = parseBedOutput(shared)


    finals = merge(blast_processed, masker_processed, shared_processed)


    ground_truth = outDir+name + "Truth.bed"

    with open(ground_truth, 'w') as outFile:
        for f in finals:
            formatted = [str(x) for x in f]
            print formatted
            start = name.rfind('/')+1
            chromName = name[start:]

            line = chromName+"\t"+"\t".join(formatted)+"\n"
            outFile.write(line)

def processChromosome(directory, repeatsFile):

    path_fa = directory+"/"+"Fa"+"/"
    path_out = directory+"/"+"Out"+"/"
    path_gt = directory+"/GT/"

    if not os.path.isdir(path_gt):
        os.mkdir(path_gt)

    repeatsFile = directory+"/"+repeatsFile 
    
    fa = os.listdir(path_fa)
    fa = [f for f in fa if f.endswith(".fa")]

    out = os.listdir(path_out)
    out = [o for o in out if o.endswith("out")]

    pairs = []

    for f in fa:
        for o in out:
            if o.startswith(f):
                
                info = (path_fa+f,path_out+o)
                pairs.append(info)
                break
    
    combined = directory+"/"+directory+"_combined.fasta"
    print combined
    repeats = buildCompleteTransposon(repeatsFile, combined)
    print "Query built"

    for p in pairs:
        
        fa_file = p[0]
        out_file = p[1]
        name =  os.path.splitext(fa_file)[0].split('/')[2]
        print "building "+name
        buildGroundTruth(repeats, out_file, fa_file,combined, name,path_gt)      

  
if __name__ == "__main__":
    

    folder = sys.argv[1]

    repbaseFile = sys.argv[2]

    processChromosome(folder, repbaseFile)






    
