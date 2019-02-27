    ### THE FOLLOWING CAN PRODUCE EIGHT ERRORS IN A ROBUST WAY ###
    """
    for i in range(0,len(upseq1)):
        codon = ""
        if len(codon) < 3:
            codon = upseq[i:i+3]
            if codon == start:
                started = True
                print("I've started!")
            if started == True:
                if (len(codon) == 3) and codon in stopcodons:
                    break
                elif codon in genetic_code.keys() and genetic_code[codon]:
                    aminos1 += genetic_code[codon]
                    codon = ""
                    print(aminos1)
            elif (started != True) and i == (len(upseq1)-1):
                return []
        i + 3

        for i in range(0,len(upseq2)):
            codon = ""
            if len(codon) < 3:
                codon = upseq2[i:i+3]
                if codon == start:
                    started = True
                if started == True:
                    if (len(codon) == 3) and genetic_code[codon] in stopcodons:
                        break
                    elif codon in genetic_code.keys() and genetic_code[codon]:
                        aminos2 += genetic_code[codon]
                        codon = ""
                        
            i + 3

        ### END EIGHT ERRORS PRODUCED ROBUSTLY ###"""

            #elif started and len(codon) == 3:
            #if codon in genetic_code.keys():
            #    aminos1 += genetic_code[codon]
            #    i + 3
            #    codon = ""
        #elif len(codon) == 3:
        #    if codon == start:
        #        started = True
        #        aminos1 += genetic_code[codon]
        #        codon = ""
        #        i + 3

    """started = False
    for letter in upseq2:
        codon = ""
        if len(codon) < 3:
            codon += letter
        elif started and len(codon) == 3:
            if codon in genetic_code.keys():
                aminos2 += genetic_code[codon]
                codon = ""
        elif len(codon) == 3:
            if codon == start:
                started = True
                aminos2 += genetic_code[codon]
                codon = ""

    started = False
    for letter in upseq3:
        codon = ""
        if len(codon) < 3:
            codon += letter
        elif started and len(codon) == 3:
            if codon in genetic_code.keys():
                aminos3 += genetic_code[codon]
                codon = ""
        elif len(codon) == 3:
            if codon == start:
                started = True
                aminos2 += genetic_code[codon]
                codon = ""
    

"""
    
    masteraminos = []
    print(aminos1)
    if len(aminos1) > 0:
        masteraminos.append(aminos1)

    print(aminos2)
    if len(aminos2) > 0:
        masteraminos.append(aminos2)

    print(aminos3)
    if len(aminos3) > 0:
        masteraminos.append(aminos3)

    return masteraminos



    

    ### CODE FROM HERE WILL PRODUCE EIGHT ERRORS ###
    """start = "AUG"
    sequence = rna_sequence.upper()
    startpos = rna_sequence.find(start)
    stopcodons = ["UGA", "UAA", "UAG"]

    newseq = sequence[startpos::]

    i=0
    codon = ""
    aminos1 = ""
    for seq in newseq:
        codon = newseq[i:i+3]
        i += 3
        if codon in stopcodons:
            return aminos1
        elif codon in genetic_code.keys():
            aminos1 += genetic_code[codon]
            codon = ""
        #print(aminos1)


    i=1
    codon = ""
    aminos2 = ""
    for seq in newseq:
        codon = newseq[i:i+3]
        i += 3
        if codon in stopcodons:
            return aminos2
        elif len(codon) > 3 and codon in genetic_code.keys():
            aminos2 += genetic_code[codon]
            codon = ""
        #print(aminos2)

    i=2
    codon = ""
    aminos3 = ""
    for seq in newseq:
        codon = newseq[i:i+3]
        i += 3
        if codon in stopcodons:
            return aminos3
        elif codon in genetic_code.keys():
            aminos3 += genetic_code[codon]
            codon = ""
        #print(aminos3)

    threeaminos = [aminos1, aminos2, aminos3]

    masteraminos = []

    for somelist in threeaminos:
        if len(somelist) > 0:
            masteraminos.append(somelist)
    
    if "*" in masteraminos:
        masteraminos.remove("*")


    print(aminos1, aminos2, aminos3)
    return masteraminos

    ### END CODE TO PRODUCE EIGHT ERRORS ###"""



