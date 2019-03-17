#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """
    rna_sequence = rna_sequence.upper()

    if len(rna_sequence) < 3:
        return ""
    elif ("UGA" or "UAA" or "UAG") in  rna_sequence[0:3]:
        return ""
    else:
        aminos = ""
        i = 0
        stopcodons = ["UGA", "UAA", "UAG"]
        #if (len(rna_sequence) % 3 == 0):
        for things in rna_sequence:
            codon = rna_sequence[i:i+3]
            if codon in stopcodons:
                return aminos
            elif codon in genetic_code.keys():
                aminos += genetic_code[codon]
                #print(aminos)
                i += 3
        return aminos
        """else:
            for codons in rna_sequence:
                codon = rna_sequence[i:i+3]
                i += 3
                if codon == ("UGA" or "UAG" or "UAA"):
                    aminos += genetic_code[codon]
                    return aminos
                if codon in genetic_code.keys():
                    aminos += genetic_code[codon]
            aminos = aminos.upper()
            return aminos"""

    """else:
        i = 0
        aminos = ""
        threebases = ""
        for letter in rna_sequence:
            if len(threebases) < 3:
                threebases += rna_sequence[i]
                i += 1
            elif len(threebases) == 3:
                if genetic_code[threebases]:
                    aminos += genetic_code[threebases]
                    threebases = ""
        return aminos"""







    



def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
    """
    upseq = rna_sequence.upper()
    upseq1 = upseq[0:]
    upseq2 = upseq[1:]
    upseq3 = upseq[2:]


    start = "AUG"
    stopcodons = ["UGA", "UAG", "UAA", "*"]
    started = False

    preaminos1 = upseq1.find(start)

    aminos1 = upseq1[preaminos1:]


    preaminos2 = aminos1[1:]

    startaminos2 = preaminos2.find(start)

    aminos2 = preaminos2[startaminos2:]

    preaminos3 = aminos2[1:]

    startaminos3 = preaminos3.find(start)

    aminos3 = preaminos3[startaminos3:]

    
    if aminos1:
        fullaminos1 = translate_sequence(aminos1, genetic_code)
    else:
        fullaminos1 = []

    if aminos2:
        fullaminos2 = translate_sequence(aminos2, genetic_code)
    else:
        fullaminos2 = []

    if aminos3:
        fullaminos3 = translate_sequence(aminos3, genetic_code)
    else:
        fullaminos3 = []

    alltogethernow = []
    if fullaminos1:
        alltogethernow.append(fullaminos1)
    if fullaminos2:
        alltogethernow.append(fullaminos2)
    if fullaminos3:
        alltogethernow.append(fullaminos3)
    #alltogethernow = [fullaminos1, fullaminos2, fullaminos3]
    #i = 0
    print(aminos1, aminos2, aminos3)

    return alltogethernow

    #print(upseq1, upseq2, upseq3)

    upseq1done = []

    #for i in range(0,len(upseq1)):
    #    if upseq[i:i+3] == start:
    #        upseq1done = translate_sequence(upseq1, genetic_code)
    #        return upseq1done
    #        break
    #    elif upseq[i:i+3] != start:
    #        i+1
    #    elif len(upseq1) == 0:
    #        upseq1done = []

    #return upseq1done




def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, and empty string is returned.
    """
    
    ### JAMIE'S WAY ###
    #sequence = sequence.upper()
    #rev_seq = ""

    #for c in sequence:
    #    rev_seq += c
    #return rev_seq

    ### END JAMIE'S WAY ###

    ### BEGIN JAMIE'S SECOND WAY ###
    sequence = sequence.upper()
    return sequence[::-1]

    ### END JAMIE'S SECOND WAY ###

    sequence = sequence.upper()

    if len(sequence) == 0:
        return ""

    rseq = sequence[::-1]
    rseq = rseq.upper()

    return rseq
    #All get_reverse tests cleared

def get_complement(sequence):
    """Get the complement of `sequence`.
    
    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, and empty string is returned.
    """
    ### BEGIN JAMIE'S WAY ###

    # Jamie's way is literally identical to my way minus the variable names

    sequence = sequence.upper()

    #comp_seq = ""
    #for c in sequence:
    #    comp_seq += comp_bases[c]

    ### END JAMIE'S WAY ###


    compdict = {
    "A" : "U",
    "U" : "A",
    "C" : "G",
    "G" : "C"
    }

    sequence = sequence.upper()

    complement = ""
    for letter in sequence:
        if letter in compdict.keys():
            complement += compdict[letter]
    return complement



def reverse_and_complement(sequence):
    """Get the reversed and complemented form of `sequence`.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, and empty string is returned.
    """

    ### JAMIE'S WAY ###

    return get_reverse(get_complement(sequence))

    ### END JAMIE'S WAY ###

    compdict = {
    "A" : "U",
    "U" : "A",
    "C" : "G",
    "G" : "C"
    }

    sequence = sequence.upper()

    if len(sequence) == 0:
        return ""

    rseq = sequence[::-1]
    rseq = rseq.upper()

    complement = ""
    for letter in rseq:
        if letter in compdict.keys():
            complement += compdict[letter]
    return complement
    



def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (three reading frames of the
    current orientation, and the reversed and complemented form) and return (as
    a string) the longest sequence of amino acids that it encodes, according to
    the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty list is returned.
    """

    ### CURRENT VERSION; ALMOST WORKING ###
    finalanswer = ""

    basetranslation = translate_sequence(rna_sequence, genetic_code)
    reverseit = get_reverse(rna_sequence)
    complement = get_complement(rna_sequence)
    revcomp = reverse_and_complement(rna_sequence)
    
    # It doesn't like the ones above, so many it'll like all the translations in different orders/arrangements?
    threeframes = get_all_translations(rna_sequence, genetic_code)
    revthreeframes = get_all_translations(reverseit, genetic_code)
    compframes = get_all_translations(complement, genetic_code)
    revcompframes = get_all_translations(revcomp, genetic_code)


    #### UGHHHHH, THE DIRECTIONS SAY "REVERSED AND COMPLEMENTED FORM" BUT NOT "REVERSED, PAUSE, AND COMPLEMENTED" (SEPARATE), WHICH IS WHY IT ISN'T WORKING ###
    everybody = threeframes+revcompframes

    lenlist = []
    thinglength = 0
    for thing in everybody:
        thinglength = len(thing)
        lenlist.append(thinglength)

    if everybody:
        thebigwinner = max(everybody, key=len)
        finalanswer = thebigwinner

    return finalanswer




    # OLD VERSIONS
    """compdict = {
    "A" : "U",
    "U" : "A",
    "C" : "G",
    "G" : "C"
    }

    sequence = rna_sequence.upper()

    caminos = ""
    for letter in sequence:
        if letter in compdict.keys():
            caminos += compdict[letter]

    
    rseq = sequence[::-1]

    raminos = ""
    for letter in sequence:
        if letter in compdict.keys():
            raminos += compdict[letter]

    aminos = []

    if len(sequence) < 3:
        return ""
    elif ("UGA" or "UAA" or "UAG") in rna_sequence[0:3]:
        return ""
    
    else:
        aminos = ""
        i = 0
        stopcodons = ["UGA", "UAA", "UAG"]
        #if (len(rna_sequence) % 3 == 0):
        for codons in raminos:
            codon = raminos[i:i+3]
            if codon in stopcodons:
                return raminos
            elif codon in genetic_code.keys():
                raminos += genetic_code[codon]
                #print(raminos)
                i += 3
        return raminos


    return raminos"""







if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
