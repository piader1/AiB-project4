#!/usr/bin/python

import sys
import numpy as np

# ===============/ SP_approx /============= #
# ==== PROGRAM FOR MSA OF N SEQUENCES ===== #


# INPUT FROM ARGUMENTS
# [1] FILE WITH LISTED FASTA SEQUENCES
# [2] SCORE MATRIX
# [3] GAP PENALTY

# OUTPUT
# MSA IN FILE NAME SPECIFIED IN ARGUMENT [4]

#
global seqs, subs, g, k, score_m  # keep identical names within functions no need to pass arguments


# ======= FUNCTIONS ======== #

def read_fasta_file(filename):
    # Reads the given FASTA file f and returns a dictionary of sequences.
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
                current_sequence_lines = []
                sequences_lines[sequence_name] = current_sequence_lines
            else:
                if current_sequence_lines is not None:
                    current_sequence_lines.append(line)
    sequences = {}
    for name, lines in sequences_lines.items():
        sequences[name] = ''.join(lines)
    return sequences


def make_matrix(m, n):
    # two dimensional matrix of size m,n
    return[[0] * n for _ in range(m)]


def make_nd_matrix(lseq):
    # input is list of each dimension lengtg
    # returns N dimensional matrix
    l_s = [i + 1 for i in lseq]
    return np.zeros(shape=l_s, dtype=float)


def get_fastas_array(d):
    # converts secuences from dicctionary to list of indeces
    return [translate_bases_to_indices(i) for i in d.values()]


def get_shape(d):
    # takes list of sequences
    # returns "shape" list  of their legth
    return [len(i) for i in d]


def read_subs_marix(filename):
    # Read score table 0
    with open(filename) as fp:
        return [[int(s) for s in line.split() if s.isdigit()] for line in fp]


def translate_bases_to_indices(obs):
    # Bases to indeces run on fasta input
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return [mapping[symbol.upper()] for symbol in obs]


def translate_indices_to_bases(indices):
    # Indeces to Bases
    mapping = {0: 'A', 1: 'C', 2: 'G', 3: 'T', '-': '-'}
    return [mapping[symbol] for symbol in indices]


def write_file(aligned, writefn):
    # Takes list of alignments and writes a file 6
    with open(writefn, 'w') as f:
        for i in range(len(aligned)):
            a = ''.join(str(i) for i in aligned[i])
            f.write(">" + 'seq' + "\n" + a + "\n")
    f.close()

def upper_bound(score_m, k):
    # finds upper bound from optimal score matrix
    return  (len(score_m)-1) * sum(score_m[k])

def lower_bound(score_m, k):
    # finds lower bound from optimal score matrix
    return 0.5 * len(seqs) * sum(score_m[k])

def cost_linear(A, B): # A,B,subs, g
    # Filling the cost table
    T = make_matrix(len(A)+1, len(B)+1)
    for i in range(len(A)+1):
        for j in range(len(B)+1):
            if (T[i][j] == 0):
                v1 = v2 = v3 = v4 = float("inf")

                if (i > 0) and (j > 0):
                    v1 = T[i-1][j-1] + subs[A[i-1]][B[j-1]]
                if (i > 0) and (j >= 0):
                    v2 = T[i-1][j] + g
                if (i >= 0) and (j > 0):
                    v3 = T[i][j-1] + g
                if (i == 0) and (j == 0):
                    v4 = 0
                T[i][j] = min(v1, v2, v3, v4)

    #print('Optimal cost: {}'.format(T[i][j]))
    return T[i][j]


def backtrack(T, A, B): #, subs, g
    # backtracking
    outputA = []
    outputB = []
    j = len(T[0]) - 1 #lenght of A and B from matrix
    i = len(T) - 1
    while i > 0 and j > 0:
        if (i > 0) and (j > 0) and T[i][j] == T[i-1][j-1] + subs[A[i-1]][B[j-1]]:
            #substitution case
            outputA.append(A[i-1])
            outputB.append(B[j-1])
            j = j-1
            i = i-1
        if (i>0) and (j>= 0) and T[i][j] == T[i-1][j] + g:
            # gap
            outputA.append(A[i-1])
            outputB.append('-')
            i = i -1
        if (i>=0) and (j> 0) and (T[i][j] == T[i][j-1] + g):
            outputA.append('-')
            outputB.append(B[j-1])
            j = j -1

    # reverse indices and translate them to bases
    colB = translate_indices_to_bases(outputB)
    colA = translate_indices_to_bases(outputA)
    colA.reverse()
    colB.reverse()
    return colA, colB


def find_s_c(s):
    # make parwise compareason of sequences and choses centroidn           # TRY TO SEPARATE FUNCTIONALITES
    # returns centroids index                                              # USE IN ANOTHER FUNCTION
    n = len(s)
    m = make_matrix(n, n + 1)
    optm = []
    min = float('inf')

    for i in range(n):
        for j in range(n):
            if i!=j:
                m[i][j] = cost_linear(s[i], s[j])
                m[i][n] = m[i][n] + m[i][j]
        optm.append(m[i][n])

    for k in range(len(optm)):
        if optm[k] < min:
            optm[k] = min
            s_c = k

    return k, m


# ============ SUM OF PAIR EXACT COST ALIGNMENT FOR 3 SEQUENCES ============  #

def exact_3_cost(A, B, C):  # A,B,C,subs,g,
    # cost of exact multiple alignment
    D = make_nd_matrix(get_shape(seqs))
    I, J, K = get_shape(seqs)
    for i in range(I + 1):
        for j in range(J + 1):
            for k in range(K + 1):
                v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float("inf")
                if i == 0 and j == 0 and k == 0:
                    v0 = 0
                if i > 0 and j > 0 and k > 0:
                    v1 = D[i - 1][j - 1][k - 1] + subs[seqs[0][i - 1]][seqs[1][j - 1]] + subs[seqs[1][j - 1]][seqs[2][k - 1]] + subs[seqs[0][i - 1]][seqs[2][k - 1]]
                if i > 0 and j > 0 and k >= 0:
                    v2 = D[i - 1][j - 1][k] + subs[seqs[0][i - 1]][seqs[1][j - 1]] + g + g
                if i > 0 and j >= 0 and k > 0:
                    v3 = D[i - 1][j][k - 1] + g + subs[seqs[0][i - 1]][seqs[2][k - 1]] + g
                if i >= 0 and j > 0 and k > 0:
                    v4 = D[i][j - 1][k - 1] + g + g + subs[seqs[1][j - 1]][seqs[2][k - 1]]
                if i > 0 and j >= 0 and k >= 0:
                    v5 = D[i - 1][j][k] + g + g
                if i >= 0 and j > 0 and k >= 0:
                    v6 = D[i][j - 1][k] + g + g
                if i >= 0 and j >= 0 and k > 0:
                    v7 = D[i][j][k - 1] + g + g
                D[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)

    return D


# ============ BACKTRACK THE ALIGNMEN =============#

def backtrack_MSA(D):
    outA = []
    outB = []
    outC = []
    i = len(D) - 1
    j = len(D[0]) - 1
    k = len(D[0][0]) - 1
    while i > 0 or j > 0 or k > 0:
        if (i > 0 and j > 0 and k > 0) and D[i][j][k] == D[i - 1][j - 1][k - 1] + subs[seqs[0][i - 1]][seqs[1][j - 1]] + \
                subs[seqs[1][j - 1]][seqs[2][k - 1]] + subs[seqs[0][i - 1]][seqs[2][k - 1]]:
            outA.append(seqs[0][i - 1])
            outB.append(seqs[1][j - 1])
            outC.append(seqs[2][k - 1])
            i -= 1
            j -= 1
            k -= 1
        if (k >= 0) and D[i][j][k] == D[i - 1][j - 1][k] + subs[seqs[0][i - 1]][seqs[1][j - 1]] + g + g:
            outA.append(seqs[0][i - 1])
            outB.append(seqs[1][j - 1])
            outC.append('-')
            i -= 1
            j -= 1
        if (j >= 0) and D[i][j][k] == D[i - 1][j][k - 1] + g + subs[seqs[0][i - 1]][seqs[2][k - 1]] + g:
            outA.append(seqs[0][i - 1])
            outB.append('-')
            outC.append(seqs[2][k - 1])
            i -= 1
            k -= 1
        if (i >= 0) and D[i][j][k] == D[i][j - 1][k - 1] + g + g + subs[seqs[1][j - 1]][seqs[2][k - 1]]:
            outA.append('-')
            outB.append(seqs[1][j - 1])
            outC.append(seqs[2][k - 1])
            j -= 1
            k -= 1
        if (j >= 0 and k >= 0) and D[i][j][k] == D[i - 1][j][k] + g + g:
            outA.append(seqs[0][i - 1])
            outB.append('-')
            outC.append('-')
            i -= 1
        if (i >= 0 and k >= 0) and D[i][j][k] == D[i][j - 1][k] + g + g:
            outA.append('-')
            outB.append(seqs[1][j - 1])
            outC.append('-')
            j -= 1
        if (i >= 0 and j >= 0) and D[i][j][k] == D[i][j][k - 1] + g + g:
            outA.append('-')
            outB.append('-')
            outC.append(seqs[2][k - 1])
            k -= 1

    return translate_indices_to_bases(outA), translate_indices_to_bases(outB), translate_indices_to_bases(outC)


def ma_approx(seqs, k):
    # Make Alignments With Sc
    M = []
    for i in range(len(seqs)-1):
        if i != k:
            a, b = backtrack(cost_linear(seqs[k]), seqs[i], seqs[k], seqs[i])


# ===================== EXECUTE =====================#

seqs = get_fastas_array(read_fasta_file(sys.argv[1]))
subs = read_subs_marix(sys.argv[2])
g = int(sys.argv[3])
k, score_m = find_s_c(seqs)
#fn = sys.argv[4]


# ma_approx(seqs, find_s_c(seqs))

print(upper_bound(score_m, k), lower_bound(score_m, k))

print(seqs[0][5])

print(exact_3_cost(seqs[0], seqs[1], seqs[2]))
































