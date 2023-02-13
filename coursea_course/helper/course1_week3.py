from typing import *
import numpy as np
import scipy.stats as sp
import math
import functools as ft
import operator
import itertools

import coursea_course.helper.course1_week2 as hf2
import coursea_course.helper.constants as consts


def motif_enumeration(dna: List[str], k: int, d: int):
    patterns = set()
    dna0 = dna[0]
    n = len(dna0)
    for i in range(n - k + 1):
        kMerPattern = dna0[i: i + k]
        neighbors = hf2.neighbors(kMerPattern, d)
        for neighbor in neighbors:
            flag = True
            for _dna in dna:
                if hf2.approximate_pattern_count(_dna, neighbor, d) < 1:
                    flag = False
                    break
            if flag:
                patterns.add(neighbor)

    return patterns


def matrix_entropy(motifs: List[str]):
    L = len(motifs[0])
    n = len(motifs)
    count = {k: [0 for _ in range(L)] for k in "ACGT"}
    for t, _dna in enumerate(motifs):
        for i in range(L):
            nucleotide = motifs[t][i]
            count[nucleotide][i] += 1

    profile = {k: (np.array(v) / n).tolist() for k, v in count.items()}

    profile_matrix = np.array(list(profile.values()))

    return np.sum(sp.entropy(profile_matrix, base=2))


# my implementation of entropy(), found out that there's already a function in scipy
def entropy(col: List[float]):
    totalEntropy = 0
    for num in col:
        if num > 0:
            totalEntropy += num * math.log(num, 2)
    return -totalEntropy


def Motifs(pattern: str, dna: List[str]) -> List[str]:
    k, t, n = len(pattern), len(dna), len(dna[0])
    motifs = []
    for _dna in dna:
        positions = [k] * (n - k + 1)
        for i in range(n - k + 1):
            positions[i] = hf2.hamming_distance(pattern, _dna[i: i + k])
        minimum_distance = min(positions)
        index = positions.index(minimum_distance)
        motifs.append(_dna[index: index + k])
    return motifs


def d(pattern: str, dna: List[str]):
    motifs = Motifs(pattern, dna)
    return sum([hf2.hamming_distance(motif, pattern) for motif in motifs])


def median_string(dna: List[str], k: int):
    distance = np.inf
    median = "A" * k
    all_patterns = hf2.neighbors(median, k)
    for pattern in all_patterns:
        pattern_distance = d(pattern, dna)
        if distance > pattern_distance:
            distance = pattern_distance
            median = pattern
    return median


def profile_most_probable_kmer(dna: str, k: int, profile: List[List[float]]):
    # print(dna, k, profile, sep="\n")
    n = len(dna)
    n_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    prob = 0
    most_probable = dna[:k]
    for i in range(n - k + 1):
        pattern = dna[i: i + k]
        n_prob = ft.reduce(operator.mul, [profile[n_map[nucleotide]][j] for j, nucleotide in enumerate(pattern)])
        # print(pattern, n_prob)
        if n_prob > prob:
            most_probable = pattern
            prob = n_prob
    return most_probable


def Profile(motifs: List[str]):
    L = len(motifs[0])
    n = len(motifs)
    count = {k: [0 for _ in range(L)] for k in "ACGT"}
    for t in range(n):
        for i in range(L):
            nucleotide = motifs[t][i]
            count[nucleotide][i] += 1

    profile = {k: (np.array(v) / n).tolist() for k, v in count.items()}

    return list(profile.values())


def score(motifs: List[str]):
    consensus = find_consensus(motifs)
    return sum([hf2.hamming_distance(motif, consensus) for motif in motifs])


def find_consensus(motifs: List[str]):
    profile = Profile(motifs)
    n = len(motifs[0])
    np_profile = np.array(profile)
    max_vals = np_profile.argmax(axis=0)
    n_map = {0: "A", 1: "C", 2: "G", 3: "T"}
    consensus = ""
    for i in range(n):
        max_val = max_vals[i]
        consensus += n_map[max_val]
    return consensus


# I have spent hours trying to find what's wrong with my algorithm (literally hours) found out that
# the issue was in the [profile_most_probable_kmer()] function. The issue was solved by initiating
# [most_probable] with the first k-mer in the dna.
def greedy_motif_search(dna: List[str], k: int, t: int):
    n = len(dna[0])
    best_motifs = [_dna[:k] for _dna in dna]
    for i in range(n - k + 1):
        motifs = [dna[0][i: i + k]]
        for j in range(1, t):
            profile = [[y + 1 for y in x] for x in Profile(motifs)]
            most_probable = profile_most_probable_kmer(dna[j], k, profile)
            motifs.append(most_probable)
        # print(motifs)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]
    return best_motifs


def distance_between_pattern_and_strings(pattern: str, dna: List[str]):
    k = len(pattern)
    n = len(dna[0])
    distance = 0
    for _dna in dna:
        hamming_distance = np.inf
        for i in range(n - k + 1):
            pattern_0 = _dna[i:i + k]
            hamming_distance_0 = hf2.hamming_distance(pattern, pattern_0)
            if hamming_distance > hamming_distance_0:
                hamming_distance = hamming_distance_0
        distance += hamming_distance
    return distance
