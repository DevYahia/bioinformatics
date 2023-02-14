import coursea_course.helper.course1_week3 as c1w3
from typing import *
import random as rnd
import math
import operator
import functools as ft


def probability_pattern_profile(pattern: str, profile: List[List[float]]):
    n_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    return ft.reduce(operator.mul, [profile[n_map[nucleotide]][j] for j, nucleotide in enumerate(pattern)])


def random_motifs(dna: List[str], k: int):
    return [(lambda text, n: text[n:n + k])(_dna, rnd.randint(0, len(_dna) - k)) for _dna in dna]


def motifs_from_profile(profile: List[List[float]], dna: List[str]):
    k = len(profile[0])
    return [c1w3.profile_most_probable_kmer(_dna, k, profile) for _dna in dna]


def randomized_motif_search(dna: List[str],
                            k: int,
                            t: int,
                            iterations: int = 1000,
                            print_min=False):
    def _inner_loop():
        motifs = random_motifs(dna, k)
        _best_motifs = motifs[:]
        while True:
            profile = c1w3.Profile(_best_motifs, laplace=True)
            motifs = motifs_from_profile(profile, dna)
            if c1w3.score(motifs) < c1w3.score(_best_motifs):
                _best_motifs = motifs[:]
            else:
                break
        return _best_motifs

    best_motifs = _inner_loop()

    for i in range(iterations):
        last_motifs = _inner_loop()
        score = c1w3.score(last_motifs)
        if score < c1w3.score(best_motifs):
            best_motifs = last_motifs[:]
            if print_min:
                consensus = c1w3.find_consensus(last_motifs)
                print(f"NEW Minimum Found! {i=} \t {score=} \t {consensus=}")

    return best_motifs


def prob_motifs(k: int, t: int, n: int, a: int):
    return math.pow(1 / (n - k + 1), a) * math.pow((n - k) / (n - k + 1), t - a) * t


def profile_randomly_generated_kmer(dna: str, k: int, profile: List[List[float]]):
    n = len(dna)
    probs = [probability_pattern_profile(dna[i:i + k], profile) for i in range(n - k + 1)]
    # probs_sum = sum(probs)
    # normalized_probs = [x / probs_sum for x in probs]
    i = rnd.choices(list(range(n - k + 1)), weights=probs, k=1)[0]
    return dna[i:i + k]


def gibbs_sampler(dna: List[str],
                  k: int,
                  t: int,
                  N: int,
                  iterations: int = 20,
                  print_min=False):
    def _inner_loop():
        motifs = random_motifs(dna, k)
        _best_motifs = motifs[:]
        for j in range(N):
            _i = rnd.randint(0, t - 1)
            motifs.pop(_i)
            profile = c1w3.Profile(motifs, laplace=True)
            motif = profile_randomly_generated_kmer(dna[_i], k, profile)
            motifs.insert(_i, motif)
            if c1w3.score(motifs) < c1w3.score(_best_motifs):
                _best_motifs = motifs[:]
        return _best_motifs

    best_motifs = _inner_loop()

    for i in range(iterations - 1):
        last_motifs = _inner_loop()
        score = c1w3.score(last_motifs)
        if score < c1w3.score(best_motifs):
            best_motifs = last_motifs[:]
            if print_min:
                consensus = c1w3.find_consensus(last_motifs)
                print(f"NEW Minimum Found! {i=} \t {score=} \t {consensus=}")

    return best_motifs
