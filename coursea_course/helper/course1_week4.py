import coursea_course.helper.course1_week1 as c1w1
import coursea_course.helper.course1_week2 as c1w2
import coursea_course.helper.course1_week3 as c1w3
from typing import *
import random as rnd


def motifs_from_profile(profile: List[List[float]], dna: List[str]):
    k = len(profile[0])
    return [c1w3.profile_most_probable_kmer(_dna, k, profile) for _dna in dna]


def randomized_motif_search(dna: List[str], k: int, t: int, iterations: int = 1000, print_min=False):
    def _inner_loop():
        motifs = [(lambda text, n: text[n:n + k])(_dna, rnd.randint(0, len(_dna) - k)) for _dna in dna]
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
