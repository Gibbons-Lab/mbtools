/*
 * Copyright 2018 Christian Diener <mail[at]cdiener.com>
 *
 * Apache License 2.0.
 */

#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <sstream>

using namespace Rcpp;

unsigned int miniter = 20;

// Convert an int vector to string representation. Used to hash the vectors.
std::string vec_to_str(std::vector<int> vec) {
    std::stringstream ss;
    for (unsigned int i=0; i<vec.size(); ++i) {
        ss << "|" << vec[i];
    }
    return ss.str();
}

// Compresses the mapping of reads to transcripts and their weights into
// equivalence classes.
//
// This means that reads that map to exactly the same transcripts get lumped
// together. This is abit of hacky data structure to optimize speed. It
// maps hash(transcripts) -> [transcripts, sum(weights) over reads] which is
// why the weights have to be integers.
std::unordered_map<std::string, std::vector<int> > equivalence_classes(
    std::vector<std::vector<int> > reads_to_txs,
    std::vector<std::vector<int> > reads_to_weights) {
    std::unordered_map<std::string, std::vector<int> > ecs;
    std::string ec;
    std::vector<int> txs;
    unsigned int start;
    for (unsigned int i=0; i<reads_to_txs.size(); ++i) {
        txs = reads_to_txs[i];
        std::sort(txs.begin(), txs.end());
        ec = vec_to_str(txs);
        if (ecs.count(ec) > 0) {
            ecs[ec][0] += 1;
            start = (ecs[ec].size() + 1)/2;
            for (unsigned int k=start; k<ecs[ec].size(); k++) {
                ecs[ec][k] += reads_to_weights[i][k-start];
            }
        } else {
            txs.insert(txs.begin(), 1);
            txs.insert(txs.end(), reads_to_weights[i].begin(),
                       reads_to_weights[i].end());
            ecs[ec] = txs;
        }
    }

    return ecs;
}

//' Calculate the effective transcript lengths. This is the mean number of
//' positions in the transcript the fragment could map to.
//'
//' @param txlengths The sequence lengths for each transcript.
//' @param rdlengths The length of all mapped fragments. Mapped length of the
//'  read after accounting for mismatches and indels.
//' @return The effective lengths.
// [[Rcpp::export]]
NumericVector effective_lengths(NumericVector txlengths,
                                NumericVector rdlengths) {
    NumericVector efflen(txlengths.size());
    NumericVector subset;
    double rdmean = mean(rdlengths);
    for (unsigned int i=0; i<txlengths.size(); ++i) {
        if (txlengths[i] >= rdmean) {
            efflen[i] = txlengths[i] - rdmean + 1;
        } else {
            subset = rdlengths[rdlengths <= txlengths[i]];
            efflen[i] = txlengths[i] - mean(subset) + 1;
        }
    }
    return round(efflen, 0);
}

//' Count transcripts using an Expectation Maximization (EM) algorithm.
//'
//' @param txreads A minimal matrix where each row corresponds to an
//'   alignment. The first column denotes transcript indices in [0, ntx-1] and
//'   the second column denotes read indices in [0, nr -1].
//' @param txlengths The sequence lengths for each transcript. Ideally those
//'   should be the effective transcript lengths meaning the overall number
//'   of possible alignment start positions in a transcript.
//' @param weights Weights for individual alignments.
//' @param ntx The total number of unique transcripts.
//' @param nr The total number of unique reads.
//' @param maxit Maximum number of EM iterations.
//' @param reltol The relative tolerance for convergence.
//' @param abstol The absolute tolerance for convergence.
//' @return A list with the following components.
//'     \describe{
//'      \item{p}{The length-normalized and read scaled transcript counts}
//'      \item{iterations}{The number of used EM iterations}
//'      \item{num_ecs}{The number of equivalence classes}
//'      \item{change}{The last osbserved absolute change in transcript counts}
//'     }
// [[Rcpp::export]]
List em_count(NumericMatrix txreads, NumericVector txlengths,
              NumericVector weights,
              int ntx, int nr, unsigned int maxit=1000,
              double reltol=0.01, double abstol=0.01) {
    std::vector<std::vector<int> > reads_to_txs(nr);
    std::vector<std::vector<int> > reads_to_weights(nr);
    NumericVector p(ntx, (double) nr / ntx);
    NumericVector pnew(ntx, 0.0);
    NumericVector change, cutoffs;
    std::vector<int> txs;
    LogicalVector is_large;
    double read_sum, count;
    std::unordered_map<std::string, std::vector<int> > ecs;
    unsigned int w_start;

    for (int i=0; i<txreads.nrow(); ++i) {
        if (txreads(i, 0) >= ntx) stop("wrong number of transcripts");
        if (txreads(i, 1) >= nr) stop("wrong number of reads");
        reads_to_txs[txreads(i, 1)].push_back(txreads(i, 0));
    }
    for (int i=0; i<txreads.nrow(); ++i) {
        reads_to_weights[txreads(i, 1)].push_back(weights[i]);
    }

    ecs = equivalence_classes(reads_to_txs, reads_to_weights);

    unsigned int k;
    for (k=0; k<maxit; ++k) {
        checkUserInterrupt();
        for (int i=0; i<ntx; ++i) {
            pnew[i] = 0.0;
        }
        for (std::pair<std::string, std::vector<int> > el : ecs) {
            read_sum = 0.0;
            txs = el.second;
            count = (double) txs[0];
            w_start = (txs.size() + 1)/2;
            for (unsigned int txi=1; txi<w_start; ++txi) {
                // This is P(txs) * mean_r(weight)
                read_sum += p[txs[txi]] * txs[txi + w_start - 1] / count;
            }
            for (unsigned int txi=1; txi<w_start; ++txi) {
                // Because P(txs) * mean_r(weight) * count / read_sum
                // = P(txs) * 1/count * sum_r(weights) * count / read_sum
                // = P(txs) * sum_r(weights) / read_sum
                pnew[txs[txi]] += p[txs[txi]] * txs[txi + w_start - 1] /
                                  read_sum;
            }
        }
        pnew = (pnew / txlengths);
        pnew = pnew * nr / sum(pnew);
        change = abs(pnew - p);
        cutoffs = reltol * pnew + abstol;
        if (is_true(all(change < cutoffs)) && k >= miniter) {
            break;
        }
        for (int i=0; i<ntx; ++i) p[i] = pnew[i];
    }

    return List::create(_["p"] = pnew,
                        _["iterations"] = k,
                        _["change"] = change,
                        _["ecs"] = ecs);
}
