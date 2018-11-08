/*
 * Copyright 2018 Christian Diener <mail[at]cdiener.com>
 *
 * Apache License 2.0.
 */

#include <Rcpp.h>
#include <vector>
#include <numeric>

using namespace Rcpp;

// [[Rcpp::export]]
List em_count(NumericMatrix txreads, NumericVector txlengths,
                       int ntx, int nr, int maxit=1000, double cutoff=0.01) {
    std::vector<std::vector<int> > reads_to_txs(nr);
    NumericVector p(ntx, 1.0 / ntx);
    NumericVector pnew(ntx, 0.0);
    double read_sum;

    for (unsigned int i=0; i<txreads.nrow(); ++i) {
        if (txreads(i, 0) >= ntx) stop("wrong number of transcripts");
        if (txreads(i, 1) >= nr) stop("wrong number of reads");
        reads_to_txs[txreads(i, 1)].push_back(txreads(i, 0));
    }

    unsigned int k;
    for (k=0; k<maxit; ++k) {
        for (unsigned int i=0; i<ntx; ++i) {
            pnew[i] = 0.0;
        }
        for (unsigned int i=0; i<nr; i++) {
            read_sum = 0.0;
            for (unsigned int txi=0; txi<reads_to_txs[i].size(); ++txi) {
                read_sum += p[reads_to_txs[i][txi]];
            }
            for (unsigned int txi=0; txi<reads_to_txs[i].size(); ++txi) {
                pnew[reads_to_txs[i][txi]] += p[reads_to_txs[i][txi]] / read_sum;
            }
        }
        pnew = pnew / nr;
        pnew = (pnew / txlengths);
        pnew = pnew / sum(pnew);
        if (max(abs(pnew - p)) < cutoff) {
            break;
        }
        for (unsigned int i=0; i<ntx; ++i) p[i] = pnew[i];
    }

    return List::create(_["p"] = pnew, _["iterations"] = k);
}
