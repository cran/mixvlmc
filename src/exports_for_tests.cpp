#include "utils.h"
using namespace Rcpp;

//[[Rcpp::export]]
double kl_crit(IntegerVector p, IntegerVector q) {
  if(p.size() != q.size()) {
    stop("Cannot use kl_crit with vectors of different lengths");
  }
  int p_sum = sum(p);
  int q_sum = sum(q);
  int nx = p.size();
  std::unordered_map<int, int>* p_counts = new std::unordered_map<int, int>();
  std::unordered_map<int, int>* q_counts = new std::unordered_map<int, int>();
  for(int i = 0; i < nx; i++) {
    if(p[i] > 0) {
      (*p_counts)[i] = p[i];
    }
    if(q[i] > 0) {
      (*q_counts)[i] = q[i];
    }
  }
  double res = kl_criterion(p_counts, p_sum, q_counts, q_sum);
  delete p_counts;
  delete q_counts;
  return res;
}

//[[Rcpp::export]]
IntegerVector mixvlmc_sample(IntegerVector p, int n) {
  int p_sum = sum(p);
  std::unordered_map<int, int>* p_counts = new std::unordered_map<int, int>();
  int nx = p.size();
  for(int i = 0; i < nx; i++) {
    if(p[i] > 0) {
      (*p_counts)[i] = p[i];
    }
  }
  RNGScope scope;
  IntegerVector res(n);
  for(int i = 0; i < n; i++) {
    res[i] = sample(p_counts, p_sum);
  }
  delete p_counts;
  return res;
}

//[[Rcpp::export]]
IntegerVector mixvlmc_sample2(IntegerVector p, int n) {
  int p_sum = sum(p);
  std::unordered_map<int, int>* p_counts = new std::unordered_map<int, int>();
  int nx = p.size();
  for(int i = 0; i < nx; i++) {
    if(p[i] > 0) {
      (*p_counts)[i] = p[i];
    }
  }
  RNGScope scope;
  IntegerVector res(n);
  for(int i = 0; i < n; i++) {
    res[i] = sample2(p_counts, nx - 1, p_sum);
  }
  delete p_counts;
  return res;
}
