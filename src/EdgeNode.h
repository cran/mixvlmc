#ifndef EDGENODE_H
#define EDGENODE_H
#include <Rcpp.h>
#include "ExtractionConditions.h"
#include "SubSequence.h"

//' @name EdgeNode
//' Suffix tree edge and node representation
//'
//' An instance of EdgeNode represents both the edge linking a node of the
//' suffix tree to another node and the target node itself.
//'
//' @noRd
class EdgeNode {
 public:
  // the parent node
  EdgeNode* parent;
  // start position in the NumericVector (included)
  int start;
  // end position in the NumericVector (excluded)
  int end;
  // nodes linked from this node, with indexed by the first value of the
  // corresponding edge
  std::unordered_map<int, EdgeNode*> children;
  // suffix link: if this node represents w and there is a in the alphabet
  // such that w=av and v is represented by the node P, that suffix points to P
  EdgeNode* suffix;

  // reverse map: if this map represents w and a is a key of the map, then
  // the value associated to a points to an EdgeNode that represents aw' where
  // w' is a maximal prefix of w such that aw' appears in the original sequence
  std::unordered_map<int, EdgeNode*>* reverse;

  // number of occurrences of the sequence represented by the path to this node
  // and the node itself. Computed post hoc.
  int total_count;

  // detailed counts of elements before the occurrences of the
  // sequence represented by this node.
  std::unordered_map<int, int>* counts;

  // positions of the sequence(s) represented by this node in
  // the original sequence
  std::vector<int>* positions;

  // Depth of the node, i.e. length of the sequence it represents
  // Computed post hoc.
  int depth;

  EdgeNode(EdgeNode* _parent, int s, int e);

  ~EdgeNode();

  EdgeNode* clone_no_relatives() const;

  EdgeNode* clone_only_counts() const;

  void setSuffix(EdgeNode* suffix_) { suffix = suffix_; }

  // compute the current edge length: most edges are open ended and thus
  // end is only defined by the length of the string processed so far. To
  // avoid repeated updates to end, we set it to the length of the NumericVector
  // to process. So it brings information only when it is smaller than current.
  int edge_length(int current) const {
    return std::min(end, current + 1) - start;
  }

  // compute the final edge length
  int edge_length() const { return end - start; }

  // compute the current edge label
  std::string edge_label(const Rcpp::IntegerVector& x, int current) const;

  void print_tree(std::string pre,
                  const Rcpp::IntegerVector& x,
                  int cend) const;

  // computes the number of occurrences of the string represented by a node.
  // This is to be called only on a full suffix tree where
  // each leaf corresponds to a suffix ended by a sentinel.
  void compute_total_count();

  // compute the counts of the values
  // that are before each instance of subsequence represented
  // by this node using an additional first term for the
  // longest suffix. Compute depth, total counts and maximal depth as a by
  // product. Keep the positions of each subsequence if keep_position is
  // true. This can lead to a quadratic memory occupation with respect
  // to the length of the original sequence.
  void compute_counts(int first,
                      const Rcpp::IntegerVector& x,
                      bool keep_position,
                      int cdepth,
                      int& mdepth);

  bool raw_contexts(const Rcpp::IntegerVector& x,
                    int nb_vals,
                              std::vector<int>& pre,
                              std::vector<const EdgeNode*>& subs,
                      std::vector<Rcpp::IntegerVector>& ctxs) const;

  // insert in subs (recursively) all subsequences of length at most
  // max_length that appear at least min_counts time in the original
  // sequence (x). It only_ctx is true, the subsequences must be contexts
  // in addition. Returns true if it is a subsequence.
  bool subsequences(const ExtractionConditions& when,
                    const ExtractionContent& what,
                    const Rcpp::IntegerVector& x,
                    int nb_vals,
                    std::vector<int>& pre,
                    std::vector<SubSequence*>& subs) const;

  // prune recursively the suffix tree based on:
  // - min_counts: a subsequence must appear at least min_counts times to
  //   be kept
  // - max_length: a subsequence cannot exceed max_length to be kept
  // - K: if a subsequence pass the previous two criteria, its preceding
  //   value distribution is tested against the one of its parent node. If
  //   the kl based criterion is smaller than K, the node is tentatively pruned.
  //   Only nodes without children are effectively pruned.
  //   This test is conducted only if K>0.
  // Updates maximal depth on the fly.
  // sentinels nodes are always removed during this pruning phase. The
  // number of contexts is computed on the fly (based on the above definition).
  bool prune(int min_counts,
             int max_length,
             double K,
             int nb_vals,
             int nx,
             int& mdepth,
             int& nb_ctx);

  // return a new collection of EdgeNodes (a tree) using the same
  // pruning rules as above. The clone is a deep one, fully independent
  // from the original tree.
  EdgeNode* clone_prune(int min_counts,
                        int max_length,
                        double K,
                        int nb_vals,
                        int nx,
                        int& mdepth,
                        int& nb_ctx) const;

  // compute significant cut off values
  double cutoff(std::set<double>& co) const;

  // build a flat representation of the tree where contexts are
  // stored in the tree_structure vector as the list of the
  // indices of there children. The tree_counts vector
  // contains the counts associated to each context.
  int flatten(const Rcpp::IntegerVector& x,
              int nb_vals,
              std::vector<Rcpp::IntegerVector>& tree_structure,
              std::vector<Rcpp::IntegerVector>& tree_counts) const;

  // make all nodes explicit
  void make_explicit(const Rcpp::IntegerVector& x);

  // compute the reverse maps recursively
  void compute_reverse(const Rcpp::IntegerVector& x,
                       const std::unordered_map<int, EdgeNode*>* parent_map);

  // counts the number of full nodes in the subtree
  int count_full_nodes(int nb_vals) const;

  // compute the log likelihood of the subtree when it is interpreted as
  // a vlmc
  double loglikelihood(int nb_vals) const;

  // clone the current node (and subtree) using less memory
  // by removing the positions
  EdgeNode* clone_trim() const ;
};

#endif
