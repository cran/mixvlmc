#include <algorithm>
#include "EdgeNode.h"
#include "utils.h"

using namespace Rcpp;

std::string to_string(const std::vector<int>& pre) {
  std::string res = "";
  for(int k : pre) {
    res += std::to_string(k);
  }
  return res;
}

EdgeNode::EdgeNode(EdgeNode* _parent, int s, int e)
    : parent(_parent),
      start(s),
      end(e),
      suffix(nullptr),
      reverse(nullptr),
      total_count(0),
      counts(nullptr),
      positions(nullptr),
      depth(0) {}

EdgeNode::~EdgeNode() {
  for(auto child : children) {
    delete child.second;
  }
  if(reverse != nullptr) {
    delete reverse;
  }
  if(counts != nullptr) {
    delete counts;
  }
  if(positions != nullptr) {
    delete positions;
  }
}

EdgeNode* EdgeNode::clone_no_relatives() const {
  EdgeNode* result = clone_only_counts();
  if(positions != nullptr) {
    result->positions =
        new std::vector<int>(positions->begin(), positions->end());
  }
  return result;
}

EdgeNode* EdgeNode::clone_only_counts() const {
  EdgeNode* result = new EdgeNode(nullptr, start, end);
  result->total_count = total_count;
  if(counts != nullptr) {
    result->counts =
        new std::unordered_map<int, int>(counts->begin(), counts->end());
  }
  result->depth = depth;
  return result;
}

// # nocov start
std::string EdgeNode::edge_label(const IntegerVector& x, int current) const {
  std::string res = "";
  int size = std::min(end, current + 1);
  for(int i = start; i < size; i++) {
    if(i < x.size()) {
      res += std::to_string(x[i]);
    } else {
      res += "$";
    }
  }
  return res;
}

void EdgeNode::print_tree(std::string pre,
                          const IntegerVector& x,
                          int cend) const {
  Rcout << pre << this << " ~ " << depth;
  Rcout << "\n";
  if(suffix != nullptr) {
    Rcout << pre << "sf " << suffix << "\n";
  }
  if(counts != nullptr) {
    Rcout << pre << counts_to_string(counts) << "\n";
  }
  if(reverse != nullptr) {
    for(auto rev : *reverse) {
      Rcout << pre << rev.first << " -> " << rev.second << "\n";
    }
  }
  if(positions != nullptr) {
    Rcout << pre << "{";
    int np = (int)positions->size() - 1;
    for(int i = 0; i < np; i++) {
      Rcout << (*positions)[i] << ", ";
    }
    Rcout << (*positions)[np] << "}\n";
  }
  for(auto child : children) {
    Rcout << pre << " [" << child.first << "] -> "
          << child.second->edge_label(x, cend);
    Rcout << " (" << (child.second->start) << " - " << (child.second->end)
          << ")\n";
    child.second->print_tree(pre + "  ", x, cend);
  }
}
// # nocov end

void EdgeNode::compute_total_count() {
  if(children.size() == 0) {
    total_count = 1;
  } else {
    total_count = 0;
    for(auto child : children) {
      child.second->compute_total_count();
      total_count += child.second->total_count;
    }
  }
}

void EdgeNode::compute_counts(int first,
                              const Rcpp::IntegerVector& x,
                              bool keep_position,
                              int cdepth,
                              int& mdepth) {
  depth = cdepth + edge_length();
  if(depth > mdepth) {
    mdepth = depth;
  }
  counts = new std::unordered_map<int, int>{};
  if(keep_position) {
    positions = new std::vector<int>{};
  }
  if(children.size() == 0) {
    // this is a leaf, therefore a suffix which has a single
    // preceding character at position nx-depth
    int pos = x.size() - depth;
    if(keep_position) {
      positions->push_back(pos + 1);
    }
    int val;
    if(pos >= 0) {
      val = x[pos];
    } else {
      val = first;
    }
    (*counts)[val] = 1;
    total_count = 1;
  } else {
    total_count = 0;
    for(auto child : children) {
      child.second->compute_counts(first, x, keep_position, depth, mdepth);
      total_count += child.second->total_count;
      if(keep_position) {
        positions->insert(positions->end(), child.second->positions->begin(),
                          child.second->positions->end());
      }
      // update counts
      for(auto count : *(child.second->counts)) {
        auto current = counts->find(count.first);
        if(current != counts->end()) {
          current->second += count.second;
        } else {
          (*counts)[count.first] = count.second;
        }
      }
    }
  }
}

bool EdgeNode::raw_contexts(const IntegerVector& x,
                            int nb_vals,
                            std::vector<int>& pre,
                            std::vector<const EdgeNode*>& subs,
                            std::vector<IntegerVector>& ctxs) const {
  bool is_sub = false;
  // let us first handle the root
  if(start < 0) {
    int nb_sub = 0;
    for(auto child : children) {
      if(child.first >= 0 &&
         child.second->raw_contexts(x, nb_vals, pre, subs, ctxs)) {
        nb_sub++;
      }
    }
    if(nb_sub < nb_vals) {
      // we consider the empty context
      subs.push_back(this);
      ctxs.push_back(IntegerVector(pre.begin(), pre.end()));
      return true;
    } else {
      return false;
    }
  } else {
    // non root
    size_t before = pre.size();
    // if the edge length is larger than one, then we are sure to keep the
    // subsequence as a context (single child)
    if(edge_length() > 1) {
      is_sub = true;
    }
    // for the raw context approach, we postpone any implicit node
    // handling
    // we do not keep the sentinel value if it is present
    int the_end = std::min((int)x.size(), end) - 1;
    for(int i = start; i < the_end; i++) {
      pre.push_back(x[i]);
      subs.push_back(this);
      ctxs.push_back(IntegerVector(pre.begin(), pre.end()));
    }
    // add the final element of the edge
    pre.push_back(x[the_end]);
    // let us handle the children
    int nb_sub = 0;
    is_sub = true;
    for(auto child : children) {
      // do not recurse in sentinel nodes
      if(child.first >= 0 &&
         child.second->raw_contexts(x, nb_vals, pre, subs, ctxs)) {
        nb_sub++;
      }
    }
    // the current subsequence is a context if some of the potential longer
    // subsequences have not
    // been included
    if(nb_sub < nb_vals) {
      subs.push_back(this);
      ctxs.push_back(IntegerVector(pre.begin(), pre.end()));
    }
    pre.resize(before);
  }
  return is_sub;
}

bool EdgeNode::subsequences(const ExtractionConditions& when,
                            const ExtractionContent& what,
                            const Rcpp::IntegerVector& x,
                            int nb_vals,
                            std::vector<int>& pre,
                            std::vector<SubSequence*>& subs) const {
  bool is_sub = false;
  if(total_count >= when.min_counts) {
    // let us first handle the root
    if(start < 0) {
      int nb_sub = 0;
      for(auto child : children) {
        if(child.first >= 0 &&
           child.second->subsequences(when, what, x, nb_vals, pre, subs)) {
          nb_sub++;
        }
      }
      if(when.only_ctx && nb_sub < nb_vals) {
        // we consider the empty context
        subs.push_back(new SubSequence(pre, this, what));
        return true;
      } else {
        return false;
      }
    } else {
      // non root
      size_t before = pre.size();
      // if the edge length is larger than one, then we are sure to keep the
      // subsequence as a context (single child)
      if(edge_length() > 1) {
        is_sub = true;
      }
      // first progress on the edge if needed (we do not include the sentinel)
      int the_end = std::min((int)x.size(), end) - 1;
      bool ml_reached = false;
      for(int i = start; i < the_end; i++) {
        pre.push_back(x[i]);
        if(pre.size() <= (size_t)when.max_length) {
          // always contexts
          subs.push_back(new SubSequence(pre, this, what));
        } else {
          ml_reached = true;
          break;
        }
      }
      if(!ml_reached) {
        // add the final element of the edge
        pre.push_back(x[the_end]);
        // conditional inclusion only, see below
        if(pre.size() < (size_t)when.max_length) {
          // we may have longer subsequences
          // we will keep this subsequence
          is_sub = true;
          int nb_sub = 0;
          for(auto child : children) {
            // do not recurse in sentinel nodes
            if(child.first >= 0 &&
               child.second->subsequences(when, what, x, nb_vals, pre, subs)) {
              nb_sub++;
            }
          }
          // the current subsequence is a context if all subsequences should be
          // added or if some of the potential longer subsequences have not
          // been included
          if((!when.only_ctx) || nb_sub < nb_vals) {
            // this is context/subsequence
            subs.push_back(new SubSequence(pre, this, what));
          }
        } else if(pre.size() == (size_t)when.max_length) {
          // no child tested
          subs.push_back(new SubSequence(pre, this, what));
          is_sub = true;
        }
      }
      pre.resize(before);
    }
  }
  // if total_count is too small, this is also the case for the children
  // so we stop recursion here
  return is_sub;
}

bool EdgeNode::prune(int min_counts,
                     int max_length,
                     double K,
                     int nb_vals,
                     int nx,
                     int& mdepth,
                     int& nb_ctx) {
  if(total_count >= min_counts) {
    if(depth > max_length) {
      // depth based pruning
      // we have to insert explicit nodes if a part of the edge is kept
      if(depth - edge_length() + 1 > max_length) {
        // nothing to keep
        return true;
      } else {
        // part of the edge should be kept
        // we do that by reducing the end counter.
        // This depends on K. If K>0, we will never keep
        // implicit nodes, so we can directly set end to start+1 (and test for
        // criterion based prunability)
        // if K<=0, we may keep implicit nodes up to the maximum depth.
        // in any case, we remove all the children
        for(auto child : children) {
          delete child.second;
        }
        children.clear();
        int allowance;
        if(K <= 0) {
          allowance = max_length - depth + edge_length();
        } else {
          if(parent != nullptr) {
            double crit = kl_criterion(counts, total_count, parent->counts,
                                       parent->total_count);
            if(crit < K) {
              // prune
              return true;
            }
          }
          allowance = 1;
        }
        depth = depth - edge_length() + allowance;
        if(depth > mdepth) {
          mdepth = depth;
        }
        end = start + allowance;
        nb_ctx += allowance;
        return false;
      }
    } else {
      // recursive processing
      int nb_sub = 0;
      for(auto child = children.begin(); child != children.end();) {
        if(child->first < 0) {
          // this a sentinel node, we can remove it safely
          delete child->second;
          child = children.erase(child);
        } else {
          bool result = child->second->prune(min_counts, max_length, K, nb_vals,
                                             nx, mdepth, nb_ctx);
          if(result) {
            delete child->second;
            child = children.erase(child);
          } else {
            ++child;
            nb_sub++;
          }
        }
      }
      if(nb_sub == 0 && K > 0 && parent != nullptr) {
        // potential pruning
        double crit = kl_criterion(counts, total_count, parent->counts,
                                   parent->total_count);
        if(crit < K) {
          // prune
          return true;
        } else {
          // even if we do not prune, when nb_sub==0, we do not keep implicit
          // nodes
          depth = depth - edge_length() + 1;
          end = start + 1;
        }
      }
      if(edge_length() > 1) {
        // count intermediate contexts
        if(end > nx) {
          nb_ctx += edge_length() - 2;
        } else {
          nb_ctx += edge_length() - 1;
        }
      }
      if(nb_sub < nb_vals) {
        // potential
        nb_ctx++;
      }
      if(depth > mdepth) {
        mdepth = depth;
      }
      return false;
    }
  } else {
    // count based pruning
    return true;
  }
}

EdgeNode* EdgeNode::clone_prune(int min_counts,
                                int max_length,
                                double K,
                                int nb_vals,
                                int nx,
                                int& mdepth,
                                int& nb_ctx) const {
  if(total_count >= min_counts) {
    if(depth > max_length) {
      // depth based pruning
      // we have to insert explicit nodes if a part of the edge is kept
      if(depth - edge_length() + 1 > max_length) {
        // nothing to keep
        return nullptr;
      } else {
        // part of the edge should be kept
        // we do that by reducing the end counter
        // This depends on K. If K>0, we will never keep
        // implicit nodes, so we can directly set end to start+1 (and test for
        // criterion based prunability)
        // if K<=0, we may keep implicit nodes up to the maximum depth.
        // in any case, the node will have no children
        int allowance;
        if(K <= 0) {
          allowance = max_length - depth + edge_length();
        } else {
          if(parent != nullptr) {
            double crit = kl_criterion(counts, total_count, parent->counts,
                                       parent->total_count);
            if(crit < K) {
              // ctiretion based prune
              return nullptr;
            }
          }
          allowance = 1;
        }
        EdgeNode* result = clone_no_relatives();
        result->end = start + allowance;
        result->depth = depth - edge_length() + allowance;
        if(result->depth > mdepth) {
          mdepth = result->depth;
        }
        nb_ctx += allowance;
        return result;
      }
    } else {
      EdgeNode* result = clone_no_relatives();
      // recursive processing
      int nb_sub = 0;
      for(auto child = children.begin(); child != children.end(); ++child) {
        if(child->first >= 0) {
          // we ignore sentinel nodes
          EdgeNode* new_child = child->second->clone_prune(
              min_counts, max_length, K, nb_vals, nx, mdepth, nb_ctx);
          if(new_child != nullptr) {
            result->children[child->first] = new_child;
            new_child->parent = result;
            nb_sub++;
          }
        }
      }
      // we may still prune the node
      if(nb_sub == 0 && K > 0 && parent != nullptr) {
        // potential pruning
        double crit = kl_criterion(counts, total_count, parent->counts,
                                   parent->total_count);
        if(crit < K) {
          // prune
          delete result;
          return nullptr;
        } else {
          // even if we do not prune, when nb_sub==0, we do not keep implicit
          // nodes
          result->depth = depth - edge_length() + 1;
          result->end = start + 1;
        }
      }
      if(result->edge_length() > 1) {
        // count intermediate contexts
        if(result->end > nx) {
          nb_ctx += result->edge_length() - 2;
        } else {
          nb_ctx += result->edge_length() - 1;
        }
      }
      if(nb_sub < nb_vals) {
        nb_ctx++;
      }
      if(result->depth > mdepth) {
        mdepth = result->depth;
      }
      return result;
    }
  } else {
    // count based pruning
    return nullptr;
  }
}

double EdgeNode::cutoff(std::set<double>& co) const {
  double local_kl = 0;
  if(parent != nullptr) {
    local_kl =
        kl_criterion(counts, total_count, parent->counts, parent->total_count);
  }
  double max_kl_child = 0;
  for(auto child : children) {
    if(child.first >= 0) {
      double child_kl = child.second->cutoff(co);
      if(child_kl > max_kl_child) {
        max_kl_child = child_kl;
      }
    }
  }
  if(local_kl > max_kl_child) {
    co.insert(local_kl);
    return local_kl;
  } else {
    return max_kl_child;
  }
}

int EdgeNode::flatten(const Rcpp::IntegerVector& x,
                      int nb_vals,
                      std::vector<Rcpp::IntegerVector>& tree_structure,
                      std::vector<Rcpp::IntegerVector>& tree_counts) const {
  IntegerVector f_by = map_to_counts(counts, nb_vals - 1);
  int pos = tree_structure.size();
  int the_end = end;
  if(the_end > x.size()) {
    // sentinel, no children
    the_end = x.size();
  }
  int sub_pos = pos;
  for(int i = start; i < the_end - 1; i++) {
    IntegerVector t_children(nb_vals, R_NaInt);
    tree_counts.push_back(f_by);
    t_children[x[i + 1]] = sub_pos + 2;    // +1 for R
    tree_structure.push_back(t_children);  // in position sub_pos
    sub_pos++;
  }
  if(children.size() > 0) {
    IntegerVector t_children(nb_vals, R_NaInt);
    tree_structure.push_back(t_children);
    tree_counts.push_back(f_by);
    int child_pos = tree_structure.size() - 1;
    for(auto child : children) {
      if(child.first >= 0) {
        t_children[child.first] =
            child.second->flatten(x, nb_vals, tree_structure, tree_counts);
      }
    }
    tree_structure[child_pos] = t_children;
  } else {
    IntegerVector empty{};
    tree_structure.push_back(empty);
    tree_counts.push_back(f_by);
  }
  return pos + 1;  // +1 for R
}

void EdgeNode::make_explicit(const Rcpp::IntegerVector& x) {
  if(edge_length() > 1) {
    // insert explicit nodes
    int el = edge_length() - 1;
    int base_depth = parent->depth + 1;
    EdgeNode* previous = parent;
    int from = x[start];
    int pos = start;
    for(int i = 0; i < el; i++) {
      EdgeNode* node = new EdgeNode(previous, pos, pos + 1);
      previous->children[from] = node;
      node->total_count = total_count;
      node->depth = base_depth + i;
      if(counts != nullptr) {
        node->counts =
            new std::unordered_map<int, int>(counts->begin(), counts->end());
      }
      if(positions != nullptr) {
        node->positions =
            new std::vector<int>(positions->begin(), positions->end());
      }
      previous = node;
      pos++;
      if(pos < x.size()) {
        from = x[pos];
      } else {
        from = -1;  // sentinel
      }
    }
    start = end - 1;
    parent = previous;
    previous->children[from] = this;
  }
  for(auto child : children) {
    if(child.first >= 0) {
      child.second->make_explicit(x);
    }
  }
}

void EdgeNode::compute_reverse(
    const Rcpp::IntegerVector& x,
    const std::unordered_map<int, EdgeNode*>* parent_map) {
  // should check for existing map
  reverse = new std::unordered_map<int, EdgeNode*>{};
  if(start < x.size()) {
    for(auto prev : *parent_map) {
      if(depth == prev.second->depth) {
        // full matching, maybe be increased
        auto child = prev.second->children.find(x[start]);
        if(child != prev.second->children.end()) {
          // success, we can extend the match
          (*reverse)[prev.first] = child->second;
        } else {
          // nope
          (*reverse)[prev.first] = prev.second;
        }
      } else {
        // partial matching, no luck
        (*reverse)[prev.first] = prev.second;
      }
    }
  }
  for(auto child : children) {
    if(child.first >= 0) {
      child.second->compute_reverse(x, reverse);
    }
  }
}

int EdgeNode::count_full_nodes(int nb_vals) const {
  int local;
  if((int)children.size() == nb_vals) {
    local = 1;
  } else {
    local = 0;
  }
  for(auto child : children) {
    if(child.first >= 0) {
      local += child.second->count_full_nodes(nb_vals);
    }
  }
  return local;
}

double EdgeNode::loglikelihood(int nb_vals) const {
  std::vector<int> lcounts(nb_vals, 0);
  for(auto count : *counts) {
    if(count.second > 0) {
      lcounts[count.first] = count.second;
    }
  }
  double ll = 0.0;
  for(auto child : children) {
    if(child.first >= 0) {
      ll += child.second->loglikelihood(nb_vals);
      for(auto count : *(child.second->counts)) {
        lcounts[count.first] -= count.second;
      }
    }
  }
  for(int i = 0; i < nb_vals; i++) {
    if(lcounts[i] > 0) {
      ll += lcounts[i] * log(((double)(*counts)[i]) / total_count);
    }
  }
  return ll;
}

EdgeNode* EdgeNode::clone_trim() const {
  EdgeNode* result = clone_only_counts();
  for(auto child : children) {
    EdgeNode* new_child = child.second->clone_trim();
    new_child->parent = result;
    result->children[child.first] = new_child;
  }
  return result;
}
