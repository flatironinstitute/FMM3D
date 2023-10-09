#ifndef _SCTL_TREE_
#define _SCTL_TREE_

#include <sctl/common.hpp>
#include SCTL_INCLUDE(comm.hpp)
#include SCTL_INCLUDE(morton.hpp)
#include SCTL_INCLUDE(vtudata.hpp)
#include SCTL_INCLUDE(ompUtils.hpp)

#include <string>
#include <vector>
#include <algorithm>

namespace SCTL_NAMESPACE {

template <Integer DIM> class Tree {
  public:

    struct NodeAttr {
      unsigned char Leaf : 1, Ghost : 1;
    };

    struct NodeLists {
      Long p2n;
      Long parent;
      Long child[1 << DIM];
      Long nbr[sctl::pow<DIM,Integer>(3)];
    };

    static constexpr Integer Dim();

    Tree(const Comm& comm_ = Comm::Self());

    ~Tree();

    const Vector<Morton<DIM>>& GetPartitionMID() const;
    const Vector<Morton<DIM>>& GetNodeMID() const;
    const Vector<NodeAttr>& GetNodeAttr() const;
    const Vector<NodeLists>& GetNodeLists() const;
    const Comm& GetComm() const;

    template <class Real> void UpdateRefinement(const Vector<Real>& coord, Long M = 1, bool balance21 = 0, bool periodic = 0);

    template <class ValueType> void AddData(const std::string& name, const Vector<ValueType>& data, const Vector<Long>& cnt);

    template <class ValueType> void GetData(Vector<ValueType>& data, Vector<Long>& cnt, const std::string& name) const;

    template <class ValueType> void ReduceBroadcast(const std::string& name);

    template <class ValueType> void Broadcast(const std::string& name);

    void DeleteData(const std::string& name);

    void WriteTreeVTK(std::string fname, bool show_ghost = false) const;

  protected:

    void GetData_(Iterator<Vector<char>>& data, Iterator<Vector<Long>>& cnt, const std::string& name);

    static void scan(Vector<Long>& dsp, const Vector<Long>& cnt);

    template <typename A, typename B> struct SortPair {
      int operator<(const SortPair<A, B> &p1) const { return key < p1.key; }
      A key;
      B data;
    };

  private:

    Vector<Morton<DIM>> mins;
    Vector<Morton<DIM>> node_mid;
    Vector<NodeAttr> node_attr;
    Vector<NodeLists> node_lst;

    std::map<std::string, Vector<char>> node_data;
    std::map<std::string, Vector<Long>> node_cnt;

    Vector<Morton<DIM>> user_mid;
    Vector<Long> user_cnt;

    Comm comm;
};

template <class Real, Integer DIM, class BaseTree = Tree<DIM>> class PtTree : public BaseTree {
  public:

    PtTree(const Comm& comm = Comm::Self());

    ~PtTree();

    void UpdateRefinement(const Vector<Real>& coord, Long M = 1, bool balance21 = 0, bool periodic = 0);

    void AddParticles(const std::string& name, const Vector<Real>& coord);

    void AddParticleData(const std::string& data_name, const std::string& particle_name, const Vector<Real>& data);

    void GetParticleData(Vector<Real>& data, const std::string& data_name) const;

    void DeleteParticleData(const std::string& data_name);

    void WriteParticleVTK(std::string fname, std::string data_name, bool show_ghost = false) const;

    static void test() {
      Long N = 100000;
      Vector<Real> X(N*DIM), f(N);
      for (Long i = 0; i < N; i++) { // Set coordinates (X), and values (f)
        f[i] = 0;
        for (Integer k = 0; k < DIM; k++) {
          X[i*DIM+k] = pow<3>(drand48()*2-1.0)*0.5+0.5;
          f[i] += X[i*DIM+k]*k;
        }
      }

      PtTree<Real,DIM> tree;
      tree.AddParticles("pt", X);
      tree.AddParticleData("pt-value", "pt", f);
      tree.UpdateRefinement(X, 1000); // refine tree with max 1000 points per box.

      { // manipulate tree node data
        const auto& node_lst = tree.GetNodeLists(); // Get interaction lists
        //const auto& node_mid = tree.GetNodeMID();
        //const auto& node_attr = tree.GetNodeAttr();

        // get point values and count for each node
        Vector<Real> value;
        Vector<Long> cnt, dsp;
        tree.GetData(value, cnt, "pt-value");

        // compute the dsp (the point offset) for each node
        dsp.ReInit(cnt.Dim()); dsp = 0;
        omp_par::scan(cnt.begin(), dsp.begin(), cnt.Dim());

        Long node_idx = 0;
        for (Long i = 0; i < cnt.Dim(); i++) { // find the tree node with maximum points
          if (cnt[node_idx] < cnt[i]) node_idx = i;
        }

        for (Long j = 0; j < cnt[node_idx]; j++) { // for this node, set all pt-value to -1
          value[dsp[node_idx]+j] = -1;
        }

        for (const Long nbr_idx : node_lst[node_idx].nbr) { // loop over the neighbors and set pt-value to 2
          if (nbr_idx >= 0 && nbr_idx != node_idx) {
            for (Long j = 0; j < cnt[nbr_idx]; j++) {
              value[dsp[nbr_idx]+j] = 2;
            }
          }
        }
      }

      // Generate visualization
      tree.WriteParticleVTK("pt", "pt-value");
      tree.WriteTreeVTK("tree");
    }

  private:

    std::map<std::string, Long> Nlocal;
    std::map<std::string, Vector<Morton<DIM>>> pt_mid;
    std::map<std::string, Vector<Long>> scatter_idx;
    std::map<std::string, std::string> data_pt_name;
};

}

#include SCTL_INCLUDE(tree.txx)

#endif //_SCTL_TREE_
