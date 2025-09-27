#pragma once

#include <CL/sycl.hpp>

#include "sycl_meta.hpp"

namespace tf {

// ----------------------------------------------------------------------------
// syclGraph class
// ----------------------------------------------------------------------------

// class: syclGraph
class syclGraph : public CustomGraphBase {

  friend class syclNode;
  friend class syclTask;
  friend class syclFlow;
  friend class Taskflow;
  friend class Executor;

  constexpr static int OFFLOADED = 0x01;
  constexpr static int TOPOLOGY_CHANGED = 0x02;

  public:

    syclGraph() = default;
    ~syclGraph() = default;

    syclGraph(const syclGraph&) = delete;
    syclGraph(syclGraph&&);

    syclGraph& operator = (const syclGraph&) = delete;
    syclGraph& operator = (syclGraph&&);

    template <typename... ArgsT>
    syclNode* emplace_back(ArgsT&&...);

    bool empty() const;

    void clear();
    void dump(std::ostream&, const void*, const std::string&) const override final;

  private:

    int _state {0};

    std::vector<std::unique_ptr<syclNode>> _nodes;
};

// ----------------------------------------------------------------------------
// syclNode definitions
// ----------------------------------------------------------------------------

// class: syclNode
class syclNode {

  friend class syclGraph;
  friend class syclTask;
  friend class syclFlow;
  friend class Taskflow;
  friend class Executor;

  struct Empty {
  };

  struct CGH {

    std::function<void(sycl::handler&)> work;

    template <typename F>
    CGH(F&& func) : work {std::forward<F>(func)} {}
  };

  using handle_t = std::variant<
    Empty,
    CGH
  >;

  public:

  // variant index
  constexpr static auto EMPTY = get_index_v<Empty, handle_t>;
  constexpr static auto COMMAND_GROUP_HANDLER = get_index_v<CGH, handle_t>;

    syclNode() = delete;

    template <typename... ArgsT>
    syclNode(syclGraph&, ArgsT&&...);

  private:

    syclGraph& _graph;

    std::string _name;

    int _level;

    sycl::event _event;

    handle_t _handle;

    SmallVector<syclNode*> _successors;
    SmallVector<syclNode*> _dependents;

    void _precede(syclNode*);
};

// ----------------------------------------------------------------------------
// syclNode definitions
// ----------------------------------------------------------------------------

// Constructor
template <typename... ArgsT>
syclNode::syclNode(syclGraph& g, ArgsT&&... args) :
  _graph  {g},
  _handle {std::forward<ArgsT>(args)...} {
}

// Procedure: _precede
inline void syclNode::_precede(syclNode* v) {
  _graph._state |= syclGraph::TOPOLOGY_CHANGED;
  _successors.push_back(v);
  v->_dependents.push_back(this);
}

// ----------------------------------------------------------------------------
// syclGraph definitions
// ----------------------------------------------------------------------------

// Move constructor
inline syclGraph::syclGraph(syclGraph&& g) :
  _nodes {std::move(g._nodes)} {

  assert(g._nodes.empty());
}

// Move assignment
inline syclGraph& syclGraph::operator = (syclGraph&& rhs) {

  // lhs
  _nodes = std::move(rhs._nodes);

  assert(rhs._nodes.empty());

  return *this;
}

// Function: empty
inline bool syclGraph::empty() const {
  return _nodes.empty();
}

// Procedure: clear
inline void syclGraph::clear() {
  _state = syclGraph::TOPOLOGY_CHANGED;
  _nodes.clear();
}

// Function: emplace_back
template <typename... ArgsT>
syclNode* syclGraph::emplace_back(ArgsT&&... args) {

  _state |= syclGraph::TOPOLOGY_CHANGED;

  auto node = std::make_unique<syclNode>(std::forward<ArgsT>(args)...);
  _nodes.emplace_back(std::move(node));
  return _nodes.back().get();

  // TODO: object pool

  //auto node = new syclNode(std::forward<ArgsT>(args)...);
  //_nodes.push_back(node);
  //return node;
}

// Procedure: dump the graph to a DOT format
inline void syclGraph::dump(
  std::ostream& os, const void* root, const std::string& root_name
) const {

  // recursive dump with stack
  std::stack<std::tuple<const syclGraph*, const syclNode*, int>> stack;
  stack.push(std::make_tuple(this, nullptr, 1));

  int pl = 0;

  while(!stack.empty()) {

    auto [graph, parent, l] = stack.top();
    stack.pop();

    for(int i=0; i<pl-l+1; i++) {
      os << "}\n";
    }

    if(parent == nullptr) {
      if(root) {
        os << "subgraph cluster_p" << root << " {\nlabel=\"syclFlow: ";
        if(root_name.empty()) os << 'p' << root;
        else os << root_name;
        os << "\";\n" << "color=\"red\"\n";
      }
      else {
        os << "digraph syclFlow {\n";
      }
    }
    else {
      os << "subgraph cluster_p" << parent << " {\nlabel=\"syclSubflow: ";
      if(parent->_name.empty()) os << 'p' << parent;
      else os << parent->_name;
      os << "\";\n" << "color=\"purple\"\n";
    }

    for(auto& v : graph->_nodes) {

      os << 'p' << v.get() << "[label=\"";
      if(v->_name.empty()) {
        os << 'p' << v.get() << "\"";
      }
      else {
        os << v->_name << "\"";
      }
      os << "];\n";

      for(const auto s : v->_successors) {
        os << 'p' << v.get() << " -> " << 'p' << s << ";\n";
      }

      if(v->_successors.size() == 0) {
        if(parent == nullptr) {
          if(root) {
            os << 'p' << v.get() << " -> p" << root << ";\n";
          }
        }
        else {
          os << 'p' << v.get() << " -> p" << parent << ";\n";
        }
      }
    }

    // set the previous level
    pl = l;
  }

  for(int i=0; i<pl; i++) {
    os << "}\n";
  }

}


}  // end of namespace tf -----------------------------------------------------


