#pragma once

/// @file KDTree.hpp
/// @author J. Frederico Carvalho
///
/// This is an adaptation of the KD-tree implementation in rosetta code
///  https://rosettacode.org/wiki/K-d_tree
/// It is a reimplementation of the C code using C++.
/// It also includes a few more queries than the original

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <vector>

/// The point type (vector of double precision floats)
using point_t = std::vector<double>;

/// Array of indices
using indexArr = std::vector<size_t>;

/// Pair of point and Index
using pointIndex = typename std::pair<std::vector<double>, size_t>;

class KDNode {
  public:
    using KDNodePtr = std::shared_ptr<KDNode>;
    size_t index;
    point_t x;
    KDNodePtr left;
    KDNodePtr right;

    // initializer
    KDNode();
    KDNode(point_t const&, size_t const&, KDNodePtr const&, KDNodePtr const&);
    KDNode(pointIndex const&, KDNodePtr const&, KDNodePtr const&);
    ~KDNode();

    // getter
    double coord(size_t const&);

    // conversions
    explicit operator bool();
    explicit operator point_t();
    explicit operator size_t();
    explicit operator pointIndex();
};

using KDNodePtr = std::shared_ptr<KDNode>;

KDNodePtr NewKDNodePtr();

// square euclidean distance
inline double dist2(point_t const&, point_t const&);
inline double dist2(KDNodePtr const&, KDNodePtr const&);

// Need for sorting
class comparer {
  public:
    size_t idx;
    explicit comparer(size_t idx_);
    inline bool compare_idx(std::pair<std::vector<double>, size_t> const&, //
                            std::pair<std::vector<double>, size_t> const&  //
    );
};

using pointIndexArr = typename std::vector<pointIndex>;

inline void sort_on_idx(pointIndexArr::iterator const&, //
                        pointIndexArr::iterator const&, //
                        size_t idx);

using pointVec = std::vector<point_t>;

class KDTree {

  public:
    KDTree() = default;

    /// Build a KDtree
    explicit KDTree(pointVec point_array);

    /// Get the point which lies closest to the input point.
    /// @param pt input point.
    point_t nearest_point(point_t const& pt);

    /// Get the index of the point which lies closest to the input point.
    ///
    /// @param pt input point.
    size_t nearest_index(point_t const& pt);

    /// Get the point and its index which lies closest to the input point.
    ///
    /// @param pt input point.
    pointIndex nearest_pointIndex(point_t const& pt);

    /// Get both the point and the index of the points closest to the input
    /// point.
    ///
    /// @param pt input point.
    /// @param num_nearest Number of nearest points to return.
    ///
    /// @returns a vector containing the points and their respective indices
    /// which are at a distance smaller than rad to the input point.
    pointIndexArr nearest_pointIndices(point_t const& pt,
                                       size_t const& num_nearest);

    /// Get the nearest set of points to the given input point.
    ///
    /// @param pt input point.
    /// @param num_nearest Number of nearest points to return.
    ///
    /// @returns a vector containing the points which are at a distance smaller
    /// than rad to the input point.
    pointVec nearest_points(point_t const& pt, size_t const& num_nearest);

    /// Get the indices of points closest to the input point.
    ///
    /// @param pt input point.
    /// @param num_nearest Number of nearest points to return.
    ///
    /// @returns a vector containing the indices of the points which are at a
    /// distance smaller than rad to the input point.
    indexArr nearest_indices(point_t const& pt, size_t const& num_nearest);

    /// Get both the point and the index of the points which are at a distance
    /// smaller than the input radius to the input point.
    ///
    /// @param pt input point.
    /// @param rad input radius.
    ///
    /// @returns a vector containing the points and their respective indices
    /// which are at a distance smaller than rad to the input point.
    pointIndexArr neighborhood(point_t const& pt, double const& rad);

    /// Get the points that are at a distance to the input point which is
    /// smaller than the input radius.
    ///
    /// @param pt input point.
    /// @param rad input radius.
    ///
    /// @returns a vector containing the points which are at a distance smaller
    /// than rad to the input point.
    pointVec neighborhood_points(point_t const& pt, double const& rad);

    /// Get the indices of points that are at a distance to the input point
    /// which is smaller than the input radius.
    ///
    /// @param pt input point.
    /// @param rad input radius.
    ///
    /// @returns a vector containing the indices of the points which are at a
    /// distance smaller than rad to the input point.
    indexArr neighborhood_indices(point_t const& pt, double const& rad);

  private:
    KDNodePtr make_tree(pointIndexArr::iterator const& begin,
                        pointIndexArr::iterator const& end,
                        size_t const& level);

    void knearest_(KDNodePtr const& branch, point_t const& pt,
                   size_t const& level, size_t const& num_nearest,
                   std::list<std::pair<KDNodePtr, double>>& k_nearest_buffer);

    void node_query_(KDNodePtr const& branch, point_t const& pt,
                     size_t const& level, size_t const& num_nearest,
                     std::list<std::pair<KDNodePtr, double>>& k_nearest_buffer);

    // default caller
    KDNodePtr nearest_(point_t const& pt);

    void neighborhood_(KDNodePtr const& branch, point_t const& pt,
                       double const& rad2, size_t const& level,
                       pointIndexArr& nbh);

    KDNodePtr root_;
    KDNodePtr leaf_;
};