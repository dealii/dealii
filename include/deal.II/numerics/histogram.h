// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__histogram_h
#define __deal2__histogram_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector.h>
#include <vector>
#include <string>

DEAL_II_NAMESPACE_OPEN


/**
 * This class provides some facilities to generate 2d and 3d histograms.
 * It is used by giving it one or several data sets and a rule how to
 * break the range of values therein into intervals (e.g. linear spacing
 * or logarithmic spacing of intervals). The values are then sorted into
 * the different intervals and the number of values in each interval is
 * stored for output later. In case only one data set was given, the
 * resulting histogram will be a 2d one, while it will be a 3d one if
 * more than one data set was given. For more than one data set, the same
 * intervals are used for each of them anyway, to make comparison easier.
 *
 *
 * <h3>Ways to generate the intervals</h3>
 *
 * At present, the following schemes for interval spacing are implemented:
 * <ul>
 * <li> Linear spacing: The intervals are distributed in constant steps
 *    between the minimum and maximum values of the data.
 * <li> Logaritmic spacing: The intervals are distributed in constant
 *    steps between the minimum and maximum values of the logs of the values.
 *    This scheme is only useful if the data has only positive values.
 *    Negative and zero values are sorted into the leftmost interval.
 * </ul>
 *
 * To keep programs extendible, you can use the two functions
 * @p get_interval_spacing_names and @p parse_interval_spacing, which always
 * give you a complete list of spacing formats presently supported and are
 * able to generate the respective value of the @p enum. If you use them,
 * you can write your program in a way such that it only needs to be
 * recompiled to take effect of newly added formats, without changing your
 * code.
 *
 *
 * <h3>Output formats</h3>
 *
 * At present, only GNUPLOT output is supported.
 *
 *
 * @ingroup textoutput
 * @author Wolfgang Bangerth, 1999
 */
class Histogram
{
public:
  /**
   * Definition of several ways to arrange
   * the spacing of intervals.
   */
  enum IntervalSpacing
  {
    linear, logarithmic
  };


  /**
   * Take several lists of values, each on
   * to produce one histogram that will
   * then be arrange one behind each other.
   *
   * Using several data sets at once allows
   * to compare them more easily, since
   * the intervals into which the data is
   * sorted is the same for all data sets.
   *
   * The histograms will be arranged such
   * that the computed intervals of the
   * <tt>values[i][j]</tt> form the x-range,
   * and the number of values in each
   * interval will be the y-range (for
   * 2d plots) or the z-range (for 3d
   * plots). For 3d plots, the @p y_values
   * parameter is used to assign each
   * data set a value in the y direction,
   * which is the depth coordinate in the
   * resulting plot. For 2d plots,
   * the @p y_values are ignored.
   *
   * If you give only one data set, i.e.
   * <tt>values.size()==1</tt>, then the
   * resulting histogram will be a 2d
   * one.
   *
   * @p n_intervals denotes the number of
   * intervals into which the data will be
   * sorted; @p interval_spacing denotes
   * the way the bounds of the intervals
   * are computed. Refer to the general
   * documentation for more information
   * on this.
   */
  template <typename number>
  void evaluate (const std::vector<Vector<number> > &values,
                 const std::vector<double>                   &y_values,
                 const unsigned int                           n_intervals,
                 const IntervalSpacing                        interval_spacing = linear);

  /**
   * This function is only a wrapper
   * to the above one in case you have
   * only one data set.
   */
  template <typename number>
  void evaluate (const Vector<number>          &values,
                 const unsigned int             n_intervals,
                 const IntervalSpacing          interval_spacing = linear);

  /**
   * Write the histogram computed by
   * the @p evaluate function to a
   * stream in a format suitable to
   * the GNUPLOT program. The function
   * generates 2d or 3d histograms.
   */
  void write_gnuplot (std::ostream &out) const;

  /**
   * Return allowed names for the interval
   * spacing as string. At present this
   * is "linear|logarithmic".
   */
  static std::string get_interval_spacing_names ();

  /**
   * Get a string containing one of the
   * names returned by the above function
   * and return the respective value of
   * @p IntervalSpacing. Throw an error
   * if the string is no valid one.
   */
  static IntervalSpacing parse_interval_spacing (const std::string &name);

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Exception.
   */
  DeclException0 (ExcEmptyData);
  /**
   * Exception.
   */
  DeclException0 (ExcInvalidIntervals);
  /**
   * Exception.
   */
  DeclException0 (ExcInvalidData);
  /**
   * Exception.
   */
  DeclException2 (ExcIncompatibleArraySize,
                  int, int,
                  << "The two array sizes " << arg1 << " and " << arg2
                  << " must match, but don't.");
  /**
   * Exception.
   */
  DeclException1 (ExcInvalidName,
                  std::string,
                  << "The given name <" << arg1
                  << "> does not match any of the known formats.");

private:
  /**
   * Structure denoting one
   * interval.
   */
  struct Interval
  {
    /**
     * Constructor. Sets the
     * bounds and sets the number of
     * values in this interval to zero.
     */
    Interval (const double left_point,
              const double right_point);

    /**
     * Determine an estimate for the
     * memory consumption (in bytes)
     * of this object.
     */
    std::size_t memory_consumption () const;

    /**
     * Left bound of the interval.
     */
    double       left_point;

    /**
     * Right bound of the interval.
     */
    double       right_point;

    /**
     * Number of values in this
     * interval.
     */
    unsigned int content;
  };

  /**
   * "Less-than" operation which
   * finds the minimal positive value
   * by sorting zero and negative value
   * to be larger than the largest positive
   * number. Used to find the lower bound
   * of the leftmost interval in the
   * logarithmic case interval spacing
   * scheme.
   *
   * Return @p true, if (<tt>n1<n2</tt>,
   * and (<tt>n1>0</tt> or <tt>n2<0</tt>)), or
   * (n2<n1 and n1>0 and n2<=0). This
   * in effect sorts all negativ
   * numbers to be larger than the
   * largest positive number.
   */
  template <typename number>
  static bool logarithmic_less (const number n1,
                                const number n2);

  /**
   * Vector holding one set of intervals
   * for each data set given to the
   * @p evaluate function.
   */
  std::vector<std::vector<Interval> > intervals;

  /**
   * Values for the depth axis of 3d
   * histograms. Stored in the @p evaluate
   * function.
   */
  std::vector<double>            y_values;
};


DEAL_II_NAMESPACE_CLOSE

#endif
