/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      

    (2006) kraftche@cae.wisc.edu    
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file QualityAssessor.hpp

Header file for the Mesquite::QualityAssessor class

  \author Thomas Leurent
  \date   2002-05-01
  \author Jason Kraftcheck
  \date   2005-03-09
 */


#ifndef MSQ_QUALITYASSESSOR_HPP
#define MSQ_QUALITYASSESSOR_HPP

#include "Mesquite.hpp"
#include "Instruction.hpp"
#include "MeshInterface.hpp"

#include <list>
#include <string>
#include <vector>
#include <iosfwd>

namespace MESQUITE_NS 
{

   class QualityMetric;
   class MsqError;
   class ParallelHelper;

  /*! \class QualityAssessor

      \brief A QualityAssessor instance can be inserted into an 
      InstructionQueue to calculate and summarize registered
      QualityMetric or QualityMetrics for the mesh.

      QualityAssessor provides a summary of the mesh quality.  An 
      instance of QualityAssessor may be inserted in the InstructionQueue
      at any point to print a summary of the mesh quality at that
      time during the optimization.  The application is expected to
      register QualityMetric instances to be used in assessing the
      mesh quality.  If no QualityMetrics are registered, the only
      assessment that will be performed is a simple count of inverted
      elements.

      The "stopping assessor" and "stopping function", if set,
      determinte the value reported to Mesquite for the overall
      run of of the QualityAssessor.
      
      All summary data except the histogram and custom power-means is 
      accumulated for all registered metrics.  QualityAssessment data
      can be obtained through three different mechanisms:
        - QualityAssessor can print a summary to std::out or a specified
          output stream after they are calculated.
        - The get_results and get_all_results methods can be used to
          obtain the sumary data programatically.
        - Quality measures for element-based or vertex-based metrics can
          be stored on the corresponding entities using "tags".
      
  */
  class QualityAssessor : public Instruction
  {
  public:
    
    /**\brief Simple consturctor.  Metrics registered separately.
     *\param print_summary_to_std_out If true, summary of mesh quality
     *                will be written to std::out.  If false, quality
     *                assessment will be available via the get_results
     *                and get_all_results methods, but will not be printed.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( bool print_summary_to_std_out = true,
                     bool free_elements_only = true,
                     const char* inverted_element_tag_name = 0,
                     std::string name = std::string() );

    /**\brief Simple consturctor.  Metrics registered separately.
     *\param output_stream IO stream to which to write a summary of the 
     *                mesh quality.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( std::ostream& output_stream,
                     bool free_elements_only = true,
                     const char* inverted_element_tag_name = 0,
                     std::string name = std::string() );
                     
    /**\brief Construct and register a QualityMetric
     *\param output_stream IO stream to which to write a summary of the 
     *                mesh quality.
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param historgram_intervals If non-zero, a histogram of quality metric
     *                values composed of the specified number of intervals
     *                will be generated.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( std::ostream& output_stream,
                     QualityMetric* metric, 
                     int histogram_intervals = 0,
                     double power_mean = 0.0,
                     bool free_elements_only = true,
                     const char* metric_value_tag_name = 0,
                     const char* inverted_element_tag_name = 0,
                     std::string name = std::string() );

                     
    /**\brief Construct and register a QualityMetric
     *\param print_summary_to_std_out If true, summary of mesh quality
     *                will be written to std::out.  If false, quality
     *                assessment will be available via the get_results
     *                and get_all_results methods, but will not be printed.
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param historgram_intervals If non-zero, a histogram of quality metric
     *                values composed of the specified number of intervals
     *                will be generated.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param free_elements_only If true, only quality values that depend
     *                on at least one free vertex will be uses.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     *\param  name    Name to include in output.  Useful if several QualityAssessors
     *                are in use at the same time.
     */
    MESQUITE_EXPORT 
    QualityAssessor( QualityMetric* metric, 
                     int histogram_intervals = 0,
                     double power_mean = 0.0,
                     bool free_elements_only = true,
                     const char* metric_value_tag_name = 0,
                     bool print_summary_to_std_out = true,
                     const char* inverted_element_tag_name = 0,
                     std::string name = std::string() );

    MESQUITE_EXPORT virtual ~QualityAssessor();
    
      //! Provides a name to the QualityAssessor (use it for default name in constructor).
    MESQUITE_EXPORT void set_name(std::string name) { qualityAssessorName = name; };
      //! Retrieves the QualityAssessor name. A default name should be set in the constructor.
    MESQUITE_EXPORT virtual std::string get_name() const { return qualityAssessorName; }
    
      /**\brief All elements or only improvable ones.
       *
       * If set to true, the quality assessment results will include
       * quality values only for those elements (or more precisely metric
       * sample points) which are influenced by at least one free vertex.
       * That is, quality for elements (or sample points) that the sovler
       * cannot improve (e.g. an element with all vertices fixed) will not
       * be included in the quality assessment.
       *
       * If set to false, quality for all elements will be assessed.
       */
    MESQUITE_EXPORT void measure_free_samples_only( bool yesno )
      { skipFixedSamples = yesno; }
      
      /**\brief All elements or only improvable ones.
       *
       * If set to true, the quality assessment results will include
       * quality values only for those elements (or more precisely metric
       * sample points) which are influenced by at least one free vertex.
       * That is, quality for elements (or sample points) that the sovler
       * cannot improve (e.g. an element with all vertices fixed) will not
       * be included in the quality assessment.
       *
       * If set to false, quality for all elements will be assessed.
       */
    MESQUITE_EXPORT bool measuring_free_samples_only() const
      { return skipFixedSamples; }
    
    
    /**\brief Register a QualityMetric for use in quality assessment.
     *
     * Add a quality metric to the list of metrics used to assess
     * the quality of the mesh.
     *
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param historgram_intervals If non-zero, a histogram of quality metric
     *                values composed of the specified number of intervals
     *                will be generated.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     */
    MESQUITE_EXPORT 
    void add_quality_assessment( QualityMetric* metric,
                                 int histogram_intervals = 0,
                                 double power_mean = 0.0,
                                 const char* metric_value_tag_name = 0,
                                 const char* metric_label = 0 );
    
    /**\brief Same as add_quality_assessment, except that the average
     *        metric value is also used as the return value from loop_over_mesh.
     *        Specify a power_mean value to control which average is used.
     */
    MESQUITE_EXPORT 
    void set_stopping_assessment( QualityMetric* metric,
                                  int histogram_intervals = 0,
                                  double power_mean = 0.0,
                                  const char* metric_value_tag_name = 0,
                                  const char* metric_label = 0 );

    /**\brief Register a QualityMetric for use in quality assessment.
     *
     * Add a quality metric to the list of metrics used to assess
     * the quality of the mesh.  Specify more parameters controlling
     * histogram.
     *
     *\param metric   QualtiyMetric to register for use in assessing mesh
     *                quality.  Will also be used as the "stopping assessment".
     *\param min      Minimum of histogram rnage.
     *\param max      Maximum of histogram range.
     *\param intervals Histogram intervals.
     *\param power_mean If non-zero, in addition to the normal summary
     *                statistics for the quality metric, an additional
     *                general power mean with the specified power will
     *                be calculated.
     *\param metric_value_tag_name  If a non-null value is specified,
     *                a tag with the specified name will be populated
     *                with quality values for individual elements or
     *                vertices if metric is an element-based or vertex-
     *                based metric.  If metric is not element-based or
     *                vertex-based, this argument has no effect.
     *\param  inverted_element_tag_name  If a non-null value is specified,
     *                an integer tag with the specified name will be used
     *                it store value of 0 for normal elements and 1 for
     *                inverted elements.
     */
    MESQUITE_EXPORT 
    void add_histogram_assessment( QualityMetric* qm, 
                                   double min, 
                                   double max,
                                   int intervals,
                                   double power_mean = 0.0,
                                   const char* metric_value_tag_name = 0,
                                   const char* metric_label = 0 );
    
    virtual MESQUITE_EXPORT
    void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                           const Settings* settings,
                           MsqError& err );
    
      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    virtual MESQUITE_EXPORT
    double loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                           const Settings* settings,
                           MsqError &err);

      //! Does one sweep over the mesh and assess the quality with the metrics previously added.
    virtual MESQUITE_EXPORT
    double loop_over_mesh( ParallelMesh* mesh,
                           MeshDomain* domain,
                           const Settings* settings,
                           MsqError &err);

      //! Do not print results of assessment.
    MESQUITE_EXPORT void disable_printing_results()
       {
         printSummary = false;
       }
      
      //! Print accumulated summary data to specified stream. 
    MESQUITE_EXPORT void print_summary( std::ostream& stream ) const;
    
      //! True if any metric evaluated to an invalid value
      //! for any element
    MESQUITE_EXPORT bool invalid_elements() const;

      //! Provides the number of inverted elements, inverted_elmes,
      //!  and the number of elements whose orientation can not be
      //!  determined, undefined_elems.
      //! Returns false if this information is not yet available.
      //! Returns true, otherwise.
    MESQUITE_EXPORT bool get_inverted_element_count(int &inverted_elems,
                                                    int &inverted_samples,
                                                    MsqError &err);
    
      //! Reset calculated data 
    MESQUITE_EXPORT void reset_data();

      //! Produces two historgrams on a single scale from a before
      //! optimization histogram and an after optimization histogram.
      //! The histogram intervals are adjusted to include the enitre
      //! range of both histograms, the horizontal interval value bars 
      //! are adjusted to be on the same scale, and the metric quality 
      //! values are placed in the correct quality value 'bin' based
      //! on the new interval scale.

    MESQUITE_EXPORT void scale_histograms(QualityAssessor* optimal);
    
    MESQUITE_EXPORT void tag_inverted_elements( std::string tagname ) 
      { invertedTagName = tagname; }
    MESQUITE_EXPORT void dont_tag_inverted_elements()
      { invertedTagName.clear(); }
    MESQUITE_EXPORT bool tagging_inverted_elements() const
      { return !invertedTagName.empty(); }
    
    MESQUITE_EXPORT void tag_fixed_elements( std::string tagname ) 
      { fixedTagName = tagname; }
    MESQUITE_EXPORT void dont_tag_fixed_elements()
      { fixedTagName.clear(); }
    MESQUITE_EXPORT bool tagging_fixed_elements() const
      { return !fixedTagName.empty(); }
                                                      
    /** \brief Per-metric QualityAssessor data
     *
     * The Assessor class holds QualityAssessor data for
     * each metric added by the calling application, including
     * a pointer to the metric instance, QAFunction flags
     * dictating what is to be calculated and output, histogram
     * parameters, and the variables used to accumulate results as
     * the QualityAssessor is running.  It also provides 
     * methods to access the calculated data once the QualityAssessor
     * pass is completed.
     */
    class Assessor
    {
      public:
      
        MESQUITE_EXPORT Assessor( QualityMetric* metric, const char* name = 0 );
        
        MESQUITE_EXPORT double get_average() const ;
        MESQUITE_EXPORT double get_maximum() const { return maximum; }
        MESQUITE_EXPORT double get_minimum() const { return minimum; }
        MESQUITE_EXPORT double get_rms()     const ;
        MESQUITE_EXPORT double get_stddev()  const ;
        MESQUITE_EXPORT double get_power_mean() const;
        MESQUITE_EXPORT double get_power() const   { return pMean; } 
        MESQUITE_EXPORT int get_count() const { return count; }
        
        MESQUITE_EXPORT bool have_power_mean() const { return 0.0 != pMean; }
        
        MESQUITE_EXPORT int get_invalid_element_count() const { return numInvalid; }
        
        /** Get historgram of data, if calculated.
         *\param lower_bound_out  The lower bound of the histogram
         *\param upper_bound_out  The upper bound of the histogram
         *\param counts_out       An array of counts of elements where
         *              the first entry is the number of elements for
         *              which the metric is below the lower bound, the
         *              last entry is the number of elements above the
         *              upper bound, and all other values are the counts
         *              for histogram intervals between the lower and
         *              upper bounds.
         */
        MESQUITE_EXPORT 
        void get_histogram( double& lower_bound_out,
                            double& upper_bound_out,
                            std::vector<int>& counts_out,
                            MsqError& err ) const;
                            
        /** Reset all calculated data */
        MESQUITE_EXPORT void reset_data();
       
        /** Print the histogram */
        MESQUITE_EXPORT void print_histogram( std::ostream&, int width = 0 ) const;

        /** Get the QualityMetric */
        MESQUITE_EXPORT QualityMetric* get_metric() const { return qualMetric; }
        
        MESQUITE_EXPORT const std::string& get_label() const { return mLabel; }
        
        /** Add a value to the running counts */
        MESQUITE_EXPORT void add_value( double metric_value );
        
        /** Add a value to the hisogram data */
        MESQUITE_EXPORT void add_hist_value( double metric_value );
        
        /** Note invalid result */
        MESQUITE_EXPORT void add_invalid_value() ;
        
        MESQUITE_EXPORT bool have_histogram() const
          { return !histogram.empty(); }
        
        /** If range of histogram has not yet been determined,
          * calculate it from the min/max values.
          */
        MESQUITE_EXPORT void calculate_histogram_range();
        
        MESQUITE_EXPORT bool write_to_tag() const { return !tagName.empty(); }
        
        MESQUITE_EXPORT void set_stopping_function( bool value )
          { stoppingFunction = value; }
        
        MESQUITE_EXPORT bool stopping_function( ) const
          { return stoppingFunction; }
        
        MESQUITE_EXPORT double stopping_function_value() const;

       
      private:
      
        friend class QualityAssessor;
        
        QualityMetric *const qualMetric; //< The quality metric
        std::string mLabel;
        
        unsigned long count;  //< The total number of times the metric was evaluated
        
        double sum;       //< The sum of the metric over all elements
        double maximum;   //< The maximum of the metric
        double minimum;   //< The minimum value of the metric
        double sqrSum;    //< The sum of the square of the metric values
        double pSum;      //< The sum of the metric values raised to the pMean power
        unsigned long numInvalid;  //< Count of invalid metric values
        
        double pMean;     //< Power for general power-mean.
            
        /** The histogram counts, where the first and last values are
         * counts of values below the lower bound and above the upper
         * bound, respectively.  The remaining values are the histogram
         * counts.
         */
        bool haveHistRange;
        double histMin;   //< Lower bound of histogram
        double histMax;   //< Upper bound of histogram
        std::vector<int> histogram;
        
        std::string tagName; //< Tag to which to write metric values.
        /** Cached tag handle */
        TagHandle tagHandle;
        
        /** Value is return value for all of QualityAssessor */
        bool stoppingFunction;
        
        int referenceCount;

        enum AssessSchemes assessScheme;

     };    
        
    typedef std::list<Assessor*> list_type;
        
    /** \brief Request summary data for a specific QualityMetric 
     * This method allows the application to request the summary
     * data for a metric it has registered with the QualityAssessor.
     * If the passed QualityMetric has not been registered with the
     * QualityAssessor instance, NULL is returned.
     */
    MESQUITE_EXPORT
    const Assessor* get_results( QualityMetric* metric ) const;
        
    /** \brief Request summary data for a specific QualityMetric 
     * This method allows the application to request the summary
     * data for a metric it has registered with the QualityAssessor.
     * If the passed QualityMetric has not been registered with the
     * QualityAssessor instance, NULL is returned.
     */
    MESQUITE_EXPORT
    const Assessor* get_results( const char* metric_name ) const;
    
    /** \brief Get list of all summary data.
     *  Return a const reference to the internal list of 
     *  calculated data.
     */
   const list_type& get_all_results() const
      { return assessList; }
   
   MESQUITE_EXPORT
   QualityAssessor( const QualityAssessor& copy );
   
   MESQUITE_EXPORT
   QualityAssessor& operator=( const QualityAssessor& copy );
   
  private:
    
    struct Data {
        /** Count of inverted elements. */
      int invertedElementCount;

        /** Count of inverted Jacobians within elements.*/
      int invertedSampleCount;

        /** Number of elements */
      size_t elementCount;

        /** Total number of element Jacobians tested for inversion */
      size_t sampleCount;

        /** Number of elements with at least one free vertex */
      size_t freeElementCount;
      
      int referenceCount;
      
      Data() : referenceCount(1) { clear(); }
      
      void clear() {
        invertedElementCount = invertedSampleCount = -1;
        elementCount = sampleCount = freeElementCount = 0;
      }
    };
  
    Data* myData;
  

      //! Common code for serial and parallel loop_over_mesh
    double loop_over_mesh_internal( MeshDomainAssoc* mesh_and_domain,
                                    const Settings* settings,
                                    ParallelHelper* helper,
                                    MsqError &err);
  
    /** Find an Assessor corresponding to the passed
     *  QualityMetric, or create it if is not found in
     *  the list.
     */
    list_type::iterator find_or_add( QualityMetric* qm, const char* label = 0 );

    /** Find an Assessor corresponding to the passed
     *  QualityMetric, or create it if is not found in
     *  the list.
     */
    list_type::iterator find_stopping_assessment();
   
    /** Find or create tag */
    TagHandle get_tag( Mesh* mesh,
                       std::string name, 
                       Mesh::TagType type, 
                       unsigned size, 
                       MsqError& err );
   
    /** Try to determine if output stream is a terminal and if so, 
        the width of that terminal.  Returns zero if width cannot
        be determined. */
    int get_terminal_width() const;

    static std::string element_name_as_string(int enum_name);
 
    static double round_to_3_significant_digits(double number);


    /** Name */
    std::string qualityAssessorName;  
    
    /** List of quality metrics and corresponding data */
    list_type assessList;
   
    /** Stream to which to write summary of metric data */
    std::ostream& outputStream;
    /** Disable printing */
    bool printSummary;
    
    std::string invertedTagName, fixedTagName;
    
    bool skipFixedSamples;

    int elementTypeCount[MIXED - POLYGON+1];

    bool invalid_values;  // set to true when a target metric contains inverted elements

  };

  
} //namespace


#endif // QualityAssessor_hpp
