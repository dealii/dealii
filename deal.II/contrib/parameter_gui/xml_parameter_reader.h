//---------------------------------------------------------------------------
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef XMLPARAMETERREADER_H
#define XMLPARAMETERREADER_H

#include <QIcon>
#include <QXmlStreamReader>
#include <QTreeWidget>
#include <QTreeWidgetItem>


namespace dealii
{
/*! @addtogroup ParameterGui
 *@{
 */
  namespace ParameterGui
  {
/**
 * The XMLParameterReader class provides an interface to parse parameters from XML files to a QTreeWidget.
 * This class makes extensive use of the QXmlStreamReader class, which implements the basic functionalities
 * for parsing XML files.
 *
 * @note This class is used in the graphical user interface for the @ref ParameterHandler class.
 *       It is not compiled into the deal.II libraries and can not be used by applications using deal.II.
 *
 * @ingroup ParameterGui
 * @author Martin Steigemann, Wolfgang Bangerth, 2010
 */
    class XMLParameterReader
    {
      public:
				     /**
				      * Constructor.
				      * The parameter values will be stored in @p tree_widget.
				      */
        XMLParameterReader (QTreeWidget *tree_widget);
				     /**
				      * This function reads the parameters from @p device into the <tt>tree_widget</tt>.
				      * We use the QXmlStreaReader class for this.
				      * There must be a start element
				      * <code>&lt;ParameterHandler&gt;</code>
				      * and an end element <code>&lt;/ParameterHandler&gt;</code>
				      * otherwise an exception is thrown.
				      */
        bool read_xml_file (QIODevice *device);
				     /**
				      * This function returns an error message.
				      */
        QString error_string () const;

      private:
				     /**
				      * This function implements a loop over the XML file
				      * and parses XML elements. It calls @ref read_subsection_element
				      * till the <code>&lt;/ParameterHandler&gt;</code> element is found
				      * or the end of the file is reached. In this case, an exception is thrown.
				      */
        void parse_parameters ();
				     /**
				      * This functions parses a <tt>subsection</tt>.
				      * and adds it as a child to @p parent.
				      * If the next element is <code>&lt;value&gt;</code>,
				      * this functions calls @ref read_parameter_element
				      * otherwise the function itself recursively.
				      */
        void read_subsection_element (QTreeWidgetItem *parent);
				     /**
				      * This function parses a <tt>parameter</tt> and
				      * and adds it as a child to @p parent.
				      * A <tt>parameter</tt> description consists of five elements:
				      * @code
				      *   <value>value</value>
				      *   <default_value>default_value</default_value>
				      *   <documentation>documentation</documentation>
				      *   <pattern>pattern</pattern>
				      *   <pattern_description>[pattern_description]</pattern_description>
				      * @endcode
				      * If a <tt>parameter</tt> description is incomplete, an exception
				      * is thrown.
				      */
        void read_parameter_element (QTreeWidgetItem *parent);
				     /**
				      * Reimplemented from the @ref ParameterHandler class.
				      * Unmangle a string @p s into its original form.
				      */
        QString  demangle (const QString &s);
				     /**
				      * This helper function creates a new child of @p item in the tree.
				      */
        QTreeWidgetItem * create_child_item(QTreeWidgetItem *item);
				     /**
				      * The QXmlStreamReader object for reading XML elements.
				      */
        QXmlStreamReader  xml;
				     /**
				      * A pointer to the tree structure.
				      */
        QTreeWidget * tree_widget;
				     /**
				      * An icon for subsections in the tree structure.
				      */
        QIcon  subsection_icon;
				     /**
				      * An icon for parameters in the tree structure.
				      */
        QIcon  parameter_icon;
    };
  }
/**@}*/
}


#endif
