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

#ifndef XMLPARAMETERWRITER_H
#define XMLPARAMETERWRITER_H

#include <QXmlStreamWriter>
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
 * The XMLParameterWriter class provides an interface to write parameters stored in a QTreeWidget to a file in XML format.
 * This class makes extensive use of the QXmlStreamWriter class, which implements the basic functionalities for writing
 * XML files.
 *
 * @note This class is used in the graphical user interface for the @ref ParameterHandler class.
 *       It is not compiled into the deal.II libraries and can not be used by applications using deal.II.
 *
 * @ingroup ParameterGui
 * @author Martin Steigemann, Wolfgang Bangerth, 2010
 */
    class XMLParameterWriter
    {
      public:
				     /**
				      * Constructor.
				      * Parameter values from @p tree_widget will be written.
				      */
        XMLParameterWriter (QTreeWidget *tree_widget);
				     /**
				      * This function writes the parameter values stored in <tt>tree_widget</tt>
				      * to @p device in XML format. We use the QXmlStreaWriter class
				      * for this. The root element is
				      * <code>&lt;ParameterHandler&gt;</code>
				      */
        bool write_xml_file (QIODevice *device);

      private:
				     /**
				      * This function writes a given @p item of <tt>tree_widget</tt>
				      * to a file in XML format. For this the QXmlStreamWriter class is used.
				      * If the @p item is a parameter, the elements that describes this parameter
				      * are written:
				      * @code
				      *   <value>value</value>
				      *   <default_value>default_value</default_value>
				      *   <documentation>documentation</documentation>
				      *   <pattern>pattern</pattern>
				      *   <pattern_description>[pattern_description]</pattern_description>
				      * @endcode
				      * If the @p item is a subsection, a start element <code>this_subsection</code> is written
				      * and <tt>write_item</tt> is called recursively to write the next <tt>item</tt>.
				      */
        void write_item (QTreeWidgetItem *item);
				     /**
				      * Reimplemented from the @ref ParameterHandler class.
				      * Mangle a string @p s so that it
				      * doesn't contain any special
				      * characters or spaces.
				      */
        QString  mangle (const QString &s);
				     /**
				      * An QXmlStreamWriter object
				      * which implements the functionalities
				      * we need for writing parameters to XML files.
				      */
        QXmlStreamWriter  xml;
				     /**
				      * A pointer to the QTreeWidget structure
				      * which stores the parameters.
				      */
        QTreeWidget * tree_widget;
    };
  }
/**@}*/
}


#endif
