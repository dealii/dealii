
#ifndef XMLPARAMETERREADER_H
#define XMLPARAMETERREADER_H

#include <QIcon>
#include <QXmlStreamReader>
#include <QTreeWidget>
#include <QTreeWidgetItem>


/**
 * The XMLParameterReader class provides an interface to parse parameters from XML files to a QTreeWidget.
 * This class makes extensive use of the QXmlStreamReader class, which implements the basic functionalities
 * for parsing XML files.
 *
 * @ingroup gui
 * @author Martin Steigemann, Wolfgang Bangerth, December 2010
 */
class XMLParameterReader
{
  public:
				     /**
				      * Constructor.
				      */
    XMLParameterReader (QTreeWidget *tree_widget);
				     /**
				      * This function reads the parameters from <tt>device</tt> into the <tt>tree_widget</tt>.
				      * We use the QXmlStreaReader class for this.
				      * There must be a start element
				      * <code><ParameterHandler></code>
				      * and an end element <code></ParameterHandler></code>
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
				      * till the <code></ParameterHandler></code> element is found
				      * or the end of the file is reached. In this case, an exception is thrown.
				      */
    void parse_parameters ();
				     /**
				      * This functions parses a <tt>subsection</tt>.
				      * If the next element is <code><value></code>,
				      * this functions calls @ref read_parameter_element
				      * otherwise the function itself recursively.
				      */
    void read_subsection_element (QTreeWidgetItem *parent);
				     /**
				      * This functions parses a <tt>parameter</tt>.
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
				      * Reimplemented from the ParameterHandler class.
				      * Unmangle a string into its original form.
				      */
    QString  demangle (const QString &s);
				     /**
				      * This helper function creates a new item in the tree.
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

#endif
