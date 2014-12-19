// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by Martin Steigemann and Wolfgang Bangerth
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


#include <QtGui>

#include "xml_parameter_reader.h"


namespace dealii
{
  namespace ParameterGui
  {
    XMLParameterReader::XMLParameterReader(QTreeWidget *tree_widget)
                      : tree_widget(tree_widget)
    {
      QStyle * style = tree_widget->style();

      subsection_icon.addPixmap(style->standardPixmap(QStyle::SP_DirClosedIcon), QIcon::Normal, QIcon::Off);
      subsection_icon.addPixmap(style->standardPixmap(QStyle::SP_DirOpenIcon), QIcon::Normal, QIcon::On);

      parameter_icon.addPixmap(style->standardPixmap(QStyle::SP_FileIcon));
    }



    bool XMLParameterReader::read_xml_file(QIODevice *device)
    {
      xml.setDevice(device);

		// We look for a StartElement "ParameterHandler"
		// and start parsing after this.
		//  <ParameterHandler>
		//   <subsection>
		//    ...
		//   </subsection>
		//  </ParameterHandler>

      while (xml.readNext() != QXmlStreamReader::Invalid)
        {
          if (xml.isStartElement())
            if (xml.name() == "ParameterHandler")
              {
                parse_parameters();

                return !xml.error();;
              };
        };

      xml.raiseError(QObject::tr("The file is not an ParameterHandler XML file."));

      return !xml.error();
    }



    QString XMLParameterReader::error_string() const
    {
      return QObject::tr("%1\nLine %2, column %3")
             .arg(xml.errorString())
             .arg(xml.lineNumber())
             .arg(xml.columnNumber());
    }



    void XMLParameterReader::parse_parameters()
    {
      Q_ASSERT(xml.isStartElement() && xml.name() == "ParameterHandler");

      while (xml.readNext() != QXmlStreamReader::Invalid)	// go to the next <start_element>
        {							// if it is the closing element of ParameterHandler,
          if (xml.isEndElement() && xml.name() == "ParameterHandler")
            break;						// break the loop

          if (xml.isStartElement())				// if it is a start element
            read_subsection_element(0);			// it must be a subsection or a parameter
        };
    }



    void XMLParameterReader::read_subsection_element(QTreeWidgetItem *parent)
    {
		// The structure of the parameter file is assumed to be of the form
		//
		//  <subsection>
		//    <subsection>
		//      ...
		//        <parameter>
		//          <value> ... </value>
		//          ...
		//          <pattern_description> ... </pattern_description>
		//        </parameter>
		//        <parameter>
		//        ...
		//        </parameter>
		//        ...
		//    </subsection>
		//    <subsection>
		//      ...
		//    </subsection>
		//    ...
		//  </subsection>
		//
		// Any subsection has a user-specified name also any parameter, but we do not know
		// the userspecified names and we can not assume anything. So, when parsing the file,
		// we do not know, if the actual <start_element> is a <subsection> or a <parameter>
		// in a subsection. To decide, if the element is a subsection- or a parameter-name,
		// we assume, that if the next <start_element> is <value>, we have a <parameter>
		// and a parameter has the entries <value>, <default_value>, <documentation>,
		// <pattern> and <pattern_description>

      Q_ASSERT(xml.isStartElement());					// the actual element is <subsection>

      QTreeWidgetItem * subsection = create_child_item(parent);		// create a new subsection in the tree

      subsection->setIcon(0, subsection_icon);				// set the icon,
      subsection->setText(0, demangle(xml.name().toString()));		// the name

      tree_widget->setItemExpanded(subsection, 0);			// and the folder is not expanded

      while (xml.readNext() != QXmlStreamReader::Invalid)		// read the next element
        {
          if (xml.isEndElement())					// if the next element is </subsection>, break the loop
            break;

          if (xml.isStartElement())					// if it is a start element
            {
              if (xml.name() == "value")				// it can be <value>, then we have found a parameter,
                {
                  subsection->setFlags(subsection->flags() | Qt::ItemIsEditable);	// values can be edited,
                  read_parameter_element (subsection);
                }
              else							// or it can be a new <subsection>
                {
                  subsection->setFlags(subsection->flags() | Qt::NoItemFlags);		// subsections can not be edited,
                  read_subsection_element (subsection);
                };
            };
        };
    }



    void XMLParameterReader::read_parameter_element(QTreeWidgetItem *parent)
    {
      Q_ASSERT(xml.isStartElement() && xml.name() == "value");		// the actual element is <value>,
									// then we have found a parameter-item
      QString value = xml.readElementText();				// read the element text
      parent->setText(1, value);					// and store as text to the item
      parent->setIcon(0, parameter_icon);				// change the icon of parent

      while (xml.readNext() != QXmlStreamReader::Invalid)				// go to the next <start_element>
        {
          if (xml.isStartElement())
            {
              if (xml.isStartElement() && xml.name() == "default_value")		// if it is <default_value>
                {
                  QString default_value = xml.readElementText();			// store it
                  parent->setText(2, default_value);
                }
              else if (xml.isStartElement() && xml.name() == "documentation")		// if it is <documentation>
                {
                  QString documentation = xml.readElementText();			// store it
                  parent->setText(3, documentation);

                  if (!documentation.isEmpty())						// if there is a documentation,
                    {
                      parent->setToolTip(0, "Documentation: " + documentation);		// set Documentation as ToolTip for both columns
                      parent->setToolTip(1, "Documentation: " + documentation);
                      parent->setStatusTip(0, "Documentation: " + documentation);	// and as StatusTip for the first column also
                    };
                }
              else if (xml.isStartElement() && xml.name() == "pattern")			// if it is <pattern>
                {
                  QString pattern = xml.readElementText();				// store it as text
                  parent->setText(4, pattern);						// we only need this value
											// for writing back to XML later
                }
              else if (xml.isStartElement() &&  xml.name() == "pattern_description")	// if it is <pattern_description>
                {
                  QString pattern_description = xml.readElementText();			// store it as text
                  parent->setText(5, pattern_description);
											// show the type and default
											// in the StatusLine
                  parent->setStatusTip(1, "Type: " + pattern_description + "   Default: " + parent->text(2));

						// in order to store values as correct data types,
						// we check the following types in the pattern_description:

                  QRegExp  rx_string("\\b(Anything|FileName|DirectoryName|Selection|List|MultipleSelection)\\b"),	
                           rx_integer("\\b(Integer)\\b"),
                           rx_double("\\b(Float|Floating|Double)\\b"),
                           rx_bool("\\b(Selection true|false)\\b");

                  if (rx_string.indexIn (pattern_description) != -1)			// the type "Anything" or "Filename"
                    {
                      QString value = parent->text(1);					// store as a QString

                      parent->setData(1, Qt::EditRole, value);				// and set the data in the item
                      parent->setData(1, Qt::DisplayRole, value);
                    }
                  else if (rx_integer.indexIn (pattern_description) != -1)		// if the tpye is "Integer"
                    {
                      QString text = parent->text(1);

                      bool ok = true;

                      int value = text.toInt(&ok);					// we convert the string to int

                      if (ok)								// and store
                        {
                          parent->setData(1, Qt::EditRole, value);
                          parent->setData(1, Qt::DisplayRole, value);
                        }
                      else								// otherwise raise an error
                        xml.raiseError(QObject::tr("Cannot convert integer type to integer!"));
                    }
                  else if (rx_double.indexIn (pattern_description) != -1)		// the same with "Float"
                    {
                      QString text = parent->text(1);

                      bool ok  = true;

                      double value = text.toDouble(&ok);

                      if (ok)
                        {
                          parent->setData(1, Qt::EditRole, value);
                          parent->setData(1, Qt::DisplayRole, value);
                        }
                      else
                        xml.raiseError(QObject::tr("Cannot convert double type to double!"));
                    };

                  if (rx_bool.indexIn (pattern_description) != -1)				// and booleans
                    {
                      QRegExp  test(parent->text(1));

                      bool value = true;

                      if (test.exactMatch("true"))
                        value = true;
                      else if (test.exactMatch("false"))
                        value = false;
                      else
                        xml.raiseError(QObject::tr("Cannot convert boolen type to boolean!"));

                      parent->setText(1, "");						// this is needed because we use
                      parent->setData(1, Qt::EditRole, value);				// for booleans the standard
                      parent->setData(1, Qt::DisplayRole, value);				// delegate
                    };

                  break;									// and break the loop
                }
              else
                {									// if there is any other element, raise an error
                  xml.raiseError(QObject::tr("Incomplete or unknown Parameter!"));
                  break;								// and break the loop, here
                };									// we assume the special structure
            };										// of the parameter-file!
        };
    }



    QTreeWidgetItem *XMLParameterReader::create_child_item(QTreeWidgetItem *item)
    {
      QTreeWidgetItem * child_item;							// create a new child-item

      if (item)
        child_item = new QTreeWidgetItem(item);						// if item is not empty,
      else										// append the new item as a child
        child_item = new QTreeWidgetItem(tree_widget);					// otherwise create a new item
											// in the tree

      child_item->setData(0, Qt::DisplayRole, xml.name().toString());			// set xml.tag_name as data
      child_item->setText(0, xml.name().toString());					// set xml.tag_name as data

      return child_item;
    }



    QString XMLParameterReader::demangle (const QString &s)
    {
      std::string  s_temp (s.toStdString()); 		// this function is copied from the ParameterHandler class

      std::string u;
      u.reserve (s_temp.size());

      for (unsigned int i=0; i<s_temp.size(); ++i)
        if (s_temp[i] != '_')
          u.push_back (s_temp[i]);
        else
          {
            Q_ASSERT (i+2 < s_temp.size());

            unsigned char c = 0;
            switch (s_temp[i+1])
              {
                case '0':  c = 0 * 16;  break;
                case '1':  c = 1 * 16;  break;
                case '2':  c = 2 * 16;  break;
                case '3':  c = 3 * 16;  break;
                case '4':  c = 4 * 16;  break;
                case '5':  c = 5 * 16;  break;
                case '6':  c = 6 * 16;  break;
                case '7':  c = 7 * 16;  break;
                case '8':  c = 8 * 16;  break;
                case '9':  c = 9 * 16;  break;
	    case 'a':  c = 10 * 16;  break;
	    case 'b':  c = 11 * 16;  break;
	    case 'c':  c = 12 * 16;  break;
	    case 'd':  c = 13 * 16;  break;
	    case 'e':  c = 14 * 16;  break;
	    case 'f':  c = 15 * 16;  break;
	    default:
		  Q_ASSERT (false);
	  }
	switch (s_temp[i+2])
	  {
	    case '0':  c += 0;  break;
	    case '1':  c += 1;  break;
	    case '2':  c += 2;  break;
	    case '3':  c += 3;  break;
	    case '4':  c += 4;  break;
	    case '5':  c += 5;  break;
	    case '6':  c += 6;  break;
	    case '7':  c += 7;  break;
	    case '8':  c += 8;  break;
	    case '9':  c += 9;  break;
	    case 'a':  c += 10;  break;
	    case 'b':  c += 11;  break;
	    case 'c':  c += 12;  break;
	    case 'd':  c += 13;  break;
	    case 'e':  c += 14;  break;
	    case 'f':  c += 15;  break;
	    default:
		  Q_ASSERT (false);
	  }

    	u.push_back (static_cast<char>(c));

					 // skip the two characters
	    i += 2;
          }

      QString  v (u.c_str());

      return v;
    }
  }
}

