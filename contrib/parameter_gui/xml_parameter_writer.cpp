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

#include "xml_parameter_writer.h"


namespace dealii
{
  namespace ParameterGui
  {
    XMLParameterWriter::XMLParameterWriter(QTreeWidget *tree_widget)
                      : tree_widget(tree_widget)
    {
      xml.setAutoFormatting(true);						// enable auto-formatting
    }



    bool XMLParameterWriter::write_xml_file(QIODevice *device)
    {
      xml.setDevice(device);							// setup the output device
      xml.writeStartDocument();							// write the head <?xml ... ?>
      xml.writeStartElement("ParameterHandler");				// write the root element <ParameterHandler>
										// loop over the elements
      for (int i = 0; i < tree_widget->topLevelItemCount(); ++i)
        write_item(tree_widget->topLevelItem(i));				// and write the items

      xml.writeEndDocument()		;					// close the first element

      return true;
    }



    void XMLParameterWriter::write_item(QTreeWidgetItem *item)
    {
      QString tag_name = mangle(item->text(0));					// store the element name

      xml.writeStartElement(tag_name);						// and write <tag_name> to the file

      if (!item->text(1).isEmpty())						// if the "value"-entry of this item is not empty
        {									// we have a parameter
          xml.writeTextElement("value", item->data(1,Qt::EditRole).toString());
          xml.writeTextElement("default_value", item->text(2));			// and we write its values
          xml.writeTextElement("documentation", item->text(3));
          xml.writeTextElement("pattern", item->text(4));
          xml.writeTextElement("pattern_description", item->text(5));
        };

      for (int i = 0; i < item->childCount(); ++i)				// go over the childrens recursively
        write_item(item->child(i));

      xml.writeEndElement();							// write closing </tag_name>
    }



    QString XMLParameterWriter::mangle (const QString &s)
    {
      std::string  s_temp (s.toStdString()); 					// this function is copied from
										// the ParameterHandler class
      std::string u;								// and adapted to mangle QString
      u.reserve (s_temp.size());

      static const std::string allowed_characters
        ("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

				   // for all parts of the string, see
				   // if it is an allowed character or
				   // not
      for (unsigned int i=0; i<s_temp.size(); ++i)
        if (allowed_characters.find (s_temp[i]) != std::string::npos)
          u.push_back (s_temp[i]);
        else
          {
	    u.push_back ('_');
	    static const char hex[16]
	      = { '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
	    u.push_back (hex[static_cast<unsigned char>(s_temp[i])/16]);
	    u.push_back (hex[static_cast<unsigned char>(s_temp[i])%16]);
          }

      QString  v (u.c_str());

      return v;
    }
  }
}
