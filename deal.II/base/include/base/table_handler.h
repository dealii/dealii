/*----------------------------   table_handler.h     ---------------------------*/
/*      $Id$                 */
#ifndef __table_handler_H
#define __table_handler_H
/*----------------------------   table_handler.h     ---------------------------*/



#include <base/exceptions.h>

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include <fstream.h>



/**
 * Abstract base class for the #TableEntry# class. See there.
 * This class is not to be used by the user.
 *
 * @author Ralf Hartmann, 1999
 */
class TableEntryBase 
{
  public:
				     /**
				      * Constructor.
				      */
    TableEntryBase();
				     /**
				      * Virtual destructor.
				      */    
    virtual ~TableEntryBase();
    
				     /**
				      * Write the table entry as text.
				      */
    virtual void write_text (ostream &) const =0;

				     /**
				      * Write the table entry in tex format.
				      */
    virtual void write_tex (ostream &) const =0;
};


/**
 * A #TableEntry# stores the value of an table entry.
 * The value type of this table entry is arbitrary. For
 * a #TableEntry<typename value_type># with not common value
 * type overload the output function.
 * This class is not to be used by the user.
 *
 * For more detail see the #TableHandler# class.
 *
 * @author Ralf Hartmann, 1999
 */
template <typename value_type>
class TableEntry : public TableEntryBase 
{
  public:
				     /**
				      * Constructor.
				      */
    TableEntry(const value_type value);

				     /**
				      * Destructor.
				      */
    virtual ~TableEntry();

				     /**
				      * Returns the value of this
				      * table entry.
				      */
    value_type value() const;

				     /**
				      * Write the table entry as text.
				      */
    virtual void write_text (ostream &out) const;

				     /**
				      * Write the table entry in tex format.
				      */
    virtual void write_tex (ostream &out) const;

  private:
				     /**
				      * Stored value.
				      */
    const value_type val;
};

/**
 * The #TableHandler# stores #TableEntries# of arbitrary value type and 
 * writes the table as text or in tex format to a tex file. The value type
 * actually may vary from column to column as from row to row. 
 *
 * \subsection{Usage}
 *
 * The most important function is the template function 
 * #add_value(const string &key, const value_type value)#, that adds a column
 * with the name #key# to the table if this column does not yet exist and adds the 
 * value of #value_type# (e.g. unsigned int, double, string, ...) to this column.
 * After the table is complete there are different possibilities of output, e.g.
 * into a tex file with #write_tex(ofstream &file)# or as text with 
 * #write_text (ostream &out)#. Two (or more) columns may be merged into a
 * 'supercolumn' by twice (or multiple) calling #add_column_to_supercolumn(...)#,
 * see there. Additionally there is a function to set for each column 
 * the precision of the output of numbers, and there are several functions to
 * prescribe the format and the captions the columns are written with in tex mode.
 *
 * \subsection{Example}
 * This is a simple example demonstrating the usage of this class. The first column
 * includes the numbers i=1..n, the second 1^2...n^2, the third sqrt(1)...sqrt(n),
 * where the second and third columns are merged into one supercolumn 
 * with the superkey `squares and roots'. Additionally the first column is
 * aligned to the right (the default was `centered') and the precision of
 * the square roots are given to be 6 (instead of 4 as default).
 *
 * \begin{verbatim}
 * TableHandler table();
 *
 * for (unsigned int i=1; i<=n; ++i)
 *   {
 *     table.add_value("numbers", i);
 *     table.add_value("squares", i*i);
 *     table.add_value("square roots", sqrt(i));
 *   }
 *                                  // merge the second and third column
 * table.add_column_to_supercolumn("squares", "squares and roots");
 * table.add_column_to_supercolumn("square roots", "squares and roots");
 *
 *                                  // additional settings
 * table.set_tex_format("numbers", "r");
 * table.set_precision("square roots", 6);
 *
 *                                  // output
 * ofstream out_file("number_table.tex");
 * table.write_tex(out_file);
 * out_file.close();
 * \end{verbatim}
 *
 * @author Ralf Hartmann, 1999
 */
class TableHandler
{
  public:

				     /**
				      * Constructor.
				      */
    TableHandler();
    
				     /**
				      * Adds a column (if not yet existent) with
				      * the key #key#
				      * and adds the value of #value_type#
				      * to the column.
				      */
    template <typename value_type>
    void add_value(const string &key, const value_type value);
    
				     /**
				      * Creates a sypercolumn (if not yet
				      * existent) and includes column to it.
				      * The keys of the column and the supercolumn
				      * are key and superkey respectively.
				      * To merge two columns c1 and c2 to
				      * a supercolumn sc, hence call
				      * #add_column_to_supercolumn(c1,sc)# and
				      * #add_column_to_supercolumn(c2,sc)#.
				      *
				      * Concerning the order of the columns,
				      * the supercolumn replaces the first
				      * column that is added to the supercolumn.
				      * Within the supercolumn the order
				      * the order of output follows the order
				      * the columns are added to the supercolumn.
				      */
    void add_column_to_supercolumn(const string &key, const string &superkey);

				     /**
				      * Change the order of columns and
				      * supercolumns in the table.
				      */
    void set_column_order(const vector<string> &new_order);
    
				     /**
				      * Sets the output format of a column.
				      */
    void set_precision(const string &key, const unsigned int precision);

				     /**
				      * Sets the caption of the column #key#
				      * for tex output.
				      */
    void set_tex_caption(const string &key, const string &tex_caption);

				     /**
				      * Sets the caption the the supercolumn
				      * #superkey# for tex output.
				      */
    void set_tex_supercaption(const string &superkey, const string &tex_supercaption);

				     /**
				      * Sets the tex output format of a column.
				      * e.g. "c", "r", "l".
				      */
    void set_tex_format(const string &key, const string &format);

				     /**
				      * Write table as formatted text, e.g.
				      * to the standard output.
				      */
    void write_text (ostream &out) const;

				     /**
				      * Write table as a tex file.
				      */
    void write_tex (ofstream &file) const;
   
    				     /**
				      * Exception
				      */
    DeclException1 (ExcColumnNotExistent,
		    string,
		    << "Column <" << arg1 << "> does not exist.");

    				     /**
				      * Exception
				      */
    DeclException1 (ExcSuperColumnNotExistent,
		    string,
		    << "Supercolumn <" << arg1 << "> does not exist.");

    				     /**
				      * Exception
				      */
    DeclException1 (ExcColumnOrSuperColumnNotExistent,
		    string,
		    << "Column or supercolumn <" << arg1 << "> does not exist.");
    
    				     /**
				      * Exception
				      */
    DeclException4 (ExcWrongNumberOfDataEntries,
		    string, int, string, int,
		    << "Column <" << arg1 << "> has got " << arg2
		    << "rows, but Column <" << arg3 << "> has got " << arg4 << ".");

    				     /**
				      * Exception
				      */
    DeclException1 (ExcUndefinedTexFormat,
		    string,
		    << "<" << arg1 << "> is not a tex column format. Use l,c,r.");

  protected:


    struct Column
    {
					 /**
					  * Constructor needed by stl_map.h
					  */
	Column();

					 /**
					  * Constructor.
					  */
	Column(const string &tex_caption);
	
					 /**
					  * Destructor.
					  */
	~Column();
	
	vector<TableEntryBase *> entries;
	
					 /**
					  * The caption of the column in tex output.
					  * By default, this is the key string that
					  * is given to the #TableHandler# by
					  * #TableHandler::add_value(...)#. This may
					  * be changed by calling
					  * #TableHandler::set_tex_caption(...)#.
					  */
	string tex_caption;

					 /**
					  * The column format in tex output.
					  * By default, this is #"c"#, meaning
					  * `centered'. This may
					  * be changed by calling
					  * #TableHandler::set_tex_format(...)#
					  * with #"c", "r", "l"# for centered,
					  * right or left.
					  */
	
	string tex_format;

					 /**
					  * Double or float entries are written with
					  * this (by the user set) precision.
					  * The default is 4.
					  */
	unsigned int precision;

					 /**
					  * Flag that may be used by derived classes.
					  */
	unsigned int flag;
    };

				     /**
				      * Help function that gives a vector of
				      * the keys of all columns that are mentioned
				      * in column_order, where each supercolumn key
				      * is replaced by its subcolumn keys.
				      *
				      * This function implicitly checks the
				      * consistency of the data.
				      */
    void get_selected_columns(vector<string> &sel_columns) const;
    
				     /**
				      * Builtin function, that checks if
				      * the number of rows is equal in every
				      * row. This function is e.g. called before
				      * writing output.
				      */
    unsigned int check_n_rows() const;

				     /**
				      * Stores the column and 
				      * supercolumn keys in the user wanted
				      * order.
				      * By default this is the order of
				      * adding the columns. This order may be
				      * changed by #set_column_order(...)#.
				      */
    vector<string> column_order;

				     /**
				      * Maps the column keys to the columns
				      * (not supercolumns).
				      */
    map<string,Column> columns;

				     /**
				      * Maps each supercolumn key to the
				      * the keys of its subcolumns in the right order.
				      * It is allowed that a supercolumn has got
				      * the same key as a column.
				      */
    map<string, vector<string> > supercolumns;

				     /**
				      * Maps the supercolumn keys to the
				      * captions of the supercolumns that
				      * are used in tex output.
				      *
				      * By default these are just the
				      * supercolumn keys but they may be changed
				      * by #set_tex_supercaptions(...)#.
				      */
    map<string, string> tex_supercaptions;
};





/*----------------------------   table_handler.h     ---------------------------*/
/* end of #ifndef __table_handler_H */
#endif
/*----------------------------   table_handler.h     ---------------------------*/
